"""
The C generator for gau2grid collocation functions
"""

import os

from . import RSH
from . import c_pragma
from . import c_util_generator as c_util
from . import codegen
from . import order
from . import utility

_grad_indices = ["x", "y", "z"]
_hess_indices = ["xx", "xy", "xz", "yy", "yz", "zz"]


def generate_c_gau2grid(max_L, path=".", cartesian_order="row", spherical_order="gaussian", inner_block="auto", do_cf=True):
    """
    Creates the C files for the gau2grid program.

    Parameters
    ----------
    max_L : int
        The maximum angular momentum compiled for.
    path : str, optional
        The path to write the files to.
    cartesian_order : str, optional
        The order of the cartesian angular momentum
    inner_block : int, optional
        The size of the inner computation block.
    do_cf : bool, option
        Apply clang-format to the generated files or not.

    Returns
    -------
    None

    """

    cartesian_order = cartesian_order.lower()
    spherical_order = spherical_order.lower()

    # Build the codegen objects for each file
    gg_header = codegen.CodeGen(cgen=True)
    gg_phi = codegen.CodeGen(cgen=True)
    gg_grad = codegen.CodeGen(cgen=True)
    gg_hess = codegen.CodeGen(cgen=True)
    gg_spherical = codegen.CodeGen(cgen=True)
    gg_helper = codegen.CodeGen(cgen=True)
    gg_pragma = codegen.CodeGen(cgen=True)

    # Add license to header only
    c_util.write_license(gg_header)

    # Add general header comments
    for cgs in [gg_header, gg_phi, gg_grad, gg_hess, gg_spherical, gg_helper, gg_pragma]:

        cgs.write("/*", endl="")
        cgs.write(" * This is a Gau2Grid automatically generated C file.", endl="")
        cgs.write(" *", endl="")
        cgs.write(" * More details can found at the following repo:", endl="")
        cgs.write(" *   https://github.com/dgasmith/gau2grid", endl="")
        cgs.write(" */", endl="")
        cgs.blankline()

    # Write out the pragma header
    c_pragma.build_pragma_header(gg_pragma)

    gg_helper.write("#include <stdio.h>")

    # Add utility headers
    for cgs in [gg_phi, gg_grad, gg_hess, gg_spherical, gg_helper]:
        cgs.write("#include <math.h>")
        cgs.write("#include <stdio.h>")
        cgs.write("#ifdef _MSC_VER")
        cgs.write("#include <malloc.h>")
        cgs.write("#else")
        cgs.write("#include <mm_malloc.h>")
        cgs.write("#endif")
        cgs.blankline()
        cgs.write('#include "gau2grid.h"')
        cgs.write('#include "gau2grid_pragma.h"')
        cgs.blankline()

    # Header guards
    gg_header.write("#ifdef __cplusplus")
    gg_header.write('extern "C" {', endl="")
    gg_header.write("#endif")
    gg_header.blankline()
    gg_header.write("#ifndef GAU2GRID_GUARD_H")
    gg_header.write("#define GAU2GRID_GUARD_H")
    gg_header.blankline()

    gg_header.write('#include "gau2grid_pragma.h"')
    gg_header.blankline()

    # Add any information needed
    gg_helper.write("// Information helpers")
    gg_header.write("// Information helpers")

    # Maximum angular momentum
    gg_helper.write("int max_L() { return %d; }" % max_L, endl="")
    gg_helper.blankline()

    gg_header.write("int max_L()")
    gg_header.blankline()

    # Cartesian ordering
    gg_helper.write('const char* cartesian_order() { return "%s"; }' % cartesian_order, endl="")
    gg_helper.blankline()

    gg_header.write("const char* cartesian_order()")
    gg_header.blankline()

    # Spherical ordering
    gg_helper.write('const char* spherical_order() { return "%s"; }' % spherical_order, endl="")
    gg_helper.blankline()

    gg_header.write("const char* spherical_order()")
    gg_header.blankline()

    # Build out the spherical transformer

    gg_header.write("// Spherical transformers")
    for L in range(max_L + 1):
        sig = RSH.transformation_c_generator(gg_spherical, L, cartesian_order, spherical_order)
        gg_header.write(sig)
    gg_header.blankline()

    # Fast transformers
    gg_header.write("// Fast transposers")
    trans_sig = c_util.naive_transpose(gg_spherical)
    gg_header.write(trans_sig)
    fast_trans_sig = c_util.fast_transpose(gg_spherical, 8)
    gg_header.write(fast_trans_sig)
    gg_header.blankline()

    # Fast copiers
    gg_header.write("// Fast segment copiers")
    block_sig = c_util.block_copy(gg_spherical)
    gg_header.write(block_sig)

    # Loop over phi, grad, hess and build blocks for each
    gg_helper.write("// Collocation selector functions")
    helper_sigs = []
    for name, grad, cg in [("Phi", 0, gg_phi), ("Phi grad", 1, gg_grad), ("Phi Hess", 2, gg_hess)]:
        cg.blankline()
        gg_header.write("// %s computers" % name)
        cg.blankline()

        # Write out the phi builders
        sig_store = []
        for L in range(max_L + 1):
            sig = shell_c_generator(cg, L, grad=grad, cartesian_order=cartesian_order, inner_block=inner_block)
            sig_store.append(sig)
            cg.blankline()

            # Write out the header data
            gg_header.write(sig)
            gg_header.blankline()

        # Write out the convenience functions
        func_name, conv_sig = sig_store[0].split("(")
        if "deriv" in func_name:
            func_name = func_name.replace("L0_", "")
        else:
            func_name = func_name.replace("_L0", "")
        func_name += "(int L, "
        func_name += conv_sig
        helper_sigs.append(_make_call(func_name).split("(")[0])

        gg_header.write(func_name)
        gg_header.blankline()

        gg_helper.start_c_block(func_name)
        gg_helper.write("// Chooses the correct function for a given L")

        # Write out if's to choose the right L
        L = 0
        gg_helper.write("if (L == 0) {", endl="")
        for sig in sig_store:
            if L != 0:
                gg_helper.write("} else if (L == %d) {" % L, endl="")

            sig = _make_call(sig)
            gg_helper.write("    " + sig)
            L += 1

        # Handle exception
        gg_helper.write("} else {", endl="")
        gg_helper.write('    printf("Requested angular momentum exceeded compiled of %d\\n")' % max_L)
        gg_helper.write('    exit(0)')
        gg_helper.write("}", endl="")
        gg_helper.close_c_block()
        # print(func_name)

    # Finish header guard
    gg_header.write("#ifdef __cplusplus")
    gg_header.write("}", endl="")
    gg_header.write("#endif")
    gg_header.write("#endif /* GAU2GRID_GUARD_H */")

    # Write out the CG's to files
    gg_header.repr(filename=os.path.join(path, "gau2grid.h"), clang_format=do_cf)
    gg_phi.repr(filename=os.path.join(path, "gau2grid_phi.c"), clang_format=do_cf)
    gg_grad.repr(filename=os.path.join(path, "gau2grid_deriv1.c"), clang_format=do_cf)
    gg_hess.repr(filename=os.path.join(path, "gau2grid_deriv2.c"), clang_format=do_cf)
    gg_spherical.repr(filename=os.path.join(path, "gau2grid_spherical.c"), clang_format=do_cf)
    gg_helper.repr(filename=os.path.join(path, "gau2grid_helper.c"), clang_format=do_cf)
    gg_pragma.repr(filename=os.path.join(path, "gau2grid_pragma.h"))


def shell_c_generator(cg, L, function_name="", grad=0, cartesian_order="row", inner_block="auto"):

    # Grab the line start
    cg_line_start = len(cg.data)
    deriv_indices = utility.get_deriv_indices(grad)

    # Parse Keywords
    if function_name == "":
        if grad == 0:
            function_name = "gg_collocation_L%d" % L
        else:
            function_name = "gg_collocation_L%d_deriv%d" % (L, grad)

    if grad > 2:
        raise TypeError("Only grad <=2 (Hessians) is supported")

    # Set a few parameters for custom loops
    L_needs_out = False

    # Precompute temps
    ncart = int((L + 1) * (L + 2) / 2)
    nspherical = L * 2 + 1

    # Do we do multiple loops for each tmp or just one at a time?
    paritioned_loops = False
    if (grad == 1) and (L >= 3):
        paritioned_loops = True
    elif (grad == 2) and (L >= 2):
        paritioned_loops = True

    # Handle inner block, everything should fit into ~50% of L1
    # L1 is roughly 32K for data so lets say 16k max or 2048 doubles
    if inner_block == "auto":
        cache_limit_doubles = 2048

        # Basic temps + grad temps
        basic_lines = 6 + grad

        if paritioned_loops:
            # If we run partitioned loops we need this many lines
            nlines = basic_lines + ncart
        else:
            # If we run a single loop we need this many lines
            nlines = basic_lines + ncart * (1 + len(deriv_indices))

        # This could be bad when we hit AVX-512 (soon)
        if nlines * 64 > cache_limit_doubles:
            inner_block = 32
        else:
            inner_block = 64

        if nlines * 32 > cache_limit_doubles:
            print(
                "WARNING: For L=%2d and grad=%d assumed 16,384B L1 cache limit will be exceeded. This may impact performance."
                % (L, grad))

    elif isinstance(inner_block, int):
        pass
    else:
        raise ValueError("Inner block of name %s not understood" % str(inner_block))

    # Build function signature
    func_sig = "const unsigned long npoints, const double* PRAGMA_RESTRICT x, const double* PRAGMA_RESTRICT y, const double* PRAGMA_RESTRICT z, const int nprim, const double* PRAGMA_RESTRICT coeffs, const double* PRAGMA_RESTRICT exponents, const double* PRAGMA_RESTRICT center, const int spherical, double* PRAGMA_RESTRICT phi_out"

    # Add extra output vals for derivs
    for deriv in deriv_indices:
        func_sig += ", double* PRAGMA_RESTRICT phi_%s_out" % deriv

    func_sig = "void %s(%s)" % (function_name, func_sig)
    cg.start_c_block(func_sig)
    cg.blankline()

    # Figure out spacing
    cg.write("// Sizing")
    cg.write("unsigned long nblocks = npoints / %d" % inner_block)
    cg.write("nblocks += (npoints %% %d) ? 1 : 0" % inner_block)
    cg.write("const unsigned long ncart = %d" % ncart)
    cg.write("const unsigned long nspherical = %d" % nspherical)
    cg.write("const unsigned long nout = spherical ? nspherical : ncart")
    cg.blankline()

    # Build temporaries
    S_cache_tmps = ["xc", "yc", "zc", "R2", "S0"]
    if grad > 0:
        S_cache_tmps.append("S1")
    if grad > 1:
        S_cache_tmps.append("S2")

    block_malloc_name = "cache_data"
    block_malloc_sizes = [(name, inner_block) for name in S_cache_tmps]
    S_tmps = [block_malloc_name]

    # Allocate as single block on heap
    cg.write("// Allocate S temporaries, single block to stay on cache")
    _block_malloc(cg, block_malloc_name, block_malloc_sizes)

    cg.blankline()

    # Hold the expn1 and expn2 arrays
    cg.write("// Allocate exponential temporaries")
    exp_tmps = ["expn1"]
    if grad > 0:
        exp_tmps += ["expn2"]
    for tname in exp_tmps:
        cg.write(_malloc(tname, "nprim"))
    S_tmps.extend(exp_tmps)
    cg.blankline()

    # Figure out powers needed
    power_tmps = []
    if (L > 1) and paritioned_loops:
        cg.write("// Allocate power temporaries")
        power_tmps = ["xc_pow", "yc_pow", "zc_pow"]

        for tname in power_tmps:
            cg.write(_malloc(tname, inner_block * (L - 1)))

        cg.blankline()

    # Determine output tmps
    inner_tmps = []
    if L >= L_needs_out:
        cg.write("// Allocate output temporaries")

        inner_tmps = ["phi_tmp"]
        if paritioned_loops is False:
            for deriv in deriv_indices:
                inner_tmps.append("phi_%s_tmp" % deriv)

        # Malloc temps
        for tname in inner_tmps:
            cg.write(_malloc(tname, inner_block * ncart))
    cg.blankline()

    # Any declerations needed
    cg.write("// Declare doubles")
    cg.write("const double center_x = center[0]")
    cg.write("const double center_y = center[1]")
    cg.write("const double center_z = center[2]")
    cg.write("double A")
    if grad > 0:
        cg.write("double " + ", ".join("A%s" % grad.upper() for grad in _grad_indices))
    if grad > 1:
        cg.write("double " + ", ".join("A%s" % hess.upper() for hess in _hess_indices))
    cg.blankline()

    cg.write("// Build negative exponents")
    cg.start_c_block("for (unsigned long i = 0; i < nprim; i++)")
    cg.write("expn1[i] = -1.0 * exponents[i]")
    if grad > 0:
        cg.write("expn2[i] = -2.0 * exponents[i]")
    cg.close_c_block()
    cg.blankline()

    # Start outer loop
    cg.write("// Start outer block loop")
    cg.start_c_block("for (unsigned long block = 0; block < nblocks; block++)")
    cg.blankline()

    # Move data into inner buffers and compute R
    cg.blankline()
    cg.write("// Copy data into inner temps")
    cg.write("const unsigned long start = block * %d" % inner_block)
    cg.write("const unsigned long remain = ((start + %d) > npoints) ? (npoints - start) : %d" % (inner_block,
                                                                                                 inner_block))
    cg.blankline()

    cg.write("PRAGMA_VECTORIZE", endl="")
    cg.start_c_block("for (unsigned long i = 0; i < remain; i++)")
    cg.write("xc[i] = x[start + i] - center_x")
    cg.write("yc[i] = y[start + i] - center_y")
    cg.write("zc[i] = z[start + i] - center_z")

    cg.blankline()
    cg.write("// Distance")
    cg.write("R2[i] = xc[i] * xc[i]")
    cg.write("R2[i] += yc[i] * yc[i]")
    cg.write("R2[i] += zc[i] * zc[i]")

    cg.blankline()
    cg.write("// Zero out S tmps")
    cg.write("S0[i] = 0.0")
    if grad > 0:
        cg.write("S1[i] = 0.0")
    if grad > 1:
        cg.write("S2[i] = 0.0")

    cg.close_c_block()
    cg.blankline()

    # Start inner loop
    cg.write("// Start exponential block loop")
    cg.start_c_block("for (unsigned long n = 0; n < nprim; n++)")

    # Build R2
    cg.write("const double coef = coeffs[n]")
    cg.write("const double alpha_n1 = expn1[n]")
    if grad > 0:
        cg.write("const double alpha_n2 = expn2[n]")

    # Build out thoese gaussian derivs
    cg.blankline()
    cg.write("PRAGMA_VECTORIZE", endl="")
    cg.start_c_block("for (unsigned long i = 0; i < remain; i++)")
    cg.write("const double T1 = coef * exp(alpha_n1 * R2[i])")
    cg.write("S0[i] += T1")
    if grad > 0:
        cg.write("const double T2 = alpha_n2 * T1")
        cg.write("S1[i] += T2")
    if grad > 1:
        cg.write("const double T3 = alpha_n2 * T2")
        cg.write("S2[i] += T3")

    cg.close_c_block()
    cg.blankline()

    # Close off
    cg.close_c_block()
    cg.blankline()

    # Grab the inner line start
    inner_line_start = len(cg.data)
    inner_line_stop = inner_line_start + 1

    # Combine blocks

    if L == 0:
        cg.write("// Combine blocks")
        cg.write("PRAGMA_VECTORIZE", endl="")
        cg.start_c_block("for (unsigned long i = 0; i < remain; i++)")

        # Build out required S
        _S_tmps(cg, L, grad, inner_block)

        # Nothing else to be done. Copy it back to outs
        cg.write("phi_out[start + i] = S0[i]")

        if grad > 0:
            cg.blankline()
            cg.write("// Gradient AM=0 Component=0")
            cg.write("phi_x_out[start + i] = SX")
            cg.write("phi_y_out[start + i] = SY")
            cg.write("phi_z_out[start + i] = SZ")

        if grad > 1:
            cg.blankline()
            cg.write("// Hessian AM=0 Component=0")
            cg.write("phi_xx_out[start + i] = SXX")
            cg.write("phi_yy_out[start + i] = SYY")
            cg.write("phi_zz_out[start + i] = SZZ")
            cg.write("phi_xy_out[start + i] = SXY")
            cg.write("phi_xz_out[start + i] = SXZ")
            cg.write("phi_yz_out[start + i] = SYZ")

        cg.close_c_block()

        # Grab the inner line stop
        inner_line_stop = len(cg.data)

    elif paritioned_loops:

        cg.write("// Build powers")
        cg.write("PRAGMA_VECTORIZE", endl="")
        cg.start_c_block("for (unsigned long i = 0; i < remain; i++)")
        _power_tmps(cg, L, inner_block, array=True)
        cg.close_c_block()

        for dind in ["A"] + deriv_indices:
            _c_am_single_build(cg, L, cartesian_order, grad, inner_block, dind, array=True)

            dind = dind.lower()
            if dind == "a":
                dind = ""
            else:
                dind = "_" + dind

            # Transform
            cg.start_c_block("if (spherical)")
            sph_fnc = "gg_cart_to_spherical_L%d" % L
            cg.write("%s(remain, phi_tmp, %d, (phi%s_out + start), npoints)" % (sph_fnc, inner_block, dind))
            cg.write("} else {", endl="")
            cg.write("block_copy(nout, remain, phi_tmp, %d, (phi%s_out + start), npoints, 0)" % (inner_block, dind))
            cg.close_c_block()

            cg.blankline()

        # Grab the inner line stop
        inner_line_stop = len(cg.data)

    else:
        cg.write("// Combine blocks")
        cg.write("PRAGMA_VECTORIZE", endl="")
        cg.start_c_block("for (unsigned long i = 0; i < remain; i++)")

        # Build out required S
        _S_tmps(cg, L, grad, inner_block)

        # Build out required power temps if needed
        _power_tmps(cg, L, inner_block)

        # Contract temps with powers
        _c_am_full_build(cg, L, cartesian_order, grad, inner_block)

        cg.blankline()

        # End inner loop
        cg.close_c_block()

        # Grab the inner line stop
        inner_line_stop = len(cg.data)

        # Spherical/Cartesian copy out
        _tmp_to_out_copy(cg, L, deriv_indices, inner_block)

    # End outer loop
    cg.close_c_block()

    # Free up those arrays
    cg.blankline()
    for name, flist in [("S", S_tmps), ("Power", power_tmps), ("inner", inner_tmps)]:
        if len(flist) == 0: continue

        cg.write("// Free %s temporaries" % name)
        for tname in flist:
            cg.write("_mm_free(%s)" % tname)
        cg.blankline()

    # End function
    cg.close_c_block()

    # Clean up data, there are a few things easier to post-process

    # Remove any "[0 + i]"
    for x in range(cg_line_start, inner_line_stop):
        cg.data[x] = cg.data[x].replace("[0 + ", "[")

    if paritioned_loops is False:
        # Remove any "A = 1" just for the inner block
        rep_data = {}
        pos = inner_line_start
        while pos < inner_line_stop:
            line = cg.data[pos]
            #    print(line)

            # If we hit a Density line its an individual angular momentum, need to reset dict
            if ("Density" in line) or ("// Combine" in line):
                rep_data = {}
                pos += 1
                continue

            # Skip comments and blanklines
            if ("=" not in line) or ("//" in line) or ("double" in line):
                pos += 1
                continue

            # Find a single
            if (" = " in line) and ("*" not in line) and ("+" not in line) and ('/' not in line):
                key, data = line.replace(";", "").split(" = ")
                rep_data[key.strip()] = data.strip()
                cg.data.pop(pos)
                inner_line_stop -= 1
                continue

            for k, v in rep_data.items():
                tmp = line.split("= ")[1]
                if k + ";" in tmp:
                    cg.data[pos] = line.replace(k + ";", v + ";")
                elif k + " " in tmp:
                    cg.data[pos] = line.replace(k + " ", v + " ")
            pos += 1

    # Remove any " * 1"
    for x in range(cg_line_start, inner_line_stop):
        cg.data[x] = cg.data[x].replace(" * 1;", ";")
        cg.data[x] = cg.data[x].replace(" * 1.0;", ";")
        cg.data[x] = cg.data[x].replace("= 1 * ", "= ")
        cg.data[x] = cg.data[x].replace("= 1.0 * ", "= ")

    return func_sig


def _make_call(string):
    for rep in ["double* ", "bool ", "int ", "unsigned long ", "void ", "PRAGMA_RESTRICT "]:
        string = string.replace("const " + rep, "")
        string = string.replace(rep, "")
    return string


def _malloc(name, size, dtype="double"):
    # return "%s*  %s = (%s*)malloc(%s * sizeof(%s))" % (dtype, name, dtype, str(size), dtype)
    return "%s* PRAGMA_RESTRICT %s = (%s*)_mm_malloc(%s * sizeof(%s), 32)" % (dtype, name, dtype, str(size), dtype)


def _block_malloc(cg, block_name, mallocs, dtype="double"):
    tot_size = sum(x[1] for x in mallocs)
    cg.write(_malloc(block_name, tot_size))
    current_shift = 0
    for name, size in mallocs:
        cg.write("%s* PRAGMA_RESTRICT %s = %s + %d" % (dtype, name, block_name, current_shift))
        current_shift += size


def _c_am_single_build(cg, L, cartesian_order, grad, shift, specific_deriv, array=True):
    """
    Builds a unrolled angular momentum function
    """

    specific_deriv = specific_deriv.upper()
    cg.write("// Combine %s blocks" % specific_deriv)
    cg.write("PRAGMA_VECTORIZE", endl="")
    cg.start_c_block("for (unsigned long i = 0; i < remain; i++)")

    if specific_deriv == "X":
        cg.write("const double SX = S1[i] * xc[i]")
    elif specific_deriv == "Y":
        cg.write("const double SY = S1[i] * yc[i]")
    elif specific_deriv == "Z":
        cg.write("const double SZ = S1[i] * zc[i]")
    elif specific_deriv == "XY":
        cg.write("const double SX = S1[i] * xc[i]")
        cg.write("const double SY = S1[i] * yc[i]")
        cg.write("const double SXY = S2[i] * xc[i] * yc[i]")
    elif specific_deriv == "XZ":
        cg.write("const double SX = S1[i] * xc[i]")
        cg.write("const double SZ = S1[i] * zc[i]")
        cg.write("const double SXZ = S2[i] * xc[i] * zc[i]")
    elif specific_deriv == "YZ":
        cg.write("const double SY = S1[i] * yc[i]")
        cg.write("const double SZ = S1[i] * zc[i]")
        cg.write("const double SYZ = S2[i] * yc[i] * zc[i]")
    elif specific_deriv == "XX":
        cg.write("const double SX = S1[i] * xc[i]")
        cg.write("const double SXX = S2[i] * xc[i] * xc[i] + S1[i]")
    elif specific_deriv == "YY":
        cg.write("const double SY = S1[i] * yc[i]")
        cg.write("const double SYY = S2[i] * yc[i] * yc[i] + S1[i]")
    elif specific_deriv == "ZZ":
        cg.write("const double SZ = S1[i] * zc[i]")
        cg.write("const double SZZ = S2[i] * zc[i] * zc[i] + S1[i]")
    elif specific_deriv == "A":
        pass
    else:
        raise KeyError("Specific deriv %s not understood." % specific_deriv)
    cg.blankline()

    # Generator
    for idx, l, m, n in order.cartesian_order_factory(L, cartesian_order):

        l = l + 2
        m = m + 2
        n = n + 2
        ld1 = l - 1
        ld2 = l - 2
        md1 = m - 1
        md2 = m - 2
        nd1 = n - 1
        nd2 = n - 2

        # Set grads back to zero
        x_grad, y_grad, z_grad = False, False, False
        shift_idx = idx * shift

        name = "X" * ld2 + "Y" * md2 + "Z" * nd2
        if name == "":
            name = "0"

        # Gradient
        AX = _build_xyz_pow("AX", ld2, ld1, m, n, shift, array=array, rhs_only=True)
        x_grad = AX is not None
        AY = _build_xyz_pow("AY", md2, l, md1, n, shift, array=array, rhs_only=True)
        y_grad = AY is not None
        AZ = _build_xyz_pow("AZ", nd2, l, m, nd1, shift, array=array, rhs_only=True)
        z_grad = AZ is not None

        if specific_deriv == "A":
            rhs = _build_xyz_pow("A", 1.0, l, m, n, shift, array=array, rhs_only=True)
            cg.write("phi_tmp[%d + i] = %s * S0[i]" % (shift_idx, rhs))

            # Keep lines together
            continue

        # Gradients
        if specific_deriv == "X":
            rhs = _build_xyz_pow("A", 1.0, l, m, n, shift, array=array, rhs_only=True)
            cg.write("phi_tmp[%d + i] = %s * SX" % (shift_idx, rhs))

            if x_grad:
                cg.write("phi_tmp[%d + i] += %s * S0[i]" % (shift_idx, AX))

        if specific_deriv == "Y":
            rhs = _build_xyz_pow("A", 1.0, l, m, n, shift, array=array, rhs_only=True)
            cg.write("phi_tmp[%d + i] = %s * SY" % (shift_idx, rhs))

            if y_grad:
                cg.write("phi_tmp[%d + i] += %s * S0[i]" % (shift_idx, AY))

        if specific_deriv == "Z":
            rhs = _build_xyz_pow("A", 1.0, l, m, n, shift, array=array, rhs_only=True)
            cg.write("phi_tmp[%d + i] = %s * SZ" % (shift_idx, rhs))

            if z_grad:
                cg.write("phi_tmp[%d + i] += %s * S0[i]" % (shift_idx, AZ))

        # Hessian
        if specific_deriv == "XX":
            rhs = _build_xyz_pow("A", 1.0, l, m, n, shift, array=array, rhs_only=True)
            cg.write("phi_tmp[%d + i] = %s * SXX" % (shift_idx, rhs))

            # Cross term, need to write it out if specific deriv
            if x_grad:
                rhs = _build_xyz_pow("AX", ld2, ld1, m, n, shift, array=array, scale=2.0, rhs_only=True)
                cg.write("phi_tmp[%d + i] += %s * SX" % (shift_idx, rhs))

            AXX = _build_xyz_pow("AXX", ld2 * (ld2 - 1), ld2, m, n, shift, array=array, rhs_only=True)
            if AXX is not None:
                cg.write("phi_tmp[%d + i] += %s * S0[i]" % (shift_idx, AXX))

        # YY
        if specific_deriv == "YY":
            rhs = _build_xyz_pow("A", 1.0, l, m, n, shift, array=array, rhs_only=True)
            cg.write("phi_tmp[%d + i] = %s * SYY" % (shift_idx, rhs))
            if y_grad:
                rhs = _build_xyz_pow("AY", md2, l, md1, n, shift, array=array, scale=2.0, rhs_only=True)
                cg.write("phi_tmp[%d + i] += %s * SY" % (shift_idx, rhs))

            AYY = _build_xyz_pow("AYY", md2 * (md2 - 1), l, md2, n, shift, array=array, rhs_only=True)
            if AYY is not None:
                cg.write("phi_tmp[%d + i] += %s * S0[i]" % (shift_idx, AYY))

        # ZZ
        if specific_deriv == "ZZ":
            rhs = _build_xyz_pow("A", 1.0, l, m, n, shift, array=array, rhs_only=True)
            cg.write("phi_tmp[%d + i] = %s * SZZ" % (shift_idx, rhs))
            if z_grad:
                rhs = _build_xyz_pow("AZ", nd2, l, m, nd1, shift, array=array, scale=2.0, rhs_only=True)
                cg.write("phi_tmp[%d + i] += %s * SZ" % (shift_idx, rhs))

            AZZ = _build_xyz_pow("AZZ", nd2 * (nd2 - 1), l, m, nd2, shift, array=array, rhs_only=True)
            if AZZ is not None:
                cg.write("phi_tmp[%d + i] += %s * S0[i]" % (shift_idx, AZZ))

        # XY
        if specific_deriv == "XY":
            rhs = _build_xyz_pow("A", 1.0, l, m, n, shift, array=array, rhs_only=True)
            cg.write("phi_tmp[%d + i] = %s * SXY" % (shift_idx, rhs))

            if y_grad:
                rhs = AY.split(" = ")[-1]
                cg.write("phi_tmp[%d + i] += %s * SX" % (shift_idx, rhs))
            if x_grad:
                rhs = AX.split(" = ")[-1]
                cg.write("phi_tmp[%d + i] += %s * SY" % (shift_idx, rhs))

            AXY = _build_xyz_pow("AXY", ld2 * md2, ld1, md1, n, shift, array=array, rhs_only=True)
            if AXY is not None:
                cg.write("phi_tmp[%d + i] += %s * S0[i]" % (shift_idx, AXY))

        # XZ
        if specific_deriv == "XZ":
            rhs = _build_xyz_pow("A", 1.0, l, m, n, shift, array=array, rhs_only=True)
            cg.write("phi_tmp[%d + i] = %s * SXZ" % (shift_idx, rhs))
            if z_grad:
                rhs = AZ.split(" = ")[-1]
                cg.write("phi_tmp[%d + i] += %s * SX" % (shift_idx, rhs))
            if x_grad:
                rhs = AX.split(" = ")[-1]
                cg.write("phi_tmp[%d + i] += %s * SZ" % (shift_idx, rhs))

            AXZ = _build_xyz_pow("AXZ", ld2 * nd2, ld1, m, nd1, shift, array=array, rhs_only=True)
            if AXZ is not None:
                cg.write("phi_tmp[%d + i] += %s * S0[i]" % (shift_idx, AXZ))

        # YZ
        if specific_deriv == "YZ":
            rhs = _build_xyz_pow("A", 1.0, l, m, n, shift, array=array, rhs_only=True)
            cg.write("phi_tmp[%d + i] = %s * SYZ" % (shift_idx, rhs))
            if z_grad:
                rhs = AZ.split(" = ")[-1]
                cg.write("phi_tmp[%d + i] += %s * SY" % (shift_idx, rhs))
            if y_grad:
                rhs = AY.split(" = ")[-1]
                cg.write("phi_tmp[%d + i] += %s * SZ" % (shift_idx, rhs))

            AYZ = _build_xyz_pow("AYZ", md2 * nd2, l, md1, nd1, shift, array=array, rhs_only=True)
            if AYZ is not None:
                cg.write("phi_tmp[%d + i] += %s * S0[i]" % (shift_idx, AYZ))

        idx += 1
        cg.blankline()

    cg.close_c_block()
    cg.blankline()


def _c_am_full_build(cg, L, cartesian_order, grad, shift):
    """
    Builds a unrolled angular momentum function
    """

    # Generator
    for idx, l, m, n in order.cartesian_order_factory(L, cartesian_order):

        l = l + 2
        m = m + 2
        n = n + 2
        ld1 = l - 1
        ld2 = l - 2
        md1 = m - 1
        md2 = m - 2
        nd1 = n - 1
        nd2 = n - 2

        # Set grads back to zero
        x_grad, y_grad, z_grad = False, False, False
        shift_idx = idx * shift

        name = "X" * ld2 + "Y" * md2 + "Z" * nd2
        if name == "":
            name = "0"

        # Density
        cg.blankline()
        cg.write("// Density AM=%d Component=%s" % (L, name))

        cg.write(_build_xyz_pow("A", 1.0, l, m, n, shift))
        cg.write("phi_tmp[%d + i] = S0[i] * A" % shift_idx)

        if grad == 0: continue
        cg.blankline()
        cg.write("// Gradient AM=%d Component=%s" % (L, name))

        # Gradient
        cg.write("phi_x_tmp[%d + i] = SX * A" % shift_idx)
        cg.write("phi_y_tmp[%d + i] = SY * A" % shift_idx)
        cg.write("phi_z_tmp[%d + i] = SZ * A" % shift_idx)

        AX = _build_xyz_pow("AX", ld2, ld1, m, n, shift)
        if AX is not None:
            x_grad = True
            cg.write(AX)
            cg.write("phi_x_tmp[%d + i] += S0[i] * AX" % shift_idx)

        AY = _build_xyz_pow("AY", md2, l, md1, n, shift)
        if AY is not None:
            y_grad = True
            cg.write(AY)
            cg.write("phi_y_tmp[%d + i] += S0[i] * AY" % shift_idx)

        AZ = _build_xyz_pow("AZ", nd2, l, m, nd1, shift)
        if AZ is not None:
            z_grad = True
            cg.write(AZ)
            cg.write("phi_z_tmp[%d + i] += S0[i] * AZ" % shift_idx)

        # Hessian temporaries
        if grad == 1: continue

        cg.blankline()
        cg.write("// Hessian AM=%d Component=%s" % (L, name))

        # S Hess
        # We will build S Hess, grad 1, grad 2, A Hess

        # XX
        cg.write("phi_xx_tmp[%d + i] = SXX * A" % shift_idx)
        if x_grad:
            cg.write("phi_xx_tmp[%d + i] += SX * AX" % shift_idx)
            cg.write("phi_xx_tmp[%d + i] += SX * AX" % shift_idx)

        AXX = _build_xyz_pow("AXX", ld2 * (ld2 - 1), ld2, m, n, shift)
        if AXX is not None:
            rhs = AXX.split(" = ")[-1]
            cg.write("phi_xx_tmp[%d + i] += %s * S0[i]" % (shift_idx, rhs))

        # YY
        cg.write("phi_yy_tmp[%d + i] = SYY * A" % shift_idx)
        if y_grad:
            cg.write("phi_yy_tmp[%d + i] += SY * AY" % shift_idx)
            cg.write("phi_yy_tmp[%d + i] += SY * AY" % shift_idx)
        AYY = _build_xyz_pow("AYY", md2 * (md2 - 1), l, md2, n, shift)
        if AYY is not None:
            rhs = AYY.split(" = ")[-1]
            cg.write("phi_yy_tmp[%d + i] += %s * S0[i]" % (shift_idx, rhs))

        # ZZ
        cg.write("phi_zz_tmp[%d + i] = SZZ * A" % shift_idx)
        if z_grad:
            cg.write("phi_zz_tmp[%d + i] += SZ * AZ" % shift_idx)
            cg.write("phi_zz_tmp[%d + i] += SZ * AZ" % shift_idx)
        AZZ = _build_xyz_pow("AZZ", nd2 * (nd2 - 1), l, m, nd2, shift)
        if AZZ is not None:
            rhs = AZZ.split(" = ")[-1]
            cg.write("phi_zz_tmp[%d + i] += %s * S0[i]" % (shift_idx, rhs))

        # XY
        cg.write("phi_xy_tmp[%d + i] = SXY * A" % shift_idx)

        if y_grad:
            cg.write("phi_xy_tmp[%d + i] += SX * AY" % shift_idx)
        if x_grad:
            cg.write("phi_xy_tmp[%d + i] += SY * AX" % shift_idx)

        AXY = _build_xyz_pow("AXY", ld2 * md2, ld1, md1, n, shift)
        if AXY is not None:
            rhs = AXY.split(" = ")[-1]
            cg.write("phi_xy_tmp[%d + i] += %s * S0[i]" % (shift_idx, rhs))

        # XZ
        cg.write("phi_xz_tmp[%d + i] = SXZ * A" % shift_idx)
        if z_grad:
            cg.write("phi_xz_tmp[%d + i] += SX * AZ" % shift_idx)
        if x_grad:
            cg.write("phi_xz_tmp[%d + i] += SZ * AX" % shift_idx)
        AXZ = _build_xyz_pow("AXZ", ld2 * nd2, ld1, m, nd1, shift)
        if AXZ is not None:
            rhs = AXZ.split(" = ")[-1]
            cg.write("phi_xz_tmp[%d + i] += %s * S0[i]" % (shift_idx, rhs))

        # YZ
        cg.write("phi_yz_tmp[%d + i] = SYZ * A" % shift_idx)
        if z_grad:
            cg.write("phi_yz_tmp[%d + i] += SY * AZ" % shift_idx)
        if y_grad:
            cg.write("phi_yz_tmp[%d + i] += SZ * AY" % shift_idx)
        AYZ = _build_xyz_pow("AYZ", md2 * nd2, l, md1, nd1, shift)
        if AYZ is not None:
            rhs = AYZ.split(" = ")[-1]
            cg.write("phi_yz_tmp[%d + i] += %s * S0[i]" % (shift_idx, rhs))

        idx += 1
        cg.blankline()


def _build_xyz_pow(name, pref, l, m, n, inner_loop, shift=2, array=False, scale=1.0, rhs_only=False):
    """
    Builds an individual row contraction line.

    name = pref * xc_pow[n] yc_pow[m] * zc_pow[n]
    """
    l = l - shift
    m = m - shift
    n = n - shift

    if (pref <= 0) or (l < 0) or (n < 0) or (m < 0):
        return None

    if rhs_only:
        ret = ""
    else:
        ret = name + " ="

    mul = " "
    if (pref * scale) != 1.0:
        # Basically always an int
        ret += " %2.1f" % (float(pref) * scale)
        mul = " * "

    # Handle x
    if l == 1:
        # If the power is one, we can just use xc
        ret += mul + "xc[i]"
        mul = " * "
    elif l > 1:
        # If the power is greater than 1 we need to use (xc_pow - 2) as we start at 2
        if array:
            ret += mul + "xc_pow[%d + i]" % ((l - 2) * inner_loop)
        else:
            ret += mul + "xc_pow%d" % l
        mul = " * "

    # Handle y
    if m == 1:
        ret += mul + "yc[i]"
        mul = " * "
    elif m > 1:
        if array:
            ret += mul + "yc_pow[%d + i]" % ((m - 2) * inner_loop)
        else:
            ret += mul + "yc_pow%d" % m
        mul = " * "

    # Handle z
    if n == 1:
        ret += mul + "zc[i]"
        mul = " * "
    elif n > 1:
        if array:
            ret += mul + "zc_pow[%d + i]" % ((n - 2) * inner_loop)
        else:
            ret += mul + "zc_pow%d" % n
        mul = " * "

    if rhs_only:
        ret = ret.strip()

    if mul == " ":
        ret += " 1"

    return ret


def _S_tmps(cg, L, grad, inner_block):
    """
    Builds out the S power temporaries if needed
    """
    if grad > 0:
        cg.write("// Gaussian derivs (gradients)")
        cg.write("const double SX = S1[i] * xc[i]")
        cg.write("const double SY = S1[i] * yc[i]")
        cg.write("const double SZ = S1[i] * zc[i]")
    if grad > 1:
        cg.blankline()
        cg.write("// Gaussians derivs (Hessians)")
        cg.write("const double SXY = S2[i] * xc[i] * yc[i]")
        cg.write("const double SXZ = S2[i] * xc[i] * zc[i]")
        cg.write("const double SYZ = S2[i] * yc[i] * zc[i]")
        cg.write("const double SXX = S2[i] * xc[i] * xc[i] + S1[i]")
        cg.write("const double SYY = S2[i] * yc[i] * yc[i] + S1[i]")
        cg.write("const double SZZ = S2[i] * zc[i] * zc[i] + S1[i]")


def _power_tmps(cg, L, inner_block, array=False):
    if L < 2:
        return

    # L == 2
    # cg.write("PRAGMA_VECTORIZE", endl="")
    # cg.start_c_block("for (unsigned long i = 0; i < remain; i++)")
    if array:
        # Build out those power derivs
        cg.blankline()
        cg.write("// Cartesian derivs")
        cg.write("xc_pow[i] = xc[i] * xc[i]")
        cg.write("yc_pow[i] = yc[i] * yc[i]")
        cg.write("zc_pow[i] = zc[i] * zc[i]")

        if L == 2:
            cg.blankline()

        for l in range(1, (L - 1)):
            cg.write("xc_pow[%d + i] = xc_pow[%d + i] * xc[i]" % (inner_block * l, inner_block * (l - 1)))
            cg.write("yc_pow[%d + i] = yc_pow[%d + i] * yc[i]" % (inner_block * l, inner_block * (l - 1)))
            cg.write("zc_pow[%d + i] = zc_pow[%d + i] * zc[i]" % (inner_block * l, inner_block * (l - 1)))

    else:
        # Build out those power derivs
        cg.blankline()
        cg.write("// Cartesian derivs")
        cg.write("const double xc_pow2 = xc[i] * xc[i]")
        cg.write("const double yc_pow2 = yc[i] * yc[i]")
        cg.write("const double zc_pow2 = zc[i] * zc[i]")

        cg.blankline()

        for l in range(2, L):
            cg.write("const double xc_pow%d = xc_pow%d * xc[i]" % (l + 1, l))
            cg.write("const double yc_pow%d = yc_pow%d * yc[i]" % (l + 1, l))
            cg.write("const double zc_pow%d = zc_pow%d * zc[i]" % (l + 1, l))
            cg.blankline()


def _tmp_to_out_copy(cg, L, deriv_indices, inner_block):

    # Start spherical switch
    cg.blankline()
    cg.write("// Copy data back into outer temps")
    cg.start_c_block("if (spherical)")
    sph_fnc = "gg_cart_to_spherical_L%d" % L
    block_fnc = "block_copy"

    cg.write("// Phi, transform data to outer temps")
    cg.write("%s(remain, phi_tmp, %d, (phi_out + start), npoints)" % (sph_fnc, inner_block))

    for num, deriv in enumerate(deriv_indices):
        # Write out pretty headers
        if num == 0:
            cg.blankline()
            cg.write("// Gradient, transform data to outer temps")
        if num == 3:
            cg.blankline()
            cg.write("// Hessian, transform data to outer temps")

        cg.write("%s(remain, phi_%s_tmp, %d, (phi_%s_out + start), npoints)" % (sph_fnc, deriv, inner_block, deriv))

    cg.write("} else {", endl="")
    # Move data into inner buffers

    # Copy over Phi
    cg.blankline()
    cg.write("// Phi, copy data to outer temps")
    cg.write("%s(nout, remain, phi_tmp, %d, (phi_out + start), npoints, 0)" % (block_fnc, inner_block))

    # Copy over grad
    for num, deriv in enumerate(deriv_indices):
        # Write out pretty headers
        if num == 0:
            cg.blankline()
            cg.write("// Gradient, copy data to outer temps")
        if num == 3:
            cg.blankline()
            cg.write("// Hessian, copy data to outer temps")

        cg.write("%s(nout, remain, phi_%s_tmp, %d, (phi_%s_out + start), npoints, 0)" % (block_fnc, deriv, inner_block,
                                                                                         deriv))

    # cg.close_c_block()
    cg.blankline()

    # End spherical switch
    cg.close_c_block()
