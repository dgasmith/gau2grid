"""
The C generator for gau2grid collocation functions
"""

import os

from . import RSH
from . import codegen
from . import order
from . import utility
from . import c_pragma

_grad_indices = ["x", "y", "z"]
_hess_indices = ["xx", "xy", "xz", "yy", "yz", "zz"]


def generate_c_gau2grid(max_L, path=".", cart_order="row", inner_block=64, do_cf=True):
    """
    Creates the C files for the gau2grid program.

    Parameters
    ----------
    max_L : int
        The maximum angular momentum compiled for.
    path : str, optional
        The path to write the files to.
    cart_order : str, optional
        The order of the cartesian angular momentum
    inner_block : int, optional
        The size of the inner computation block.
    do_cf : bool, option
        Apply clang-format to the generated files or not.

    Returns
    -------
    None

    """

    # Build the codegen objects for each file
    gg_header = codegen.CodeGen(cgen=True)
    gg_phi = codegen.CodeGen(cgen=True)
    gg_grad = codegen.CodeGen(cgen=True)
    gg_hess = codegen.CodeGen(cgen=True)
    gg_spherical = codegen.CodeGen(cgen=True)
    gg_pybind = codegen.CodeGen(cgen=True)
    gg_helper = codegen.CodeGen(cgen=True)
    gg_pragma = codegen.CodeGen(cgen=True)

    # Add general header comments
    for cgs in [gg_header, gg_phi, gg_grad, gg_hess, gg_spherical, gg_pybind, gg_helper, gg_pragma]:
        cgs.write("// This is an automtically generated file from ...")
        cgs.write("// Blah blah blah")
        cgs.blankline()

    # Write out the pragma header
    c_pragma.build_pragma_header(gg_pragma)

    gg_helper.write("#include <stdio.h>")

    # Add utility headers
    for cgs in [gg_phi, gg_grad, gg_hess, gg_spherical, gg_helper]:
        cgs.write("#include <math.h>")
        cgs.write("#include <stdio.h>")
        cgs.write("#include <mm_malloc.h>")
        cgs.blankline()
        cgs.write('#include "gau2grid.h"')
        cgs.write('#include "gau2grid_pragma.h"')
        cgs.blankline()
        cgs.write("// Adds a few typedefs to make the world easier")
        cgs.write("typedef unsigned long size_t")
        cgs.blankline()

    # Header guards
    gg_header.write("#ifndef GAU2GRID_GUARD_H")
    gg_header.write("#define GAU2GRID_GUARD_H")
    gg_header.blankline()

    # Add any information needed
    gg_helper.write("// Information helpers")
    gg_helper.write("int max_L() { return %d; }" % max_L, endl="")
    gg_helper.blankline()
    gg_helper.write("// Collocation selector functions")

    gg_header.write("// Information helpers")
    gg_header.write("int max_L()")
    gg_header.blankline()

    # Build out the spherical transformer
    gg_header.write("// Spherical transformers")
    for L in range(max_L + 1):
        sig = RSH.transformation_c_generator(gg_spherical, L, cart_order)
        gg_header.write(sig)
    gg_header.blankline()

    # Fast transformers
    trans_sig = _naive_transpose(gg_spherical)
    gg_header.write(trans_sig)
    fast_trans_sig = _fast_transpose(gg_spherical, 8)
    gg_header.write(fast_trans_sig)
    gg_header.blankline()


    # Loop over phi, grad, hess and build blocks for each
    helper_sigs = []
    for name, grad, cg in [("Phi", 0, gg_phi), ("Phi grad", 1, gg_grad), ("Phi Hess", 2, gg_hess)]:
        cg.blankline()
        gg_header.write("// %s computers" % name)
        cg.blankline()

        # Write out the phi builders
        sig_store = []
        for L in range(max_L + 1):
            sig = shell_c_generator(cg, L, grad=grad, cart_order=cart_order, inner_block=inner_block)
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
    gg_header.write("#endif /* GAU2GRID_GUARD_H */")

    # Write out the CG's to files
    gg_header.repr(filename=os.path.join(path, "gau2grid.h"), clang_format=do_cf)
    gg_phi.repr(filename=os.path.join(path, "gau2grid_phi.cc"), clang_format=do_cf)
    gg_grad.repr(filename=os.path.join(path, "gau2grid_deriv1.cc"), clang_format=do_cf)
    gg_hess.repr(filename=os.path.join(path, "gau2grid_deriv2.cc"), clang_format=do_cf)
    gg_spherical.repr(filename=os.path.join(path, "gau2grid_spherical.cc"), clang_format=do_cf)
    gg_helper.repr(filename=os.path.join(path, "gau2grid_helper.cc"), clang_format=do_cf)
    gg_pragma.repr(filename=os.path.join(path, "gau2grid_pragma.h"))

    ### Build out the PyBind11 plugin
    gg_pybind.blankline()
    gg_pybind.write("#include <pybind11/pybind11.h>")
    gg_pybind.write("#include <pybind11/numpy.h>")
    gg_pybind.write('#include "gau2grid.h"')
    gg_pybind.blankline()

    gg_pybind.write("namespace py = pybind11")
    gg_pybind.blankline()

    # Name space CXX files
    gg_pybind.write("namespace gau2grid {", endl="")
    gg_pybind.blankline()

    # Call the execute functions
    _pybind11_func(gg_pybind, "collocation_wrapper", 0, helper_sigs[0], max_L)
    _pybind11_func(gg_pybind, "collocation_deriv1_wrapper", 1, helper_sigs[1], max_L)
    _pybind11_func(gg_pybind, "collocation_deriv2_wrapper", 2, helper_sigs[2], max_L)

    tsig = _pybind11_transpose(gg_pybind, "gg_naive_transpose", "naive_transpose_wrapper")
    tsig = _pybind11_transpose(gg_pybind, "gg_fast_transpose", "fast_transpose_wrapper")

    # Open up the pybind module
    gg_pybind.start_c_block("PYBIND11_MODULE(pygg_core, m)")
    gg_pybind.write('m.doc() = "A Python wrapper to the Gau2Grid library."')
    gg_pybind.write('m.def("max_L", &max_L)')
    gg_pybind.write('m.def("collocation", &collocation_wrapper)')
    gg_pybind.write('m.def("collocation_deriv1", &collocation_deriv1_wrapper)')
    gg_pybind.write('m.def("collocation_deriv2", &collocation_deriv2_wrapper)')
    gg_pybind.write('m.def("naive_transpose", &naive_transpose_wrapper)')
    gg_pybind.write('m.def("fast_transpose", &fast_transpose_wrapper)')
    gg_pybind.blankline()

    # Close out the pybind module
    gg_pybind.close_c_block()

    gg_pybind.blankline()

    gg_pybind.write("} // End gau2grid namespace", endl="")

    gg_pybind.repr(filename=os.path.join(path, "pygg_core.cc"), clang_format=do_cf)
    # print(gg_pybind.repr(clang_format=do_cf))


def shell_c_generator(cg, L, function_name="", grad=0, cart_order="row", inner_block=64):

    # Set a few parameters for custom loops
    L_needs_out = 1

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

    # Precompute temps
    ncart = int((L + 1) * (L + 2) / 2)
    nspherical = L * 2 + 1

    # Build function signature
    func_sig = "const size_t npoints, const double* __restrict__ x, const double* __restrict__ y, const double* __restrict__ z, const int nprim, const double* __restrict__ coeffs, const double* __restrict__ exponents, const double* __restrict__ center, const int spherical, double* __restrict__ phi_out"

    # Add extra output vals for derivs
    for deriv in deriv_indices:
        func_sig += ", double* __restrict__ phi_%s_out" % deriv

    func_sig = "void %s(%s)" % (function_name, func_sig)
    cg.start_c_block(func_sig)
    cg.blankline()

    # Figure out spacing
    cg.write("// Sizing")
    cg.write("size_t nblocks = npoints / %d" % inner_block)
    cg.write("nblocks += (npoints %% %d) ? 1 : 0" % inner_block)
    cg.write("const size_t ncart = %d" % ncart)
    cg.write("const size_t nspherical = %d" % nspherical)
    cg.write("const size_t nout = spherical ? nspherical : ncart")
    cg.blankline()

    # Build temporaries
    cg.write("// Allocate S temporaries, single block to stay on cache")
    S_cache_tmps = ["xc", "yc", "zc", "R2", "S0"]
    if grad > 0:
        S_cache_tmps.append("S1")
    if grad > 1:
        S_cache_tmps.append("S2")

    if True:
        # Allocate as single block on heap
        _block_malloc(cg, "cache_data", S_cache_tmps, inner_block)
        S_tmps = ["cache_data"]
    else:
        # Allocate as single block on stack
        for tmp in S_cache_tmps:
             cg.write("double %s[%d] __attribute__((aligned(64)))" % (tmp, inner_block))
        S_tmps = []

    cg.blankline()

    # Hold the expn1 and expn2 arrays
    exp_tmps = ["expn1"]
    if grad > 0:
        exp_tmps += ["expn2"]
    for tname in exp_tmps:
        cg.write(_malloc(tname, "nprim"))
    S_tmps.extend(exp_tmps)
    cg.blankline()

    # Figure out powers needed
    power_tmps = []
    if L > 1:
        cg.write("// Allocate power temporaries")
        power_tmps = ["xc_pow", "yc_pow", "zc_pow"]
        for tname in power_tmps:
            cg.write(_malloc(tname, inner_block * (L - 1)))
        cg.blankline()

    # Determine tmps
    inner_tmps = []
    if L >= L_needs_out:
        cg.write("// Allocate output temporaries")

        inner_tmps = ["phi_tmp"]
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
    cg.start_c_block("for (size_t i = 0; i < nprim; i++)")
    cg.write("expn1[i] = -1.0 * exponents[i]")
    if grad > 0:
        cg.write("expn2[i] = -2.0 * exponents[i]")
    cg.close_c_block()
    cg.blankline()

    # Start outer loop
    cg.write("// Start outer block loop")
    cg.start_c_block("for (size_t block = 0; block < nblocks; block++)")
    cg.blankline()

    # Move data into inner buffers
    cg.blankline()
    cg.write("// Copy data into inner temps")
    cg.write("const size_t start = block * %d" % inner_block)
    cg.write("const size_t remain = ((start + %d) > npoints) ? (npoints - start) : %d" % (inner_block, inner_block))

    cg.write("PRAGMA_VECTORIZE", endl="")
    cg.start_c_block("for (size_t i = 0; i < remain; i++)")
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
    cg.start_c_block("for (size_t n = 0; n < nprim; n++)")

    # Build R2
    cg.write("const double coef = coeffs[n]")
    cg.write("const double alpha_n1 = expn1[n]")
    if grad > 0:
        cg.write("const double alpha_n2 = expn2[n]")

    # Build out thoese gaussian derivs
    cg.blankline()
    cg.write("PRAGMA_VECTORIZE", endl="")
    cg.start_c_block("for (size_t i = 0; i < remain; i++)")
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

    # Combine blocks
    cg.write("// Combine blocks")
    cg.write("PRAGMA_VECTORIZE", endl="")
    cg.start_c_block("for (size_t i = 0; i < remain; i++)")
    _S_pow_tmps(cg, L, grad, inner_block)

    if L == 0:
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

    else:

        # Contract temps with powers
        _c_am_build(cg, L, cart_order, grad, inner_block)

        cg.blankline()

    # End inner loop
    cg.close_c_block()

    # Grab the inner line stop
    inner_line_stop = len(cg.data)

    # Spherical/Cartesian copy out
    if L >= L_needs_out:
        _tmp_to_out_copy(cg, L, deriv_indices, inner_block)

    # End outer loop
    cg.close_c_block()

    # Free up those arrays
    cg.blankline()
    for name, flist in [("S", S_tmps), ("power", power_tmps), ("inner", inner_tmps)]:
        if len(flist) == 0: continue

        cg.write("// Free %s temporaries" % name)
        for tname in flist:
            cg.write("_mm_free(%s)" % tname)
        cg.blankline()

    # End function
    cg.close_c_block()

    # Clean up data, there are a few things easier to post-process

    # Remove any "[0 + i]"
    for x in range(cg_line_start, len(cg.data)):
        cg.data[x] = cg.data[x].replace("[0 + ", "[")

    # return func_sig
    # Remove any "A = 1" just for the inner block
    rep_data = {}
    pos = inner_line_start
    while pos < inner_line_stop:
        line = cg.data[pos]

        # If we hit a Density line its an individual angular momentum, need to reset dict
        if "Density" in line:
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
            continue

        for k, v in rep_data.items():
            tmp = line.split("= ")[1]
            if k + ";" in tmp:
                cg.data[pos] = line.replace(k + ";", v + ";")
            elif k + " " in tmp:
                cg.data[pos] = line.replace(k + " ", v + " ")
        pos += 1

    # Remove any " * 1"
    for x in range(cg_line_start, len(cg.data)):
        cg.data[x] = cg.data[x].replace(" * 1;", ";")
        cg.data[x] = cg.data[x].replace(" * 1.0;", ";")
        cg.data[x] = cg.data[x].replace("= 1 * ", "= ")
        cg.data[x] = cg.data[x].replace("= 1.0 * ", "= ")

    return func_sig


def _make_call(string):
    for rep in ["double* ", "bool ", "int ", "size_t ", "void ", "__restrict__ "]:
        string = string.replace("const " + rep, "")
        string = string.replace(rep, "")
    return string


def _malloc(name, size, dtype="double"):
    # return "%s*  %s = (%s*)malloc(%s * sizeof(%s))" % (dtype, name, dtype, str(size), dtype)
    return "%s* __restrict__ %s = (%s*)_mm_malloc(%s * sizeof(%s), 32)" % (dtype, name, dtype, str(size), dtype)

def _block_malloc(cg, block_name, names, size, dtype="double"):
    tot_size = len(names) * size
    cg.write(_malloc(block_name, tot_size))
    for num, name in enumerate(names):
        cg.write("%s* __restrict__ %s = %s + %d" % (dtype, name, block_name, num * size))


def _c_am_build(cg, L, cart_order, grad, shift):
    """
    Builds a unrolled angular momentum function
    """

    # Generator
    for idx, l, m, n in order.cartesian_order_factory(L, cart_order):

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
        cg.write("// Starting AM=%d Component=%s" % (L, name))
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
            # cg.write(AYZ)
            rhs = AYZ.split(" = ")[-1]
            cg.write("phi_yz_tmp[%d + i] += %s * S0[i]" % (shift_idx, rhs))

        idx += 1
        cg.blankline()


def _build_xyz_pow(name, pref, l, m, n, inner_loop, shift=2):
    """
    Builds an individual row contraction line.

    name = pref * xc_pow[n] yc_pow[m] * zc_pow[n]
    """
    l = l - shift
    m = m - shift
    n = n - shift

    if (pref <= 0) or (l < 0) or (n < 0) or (m < 0):
        return None

    mul = " "
    if pref == 1:
        ret = name + " ="
    else:
        # Basically always an int
        ret = name + " = %2.1f" % float(pref)
        mul = " * "

    # Handle x
    if l == 1:
        # If the power is one, we can just use xc
        ret += mul + "xc[i]"
        mul = " * "
    elif l > 1:
        # If the power is greater than 1 we need to use (xc_pow - 2) as we start at 2
        ret += mul + "xc_pow[%d + i]" % ((l - 2) * inner_loop)
        mul = " * "

    # Handle y
    if m == 1:
        ret += mul + "yc[i]"
        mul = " * "
    elif m > 1:
        ret += mul + "yc_pow[%d + i]" % ((m - 2) * inner_loop)
        mul = " * "

    # Handle z
    if n == 1:
        ret += mul + "zc[i]"
        mul = " * "
    elif n > 1:
        ret += mul + "zc_pow[%d + i]" % ((n - 2) * inner_loop)
        mul = " * "

    if mul == " ":
        ret += " 1"

    return ret


def _S_pow_tmps(cg, L, grad, inner_block):
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

    # Build out those power derivs
    cg.blankline()
    if L > 1:
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
        cg.blankline()


def _pybind11_func(cg, name, grad, call_name, max_L):
    """
    A function that builds the PyBind11 wrapper functions
    """

    # Figure out what we need to add per deriv
    deriv_indices = utility.get_deriv_indices(grad)

    # Write out wrapper functions
    sig = """void %s(int L, py::array_t<double> arr_xyz, py::array_t<double> arr_coeffs,
py::array_t<double> arr_exponents, py::array_t<double> arr_center, bool spherical,
py::array_t<double> arr_out""" % name

    # Pad out deriv outputs
    for cart in deriv_indices:
        sig += ", py::array_t<double> arr_%s_out" % cart

    sig += ")"
    cg.start_c_block(sig)
    cg.blankline()

    # Grab the pointers
    cg.write('// Grab array pointers')
    cg.write('auto xyz = arr_xyz.unchecked<2>()')
    cg.write('auto coeffs = arr_coeffs.unchecked<1>()')
    cg.write('auto exponents = arr_exponents.unchecked<1>()')
    cg.write('auto center = arr_center.unchecked<1>()')
    cg.write('auto out = arr_out.mutable_unchecked<2>()')

    # Pad out deriv pointers
    for cart in deriv_indices:
        cg.write("auto out_%s = arr_%s_out.mutable_unchecked<2>()" % (cart, cart))

    cg.blankline()

    # Run through checks
    cg.write('// XYZ is of size 3')
    cg.start_c_block('if (L > %d)' % max_L)
    cg.write(
        '    throw std::invalid_argument("Exceeded compiled angular momentum of %d. Please recompile with a higher angular momentum.\\n")'
        % max_L)
    cg.close_c_block()

    cg.write('// XYZ is of size 3')
    cg.start_c_block('if (arr_xyz.shape(0) != 3)')
    cg.write('    throw std::length_error("Length of XYZ array must be (3, n).\\n")')
    cg.close_c_block()
    cg.blankline()

    cg.write('// Coeff matches exponent shape')
    cg.start_c_block('if (coeffs.shape(0) != exponents.shape(0))')
    cg.write('    throw std::length_error("Length of coefficients and exponents must match.\\n")')
    cg.close_c_block()
    cg.blankline()

    cg.write('// Center is of size 3')
    cg.start_c_block('if (center.shape(0) != 3)')
    cg.write('    throw std::length_error("Length of center vector must be 3 (X, Y, Z).\\n")')
    cg.close_c_block()
    cg.blankline()

    cg.write('// Make sure output length matches')
    cg.write('size_t nsize')
    cg.start_c_block('if (spherical)')
    cg.write('    nsize = 2 * L + 1')
    cg.write('} else {', endl="")
    cg.write('    nsize = ((L + 2) * (L + 1)) / 2')
    cg.close_c_block()
    cg.blankline()

    cg.start_c_block('if (out.shape(0) != nsize)')
    cg.write('    throw std::length_error("Size of the output array does not match the angular momentum.\\n")')
    cg.close_c_block()

    for cart in deriv_indices:
        cg.start_c_block('if (out_%s.shape(0) != nsize)' % cart)
        cg.write('    throw std::length_error("Size of the output %s array does not match the angular momentum.\\n")' %
                 cart.upper())
        cg.close_c_block()
    cg.blankline()

    cg.write('// Ensure lengths match')
    cg.start_c_block('if (out.shape(1) != arr_xyz.shape(1))')
    cg.write('    throw std::length_error("Size of the output array and XYZ array must be the same.\\n")')
    cg.close_c_block()

    # Pad out deriv length checkers
    for cart in deriv_indices:
        cg.start_c_block('if (out_%s.shape(1) != arr_xyz.shape(1))' % cart)
        cg.write('    throw std::length_error("Size of the output %s array and XYZ array must be the same.\\n")' %
                 cart.upper())
        cg.close_c_block()
    cg.blankline()

    cg.write("// Call the GG helper function")
    call_func = call_name + "(L, xyz.shape(1)"
    call_func += ", xyz.data(0, 0), xyz.data(1, 0), xyz.data(2, 0)"
    call_func += ", coeffs.shape(0), coeffs.data(0), exponents.data(0)"
    call_func += ", center.data(0)"
    call_func += ", spherical"
    call_func += ", out.mutable_data(0, 0)"
    for cart in deriv_indices:
        call_func += ", out_%s.mutable_data(0, 0)" % cart
    call_func += ")"

    cg.write(call_func)

    cg.close_c_block()


def _tmp_to_out_copy(cg, L, deriv_indices, inner_block):

    # Start spherical switch
    cg.blankline()
    cg.write("// Copy data back into outer temps")
    cg.start_c_block("if (spherical)")
    sph_fnc = "gg_cart_to_spherical_L%d" % L

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
    cg.blankline()
    cg.start_c_block("for (size_t n = 0; n < nout; n++)")
    cg.write("const size_t out_shift = start + n * npoints")
    cg.write("const size_t tmp_shift = n * %d" % inner_block)

    # Copy over Phi
    cg.blankline()
    cg.write("// Phi, copy data to outer temps")
    cg.start_c_block("for (size_t i = 0; i < remain; i++)")
    cg.write("phi_out[out_shift + i] = phi_tmp[tmp_shift + i]")
    cg.close_c_block()

    # Copy over grad
    for num, deriv in enumerate(deriv_indices):
        # Write out pretty headers
        if num == 0:
            cg.blankline()
            cg.write("// Gradient, copy data to outer temps")
        if num == 3:
            cg.blankline()
            cg.write("// Hessian, copy data to outer temps")

        cg.start_c_block("for (size_t i = 0; i < remain; i++)")
        cg.write("phi_%s_out[out_shift + i] = phi_%s_tmp[tmp_shift + i]" % (deriv, deriv))
        cg.close_c_block()

    cg.close_c_block()
    cg.blankline()

    # End spherical switch
    cg.close_c_block()

def _pybind11_transpose(cg, func_name, wrapper_name):

    sig = "void %s(py::array_t<double> arr_input" % wrapper_name
    sig += ", py::array_t<double> arr_output)"

    cg.start_c_block(sig)
    cg.write("auto input = arr_input.unchecked<2>()")
    cg.write("auto output = arr_output.mutable_unchecked<2>()")
    cg.write("size_t n = input.shape(0)")
    cg.write("size_t m = input.shape(1)")
    cg.blankline()

    cg.write('// Check shapes')
    cg.start_c_block('if (input.shape(0) != output.shape(1))')
    cg.write('    throw std::length_error("Input tranpose shape 0 does not match output transpose shape 1.\\n")')
    cg.close_c_block()
    cg.blankline()

    cg.start_c_block('if (input.shape(1) != output.shape(0))')
    cg.write('    throw std::length_error("Input tranpose shape 1 does not match output transpose shape 0.\\n")')
    cg.close_c_block()
    cg.blankline()

    cg.write("%s(n, m, input.data(0, 0), output.mutable_data(0, 0))" % func_name)

    cg.close_c_block()

def _naive_transpose(cg):

    sig = "void gg_naive_transpose(size_t n, size_t m, const double* __restrict__ input, double* __restrict__ output)"
    cg.start_c_block(sig)

    cg.start_c_block("for (size_t i = 0; i < n; i++)")

    # Inner block
    cg.start_c_block("for (size_t j = 0; j < m; j++)")
    cg.write("output[j * n + i] = input[i * m + j]")
    cg.close_c_block()

    # Outer block
    cg.close_c_block()

    cg.close_c_block()
    return sig


def _fast_transpose(cg, inner_block):
    """
    Builds a fast transpose
    """

    sig = "void gg_fast_transpose(size_t n, size_t m, const double* __restrict__ input, double* __restrict__ output)"
    cg.start_c_block(sig)
    cg.blankline()

    cg.write("// Temps")
    cg.write("double tmp[%d]  __attribute__((aligned(64)))" % (inner_block * inner_block))

    cg.write("// Sizing")
    cg.write("size_t nblocks = n / %d" % inner_block)
    cg.write("nblocks += (n %% %d) ? 1 : 0" % inner_block)

    cg.write("size_t mblocks = m / %d" % inner_block)
    cg.write("mblocks += (m %% %d) ? 1 : 0" % inner_block)
    # cg.write('printf("Blocks: %ld %ld\\n", nblocks, mblocks)')


    cg.write("// Outer blocks")
    cg.start_c_block("for (size_t nb = 0; nb < nblocks; nb++)")
    cg.write("const size_t nstart = nb * %d" % inner_block)
    cg.write("size_t nremain = ((nstart + %d) > n) ? (n - nstart) : %d" % (inner_block, inner_block))

    cg.start_c_block("for (size_t mb = 0; mb < mblocks; mb++)")
    cg.write("const size_t mstart = mb * %d" % inner_block)
    cg.write("size_t mremain = ((mstart + %d) > m) ? (m - mstart) : %d" % (inner_block, inner_block))

    # cg.start_c_block("if ((nremain == 0) & (mremain > 0))")
    # cg.write("nremain++;")
    # cg.close_c_block()

    # cg.start_c_block("if ((mremain == 0) & (nremain > 0))")
    # cg.write("mremain++;")
    # cg.close_c_block()
    # cg.write('printf("(n,m)%ld %ld | %ld %ld\\n", nb, mb, nremain, mremain)')


    # Pull block
    cg.write("// Copy data to inner block")
    # cg.write('printf("%ld %ld | %ld\\n ", mstart, nstart, start)')
    cg.start_c_block("for (size_t l = 0; l < nremain; l++)")
    cg.write("const size_t start = (nstart + l) * m + mstart")
    cg.write("PRAGMA_VECTORIZE", endl="")
    cg.start_c_block("for (size_t k = 0; k < mremain; k++)")

    # cg.write("tmp[l * %d + k] = input[start + k]" % inner_block)
    cg.write("tmp[k * %d + l] = input[start + k]" % inner_block)

    # cg.write('printf("(%ld %ld %lf) ", l * 2+ k, start +k, input[start + k])')
    # cg.write('printf("%%lf ", tmp[k * %d + l])' % inner_block)
    cg.close_c_block()
    cg.close_c_block()
    # cg.write('printf("\\n--\\n")')
    # cg.start_c_block("for (size_t k = 0; k < 4; k++)")
    # cg.write('printf("%lf ", tmp[k])')
    # cg.close_c_block()
    # cg.write('printf("\\n--\\n")')

    # Tranpose block
    # cg.write("// Transpose inner block")
    # cg.start_c_block("for (size_t k = 0; k < %d; k++)" % inner_block)
    # cg.start_c_block("for (size_t l = k; l < %d; l++)" % inner_block)
    # # cg.write('printf("%ld %ld \\n", k, l)')
    # cg.write("const double itmp = tmp[l * %d + k]" % inner_block)
    # cg.write("tmp[l * %d + k] = tmp[k * %d + l]" % (inner_block, inner_block))
    # cg.write("tmp[k * %d + l] = itmp" % (inner_block))
    # cg.close_c_block()
    # cg.close_c_block()
    # cg.write('printf("--\\n")')
    # cg.start_c_block("for (size_t k = 0; k < 4; k++)")
    # cg.write('printf("%lf ", tmp[k])')
    # cg.close_c_block()
    # cg.write('printf("\\n--\\n")')


    # Push block
    cg.write("// Copy data to inner block")
    cg.start_c_block("for (size_t k = 0; k < mremain; k++)")
    cg.write("const size_t start = (mstart + k) * n + nstart")
    cg.write("PRAGMA_VECTORIZE", endl="")
    cg.start_c_block("for (size_t l = 0; l < nremain; l++)")
    # cg.write('printf("(k,l) %ld %ld | %ld\\n", k, l, start+l)')

    cg.write("output[start + l] = tmp[k * %d + l]" % inner_block)
    cg.close_c_block()
    cg.close_c_block()
    # cg.write('printf("--------\\n")')

    # cg.start_c_block("for (size_t k = 0; k < %d; k++)" % inner_block)
    # cg.start_c_block("for (size_t l = 0; l < %d; l++)" % inner_block)
    # cg.write("tmp[k * %d + l] = 0.0" % inner_block)
    # cg.close_c_block()
    # cg.close_c_block()


    # Outer block
    cg.close_c_block()
    cg.close_c_block()

    cg.close_c_block()

    return sig



def generate_hello(path='.'):
    print(path)
    with open(path + '/hello.c', 'w') as fl:
        fl.write("""
/* Hello World program */

#include<stdio.h>

int main()
{
    printf("Hello World");
}
""")
