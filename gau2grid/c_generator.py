"""
The C generator for gau2grid collocation functions
"""

import os
from . import order
from . import RSH
from . import codegen

_grads = ["x", "y", "z"]
_hessians = ["xx", "xy", "xz", "yy", "yz", "zz"]


def generate_c_gau2grid(max_L, path=".", cart_order="row", inner_block=64):
    print(path)

    gg_header = codegen.CodeGen(cgen=True)
    gg_phi = codegen.CodeGen(cgen=True)
    gg_grad = codegen.CodeGen(cgen=True)
    gg_hess = codegen.CodeGen(cgen=True)

    # Add general header comments
    for cgs in [gg_header, gg_phi, gg_grad, gg_hess]:
        cgs.write("// This is an automtically generated file from ...")
        cgs.write("// Blah blah blah")
        cgs.blankline()

    # Add utility headers
    for cgs in [gg_phi, gg_grad, gg_hess]:
        cgs.write("#include <math.h>")
        cgs.write("#include <stdbool.h>")
        cgs.write("#include <stdlib.h>")
        cgs.blankline()
        cgs.write("// Adds a few typedefs to make the world easier")
        cgs.write("typedef unsigned long size_t;")
        cgs.blankline()

    # Loop over phi, grad, hess
    for name, grad, cg in [("Phi", 0, gg_phi), ("Phi grad", 1, gg_grad), ("Phi Hess", 2, gg_hess)]:
        gg_header.write("// %s computers")
        cgs.blankline()
        # Write out the phi builders
        for L in range(max_L + 1):
            sig = shell_c_generator(cg, L, grad=grad, cart_order=cart_order, inner_block=inner_block)
            cg.blankline()

            # Write out the header data
            gg_header.write(sig)
            gg_header.blankline()

    gg_header.repr(filename=os.path.join(path, "gau2grid.h"), clang_format=True)
    gg_phi.repr(filename=os.path.join(path, "gau2grid_phi.cc"), clang_format=True)
    gg_grad.repr(filename=os.path.join(path, "gau2grid_phi_grad.cc"), clang_format=True)
    gg_hess.repr(filename=os.path.join(path, "gau2grid_phi_hess.cc"), clang_format=True)


def shell_c_generator(cg, L, function_name="", grad=0, cart_order="row", inner_block=64):

    # Grab the line start
    cg_line_start = len(cg.data)

    # Parse Keywords
    if function_name == "":
        function_name = "coll_L%d_deriv%d" % (L, grad)

    if grad > 2:
        raise TypeError("Only grad <=2 (Hessians) is supported")

    # Precompute temps
    ncart = int((L + 1) * (L + 2) / 2)
    nspherical = L * 2 + 1
    grad_indices = ["x", "y", "z"]
    hess_indices = ["xx", "xy", "xz", "yy", "yz", "zz"]

    # Build function signature
    func_sig = "size_t npoints, double* x, double* y, double* z, int nprim, double* coeffs, double* exponents, double* center, bool spherical, double* phi_out"
    if grad > 0:
        func_sig += ", "
        func_sig += ", ".join("double* phi_%s_out" % grad for grad in grad_indices)
    if grad > 1:
        func_sig += ", "
        func_sig += ", ".join("double* phi_%s_out" % hess for hess in hess_indices)

    func_sig = "void %s(%s)" % (function_name, func_sig)
    cg.start_c_block(func_sig)
    cg.blankline()

    # Figure out spacing
    cg.write("// Sizing")
    cg.write("size_t nblocks = npoints / %d" % inner_block)
    cg.write("nblocks += (nblocks %% %d) ? 0 : 1" % inner_block)
    cg.write("size_t ncart = %d" % ncart)
    cg.write("size_t nspherical = %d" % nspherical)
    cg.write("size_t nout = spherical ? nspherical : ncart")
    cg.blankline()

    # Build temporaries
    cg.write("// Allocate S temporaries")
    S_tmps = ["xc", "yc", "zc", "R2", "S0"]
    if grad > 0:
        S_tmps += ["S1"] + ["S%s" % grad.upper() for grad in grad_indices]
    if grad > 1:
        S_tmps += ["S2"] + ["S%s" % hess.upper() for hess in hess_indices]
    for tname in S_tmps:
        cg.write("double*  %s = (double*)malloc(%d * sizeof(double))" % (tname, inner_block))
    cg.blankline()

    power_tmps = []
    if L > 1:
        cg.write("// Allocate power temporaries")
        power_tmps = ["xc_pow", "yc_pow", "zc_pow"]
        for tname in power_tmps:
            cg.write("double*  %s = (double*)malloc(%d * sizeof(double))" % (tname, inner_block * L))
        cg.blankline()

    cg.write("// Allocate output temporaries")
    inner_tmps = ["phi_tmp"]
    if grad > 0:
        inner_tmps += ["phi_%s_tmp" % grad for grad in grad_indices]
    if grad > 1:
        inner_tmps += ["phi_%s_tmp" % hess for hess in hess_indices]
    for tname in inner_tmps:
        cg.write("double*  %s = (double*)malloc(%d * sizeof(double))" % (tname, inner_block * ncart))
    cg.blankline()

    # Any declerations needed
    cg.write("// Declare doubles")
    cg.write("double A")
    if grad > 0:
        cg.write("double " + ", ".join("A%s" % grad.upper() for grad in grad_indices))
    if grad > 1:
        cg.write("double " + ", ".join("A%s" % hess.upper() for hess in hess_indices))
    cg.blankline()

    # Start outer loop
    cg.write("// Start outer block loop")
    cg.start_c_block("for (size_t block = 0; block < nblocks; block++)")
    cg.blankline()

    # Move data into inner buffers
    cg.blankline()
    cg.write("// Copy data into inner temps")
    cg.write("size_t start = block * %d" % inner_block)
    cg.write("size_t remain = ((start + %d) > npoints) ? (npoints - start) : %d" % (inner_block, inner_block))

    cg.start_c_block("for (size_t i = 0; i < remain; i++)")
    cg.write("xc[i] = x[start + i] - center[0]")
    cg.write("yc[i] = y[start + i] - center[1]")
    cg.write("zc[i] = z[start + i] - center[2]")
    cg.close_c_block()
    cg.blankline()

    # Grab the inner line start
    inner_line_start = len(cg.data)

    # Start inner loop
    cg.write("// Start inner block loop")
    cg.start_c_block("for (size_t i = 0; i < %d; i++)" % inner_block)

    # Build R2
    cg.blankline()
    cg.write("// Position temps")
    cg.write("R2[i] = xc[i] * xc[i] + yc[i] * yc[i] + zc[i] * zc[i]")
    cg.blankline()

    # Build out thoese gaussian derivs
    cg.blankline()
    cg.write("// Deriv tmps")
    cg.start_c_block("for (size_t n = 0; n < nprim; n++)")
    cg.write("double T1 = coeffs[n] * exp(-1.0 * exponents[n] * R2[i])")
    cg.write("S0[i] += T1")
    if grad > 0:
        cg.write("double T2 = -2.0 * exponents[n] * T1")
        cg.write("S1[i] += T2")
    if grad > 1:
        cg.write("double T3 = -2.0 * exponents[n] * T2")
        cg.write("S2[i] += T3")

    cg.close_c_block()
    cg.blankline()

    if grad > 0:
        cg.write("// S derivs")
        cg.write("SX[i] = S1[i] * xc[i]")
        cg.write("SY[i] = S1[i] * yc[i]")
        cg.write("SZ[i] = S1[i] * zc[i]")
    if grad > 1:
        cg.blankline()
        cg.write("// S Hessians")
        cg.write("SXY[i] = S2[i] * xc[i] * yc[i]")
        cg.write("SXZ[i] = S2[i] * xc[i] * zc[i]")
        cg.write("SYZ[i] = S2[i] * yc[i] * zc[i]")
        cg.write("SXX[i] = S2[i] * xc[i] * xc[i] + S1[i]")
        cg.write("SYY[i] = S2[i] * yc[i] * yc[i] + S1[i]")
        cg.write("SZZ[i] = S2[i] * zc[i] * zc[i] + S1[i]")

    # Build out those power derivs
    cg.blankline()
    if L > 1:
        cg.write("// Power tmps")
        cg.write("xc_pow[i] = xc[i] * xc[i]")
        cg.write("yc_pow[i] = yc[i] * yc[i]")
        cg.write("zc_pow[i] = zc[i] * zc[i]")
    if L == 2:
        cg.blankline()

    for l in range(2, L):
        cg.write("xc_pow[%d + i] = xc_pow[%d + i] * xc[i]" % (inner_block * l, inner_block * (l - 1)))
        cg.write("yc_pow[%d + i] = yc_pow[%d + i] * yc[i]" % (inner_block * l, inner_block * (l - 1)))
        cg.write("zc_pow[%d + i] = zc_pow[%d + i] * zc[i]" % (inner_block * l, inner_block * (l - 1)))
        cg.blankline()

    # Contract temps with powers
    _c_am_build(cg, L, cart_order, grad, inner_block)

    cg.blankline()

    # End inner loop
    cg.close_c_block()

    # Grab the inner line stop
    inner_line_stop = len(cg.data)

    # Move data into inner buffers
    cg.blankline()
    cg.write("// Copy data back into outer temps")
    cg.start_c_block("for (size_t n = 0; n < nout; n++)")

    # Copy over Phi
    cg.write("// Phi, copy data to outer temps")
    cg.start_c_block("for (size_t i = 0; i < remain; i++)")
    cg.write("phi_out[start * nout + i] = phi_tmp[%d * n + i]" % (inner_block))
    cg.close_c_block()

    # Copy over grad
    if grad > 0:
        cg.blankline()
        cg.write("// Grad, copy data to outer temps")
        for ind in grad_indices:
            cg.start_c_block("for (size_t i = 0; i < remain; i++)")
            cg.write("phi_%s_out[start * nout + i] = phi_%s_tmp[%d * n + i]" % (ind, ind, inner_block))
            cg.close_c_block()
        cg.blankline()

    if grad > 1:
        cg.write("// Hess, copy data to outer temps")
        for ind in hess_indices:
            cg.start_c_block("for (size_t i = 0; i < remain; i++)")
            cg.write("phi_%s_out[start * nout + i] = phi_%s_tmp[%d * n + i]" % (ind, ind, inner_block))
            cg.close_c_block()
        cg.blankline()

    cg.close_c_block()
    cg.blankline()

    # End outer loop
    cg.close_c_block()

    # Free up those arrays
    cg.blankline()
    for name, flist in [("S", S_tmps), ("power", power_tmps), ("inner", inner_tmps)]:
        if (len(flist) == 0): continue

        cg.write("// Free %s temporaries" % name)
        for tname in flist:
            cg.write("free(%s)" % tname)
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
        if ("=" not in line) or ("//" in line):
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
            if (k + ";" in tmp):
                cg.data[pos] = line.replace(k + ";", v + ";")
            elif (k + " " in tmp):
                cg.data[pos] = line.replace(k + " ", v + " ")
        pos += 1

    # Remove any " * 1"
    for x in range(cg_line_start, len(cg.data)):
        cg.data[x] = cg.data[x].replace(" * 1;", ";")
        cg.data[x] = cg.data[x].replace(" * 1.0;", ";")
        cg.data[x] = cg.data[x].replace("= 1 * ", "= ")
        cg.data[x] = cg.data[x].replace("= 1.0 * ", "= ")

    return func_sig


def _c_am_build(cg, L, cart_order, grad, shift):
    """
    Builds a unrolled angular momentum function
    """
    names = ["X", "Y", "Z"]

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
        tmp_ret = []

        # Set grads back to zero
        x_grad, y_grad, z_grad = False, False, False
        shift_idx = idx * shift

        name = "X" * ld2 + "Y" * md2 + "Z" * nd2
        if name == "":
            name = "0"

        # Density
        cg.write("// Density AM=%d Component=%s" % (L, name))

        cg.write(_build_xyz_pow("A", 1.0, l, m, n, shift))
        cg.write("phi_tmp[%d + i] = S0[i] * A" % shift_idx)

        if grad == 0: continue
        cg.write("// Gradient AM=%d Component=%s" % (L, name))

        # Gradient
        cg.write("phi_x_tmp[%d + i] = SX[i] * A" % shift_idx)
        cg.write("phi_y_tmp[%d + i] = SY[i] * A" % shift_idx)
        cg.write("phi_z_tmp[%d + i] = SZ[i] * A" % shift_idx)

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
        cg.write("// Hessian AM=%d Component=%s" % (L, name))

        # S Hess
        # We will build S Hess, grad 1, grad 2, A Hess

        # XX
        cg.write("phi_xx_tmp[%d + i] = SXX[i] * A" % shift_idx)
        if x_grad:
            cg.write("phi_xx_tmp[%d + i] += SX[i] * AX" % shift_idx)
            cg.write("phi_xx_tmp[%d + i] += SX[i] * AX" % shift_idx)

        AXX = _build_xyz_pow("AXX", ld2 * (ld2 - 1), ld2, m, n, shift)
        if AXX is not None:
            rhs = AXX.split(" = ")[-1]
            cg.write("phi_xx_tmp[%d + i] += %s * S0[i]" % (shift_idx, rhs))

        # YY
        cg.write("phi_yy_tmp[%d + i] = SYY[i] * A" % shift_idx)
        if y_grad:
            cg.write("phi_yy_tmp[%d + i] += SY[i] * AY" % shift_idx)
            cg.write("phi_yy_tmp[%d + i] += SY[i] * AY" % shift_idx)
        AYY = _build_xyz_pow("AYY", md2 * (md2 - 1), l, md2, n, shift)
        if AYY is not None:
            rhs = AYY.split(" = ")[-1]
            cg.write("phi_yy_tmp[%d + i] += %s * S0[i]" % (shift_idx, rhs))

        # ZZ
        cg.write("phi_zz_tmp[%d + i] = SZZ[i] * A" % shift_idx)
        if z_grad:
            cg.write("phi_zz_tmp[%d + i] += SZ[i] * AZ" % shift_idx)
            cg.write("phi_zz_tmp[%d + i] += SZ[i] * AZ" % shift_idx)
        AZZ = _build_xyz_pow("AZZ", nd2 * (nd2 - 1), l, m, nd2, shift)
        if AZZ is not None:
            rhs = AZZ.split(" = ")[-1]
            cg.write("phi_zz_tmp[%d + i] += %s * S0[i]" % (shift_idx, rhs))

        # XY
        cg.write("phi_xy_tmp[%d + i] = SXY[i] * A" % shift_idx)

        if y_grad:
            cg.write("phi_xy_tmp[%d + i] += SX[i] * AY" % shift_idx)
        if x_grad:
            cg.write("phi_xy_tmp[%d + i] += SY[i] * AX" % shift_idx)

        AXY = _build_xyz_pow("AXY", ld2 * md2, ld1, md1, n, shift)
        if AXY is not None:
            rhs = AXY.split(" = ")[-1]
            cg.write("phi_xy_tmp[%d + i] += %s * S0[i]" % (shift_idx, rhs))

        # XZ
        cg.write("phi_xz_tmp[%d + i] = SXZ[i] * A" % shift_idx)
        if z_grad:
            cg.write("phi_xz_tmp[%d + i] += SX[i] * AZ" % shift_idx)
        if x_grad:
            cg.write("phi_xz_tmp[%d + i] += SZ[i] * AX" % shift_idx)
        AXZ = _build_xyz_pow("AXZ", ld2 * nd2, ld1, m, nd1, shift)
        if AXZ is not None:
            rhs = AXZ.split(" = ")[-1]
            cg.write("phi_xz_tmp[%d + i] += %s * S0[i]" % (shift_idx, rhs))

        # YZ
        cg.write("phi_yz_tmp[%d + i] = SYZ[i] * A" % shift_idx)
        if z_grad:
            cg.write("phi_yz_tmp[%d + i] += SY[i] * AZ" % shift_idx)
        if y_grad:
            cg.write("phi_yz_tmp[%d + i] += SZ[i] * AY" % shift_idx)
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
        # If the power is greater than 1 we need to use (xc_pow - 1) as we start at 2
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
