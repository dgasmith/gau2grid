"""
The C generator for gau2grid collocation functions
"""

from . import order
from . import RSH
from . import codegen

_grads = ["x", "y", "z"]
_hessians = ["xx", "xy", "xz", "yy", "yz", "zz"]
def c_generator(L, function_name="", grad=0, cart_order="row", inner_block=64):

    if function_name == "":
        function_name = "coll_%d_%d" % (L, grad)

    cg = codegen.CodeGen(cgen=True)

    # Precompute temps
    ncart = int((L + 1) * (L + 2) / 2)
    nspherical = L * 2 + 1

    # Build function signature
    if grad == 0:
        func_sig = "size_t npoints, double* x, double* y, double* z, int nprim, double* coeffs, double* exponents, double* center, bool spherical, double* ret"
    else:
        raise KeyError("Grad larger than 2 is not yet implemented.")
    cg.start_c_block("void %s(%s)" % (function_name, func_sig))
    cg.blankline()

    # Figure out spacing
    cg.write("// Sizing")
    cg.write("size_t nblocks = npoints / %d" % inner_block)
    cg.write("nblocks += (nblocks %% %d) ? 0 : 1" % inner_block)
    cg.write("size_t ncart = %d" % ncart)
    cg.blankline()

    # Build temporaries
    cg.write("// Allocate temporaries")
    S_tmps = ["xc", "yz", "zc", "R2", "S0"]
    for tname in S_tmps:
        cg.write("double*  %s = double[%d]" % (tname, inner_block))

    L_tmps = ["xc_pow", "yz_pow", "zc_pow"]
    for tname in L_tmps:
        cg.write("double*  %s = double[%d]" % (tname, inner_block * (L + 1)))

    inner_tmps = ["rho_tmp"]
    for tname in inner_tmps:
        cg.write("double*  %s = double[%d]" % (tname, inner_block * ncart))
    cg.blankline()

    # Start outer loop
    cg.write("// Start outer block loop")
    cg.start_c_block("for (size_t block = 0; block < nblocks; block++)")

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

    # Start inner loop
    cg.write("// Start inner block loop")
    cg.start_c_block("for (size_t i = 0; i < %d; i++)" % inner_block)

    cg.blankline()
    cg.write("// Position temps")
    cg.write("R2[i] = xc[i] * xc[i] + yc[i] * yc[i] + zc[i] * zc[i]")
    cg.blankline()


    cg.blankline()
    cg.write("// Deriv tmps")
    cg.start_c_block("for (size_t n = 0; n < nprim; n++)")
    cg.write("double T1 = coeffs[n] * exp(-1.0 * exponents[n] * R2[i])")
    cg.write("S[i] += T1")
    cg.blankline()

    cg.blankline()
    cg.write("// Power tmps")
    cg.write("xc_pow[i] = xc[i]")
    cg.write("yc_pow[i] = yc[i]")
    cg.write("zc_pow[i] = zc[i]")

    for l in range(1, L):
        cg.write("xc_pow[%d + i] = xc[%d + i]" % (inner_block * l, inner_block * (l - 1)))
        cg.write("xc_pow[%d + i] = xc[%d + i]" % (inner_block * l, inner_block * (l - 1)))
        cg.write("xc_pow[%d + i] = xc[%d + i]" % (inner_block * l, inner_block * (l - 1)))
    cg.blankline()


    cg.close_c_block()
    cg.blankline()

    # End inner loop
    cg.close_c_block()

    # End outer loop
    cg.close_c_block()


    cg.write("// Free temporaries")
    for tname in (S_tmps + L_tmps + inner_tmps):
        cg.write("free %s[]" % tname)
    cg.blankline()

    cg.close_c_block()

    # return cg.repr()
    return cg.repr(clang_format=True)
