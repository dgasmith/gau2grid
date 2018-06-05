"""
This is a Python-based automatic generator.
"""

from . import RSH
from . import codegen
from . import docs_generator
from . import order
from . import utility

__built_npcoll_functions = {}


def collocation_basis(xyz, basis, grad=0, spherical=True, out=None, cartesian_order="row", spherical_order="gaussian"):
    return utility.wrap_basis_collocation(collocation, xyz, basis, grad, spherical, out, cartesian_order,
                                          spherical_order)


collocation_basis.__doc__ = docs_generator.build_collocation_basis_docs(
    "This function uses optimized NumPy expressions as a backend.")


def collocation(xyz,
                L,
                coeffs,
                exponents,
                center,
                grad=0,
                spherical=True,
                cartesian_order="row",
                spherical_order="gaussian",
                out=None):

    if grad > 2:
        raise ValueError("Only up to Hessians's of the points (grad = 2) is supported.")

    lookup = "compute_shell_%d_%s_%s" % (L, cartesian_order, spherical_order)
    if lookup not in __built_npcoll_functions:
        exec(numpy_generator(L, lookup, cartesian_order, spherical_order), globals(), __built_npcoll_functions)

    func = __built_npcoll_functions[lookup]
    return func(xyz, L, coeffs, exponents, center, grad=grad, spherical=spherical, out=out)


collocation.__doc__ = docs_generator.build_collocation_docs(
    "This function uses optimized NumPy expressions as a backend.")


def numpy_generator(L, function_name="generated_compute_numpy_shells", cartesian_order="row", spherical_order="gaussian"):
    """
    """

    # Function definition
    cg = codegen.CodeGen()
    cg.write("def %s(xyz, L, coeffs, exponents, center, grad=2, spherical=True, out=None):" % function_name)
    cg.indent()

    cg.blankline()
    cg.write("# Make sure NumPy is in locals")
    cg.write("import numpy as np")
    cg.blankline()

    cg.write("# Unpack shell data")
    cg.write("nprim = len(coeffs)")
    cg.write("npoints = xyz.shape[1]")
    cg.blankline()

    cg.write("# First compute the diff distance in each cartesian")
    cg.write("xc = xyz[0] - center[0]")
    cg.write("yc = xyz[1] - center[1]")
    cg.write("zc = xyz[2] - center[2]")
    cg.write("R2 = xc * xc + yc * yc + zc * zc")
    cg.blankline()

    # All gaussian derivatives
    cg.write("# Build up the derivates in each direction")
    cg.write("V1 = np.zeros((npoints))")
    cg.write("V2 = np.zeros((npoints))")
    cg.write("V3 = np.zeros((npoints))")
    cg.write("for K in range(nprim):")
    cg.write("    T1 = coeffs[K] * np.exp(-exponents[K] * R2)")
    cg.write("    V1 += T1")
    cg.write("    if grad > 0:")
    cg.write("       T2 = -2.0 * exponents[K] * T1")
    cg.write("       V2 += T2")
    cg.write("    if grad > 1:")
    cg.write("       T3 = -2.0 * exponents[K] * T2")
    cg.write("       V3 += T3")

    cg.blankline()
    cg.write("S0 = V1.copy()")
    cg.write("if grad > 0:")
    cg.write("   SX = V2 * xc")
    cg.write("   SY = V2 * yc")
    cg.write("   SZ = V2 * zc")
    cg.write("if grad > 1:")
    cg.write("   SXY = V3 * xc * yc")
    cg.write("   SXZ = V3 * xc * zc")
    cg.write("   SYZ = V3 * yc * zc")
    cg.write("   SXX = V3 * xc * xc + V2")
    cg.write("   SYY = V3 * yc * yc + V2")
    cg.write("   SZZ = V3 * zc * zc + V2")
    cg.blankline()

    # Directional power derivs for angular momenta > 0
    cg.write("# Power matrix for higher angular momenta")
    cg.write("xc_pow = np.zeros((L + 1, npoints))")
    cg.write("yc_pow = np.zeros((L + 1, npoints))")
    cg.write("zc_pow = np.zeros((L + 1, npoints))")
    cg.blankline()
    cg.write("xc_pow[0] = xc")
    cg.write("yc_pow[0] = yc")
    cg.write("zc_pow[0] = zc")

    cg.write("for LL in range(1, L):")
    cg.write("    xc_pow[LL] = xc_pow[LL - 1] * xc")
    cg.write("    yc_pow[LL] = yc_pow[LL - 1] * yc")
    cg.write("    zc_pow[LL] = zc_pow[LL - 1] * zc")
    cg.blankline()

    # Build output data
    cg.write("# Allocate data")
    cg.write("ncart = utility.ncartesian(L)")
    cg.write("nsph = utility.nspherical(L)")
    cg.write("if spherical:")
    cg.write("    keys = utility.get_output_keys(grad)")
    cg.write("    out = utility.validate_coll_output(grad, (nsph, npoints), out)")
    cg.write("    tmps = {k: np.zeros((ncart, npoints)) for k in keys}")
    cg.write("else:")
    cg.write("    out = utility.validate_coll_output(grad, (ncart, npoints), out)")
    cg.write("    tmps = out")
    cg.blankline()

    # Build individual angular moment
    cg.write("# Angular momentum loops")
    for l in range(L + 1):
        cg.write("if L == %d:" % l)
        cg.indent()
        _numpy_am_build(cg, l, cartesian_order)
        cg.dedent()

    cg.write("# If Cartesian were done, return")
    cg.write("if spherical is False:")
    cg.write("    return out")
    cg.blankline()

    # Now spherical transformers
    spherical_func = "spherical_trans"
    for l in range(L + 1):
        RSH.transformation_np_generator(cg, l, cartesian_order, spherical_order, function_name=spherical_func)
        cg.blankline()

    for l in range(L + 1):
        cg.write("if L == %d:" % l)
        cg.indent()
        cg.write("for k, v in out.items():")
        cg.indent()
        cg.write("%s_%d(tmps[k], out=out[k])" % (spherical_func, l))
        cg.dedent(2)
        cg.write("")

    # cg.write("print(np.linalg.norm(out['PHI']))")
    cg.write("return out")

    return cg.repr()


def _numpy_am_build(cg, L, cartesian_order):
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

        name = "X" * ld2 + "Y" * md2 + "Z" * nd2
        if name == "":
            name = "0"

        # Density
        cg.write("# Density AM=%d Component=%s" % (L, name))

        cg.write(_build_xyz_pow("A", 1.0, l, m, n))
        cg.write("tmps['PHI'][%d] = S0 * A" % idx)

        cg.write("# Gradient AM=%d Component=%s" % (L, name))
        cg.write("if grad > 0:")
        cg.indent()

        # Gradient
        cg.write("tmps['PHI_X'][%d] = SX * A" % idx)
        cg.write("tmps['PHI_Y'][%d] = SY * A" % idx)
        cg.write("tmps['PHI_Z'][%d] = SZ * A" % idx)

        AX = _build_xyz_pow("AX", ld2, ld1, m, n)
        if AX is not None:
            x_grad = True
            cg.write(AX)
            cg.write("tmps['PHI_X'][%d] += S0 * AX" % idx)

        AY = _build_xyz_pow("AY", md2, l, md1, n)
        if AY is not None:
            y_grad = True
            cg.write(AY)
            cg.write("tmps['PHI_Y'][%d] += S0 * AY" % idx)

        AZ = _build_xyz_pow("AZ", nd2, l, m, nd1)
        if AZ is not None:
            z_grad = True
            cg.write(AZ)
            cg.write("tmps['PHI_Z'][%d] += S0 * AZ" % idx)
        cg.dedent()

        # Hessian temporaries
        cg.write("# Hessian AM=%d Component=%s" % (L, name))
        cg.write("if grad > 1:")
        cg.indent()

        # S Hess
        # We will build S Hess, grad 1, grad 2, A Hess

        # XX
        cg.write("tmps['PHI_XX'][%d] = SXX * A" % idx)
        if x_grad:
            cg.write("tmps['PHI_XX'][%d] += SX * AX" % idx)
            cg.write("tmps['PHI_XX'][%d] += SX * AX" % idx)

        AXX = _build_xyz_pow("AXX", ld2 * (ld2 - 1), ld2, m, n)
        if AXX is not None:
            rhs = AXX.split(" = ")[-1]
            cg.write("tmps['PHI_XX'][%d] += %s * S0" % (idx, rhs))

        # YY
        cg.write("tmps['PHI_YY'][%d] = SYY * A" % idx)
        if y_grad:
            cg.write("tmps['PHI_YY'][%d] += SY * AY" % idx)
            cg.write("tmps['PHI_YY'][%d] += SY * AY" % idx)
        AYY = _build_xyz_pow("AYY", md2 * (md2 - 1), l, md2, n)
        if AYY is not None:
            rhs = AYY.split(" = ")[-1]
            cg.write("tmps['PHI_YY'][%d] += %s * S0" % (idx, rhs))

        # ZZ
        cg.write("tmps['PHI_ZZ'][%d] = SZZ * A" % idx)
        if z_grad:
            cg.write("tmps['PHI_ZZ'][%d] += SZ * AZ" % idx)
            cg.write("tmps['PHI_ZZ'][%d] += SZ * AZ" % idx)
        AZZ = _build_xyz_pow("AZZ", nd2 * (nd2 - 1), l, m, nd2)
        if AZZ is not None:
            rhs = AZZ.split(" = ")[-1]
            cg.write("tmps['PHI_ZZ'][%d] += %s * S0" % (idx, rhs))

        # XY
        cg.write("tmps['PHI_XY'][%d] = SXY * A" % idx)

        if y_grad:
            cg.write("tmps['PHI_XY'][%d] += SX * AY" % idx)
        if x_grad:
            cg.write("tmps['PHI_XY'][%d] += SY * AX" % idx)

        AXY = _build_xyz_pow("AXY", ld2 * md2, ld1, md1, n)
        if AXY is not None:
            rhs = AXY.split(" = ")[-1]
            cg.write("tmps['PHI_XY'][%d] += %s * S0" % (idx, rhs))

        # XZ
        cg.write("tmps['PHI_XZ'][%d] = SXZ * A" % idx)
        if z_grad:
            cg.write("tmps['PHI_XZ'][%d] += SX * AZ" % idx)
        if x_grad:
            cg.write("tmps['PHI_XZ'][%d] += SZ * AX" % idx)
        AXZ = _build_xyz_pow("AXZ", ld2 * nd2, ld1, m, nd1)
        if AXZ is not None:
            rhs = AXZ.split(" = ")[-1]
            cg.write("tmps['PHI_XZ'][%d] += %s * S0" % (idx, rhs))

        # YZ
        cg.write("tmps['PHI_YZ'][%d] = SYZ * A" % idx)
        if z_grad:
            cg.write("tmps['PHI_YZ'][%d] += SY * AZ" % idx)
        if y_grad:
            cg.write("tmps['PHI_YZ'][%d] += SZ * AY" % idx)
        AYZ = _build_xyz_pow("AYZ", md2 * nd2, l, md1, nd1)
        if AYZ is not None:
            # cg.write(AYZ)
            rhs = AYZ.split(" = ")[-1]
            cg.write("tmps['PHI_YZ'][%d] += %s * S0" % (idx, rhs))
        cg.dedent()

        idx += 1
        cg.write(" ")


def _build_xyz_pow(name, pref, l, m, n, shift=2):
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

    if l > 0:
        ret += mul + "xc_pow[%d]" % (l - 1)
        mul = " * "

    if m > 0:
        ret += mul + "yc_pow[%d]" % (m - 1)
        mul = " * "

    if n > 0:
        ret += mul + "zc_pow[%d]" % (n - 1)
        mul = " * "

    if mul == " ":
        ret += " 1"

    return ret
