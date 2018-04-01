"""
Python reference for the collocation matrix and transformation code.
"""
import numpy as np

from . import RSH
from . import docs
from . import order
from . import utility


def collocation_basis(xyz, basis, grad=0, spherical=True, out=None):
    return utility.wrap_basis_collocation(collocation, xyz, basis, grad, spherical, out)


collocation_basis.__doc__ = docs.build_collocation_basis_docs(
    "This function uses a reference NumPy expression as a backed.")


def collocation(xyz, L, coeffs, exponents, center, grad=0, spherical=True, cart_order="row", out=None):

    if grad > 2:
        raise ValueError("Only up to Hessians's of the points (grad = 2) is supported.")

    # Unpack the shell data
    nprim = len(coeffs)
    npoints = xyz.shape[1]

    # First compute the diff distance in each cartesian
    xc = xyz[0] - center[0]
    yc = xyz[1] - center[1]
    zc = xyz[2] - center[2]
    R2 = xc * xc + yc * yc + zc * zc

    # Build up the derivates in each direction
    V1 = np.zeros(npoints)
    V2 = np.zeros(npoints)
    V3 = np.zeros(npoints)
    for K in range(nprim):
        T1 = coeffs[K] * np.exp(-exponents[K] * R2)
        T2 = -2.0 * exponents[K] * T1
        T3 = -2.0 * exponents[K] * T2
        V1 += T1
        V2 += T2
        V3 += T3

    S = V1.copy()
    SX = V2 * xc
    SY = V2 * yc
    SZ = V2 * zc
    SXY = V3 * xc * yc
    SXZ = V3 * xc * zc
    SYZ = V3 * yc * zc
    SXX = V3 * xc * xc + V2
    SYY = V3 * yc * yc + V2
    SZZ = V3 * zc * zc + V2

    # SX, SY, SZ, SXX, SXZ, SXZ, SYY, SYZ, SZZ

    # Power matrix for higher angular momenta
    xc_pow = np.zeros((L + 3, npoints))
    yc_pow = np.zeros((L + 3, npoints))
    zc_pow = np.zeros((L + 3, npoints))

    xc_pow[0] = 0.0
    yc_pow[0] = 0.0
    zc_pow[0] = 0.0
    xc_pow[1] = 0.0
    yc_pow[1] = 0.0
    zc_pow[1] = 0.0
    xc_pow[2] = 1.0
    yc_pow[2] = 1.0
    zc_pow[2] = 1.0

    for LL in range(3, L + 3):
        xc_pow[LL] = xc_pow[LL - 1] * xc
        yc_pow[LL] = yc_pow[LL - 1] * yc
        zc_pow[LL] = zc_pow[LL - 1] * zc

    # Allocate data
    ncart = utility.ncartesian(L)
    nsph = utility.nspherical(L)
    if spherical:
        keys = utility.get_output_keys(grad)
        out = utility.validate_coll_output(grad, (nsph, npoints), out)
        tmps = {k: np.zeros((ncart, npoints)) for k in keys}
    else:
        out = utility.validate_coll_output(grad, (ncart, npoints), out)
        tmps = out

    # Loop over grid ordering data and compute by row
    for idx, l, m, n in order.cartesian_order_factory(L, cart_order):

        # build a few indices
        l = l + 2
        m = m + 2
        n = n + 2

        ld1 = l - 1
        ld2 = l - 2
        md1 = m - 1
        md2 = m - 2
        nd1 = n - 1
        nd2 = n - 2

        A = xc_pow[l] * yc_pow[m] * zc_pow[n]
        AX = ld2 * xc_pow[ld1] * yc_pow[m] * zc_pow[n]
        AY = md2 * xc_pow[l] * yc_pow[md1] * zc_pow[n]
        AZ = nd2 * xc_pow[l] * yc_pow[m] * zc_pow[nd1]

        tmps["PHI"][idx] = S * A
        if grad > 0:
            tmps["PHI_X"][idx] = S * AX + SX * A
            tmps["PHI_Y"][idx] = S * AY + SY * A
            tmps["PHI_Z"][idx] = S * AZ + SZ * A
        if grad > 1:
            AXY = ld2 * md2 * xc_pow[ld1] * yc_pow[md1] * zc_pow[n]
            AXZ = ld2 * nd2 * xc_pow[ld1] * yc_pow[m] * zc_pow[nd1]
            AYZ = md2 * nd2 * xc_pow[l] * yc_pow[md1] * zc_pow[nd1]
            AXX = ld2 * (ld2 - 1) * xc_pow[ld2] * yc_pow[m] * zc_pow[n]
            AYY = md2 * (md2 - 1) * xc_pow[l] * yc_pow[md2] * zc_pow[n]
            AZZ = nd2 * (nd2 - 1) * xc_pow[l] * yc_pow[m] * zc_pow[nd2]
            tmps["PHI_XX"][idx] = SXX * A + SX * AX + SX * AX + S * AXX
            tmps["PHI_YY"][idx] = SYY * A + SY * AY + SY * AY + S * AYY
            tmps["PHI_ZZ"][idx] = SZZ * A + SZ * AZ + SZ * AZ + S * AZZ
            tmps["PHI_XY"][idx] = SXY * A + SX * AY + SY * AX + S * AXY
            tmps["PHI_XZ"][idx] = SXZ * A + SX * AZ + SZ * AX + S * AXZ
            tmps["PHI_YZ"][idx] = SYZ * A + SY * AZ + SZ * AY + S * AYZ

    # Transform results back to spherical
    if spherical:
        for k, v in out.items():
            out[k][:] = RSH.cart_to_spherical_transform(tmps[k], L, cart_order)

    return out


collocation.__doc__ = docs.build_collocation_docs("This function uses a reference NumPy expression as a backed.")
