"""
Python reference for the collocation matrix and transformation code.
"""
import numpy as np

from . import order
from . import RSH


def compute_collocation(xyz, L, coeffs, exponents, center, grad=0, spherical=True, cart_order="row"):
    """
    Computes the collocation matrix for a given set of cartesian points and a contracted gaussian of the form:
        \sum_i coeff_i e^(exponent_i * R^2)


    Parameters
    ----------
    xyz : array_like
        The (N, 3) cartesian points to compute the grid on
    L : int
        The angular momentum of the gaussian
    coeffs : array_like
        The coefficients of the gaussian
    exponents : array_like
        The exponents of the gaussian
    center : array_like
        The cartesian center of the gaussian
    grad : int
        Can return cartesian gradient and Hessian per point if requested.
    spherical : bool
        Transform the resulting cartesian gaussian to spherical
    cart_order : str
        The order of the resulting cartesian basis, no effect if spherical=True

    Returns
    -------
    ret : dict of array_like
        Returns a dictionary containing the requested arrays (PHI, PHI_X, PHI_XX, etc).
        Where each matrix is of shape (ngaussian_basis x npoints)


    """

    if grad > 2:
        raise ValueError("Only up to Hessians's of the points (grad = 2) is supported.")

    # Unpack the shell data
    nprim = len(coeffs)
    npoints = xyz.shape[0]

    # First compute the diff distance in each cartesian
    xc = xyz[:, 0] - center[0]
    yc = xyz[:, 1] - center[1]
    zc = xyz[:, 2] - center[2]
    R2 = xc * xc + yc * yc + zc * zc

    # Build up the derivates in each direction
    V1 = np.zeros((npoints))
    V2 = np.zeros((npoints))
    V3 = np.zeros((npoints))
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
    ncart = int((L + 1) * (L + 2) / 2)
    output = {}
    output["PHI"] = np.zeros((ncart, npoints))
    if grad > 0:
        output["PHI_X"] = np.zeros((ncart, npoints))
        output["PHI_Y"] = np.zeros((ncart, npoints))
        output["PHI_Z"] = np.zeros((ncart, npoints))
    if grad > 1:
        output["PHI_XX"] = np.zeros((ncart, npoints))
        output["PHI_YY"] = np.zeros((ncart, npoints))
        output["PHI_ZZ"] = np.zeros((ncart, npoints))
        output["PHI_XY"] = np.zeros((ncart, npoints))
        output["PHI_XZ"] = np.zeros((ncart, npoints))
        output["PHI_YZ"] = np.zeros((ncart, npoints))
    if grad > 2:
        raise ValueError("Only grid derivatives through Hessians (grad = 2) has been implemented")

    # Loop over grid ordering data
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

        A = xc_pow[l] * yc_pow[m] * zc_pow[n]
        AX = ld2 * xc_pow[ld1] * yc_pow[m] * zc_pow[n]
        AY = md2 * xc_pow[l] * yc_pow[md1] * zc_pow[n]
        AZ = nd2 * xc_pow[l] * yc_pow[m] * zc_pow[nd1]

        output["PHI"][idx] = S * A
        if grad > 0:
            output["PHI_X"][idx] = S * AX + SX * A
            output["PHI_Y"][idx] = S * AY + SY * A
            output["PHI_Z"][idx] = S * AZ + SZ * A
        if grad > 1:
            AXY = ld2 * md2 * xc_pow[ld1] * yc_pow[md1] * zc_pow[n]
            AXZ = ld2 * nd2 * xc_pow[ld1] * yc_pow[m] * zc_pow[nd1]
            AYZ = md2 * nd2 * xc_pow[l] * yc_pow[md1] * zc_pow[nd1]
            AXX = ld2 * (ld2 - 1) * xc_pow[ld2] * yc_pow[m] * zc_pow[n]
            AYY = md2 * (md2 - 1) * xc_pow[l] * yc_pow[md2] * zc_pow[n]
            AZZ = nd2 * (nd2 - 1) * xc_pow[l] * yc_pow[m] * zc_pow[nd2]
            output["PHI_XX"][idx] = SXX * A + SX * AX + SX * AX + S * AXX
            output["PHI_YY"][idx] = SYY * A + SY * AY + SY * AY + S * AYY
            output["PHI_ZZ"][idx] = SZZ * A + SZ * AZ + SZ * AZ + S * AZZ
            output["PHI_XY"][idx] = SXY * A + SX * AY + SY * AX + S * AXY
            output["PHI_XZ"][idx] = SXZ * A + SX * AZ + SZ * AX + S * AXZ
            output["PHI_YZ"][idx] = SYZ * A + SY * AZ + SZ * AY + S * AYZ

    if spherical:
        for k, v in output.items():
            tmp = output[k].shape
            output[k] = RSH.cart_to_spherical_transform(output[k], L, cart_order)

    return output

