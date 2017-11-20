"""
A Python wrapper for the compiled GG functions.
"""

from . import utility

# Attempt to load the compiled C code
try:
    import pygg_core as ggc
    __core_found = True
except ImportError:
    __core_found = False


def c_compiled():
    """
    Checks if the c code has been compiled or not.
    """
    return __core_found


def _validate_c_import():
    if __core_found is False:
        raise ImportError("Compiled pygg_core not found. Please compile gau2grid before calling these functions.")


def collocation_basis(xyz, basis, grad=0, spherical=True, out=None):

    return utility.wrap_basis_collocation(collocation, xyz, basis, grad, spherical, out)


def collocation(xyz, L, coeffs, exponents, center, grad=0, spherical=True, out=None):
    """
    Computes the collocation matrix for a given set of cartesian points and a contracted gaussian of the form:
        \sum_i coeff_i e^(exponent_i * R^2)

    This function builds a optimized NumPy version on the fly and caches it for future use.

    Parameters
    ----------
    xyz : array_like
        The (3, N) cartesian points to compute the grid on
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
    out : dict, optional
        A dictionary of output NumPy arrays to write the data to.

    Returns
    -------
    ret : dict of array_like
        Returns a dictionary containing the requested arrays (PHI, PHI_X, PHI_XX, etc).
        Where each matrix is of shape (ngaussian_basis x npoints)
        The cartesian center of the gaussian

    """

    # Validates we loaded correctly
    _validate_c_import()

    # Find the output shape
    npoints = xyz.shape[1]
    if spherical:
        nvals = utility.nspherical(L)
    else:
        nvals = utility.ncartesian(L)

    # Build the outputs
    out = utility.validate_coll_output(grad, (nvals, npoints), out)

    # Select the correct function
    if grad == 0:
        ggc.collocation(L, xyz, coeffs, exponents, center, spherical, out["PHI"])
    elif grad == 1:
        ggc.collocation_deriv1(L, xyz, coeffs, exponents, center, spherical, out["PHI"], out["PHI_X"], out["PHI_Y"],
                               out["PHI_Z"])
    elif grad == 2:
        ggc.collocation_deriv2(L, xyz, coeffs, exponents, center, spherical, out["PHI"], out["PHI_X"], out["PHI_Y"],
                               out["PHI_Z"], out["PHI_XX"], out["PHI_XY"], out["PHI_XZ"], out["PHI_YY"], out["PHI_YZ"],
                               out["PHI_ZZ"])
    else:
        raise ValueError("Only up to Hessians's of the points (grad = 2) is supported.")

    return out