"""
A Python wrapper for the compiled GG functions.
"""

import ctypes
import ctypes.util
import numpy as np

from . import utility

# Attempt to load the compiled C code
lib_path = ctypes.util.find_library('libgg')
__lib_found = False
print(lib_path)


def _build_collocation_ctype(nout):
    """
    Builds the ctypes signatures for the libgg C lib
    """
    ret = [
        # L, npoints
        ctypes.c_int,
        ctypes.c_ulong,

        # XYZ
        np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("C", "A")),
        np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("C", "A")),
        np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("C", "A")),

        # Gaussian
        ctypes.c_int,
        np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("C", "A")),  # coef
        np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags=("C", "A")),  # exp
        np.ctypeslib.ndpointer(dtype=np.double, shape=(3, ), ndim=1, flags=("C", "A")),  # center

        # Spherical
        ctypes.c_int,
    ]

    # Pushback output
    for n in range(nout):
        ret.append(np.ctypeslib.ndpointer(dtype=np.double, ndim=2, flags=("W", "C", "A")))

    return tuple(ret)


# Bind the C object
if lib_path is not None:
    __lib_found = True

    # Bind the obect
    cgg = ctypes.CDLL(lib_path)

    # Transposes
    cgg.gg_naive_transpose.restype = None
    cgg.gg_naive_transpose.argtypes = (ctypes.c_ulong, ctypes.c_ulong, np.ctypeslib.ndpointer(),
                                       np.ctypeslib.ndpointer())

    cgg.gg_fast_transpose.restype = None
    cgg.gg_fast_transpose.argtypes = (ctypes.c_ulong, ctypes.c_ulong, np.ctypeslib.ndpointer(),
                                      np.ctypeslib.ndpointer())

    # Collocation generators
    cgg.gg_collocation.restype = None
    cgg.gg_collocation.argtypes = _build_collocation_ctype(1)

    cgg.gg_collocation.restype = None
    cgg.gg_collocation_deriv1.argtypes = _build_collocation_ctype(4)

    cgg.gg_collocation.restype = None
    cgg.gg_collocation_deriv2.argtypes = _build_collocation_ctype(10)


def c_compiled():
    """
    Checks if the c code has been compiled or not.
    """
    return __lib_found


def _validate_c_import():
    if __lib_found is False:
        raise ImportError("Compiled libgg not found. Please compile gau2grid before calling these functions.")


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

    if L > cgg.max_L():
        raise ValueError("LibGG was only compiled to AM=%d, requested AM=%d." % (cgg.max_L(), L))

    # Check XYZ
    if xyz.shape[0] != 3:
        raise ValueError("XYZ array must be of shape (3, N), found %s" % str(xyz.shape))

    # Check gaussian
    coeffs = np.asarray(coeffs, dtype=np.double)
    exponents = np.asarray(exponents, dtype=np.double)
    center = np.asarray(center, dtype=np.double)
    if coeffs.shape[0] != exponents.shape[0]:
        raise ValueError("Coefficients (N=%d) and exponents (N=%d) must have the same shape." % (coeffs.shape,
                                                                                                 exponents.shape[0]))

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
        cgg.gg_collocation(L, xyz.shape[1], xyz[0], xyz[1], xyz[2], coeffs.shape[0], coeffs, exponents, center,
                           spherical, out["PHI"])
    elif grad == 1:
        cgg.gg_collocation_deriv1(L, xyz.shape[1], xyz[0], xyz[1], xyz[2], coeffs.shape[0], coeffs, exponents, center,
                                  spherical, out["PHI"], out["PHI_X"], out["PHI_Y"], out["PHI_Z"])
    elif grad == 2:
        cgg.gg_collocation_deriv2(L, xyz.shape[1], xyz[0], xyz[1], xyz[2], coeffs.shape[0], coeffs, exponents, center,
                                  spherical, out["PHI"], out["PHI_X"], out["PHI_Y"], out["PHI_Z"], out["PHI_XX"],
                                  out["PHI_XY"], out["PHI_XZ"], out["PHI_YY"], out["PHI_YZ"], out["PHI_ZZ"])
    else:
        raise ValueError("Only up to Hessians's of the points (grad = 2) is supported.")

    return out