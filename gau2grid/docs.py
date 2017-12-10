"""
Contains several docstrings as there are several duplicate functions
"""

def build_collocation_docs(insert=""):

    doc_header =  """
    Computes the collocation matrix for a given set of cartesian points and a contracted gaussian of the form:
        \sum_i coeff_i e^(exponent_i * R^2)

    """

    param_data = """

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
    grad : int, optional (default: 0)
        Can return cartesian gradient and Hessian per point if requested.
    spherical : bool, optional (default: True)
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

    ret = doc_header
    if insert == "":
        ret += "\n"
    else:
        ret += insert

    ret += param_data
    return ret


def build_collocation_basis_docs(insert=""):
    doc_header = """
    Computes the collocation matrix for a given gaussian basis of the form:
        x^l y^m z^n \sum_i coeff_i e^(exponent_i * R^2)
    Where for a given angular momentum all combinations of l + m + n = L are explored.

    """

    param_data = """

    Expects the basis to take the form of:
        [L, coeffs, exponents, center]
    xyz : array_like
        The (3, N) cartesian points to compute the grid on
    basis : list of tuples
        Each tuple should represent a basis in the form of (L, coeffs, exponents, center).
        L : int
            The angular momentum of the gaussian
        coeffs : array_like
            The coefficients of the gaussian
        exponents : array_like
            The exponents of the gaussian
        center : array_like
            The cartesian center of the gaussian
    grad : int, default=0
        Can return cartesian gradient and Hessian per point if requested.
    spherical : bool, default=True
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

    ret = doc_header
    if insert == "":
        ret += "\n"
    else:
        ret += insert

    ret += param_data
    return ret


