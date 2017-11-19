"""
Compares the Psi4 collocation grids against the NumPy reference code.
"""

import time
import pytest
import numpy as np

import gau2grid as gg

import test_helper as th

np.set_printoptions(linewidth=120, suppress=True)

# Build the test molecule
HeC_mol = """
He 0 0 0
C 0 0 2
no_com
no_reorient
"""

# Tweakers
npoints = int(1.e2)

# Global points
np.random.seed(0)
xyzw = np.random.rand(4, npoints)


def _build_psi4_basis(mol, basis, puream=False):
    """
    Builds a Psi4 basis and its Python interpretation from a molecule string and basis string.
    """

    import psi4
    psi4.set_output_file("output.dat")

    # Can only handle cartesian data
    psi4.set_options({"PUREAM": puream})

    mol = psi4.geometry(mol)
    mol.update_geometry()
    geom = np.array(mol.geometry())

    # Convert the basis to a dictionary
    basis = psi4.core.BasisSet.build(mol, "orbital", basis, puream=puream)
    py_basis = []
    for x in range(basis.nshell()):
        shell = basis.shell(x)
        tmp = {}
        tmp["center"] = geom[shell.ncenter]

        tmp["exp"] = [shell.exp(n) for n in range(shell.nprimitive)]
        tmp["coef"] = [shell.coef(n) for n in range(shell.nprimitive)]
        tmp["am"] = shell.am
        py_basis.append(tmp)

    return basis, py_basis


def _compute_psi4_points(xyzw, basis, grad=2, puream=False):
    """
    Computes the Psi4 collocation matrices
    """

    import psi4

    # Build up an exact block
    extents = psi4.core.BasisExtents(basis, 1.e-100)
    block = psi4.core.BlockOPoints(
        psi4.core.Vector.from_array(xyzw[0]),
        psi4.core.Vector.from_array(xyzw[1]),
        psi4.core.Vector.from_array(xyzw[2]), psi4.core.Vector.from_array(xyzw[3]), extents)

    # Compute BasisFunctions on a grid
    p4_points = psi4.core.BasisFunctions(basis, npoints, basis.nbf())
    p4_points.set_deriv(grad)
    p4_points.compute_functions(block)

    psi_results = {k: np.array(v).T for k, v in p4_points.basis_values().items()}

    return psi_results


# Build up a list of tests
psi_tests = []
for basis in ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "cc-pV5Z", "cc-pV6Z"]:
    for spherical in ["cartesian", "spherical"]:
        psi_tests.append((basis, spherical))


@th.using_psi4_libxc
@pytest.mark.parametrize("basis,spherical", psi_tests)
def test_psi_collocation(basis, spherical):
    trans = "spherical" == spherical

    psi_basis, py_basis = _build_psi4_basis(HeC_mol, basis, puream=trans)

    t = time.time()
    psi_results = _compute_psi4_points(xyzw, psi_basis, puream=trans)
    psi_time = time.time() - t

    t = time.time()
    gg_results = th.compute_points_block(gg.ref.compute_collocation, xyzw, py_basis, spherical=trans)
    gg_time = time.time() - t

    # Print time with py.test -s flags
    print("")
    print("%s-%s time PSI: %8.4f GG: %8.4f" % (basis, spherical, psi_time, gg_time))

    # Print test blocks
    # for x in py_basis:
    #     print(x["am"])
    # print(psi_results["PHI"][-5:])
    # print(gg_results["PHI"][-5:])
    # print((psi_results["PHI"] / gg_results["PHI"])[-5:])
    # print(np.abs(psi_results["PHI"] - gg_results["PHI"])[-5:] < 1.e-10)

    th.compare_collocation_results(gg_results, psi_results)


@th.using_psi4_libxc
@pytest.mark.parametrize("grad", [0, 1, 2])
def test_psi_derivs(grad):
    psi_basis, py_basis = _build_psi4_basis(HeC_mol, "cc-pV6Z", puream=False)

    psi_results = _compute_psi4_points(xyzw, psi_basis, puream=False)
    gg_results = th.compute_points_block(gg.ref.compute_collocation, xyzw, py_basis, spherical=False)

    th.compare_collocation_results(gg_results, psi_results)


@th.using_psi4_libxc
@pytest.mark.parametrize("grad", [0, 1, 2])
def test_psi_derivs_spherical(grad):
    psi_basis, py_basis = _build_psi4_basis(HeC_mol, "cc-pV6Z", puream=True)

    psi_results = _compute_psi4_points(xyzw, psi_basis, puream=True)
    gg_results = th.compute_points_block(gg.ref.compute_collocation, xyzw, py_basis, spherical=True)

    th.compare_collocation_results(gg_results, psi_results)