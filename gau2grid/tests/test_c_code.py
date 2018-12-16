"""
Compare the generated C code against the NumPy reference code.
"""

import os
import sys
import time

import numpy as np
import pytest

import gau2grid as gg

np.set_printoptions(linewidth=120, suppress=True)

# Import locals
from . import ref_basis
from . import test_helper as th

# Tweakers
npoints = int(1.e3)
npoints2 = int(npoints / 2)

# Global points
np.random.seed(0)
xyzw = np.random.rand(3, npoints)

# LR points
xyzw[:, npoints2:] += 5 * np.random.rand(3, npoints2)

# Make sure the C-side has been compiled
if "GAU2GRID_FORCE_C_TEST" in os.environ:
    skip_c_test = False
else:
    skip_c_test = gg.c_compiled() is False
check_compile = pytest.mark.skipif(skip_c_test, reason="Could not find the C compiled SO for gau2grid")


@check_compile
@pytest.mark.parametrize("basis_name", ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "cc-pV5Z", "cc-pV6Z"])
@pytest.mark.parametrize("spherical", ["cartesian", "spherical"])
def test_generator_collocation(basis_name, spherical):

    trans = "spherical" == spherical
    kwargs = {
        "spherical_order": gg.spherical_order(),
        "cartesian_order": gg.cartesian_order(),
        "spherical": trans,
        "grad": 2
    }

    basis = ref_basis.test_basis[basis_name]

    t = time.time()
    gen_results = gg.collocation_basis(xyzw, basis, **kwargs)
    gg_time = time.time() - t

    t = time.time()
    ref_results = gg.ref.collocation_basis(xyzw, basis, **kwargs)
    ref_time = time.time() - t

    # Print time with py.test -s flags
    print("")
    print("%s-%s time REF: %8.4f GG: %8.4f" % (basis_name, spherical, ref_time, gg_time))

    th.compare_collocation_results(gen_results, ref_results)


@check_compile
@pytest.mark.parametrize("basis_name", ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "cc-pV5Z", "cc-pV6Z"])
@pytest.mark.parametrize("spherical", ["cartesian", "spherical"])
def test_generator_orbital(basis_name, spherical):

    trans = "spherical" == spherical
    kwargs = {
        "spherical_order": gg.spherical_order(),
        "cartesian_order": gg.cartesian_order(),
        "spherical": trans,
        "grad": 0
    }

    basis = ref_basis.test_basis[basis_name]

    t = time.time()
    phi = gg.collocation_basis(xyzw, basis, **kwargs)["PHI"]
    orbs = np.random.rand(5, phi.shape[0])
    benchmark = np.dot(orbs, phi)
    ref_time = time.time() - t

    t = time.time()
    del kwargs["grad"]
    ref_results = gg.orbital_basis(orbs, xyzw, basis, **kwargs)
    gg_time = time.time() - t

    # Print time with py.test -s flags
    print("")
    print("%s-%s time REF: %8.4f GG: %8.4f" % (basis_name, spherical, ref_time, gg_time))
    # print(benchmark)
    # print(ref_results)

    th.compare_collocation_results({"ORBITALS": benchmark}, {"ORBITALS": ref_results})


@check_compile
@pytest.mark.parametrize("grad", [0, 1, 2])
def test_generator_derivs(grad):
    kwargs = {
        "spherical_order": gg.spherical_order(),
        "cartesian_order": gg.cartesian_order(),
        "spherical": False,
        "grad": grad
    }

    basis = ref_basis.test_basis["cc-pVDZ"]

    gen_results = gg.collocation_basis(xyzw, basis, **kwargs)
    ref_results = gg.ref.collocation_basis(xyzw, basis, **kwargs)

    th.compare_collocation_results(gen_results, ref_results)


@check_compile
@pytest.mark.parametrize("grad", [0, 1, 2])
def test_generator_derivs_spherical(grad):
    kwargs = {
        "spherical_order": gg.spherical_order(),
        "cartesian_order": gg.cartesian_order(),
        "spherical": True,
        "grad": grad
    }

    basis = ref_basis.test_basis["cc-pVDZ"]

    gen_results = gg.collocation_basis(xyzw, basis, **kwargs)
    ref_results = gg.ref.collocation_basis(xyzw, basis, **kwargs)

    th.compare_collocation_results(gen_results, ref_results)


@check_compile
def test_libgg_path():
    assert "gg" in gg.cgg_path()


@check_compile
def test_spherical_order():
    assert gg.spherical_order() in ["gaussian", "cca"]


@check_compile
def test_cartesian_order():
    assert gg.cartesian_order() in ["row"]

@check_compile
@pytest.mark.parametrize("am,spherical,result", [
        (0, True, 1),
        (0, False, 1),
        (1, True, 3),
        (1, False, 3),
        (2, True, 5),
        (2, False, 6),
        (3, True, 7),
        (3, False, 10),
    ])
def test_ncomponents(am, spherical, result):
    assert gg.ncomponents(am, spherical) == result


@check_compile
@pytest.mark.parametrize("spherical", [True, False])
@pytest.mark.parametrize("am", [0, 1, 2, 3, 4])
def test_generator_orbitals_am(spherical, am):
    kwargs = {
        "spherical_order": gg.spherical_order(),
        "cartesian_order": gg.cartesian_order(),
        "spherical": spherical,
        "grad": 0
    }

    # Build a single orbital
    coeffs = [0.44135347600549724, 0.6934968471367846, 0.6641842253258472, 0.0001]
    exponents = [38.36, 5.77, 1.24, 1.e-2]
    center = [0., 0., 0.]
    L = am

    ret = gg.collocation(xyzw, L, coeffs, exponents, center, **kwargs)["PHI"]
    orb = np.random.rand(3, ret.shape[0])
    bench = np.dot(orb, ret)

    del kwargs["grad"]
    ret = gg.orbital(orb, xyzw, L, coeffs, exponents, center, **kwargs)

    # Compare the results
    th.compare_collocation_results({"ORBITALS": bench}, {"ORBITALS": ret})