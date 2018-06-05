"""
Compare the generated NumPy code against the NumPy reference code.
"""

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
xyzw = np.random.rand(4, npoints)

# LR points
xyzw[:, npoints2:] += 20 * np.random.rand(4, npoints2)


@pytest.mark.parametrize("basis_name", ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "cc-pV5Z", "cc-pV6Z"])
@pytest.mark.parametrize("spherical", ["cartesian", "spherical"])
def test_generator_collocation(basis_name, spherical):

    trans = "spherical" == spherical
    basis = ref_basis.test_basis[basis_name]

    t = time.time()
    gen_results = gg.np_gen.collocation_basis(xyzw, basis, spherical=trans, grad=2)
    gg_time = time.time() - t

    t = time.time()
    ref_results = gg.ref.collocation_basis(xyzw, basis, spherical=trans, grad=2)
    ref_time = time.time() - t

    # Print time with py.test -s flags
    print("")
    print("%s-%s time REF: %8.4f GG: %8.4f" % (basis_name, spherical, ref_time, gg_time))

    th.compare_collocation_results(gen_results, ref_results)


@pytest.mark.parametrize("grad", [0, 1, 2])
def test_generator_derivs(grad):

    basis = ref_basis.test_basis["cc-pVDZ"]

    gen_results = gg.np_gen.collocation_basis(xyzw, basis, spherical=False, grad=grad)
    ref_results = gg.ref.collocation_basis(xyzw, basis, spherical=False, grad=grad)

    th.compare_collocation_results(gen_results, ref_results)


@pytest.mark.parametrize("grad", [0, 1, 2])
def test_generator_derivs_spherical(grad):

    basis = ref_basis.test_basis["cc-pVDZ"]

    gen_results = gg.np_gen.collocation_basis(xyzw, basis, spherical=True, grad=grad)
    ref_results = gg.ref.collocation_basis(xyzw, basis, spherical=True, grad=grad)

    th.compare_collocation_results(gen_results, ref_results)

@pytest.mark.parametrize("L", [0, 1, 2])
@pytest.mark.parametrize("s_order", ["cca", "gaussian"])
def test_generator_spherical_order(L, s_order):

    coeffs = [4, 1, 0.25]
    exps = [3, 2, 1]
    center = [1, 1, 1]

    gen_results = gg.np_gen.collocation(xyzw, L, coeffs, exps, center, spherical=True, spherical_order=s_order)
    ref_results = gg.ref.collocation(xyzw, L, coeffs, exps, center, spherical=True, spherical_order=s_order)

    th.compare_collocation_results(gen_results, ref_results)