"""
Compare the generated NumPy code against the NumPy reference code.
"""

import numpy as np
import gau2grid as gg
import pytest
import time
np.set_printoptions(linewidth=120, suppress=True)

# Import locals
import ref_basis
import test_helper as th

# Tweakers
npoints = int(1.e2)

# Global points
np.random.seed(0)
xyzw = np.random.rand(4, npoints)

# Build up a list of tests
gg_tests = []
for basis in ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "cc-pV5Z", "cc-pV6Z"]:
    for spherical in ["cartesian", "spherical"]:
        gg_tests.append((basis, spherical))


@pytest.mark.parametrize("basis_name,spherical", gg_tests)
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
