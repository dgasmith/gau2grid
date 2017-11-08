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

# Tweakers
npoints = 500

# Global points
np.random.seed(0)
xyzw = np.random.rand(npoints, 4)


def _compute_points_block(func, xyzw, basis, grad=2, spherical=False):
    """
    Computes the reference collocation matrices and stitches them together
    """

    # Sum up g2g points
    tmp = []
    for shell in basis:
        shell_collocation = func(
            xyzw, shell["am"], shell["coef"], shell["exp"], shell["center"], grad=2, spherical=spherical)
        tmp.append(shell_collocation)

    g2g_results = {k: [] for k in tmp[0].keys()}
    for coll in tmp:
        for k, v in coll.items():
            g2g_results[k].append(v)

    g2g_results = {k: np.vstack(v) for k, v in g2g_results.items()}
    return g2g_results


# Build up a list of tests
gg_tests = []
for basis in ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "cc-pV5Z", "cc-pV6Z"]:
    for spherical in ["cart", "spherical"]:
        gg_tests.append((basis, spherical))


@pytest.mark.parametrize("basis_name,spherical", gg_tests)
def test_generator_collocation(basis_name, spherical):

    trans = "spherical" == spherical
    basis = ref_basis.test_basis[basis_name]

    max_am = max(shell["am"] for shell in basis)
    code = gg.generator.numpy_generator(max_am, function_name="tmp_np_gen")

    # Exec the code into a namespace
    test_namespace = {}
    exec(code, test_namespace)

    t = time.time()
    gen_results = _compute_points_block(test_namespace["tmp_np_gen"], xyzw, basis, spherical=trans)
    gg_time = time.time() - t

    t = time.time()
    ref_results = _compute_points_block(gg.ref.compute_collocation, xyzw, basis, spherical=trans)
    ref_time = time.time() - t

    print("")
    print("%s-%s time REF: %8.4f GG: %8.4f" % (basis_name, spherical, ref_time, gg_time))

    if set(ref_results) != set(ref_results):
        raise KeyError("Psi4 and GG results dicts do not match")

    for k in ref_results.keys():
        match = np.allclose(gen_results[k], ref_results[k])
        diff = np.linalg.norm(gen_results[k] - ref_results[k])
        if not match:
            raise ValueError("NumPy generator results do not match reference for %s" % k)
