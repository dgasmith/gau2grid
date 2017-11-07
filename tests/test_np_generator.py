"""
Compare the generated NumPy code against the NumPy reference code.
"""

import numpy as np
import gau2grid as gg
import pytest
np.set_printoptions(linewidth=120, suppress=True)

# Import locals
import ref_basis

# Tweakers
npoints = 100

# Global points
np.random.seed(0)
xyzw = np.random.rand(npoints, 4)


def _compute_points_block(func, xyzw, basis, grad=2, puream=False):
    """
    Computes the reference collocation matrices and stitches them together
    """

    # Sum up g2g points
    tmp = []
    for shell in basis:
        shell_collocation = func(xyzw, shell["am"], shell["coef"], shell["exp"], shell["center"], grad=2)
        tmp.append(shell_collocation)

    g2g_results = {k: [] for k in tmp[0].keys()}
    for coll in tmp:
        for k, v in coll.items():
            g2g_results[k].append(v)

    g2g_results = {k: np.vstack(v) for k, v in g2g_results.items()}
    return g2g_results


@pytest.mark.parametrize("basis", ["cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "cc-pV5Z", "cc-pV6Z"])
def test_generator_collocation(basis):
    basis = ref_basis.test_basis[basis]

    max_am = max(shell["am"] for shell in basis)
    code = gg.generator.numpy_generator(max_am, function_name="tmp_np_gen")

    # Exec the code into a namespace
    test_namespace = {}
    exec(code, test_namespace)

    gen_results = _compute_points_block(test_namespace["tmp_np_gen"], xyzw, basis)
    ref_results = _compute_points_block(gg.ref.compute_collocation, xyzw, basis)

    if set(ref_results) != set(ref_results):
        raise KeyError("Psi4 and GG results dicts do not match")

    for k in ref_results.keys():
        match = np.allclose(gen_results[k], ref_results[k])
        diff = np.linalg.norm(gen_results[k] - ref_results[k])
        if not match:
            raise ValueError("NumPy generator results do not match reference for %s" % k)
