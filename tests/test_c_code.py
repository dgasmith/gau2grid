"""
Compare the generated C code against the NumPy reference code.
"""

import numpy as np
import gau2grid as gg
import pytest
import sys
import time
np.set_printoptions(linewidth=120, suppress=True)

# Import locals
import ref_basis
import test_helper as th

c_tests = [("cc-pVDZ", "cartesian")]

# pytest.mark.skipif(th.find_pygau2grid(), reason="Could not find a compiled pygau2grid")
# @pytest.mark.parametrize("basis_name,spherical", c_tests)
# def test_generator_collocation(basis_name, spherical):

#     sys.path.insert(1, th.find_pygau2grid())
#     import pygau2grid as pgg

#     gg.np_gen.collocation
#     trans = "spherical" == spherical
#     basis = ref_basis.test_basis[basis_name]
#     basis = [basis[0]]

#     t = time.time()
#     gen_results = th.compute_points_block(gg.np_gen.collocation, xyzw, basis, spherical=trans)
#     gg_time = time.time() - t

#     t = time.time()
#     ref_results = th.compute_points_block(gg.ref.collocation, xyzw, basis, spherical=trans)
#     ref_time = time.time() - t

#     print("")
#     print("%s-%s time NP: %8.4f CGG: %8.4f" % (basis_name, spherical, ref_time, gg_time))

#     th.compare_collocation_results(gen_results, ref_results)


