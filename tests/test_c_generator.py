"""
Compare the generated C code against the NumPy reference code.
"""

import numpy as np
import gau2grid as gg
import pytest
import time
np.set_printoptions(linewidth=120, suppress=True)

# Import locals
import ref_basis
import test_helper as th

# Simply test that the code runs for now
c_gen_tests = []
for AM in range(4):
    for grad in range(3):
        c_gen_tests.append((AM, grad))


@pytest.mark.parametrize("AM,grad", c_gen_tests)
def test_c_collocation_codgen(AM, grad):
    cg = gg.codegen.CodeGen(cgen=True)
    gg.c_gen.shell_c_generator(cg, AM, grad=grad)


@pytest.mark.parametrize("AM", list(range(3)))
def test_c_spherical_trans_codgen(AM):
    cg = gg.codegen.CodeGen(cgen=True)
    gg.RSH.transformation_c_generator(cg, AM, "row")


def test_library_gen():
    gg.c_gen.generate_c_gau2grid(4, path="/tmp")