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


@pytest.mark.parametrize("AM", list(range(6)))
def test_c_collocation_codgen(AM):
    # Simply test that it runs for now
    cg = gg.codegen.CodeGen(cgen=True)
    gg.c_gen.shell_c_generator(cg, AM)
