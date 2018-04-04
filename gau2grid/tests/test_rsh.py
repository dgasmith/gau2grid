"""
Compare the generated NumPy code against the NumPy reference code.
"""

import pytest
import numpy as np

import gau2grid as gg

try:
    import mpmath
    _has_mpmath = True
except ImportError:
    _has_mpmath = False


@pytest.mark.skipif(_has_mpmath is False, reason="Did not detect mpmath module")
@pytest.mark.parametrize("AM", range(17))
def test_RSH(AM):

    pkl_data = gg.RSH.cart_to_RSH_coeffs(AM)
    bench_data = gg.RSH.cart_to_RSH_coeffs(AM, gen=True)

    assert len(pkl_data) == len(bench_data)
    for sph in range(len(pkl_data)):
        pkl_line = sorted(pkl_data[sph])
        bench_line = sorted(bench_data[sph])
        assert len(pkl_line) == len(bench_line)

        for cart in range(len(pkl_line)):
            pkl_coeff = pkl_line[cart]
            bench_coeff = bench_line[cart]

            # Make sure cartesian alignment
            assert pkl_coeff[0] == bench_coeff[0]

            # Check coefficient
            assert pytest.approx(pkl_coeff[1], abs=1.e-15) == np.float128(bench_coeff[1])


def test_RSH_max_am():
    with pytest.raises(ValueError):
        pkl_data = gg.RSH.cart_to_RSH_coeffs(17, gen=False)
