"""
Compare the generated NumPy code against the NumPy reference code.
"""

import pytest
import platform
import numpy as np
np.set_printoptions(precision=30)

import gau2grid as gg

try:
    import mpmath
    mpmath.mp.pretty = False
    _has_mpmath = True
except ImportError:
    _has_mpmath = False


def _test_shell(bench, comp):
    comp_line = sorted(comp)
    bench_line = sorted(bench)
    assert len(comp_line) == len(bench_line)

    for cart in range(len(comp_line)):
        comp_coeff = comp_line[cart]
        bench_coeff = bench_line[cart]

        # Make sure cartesian alignment
        assert comp_coeff[0] == bench_coeff[0]

        # Check coefficient
        # print(type(comp_coeff[1]), np.array([comp_coeff[1]]), np.float128(bench_coeff[1]), bench_coeff[1])
        # numpy.float128 is missing on Windows
        if platform.system() == "Windows":
            assert pytest.approx(comp_coeff[1], rel=1.e-16) == np.float64(bench_coeff[1])
        else:
            assert pytest.approx(comp_coeff[1], rel=1.e-16) == np.float128(bench_coeff[1])

    return True


@pytest.mark.skipif(platform.system() == "Windows", reason="Windows does not support caching due to missing numpy.float128")
@pytest.mark.skipif(_has_mpmath is False, reason="Did not detect mpmath module")
#@pytest.mark.parametrize("AM", range(6))
@pytest.mark.parametrize("AM", range(17))
def test_RSH(AM):
    # print("AM %d" % AM)

    pkl_data = gg.RSH.cart_to_RSH_coeffs(AM)
    bench_data = gg.RSH.cart_to_RSH_coeffs(AM, gen=True)

    assert len(pkl_data) == len(bench_data)
    for sph in range(len(pkl_data)):

        assert _test_shell(bench_data[sph], pkl_data[sph])


@pytest.mark.skipif(platform.system() == "Windows", reason="Windows does not support caching due to missing numpy.float128")
def test_RSH_max_am():
    with pytest.raises(ValueError):
        pkl_data = gg.RSH.cart_to_RSH_coeffs(17, gen=False)

def test_RSH_order_p():
    gaus = gg.RSH.cart_to_RSH_coeffs(1, order="gaussian")
    cca = gg.RSH.cart_to_RSH_coeffs(1, order="cca")

    assert _test_shell(gaus[0], cca[1])
    assert _test_shell(gaus[1], cca[2])
    assert _test_shell(gaus[2], cca[0])

def test_RSH_order_d():
    gaus = gg.RSH.cart_to_RSH_coeffs(2, order="gaussian")
    cca = gg.RSH.cart_to_RSH_coeffs(2, order="cca")

    assert _test_shell(gaus[0], cca[2])
    assert _test_shell(gaus[1], cca[3])
    assert _test_shell(gaus[2], cca[1])
    assert _test_shell(gaus[3], cca[4])
    assert _test_shell(gaus[4], cca[0])
