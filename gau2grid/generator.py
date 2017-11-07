"""
This is a Python-based automatic generator.
"""

import numpy as np

def numpy_generator(L, function_name="generated_compute_numpy_shells"):

    # Builds a few tmps
    s1 = "    "
    s2 = "    " * 2
    s3 = "    " * 3


    ret = []
    ret.append("def %s(xyz, L, coeffs, exponents, center, grad=2):" % function_name)

    ret.append("")
    ret.append(s1 + "# Make sure NumPy is in locals")
    ret.append(s1 + "import numpy as np")
    ret.append("")

    ret.append(s1 + "# Unpack shell data")
    ret.append(s1 + "nprim = len(coeffs)")
    ret.append(s1 + "npoints = xyz.shape[0]")
    ret.append("")

    ret.append(s1 + "# First compute the diff distance in each cartesian")
    ret.append(s1 + "xc = xyz[:, 0] - center[0]")
    ret.append(s1 + "yc = xyz[:, 1] - center[1]")
    ret.append(s1 + "zc = xyz[:, 2] - center[2]")
    ret.append(s1 + "R2 = xc * xc + yc * yc + zc * zc")
    ret.append("")


    ret.append(s1 + "# Build up the derivates in each direction")
    ret.append(s1 + "V1 = np.zeros((npoints))")
    ret.append(s1 + "V2 = np.zeros((npoints))")
    ret.append(s1 + "V3 = np.zeros((npoints))")
    ret.append(s1 + "for K in range(nprim):")
    ret.append(s1 + "    T1 = coeffs[K] * np.exp(-exponents[K] * R2)")
    ret.append(s1 + "    T2 = -2.0 * exponents[K] * T1")
    ret.append(s1 + "    T3 = -2.0 * exponents[K] * T2")
    ret.append(s1 + "    V1 += T1")
    ret.append(s1 + "    V2 += T2")
    ret.append(s1 + "    V3 += T3")
    ret.append("")
    ret.append(s1 + "S0 = V1.copy()")
    ret.append(s1 + "SX = V2 * xc")
    ret.append(s1 + "SY = V2 * yc")
    ret.append(s1 + "SZ = V2 * zc")
    ret.append(s1 + "SXY = V3 * xc * yc")
    ret.append(s1 + "SXZ = V3 * xc * zc")
    ret.append(s1 + "SYZ = V3 * yc * zc")
    ret.append(s1 + "SXX = V3 * xc * xc + V2")
    ret.append(s1 + "SYY = V3 * yc * yc + V2")
    ret.append(s1 + "SZZ = V3 * zc * zc + V2")
    ret.append("")

    ret.append(s1 + "# Power matrix for higher angular momenta")
    ret.append(s1 + "xc_pow = np.zeros((L + 3, npoints))")
    ret.append(s1 + "yc_pow = np.zeros((L + 3, npoints))")
    ret.append(s1 + "zc_pow = np.zeros((L + 3, npoints))")
    ret.append("")
    ret.append(s1 + "xc_pow[0] = 0.0")
    ret.append(s1 + "yc_pow[0] = 0.0")
    ret.append(s1 + "zc_pow[0] = 0.0")
    ret.append(s1 + "xc_pow[1] = 0.0")
    ret.append(s1 + "yc_pow[1] = 0.0")
    ret.append(s1 + "zc_pow[1] = 0.0")
    ret.append(s1 + "xc_pow[2] = 1.0")
    ret.append(s1 + "yc_pow[2] = 1.0")
    ret.append(s1 + "zc_pow[2] = 1.0")

    ret.append(s1 + "for LL in range(3, L + 3):")
    ret.append(s1 + "    xc_pow[LL] = xc_pow[LL - 1] * xc")
    ret.append(s1 + "    yc_pow[LL] = yc_pow[LL - 1] * yc")
    ret.append(s1 + "    zc_pow[LL] = zc_pow[LL - 1] * zc")
    ret.append("")

    ret.append(s1 + "# Allocate data")
    ret.append(s1 + "ncart = int((L + 1) * (L + 2) / 2)")
    ret.append("")

    ret.append(s1 + "output = {}")
    ret.append(s1 + "output['PHI'] = np.zeros((ncart, npoints))")
    ret.append(s1 + "if grad > 0:")
    ret.append(s1 + "    output['PHI_X'] = np.zeros((ncart, npoints))")
    ret.append(s1 + "    output['PHI_Y'] = np.zeros((ncart, npoints))")
    ret.append(s1 + "    output['PHI_Z'] = np.zeros((ncart, npoints))")
    ret.append(s1 + "if grad > 1:")
    ret.append(s1 + "    output['PHI_XX'] = np.zeros((ncart, npoints))")
    ret.append(s1 + "    output['PHI_YY'] = np.zeros((ncart, npoints))")
    ret.append(s1 + "    output['PHI_ZZ'] = np.zeros((ncart, npoints))")
    ret.append(s1 + "    output['PHI_XY'] = np.zeros((ncart, npoints))")
    ret.append(s1 + "    output['PHI_XZ'] = np.zeros((ncart, npoints))")
    ret.append(s1 + "    output['PHI_YZ'] = np.zeros((ncart, npoints))")
    ret.append(s1 + "if grad > 2:")
    ret.append(s1 + "    raise ValueError('Only grid derivatives through Hessians (grad = 2) has been implemented')")
    ret.append("")

    ret.append("# Angular momentum loops")
    for l in range(L + 1):
        ret.append(s1 + "if L == %d:" % l)
        ret.extend(numpy_am_build(l, s2))

    ret.append(s1 + "return output")

    return "\n".join(ret)

def numpy_am_build(L, spacer=""):
    ret = []
    names = ["X", "Y", "Z"]
    # Generator
    idx = 0
    for i in range(L + 1):
        l = L - i + 2
        for j in range(i + 1):
            m = i - j + 2
            n = j + 2

            ld1 = l - 1
            ld2 = l - 2
            md1 = m - 1
            md2 = m - 2
            nd1 = n - 1
            nd2 = n - 2
            name = "X" * ld2 + "Y" * md2 + "Z" * nd2
            if name == "":
                name = "0"

            # Density
            ret.append("# Density AM=%d Component=%s" % (L, name))

            ret.append(spacer + "A = xc_pow[%d] * yc_pow[%d] * zc_pow[%d]" % (l, m, n))
            ret.append(spacer + "output['PHI'][%d] = S0 * A" % idx)

            ret.append("# Gradient AM=%d Component=%s" % (L, name))
            ret.append(spacer + "AX = %d * xc_pow[%d] * yc_pow[%d] * zc_pow[%d]" % (ld2, ld1, m, n))
            ret.append(spacer + "AY = %d * xc_pow[%d] * yc_pow[%d] * zc_pow[%d]" % (md2, l, md1, n))
            ret.append(spacer + "AZ = %d * xc_pow[%d] * yc_pow[%d] * zc_pow[%d]" % (nd2, l, m, nd1))
            ret.append(spacer + "output['PHI_X'][%d] = S0 * AX + SX * A" % idx)
            ret.append(spacer + "output['PHI_Y'][%d] = S0 * AY + SY * A" % idx)
            ret.append(spacer + "output['PHI_Z'][%d] = S0 * AZ + SZ * A" % idx)

            ret.append("# Hessian AM=%d Component=%s" % (L, name))
            ret.append(spacer + "AXY = %d * %d * xc_pow[%d] * yc_pow[%d] * zc_pow[%d]" % (ld2, md2, ld1, md1, n))
            ret.append(spacer + "AXZ = %d * %d * xc_pow[%d] * yc_pow[%d] * zc_pow[%d]" % (ld2, nd2, ld1, m, nd1))
            ret.append(spacer + "AYZ = %d * %d * xc_pow[%d] * yc_pow[%d] * zc_pow[%d]" % (md2, nd2, l, md1, nd1))
            ret.append(spacer + "AXX = %d * %d * xc_pow[%d] * yc_pow[%d] * zc_pow[%d]" % (ld2, ld2 - 1, ld2, m, n))
            ret.append(spacer + "AYY = %d * %d * xc_pow[%d] * yc_pow[%d] * zc_pow[%d]" % (md2, md2 - 1, l, md2, n))
            ret.append(spacer + "AZZ = %d * %d * xc_pow[%d] * yc_pow[%d] * zc_pow[%d]" % (nd2, nd2 - 1, l, m, nd2))
            ret.append(spacer + "output['PHI_XX'][%d] = SXX * A + SX * AX + SX * AX + S0 * AXX" % idx)
            ret.append(spacer + "output['PHI_YY'][%d] = SYY * A + SY * AY + SY * AY + S0 * AYY" % idx)
            ret.append(spacer + "output['PHI_ZZ'][%d] = SZZ * A + SZ * AZ + SZ * AZ + S0 * AZZ" % idx)
            ret.append(spacer + "output['PHI_XY'][%d] = SXY * A + SX * AY + SY * AX + S0 * AXY" % idx)
            ret.append(spacer + "output['PHI_XZ'][%d] = SXZ * A + SX * AZ + SZ * AX + S0 * AXZ" % idx)
            ret.append(spacer + "output['PHI_YZ'][%d] = SYZ * A + SY * AZ + SZ * AY + S0 * AYZ" % idx)

            idx += 1
            ret.append(" ")

    return ret