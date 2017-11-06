"""
This is a Python-based automatic generator.
"""

import numpy as np


def numpy_generator(L, spacer="        "):
    ret = []
    names = ["X", "Y", "Z"]
    # Generator
    ret.append("# Order generator")
    L = 2
    idx = 0
    for i in range(L + 1):
        l = L - i + 1
        for j in range(i + 1):
            m = i - j + 1
            n = j + 1

            lp = l - 1
            mp = m - 1
            np = n - 1
            name = "X" * lp + "Y" * mp + "Z" * np

            # Density
            ret.append("# Basis %s" % name)
            ret.append("# Density")

            ret.append(spacer + "xyz = np.einsum('p,p,p->p', xc_pow[%d], yc_pow[%d], zc_pow[%d])" % (l, m, n))

            tmp = spacer
            tmp += "np.einsum('p,p->p', S0, xyz, "
            tmp += "out=phi_out[%d])" % idx
            ret.append(tmp)

            ret.append("# Gradient components")
            if lp > 0:
                tmp = spacer
                tmp += "np.einsum(',p,p,p,p->p', %d, S0, " % lp
                tmp += "xc_pow[%d], " % lp
                tmp += "yc_pow[%d], " % m
                tmp += "zc_pow[%d], " % n
                tmp += "out=phi_out_x[%d])" % idx
                ret.append(tmp)
            ret.append(spacer + "phi_out_x[%d] += SX * xyz" % idx)

            if mp > 0:
                tmp = spacer
                tmp += "np.einsum(',p,p,p,p->p', %d, S0, " % mp
                tmp += "xc_pow[%d], " % mp
                tmp += "yc_pow[%d], " % m
                tmp += "zc_pow[%d], " % n
                tmp += "out=phi_out_y[%d])" % idx
                ret.append(tmp)
            ret.append(spacer + "phi_out_y[%d] += SY * xyz" % idx)

            if np > 0:
                tmp = spacer
                tmp += "np.einsum(',p,p,p,p->p', %d, S0, " % np
                tmp += "xc_pow[%d], " % np
                tmp += "yc_pow[%d], " % m
                tmp += "zc_pow[%d], " % n
                tmp += "out=phi_out_z[%d])" % idx
                ret.append(tmp)
            ret.append(spacer + "phi_out_z[%d] += SZ * xyz" % idx)

            ret.append("# Hessian components")
            # ret.append(spacer + "AX = lp * xc_pow[l - 1] * yc_pow[m] * zc_pow[n]")
            # ret.append(spacer + "AY = mp * xc_pow[l] * yc_pow[m - 1] * zc_pow[n]")
            # ret.append(spacer + "AZ = np * xc_pow[l] * yc_pow[m] * zc_pow[n - 1]")

            if ((lp - 1) > 0):
                ret.append(spacer + "AXX = lp * (lp - 1) * xc_pow[l - 2] * yc_pow[m] * zc_pow[n]")
                ret.append(spacer + "phi_out_xx[%d] = SXX * xyz + SX * AX + SX * AX + S0 * AXX" % idx)

            if (lp > 0) and (mp > 0):
                ret.append(spacer + "AXY = lp * mp * xc_pow[l - 1] * yc_pow[m - 1] * zc_pow[n]")
                ret.append(spacer + "phi_out_xy[%d] = SXY * xyz + SX * AY + SY * AX + S0 * AXY" % idx)

            if (lp > 0) and (mp > 0):
                ret.append(spacer + "AXZ = lp * np * xc_pow[l - 1] * yc_pow[m] * zc_pow[n - 1]")
                ret.append(spacer + "phi_out_xz[%d] = SXZ * xyz + SX * AZ + SZ * AX + S0 * AXZ" % idx)

            if (np > 0) and (mp > 0):
                ret.append(spacer + "AYZ = mp * np * xc_pow[l] * yc_pow[m - 1] * zc_pow[n - 1]")
                ret.append(spacer + "phi_out_yz[%d] = SYZ * xyz + SY * AZ + SZ * AY + S0 * AYZ" % idx)

            if (mp > 0):
                ret.append(spacer + "AYY = mp * (mp - 1) * xc_pow[l] * yc_pow[m - 2] * zc_pow[n]")
                ret.append(spacer + "phi_out_yy[%d] = SYY * xyz + SY * AY + SY * AY + S0 * AYY" % idx)

            if (np > 0):
                ret.append(spacer + "AZZ = np * (np - 1) * xc_pow[l] * yc_pow[m] * zc_pow[n - 2]")
                ret.append(spacer + "phi_out_zz[%d] = SZZ * xyz + SZ * AZ + SZ * AZ + S0 * AZZ" % idx)

            idx += 1
            ret.append(" ")

    return ret