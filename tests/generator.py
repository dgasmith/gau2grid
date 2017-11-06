"""
NumPy based example project.
"""

import numpy as np

names = ["X", "Y", "Z"]
# Generator
print("Order generator")
spacer = "        "
L = 2
idx = 0
for i in range(L + 1):
    l = L - i + 1
    for j in range(i + 1):
        m = i - j + 1
        n = j + 1


        lp = l - 1;
        mp = m - 1;
        np = n - 1;
        name = "X" * lp + "Y" * mp + "Z" * np

        # Density 
        print("# Basis %s" % name)
        print("# Density")

        print(spacer + "xyz = np.einsum('p,p,p->p', xc_pow[%d], yc_pow[%d], zc_pow[%d])" % (l, m, n))

        tmp = spacer
        tmp += "np.einsum('p,p->p', S0, xyz, "
        tmp += "out=phi_out[%d])" % idx
        print(tmp)

        print("# Gradient components")
        if lp > 0:
            tmp = spacer
            tmp += "np.einsum(',p,p,p,p->p', %d, S0, " % lp
            tmp += "xc_pow[%d], " % lp
            tmp += "yc_pow[%d], " % m
            tmp += "zc_pow[%d], " % n
            tmp += "out=phi_out_x[%d])" % idx
            print(tmp)
        print(spacer + "phi_out_x[%d] += SX * xyz" % idx)

        if mp > 0:
            tmp = spacer
            tmp += "np.einsum(',p,p,p,p->p', %d, S0, " % mp
            tmp += "xc_pow[%d], " % mp
            tmp += "yc_pow[%d], " % m
            tmp += "zc_pow[%d], " % n
            tmp += "out=phi_out_y[%d])" % idx
            print(tmp)
        print(spacer + "phi_out_y[%d] += SY * xyz" % idx)
    
        if np > 0:
            tmp = spacer
            tmp += "np.einsum(',p,p,p,p->p', %d, S0, " % np
            tmp += "xc_pow[%d], " % np
            tmp += "yc_pow[%d], " % m
            tmp += "zc_pow[%d], " % n
            tmp += "out=phi_out_z[%d])" % idx
            print(tmp)
        print(spacer + "phi_out_z[%d] += SZ * xyz" % idx)


        print("# Hessian components")
        # print(spacer + "AX = lp * xc_pow[l - 1] * yc_pow[m] * zc_pow[n]")
        # print(spacer + "AY = mp * xc_pow[l] * yc_pow[m - 1] * zc_pow[n]")
        # print(spacer + "AZ = np * xc_pow[l] * yc_pow[m] * zc_pow[n - 1]")

        if ((lp - 1) > 0):
            print(spacer + "AXX = lp * (lp - 1) * xc_pow[l - 2] * yc_pow[m] * zc_pow[n]")
            print(spacer + "phi_out_xx[%d] = SXX * xyz + SX * AX + SX * AX + S0 * AXX" % idx)

        if (lp > 0) and (mp > 0):
            print(spacer + "AXY = lp * mp * xc_pow[l - 1] * yc_pow[m - 1] * zc_pow[n]")
            print(spacer + "phi_out_xy[%d] = SXY * xyz + SX * AY + SY * AX + S0 * AXY" % idx)

        if (lp > 0) and (mp > 0):
            print(spacer + "AXZ = lp * np * xc_pow[l - 1] * yc_pow[m] * zc_pow[n - 1]")
            print(spacer + "phi_out_xz[%d] = SXZ * xyz + SX * AZ + SZ * AX + S0 * AXZ" % idx)

        if (np > 0) and (mp > 0):
            print(spacer + "AYZ = mp * np * xc_pow[l] * yc_pow[m - 1] * zc_pow[n - 1]")
            print(spacer + "phi_out_yz[%d] = SYZ * xyz + SY * AZ + SZ * AY + S0 * AYZ" % idx)

        if (mp > 0):
            print(spacer + "AYY = mp * (mp - 1) * xc_pow[l] * yc_pow[m - 2] * zc_pow[n]")
            print(spacer + "phi_out_yy[%d] = SYY * xyz + SY * AY + SY * AY + S0 * AYY" % idx)

        if (np > 0):
            print(spacer + "AZZ = np * (np - 1) * xc_pow[l] * yc_pow[m] * zc_pow[n - 2]")
            print(spacer + "phi_out_zz[%d] = SZZ * xyz + SZ * AZ + SZ * AZ + S0 * AZZ" % idx)

        idx += 1
        print(" ")

