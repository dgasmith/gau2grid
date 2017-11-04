"""
NumPy based example project.
"""

import numpy as np

names = ["X", "Y", "Z"]
# Generator
print("Order generator")
spacer = "        "
L = 1
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

        print("# Derivatives")
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
        print(" ")

        idx += 1

