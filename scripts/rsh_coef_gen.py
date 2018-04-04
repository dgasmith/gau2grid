"""
Builds a new RSH dictionary for gau2grid.
"""
import gau2grid
import numpy as np

np.set_printoptions(precision=30)

do_print = True

rsh_dict = {}
for AM in range(17):
    data = gau2grid.RSH.cart_to_RSH_coeffs(AM)

    # Loop over spherical components
    for nsph, spherical in enumerate(data):
   

        title = "AM_%d_%d_" % (AM, nsph) 
        # Loop over cat components per spherical
        xyz_order = []
        coefs = []
        for cart in spherical:
            xyz, coef = cart
            xyz_order.append(xyz)
            coefs.append(coef)

        xyz_order = np.array(xyz_order, dtype=np.int)
        coefs = np.array(coefs, dtype=np.float128)

        if do_print:
            print(xyz_order)
            print(coefs)

        rsh_dict[title + "xyz"] = xyz_order
        rsh_dict[title + "coeffs"] = coefs

    if do_print:
        print(" ")

np.savez("rsh_coeffs.npz", **rsh_dict)
#print(rsh_dict)

