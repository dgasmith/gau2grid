"""
Builds a new RSH dictionary for gau2grid.
"""
import gau2grid
import numpy as np

np.set_printoptions(precision=60)
import mpmath
mpmath.mp.pretty = False

do_print = True

rsh_dict = {}
for AM in range(17):
    print("AM %d" % AM)
    data = gau2grid.RSH.cart_to_RSH_coeffs(AM, gen=True)

    # Loop over spherical components
    for nsph, spherical in enumerate(data):
        print("nsph %d" % nsph)
   

        title = "AM_%d_%d_" % (AM, nsph) 
        # Loop over cart components per spherical
        xyz_order = []
        coefs = []
        for cart in spherical:
            xyz, coef = cart
            xyz_order.append(xyz)
            #coefs.append(coef)
            coefs.append(str(coef))

        xyz_order = np.array(xyz_order, dtype=np.int)
        coefs = np.array(coefs, dtype=np.float128)
        #coefs = np.array(coefs, dtype=np.longdouble)

        if do_print:
            print(abs(spherical[-1][1]))
            print(str(np.array([str(spherical[-1][1])], dtype=np.float128))[2:-1])
            print(str(np.array([str(spherical[-1][1])], dtype=np.float64))[2:-1])
            print(xyz_order)
            print(coefs)

        rsh_dict[title + "xyz"] = xyz_order
        rsh_dict[title + "coeffs"] = coefs

    if do_print:
        print(" ")

np.savez("rsh_coeffs.npz", **rsh_dict)
#print(rsh_dict)

print("""
0.7015607600201140097969528873093529940792537139585113024319522751225065918075510677456518553489461697
0.70156076002011402703573139660875312983989715576171875

8.464890209258860703766777652251376565148731174244082586689222950404896579281534196672844804179245696
8.46489020925886070413624029384891400695778429508209228515625]
""")
print(mpmath.mpf("0.7015607600201140097969528873093529940792537139585113024319522751225065918075510677456518553489461697"))
