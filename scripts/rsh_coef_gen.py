"""
Builds a new RSH dictionary for gau2grid.
"""
import gau2grid
import numpy as np
import decimal
import pickle
np.set_printoptions(precision=60)

rsh_dict = {}
for AM in range(17):
    print("AM %d" % AM)
    data = gau2grid.RSH.cart_to_RSH_coeffs(AM, gen=True)
    rsh_dict[AM] = data

with open("rsh_coeffs.pkl", "wb") as handle:
    pickle.dump(rsh_dict, handle, protocol=2)


with open("rsh_coeffs.pkl", "rb") as handle:
    data = pickle.load(handle)

print(data[3])
