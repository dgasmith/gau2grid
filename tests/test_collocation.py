"""
NumPy based example project.
"""

import numpy as np
import psi4
np.set_printoptions(linewidth=120)

mol = psi4.geometry("""
#He 0 0 0
He 0 0 2
no_com
no_reorient
""")

# Tweakers
npoints = 2
basis = "cc-pVDZ"

# Can only handle cartesian data
psi4.set_options({"PUREAM": False})

mol.update_geometry()
geom = np.array(mol.geometry())

cart_x = np.random.rand(npoints)
cart_y = np.random.rand(npoints)
cart_z = np.random.rand(npoints)
weights = np.random.rand(npoints)

basis = psi4.core.BasisSet.build(mol, "orbital", basis)
py_basis = []
for x in range(basis.nshell()):
    shell = basis.shell(x)
    tmp = {}
    tmp["center"] = geom[shell.ncenter]

    tmp["exp"] = [shell.exp(n) for n in range(shell.nprimitive)]
    tmp["coef"] = [shell.coef(n) for n in range(shell.nprimitive)]
    tmp["am"] = shell.am
    py_basis.append(tmp)


extents = psi4.core.BasisExtents(basis, 1.e-50)
block = psi4.core.BlockOPoints(psi4.core.Vector.from_array(cart_x),
                               psi4.core.Vector.from_array(cart_y),
                               psi4.core.Vector.from_array(cart_z),
                               psi4.core.Vector.from_array(weights),
                               extents)
p4_points = psi4.core.BasisFunctions(basis, npoints, basis.nbf())
p4_points.compute_functions(block)
points = p4_points.basis_values()


def compute_shell(x, y, z, shell):

    alpha = shell["exp"]
    L = shell["am"]
    norm = shell["coef"]
    center = shell["center"]
    nprim = len(norm)
    npoints = x.size
    

    xc = x - center[0]
    yc = y - center[1]
    zc = z - center[2]
    R2 = xc * xc + yc * yc + zc * zc
    
    S0 = np.zeros(npoints)
    for K in range(nprim):
        S0 += norm[K] * np.exp(- alpha[K] * R2)
    
    # SX, SY, SZ, SXX, SXZ, SXZ, SYY, SYZ, SZZ
    
    xc_pow = np.zeros((L + 1, npoints))
    yc_pow = np.zeros((L + 1, npoints))
    zc_pow = np.zeros((L + 1, npoints))
    
    xc_pow[0] = 1.0
    yc_pow[0] = 1.0
    zc_pow[0] = 1.0
    
    for LL in range(1, L + 1):
        xc_pow[LL] = xc_pow[LL - 1] * xc
        yc_pow[LL] = xc_pow[LL - 1] * yc
        zc_pow[LL] = xc_pow[LL - 1] * zc
    
    output = np.zeros((L * 2 + 1, npoints))
    
    # S 
    if L == 0:
        output[0] = S0 * xc_pow[0] * yc_pow[0] * zc_pow[0]
    
    # P
    elif L == 1:
        output[0] = S0 * xc_pow[1] * yc_pow[0] * zc_pow[0]
        output[1] = S0 * xc_pow[0] * yc_pow[1] * zc_pow[0]
        output[2] = S0 * xc_pow[0] * yc_pow[0] * zc_pow[1]
    
    # D
    elif L == 2: 
        output[0] = S0 * xc_pow[2] * xc_pow[0] * xc_pow[0]
        output[1] = S0 * xc_pow[1] * xc_pow[1] * xc_pow[0]
        output[2] = S0 * xc_pow[1] * xc_pow[0] * xc_pow[1]
        output[3] = S0 * xc_pow[0] * xc_pow[2] * xc_pow[0]
        output[4] = S0 * xc_pow[0] * xc_pow[1] * xc_pow[1]
        output[5] = S0 * xc_pow[0] * xc_pow[0] * xc_pow[2]

    else:
        raise ValueError("AM of %d NYI." % L)

    return output

phi_tmp = []
for shell in py_basis:
    tmp = compute_shell(cart_x, cart_y, cart_z, shell)
    phi_tmp.append(tmp)

phi = np.vstack(phi_tmp)
psi_phi = points["PHI"].np.T
print(phi)
print(psi_phi)

print(np.allclose(phi, psi_phi))



# ... 
