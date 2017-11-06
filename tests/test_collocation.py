"""
NumPy based example project.
"""

import psi4
import numpy as np
np.set_printoptions(linewidth=120)

mol = psi4.geometry("""
#He 0 0 0
He 0 0 2
no_com
no_reorient
""")

# Tweakers
npoints = 3
basis = "sto-3g"
#basis = "cc-pVQZ"


# Can only handle cartesian data
psi4.set_options({"PUREAM": False})

mol.update_geometry()
geom = np.array(mol.geometry())

cart_x = np.random.rand(npoints)
cart_y = np.random.rand(npoints)
cart_z = np.random.rand(npoints)
weights = np.random.rand(npoints)
print(cart_x, cart_y, cart_z)

basis = psi4.core.BasisSet.build(mol, "orbital", basis, puream=False)
py_basis = []
for x in range(basis.nshell()):
    shell = basis.shell(x)
    tmp = {}
    tmp["center"] = geom[shell.ncenter]

    tmp["exp"] = [shell.exp(n) for n in range(shell.nprimitive)]
    tmp["coef"] = [shell.coef(n) for n in range(shell.nprimitive)]
    tmp["am"] = shell.am
    py_basis.append(tmp)
#    print(tmp)



extents = psi4.core.BasisExtents(basis, 1.e-50)
block = psi4.core.BlockOPoints(psi4.core.Vector.from_array(cart_x),
                               psi4.core.Vector.from_array(cart_y),
                               psi4.core.Vector.from_array(cart_z),
                               psi4.core.Vector.from_array(weights),
                               extents)
p4_points = psi4.core.BasisFunctions(basis, npoints, basis.nbf())
p4_points.set_deriv(1)
p4_points.compute_functions(block)
points = p4_points.basis_values()

def compute_shell(x, y, z, shell, grad=0):

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
    
    V1 = np.zeros((npoints))
    V2 = np.zeros((npoints))
    for K in range(nprim):
        T1 = norm[K] * np.exp(- alpha[K] * R2)
        T2 = -2.0 * alpha[K] * T1
        V1 += T1
        V2 += T2

    S0 = V1.copy()
    SX = V2 * xc
    SY = V2 * yc
    SZ = V2 * zc
    
    # SX, SY, SZ, SXX, SXZ, SXZ, SYY, SYZ, SZZ
    
    xc_pow = np.zeros((L + 2, npoints))
    yc_pow = np.zeros((L + 2, npoints))
    zc_pow = np.zeros((L + 2, npoints))
    
    xc_pow[0] = 0.0
    yc_pow[0] = 0.0
    zc_pow[0] = 0.0
    xc_pow[1] = 1.0
    yc_pow[1] = 1.0
    zc_pow[1] = 1.0
    
    for LL in range(2, L + 2):
        xc_pow[LL] = xc_pow[LL - 1] * xc
        yc_pow[LL] = yc_pow[LL - 1] * yc
        zc_pow[LL] = zc_pow[LL - 1] * zc
   
    ncart = int((L + 1) * (L + 2) / 2) 
    output = {}
    output["PHI"] = np.zeros((ncart, npoints))
    output["PHI_X"] = np.zeros((ncart, npoints))
    output["PHI_Y"] = np.zeros((ncart, npoints))
    output["PHI_Z"] = np.zeros((ncart, npoints))

    idx = 0
    for i in range(L + 1):
        l = L - i + 1
        for j in range(i + 1):
            m = i - j + 1
            n = j + 1

            a_lp = l - 1
            a_mp = m - 1
            a_np = n - 1
            
            output["PHI"][idx] = S0 * xc_pow[l] * yc_pow[m] * zc_pow[n]

            xyz = xc_pow[l] * yc_pow[m] * zc_pow[n];
            output["PHI"][idx] = S0 * xyz;
            output["PHI_X"][idx] = S0 * a_lp * xc_pow[l-1] * yc_pow[m] * zc_pow[n] + SX * xyz;
            output["PHI_Y"][idx] = S0 * a_mp * xc_pow[l] * yc_pow[m-1] * zc_pow[n] + SY * xyz;
            output["PHI_Z"][idx] = S0 * a_np * xc_pow[l] * yc_pow[m] * zc_pow[n-1] + SZ * xyz;
            idx += 1
        

    return output

g2g_results = {k : [] for k in list(points)}
for shell in py_basis:
    for k, v in compute_shell(cart_x, cart_y, cart_z, shell).items():
        g2g_results[k].append(v)

PHI = np.vstack(g2g_results["PHI"])
PHI_X = np.vstack(g2g_results["PHI_X"])
PHI_Y = np.vstack(g2g_results["PHI_Y"])
PHI_Z = np.vstack(g2g_results["PHI_Z"])

psi_PHI = points["PHI"].np.T
psi_PHI_X = points["PHI_X"].np.T
psi_PHI_Y = points["PHI_Y"].np.T
psi_PHI_Z = points["PHI_Z"].np.T
#print(PHI)
#print(psi_PHI)

#print(np.where(np.any(np.abs(PHI - psi_PHI) > 1.e-14, axis=1)))
wrong = np.any(np.abs(PHI - psi_PHI) > 1.e-14, axis=1)
#print(PHI[wrong])
#print(psi_PHI[wrong])

print("Psi4 and Numpy PHI   same: %s" % np.allclose(PHI, psi_PHI))
print("Psi4 and Numpy PHI_X same: %s" % np.allclose(PHI, psi_PHI))
print("Psi4 and Numpy PHI_Y same: %s" % np.allclose(PHI, psi_PHI))
print("Psi4 and Numpy PHI_Z same: %s" % np.allclose(PHI, psi_PHI))

# ... 
