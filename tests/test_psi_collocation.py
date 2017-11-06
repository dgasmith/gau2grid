"""
NumPy based example project.
"""

import psi4
import numpy as np
np.set_printoptions(linewidth=120)

mol = psi4.geometry("""
He 0 0 0
He 0 0 2
no_com
no_reorient
""")

# Tweakers
npoints = 50
basis = "sto-3g"
basis = "cc-pVQZ"


# Can only handle cartesian data
psi4.set_options({"PUREAM": False})

mol.update_geometry()
geom = np.array(mol.geometry())

cart_x = np.random.rand(npoints)
cart_y = np.random.rand(npoints)
cart_z = np.random.rand(npoints)
weights = np.random.rand(npoints)
#print(cart_x, cart_y, cart_z)

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
p4_points.set_deriv(2)
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
    V3 = np.zeros((npoints))
    for K in range(nprim):
        T1 = norm[K] * np.exp(- alpha[K] * R2)
        T2 = -2.0 * alpha[K] * T1
        T3 = -2.0 * alpha[K] * T2
        V1 += T1
        V2 += T2
        V3 += T3

    S = V1.copy()
    SX = V2 * xc
    SY = V2 * yc
    SZ = V2 * zc
    SXY = V3 * xc * yc
    SXZ = V3 * xc * zc
    SYZ = V3 * yc * zc
    SXX = V3 * xc * xc + V2
    SYY = V3 * yc * yc + V2
    SZZ = V3 * zc * zc + V2
    
    # SX, SY, SZ, SXX, SXZ, SXZ, SYY, SYZ, SZZ
    
    xc_pow = np.zeros((L + 3, npoints))
    yc_pow = np.zeros((L + 3, npoints))
    zc_pow = np.zeros((L + 3, npoints))
    
    xc_pow[0] = 0.0
    yc_pow[0] = 0.0
    zc_pow[0] = 0.0
    xc_pow[1] = 0.0
    yc_pow[1] = 0.0
    zc_pow[1] = 0.0
    xc_pow[2] = 1.0
    yc_pow[2] = 1.0
    zc_pow[2] = 1.0
    
    for LL in range(3, L + 3):
        xc_pow[LL] = xc_pow[LL - 1] * xc
        yc_pow[LL] = yc_pow[LL - 1] * yc
        zc_pow[LL] = zc_pow[LL - 1] * zc
   
    ncart = int((L + 1) * (L + 2) / 2) 
    output = {}
    output["PHI"] = np.zeros((ncart, npoints))
    output["PHI_X"] = np.zeros((ncart, npoints))
    output["PHI_Y"] = np.zeros((ncart, npoints))
    output["PHI_Z"] = np.zeros((ncart, npoints))
    output["PHI_XX"] = np.zeros((ncart, npoints))
    output["PHI_YY"] = np.zeros((ncart, npoints))
    output["PHI_ZZ"] = np.zeros((ncart, npoints))
    output["PHI_XY"] = np.zeros((ncart, npoints))
    output["PHI_XZ"] = np.zeros((ncart, npoints))
    output["PHI_YZ"] = np.zeros((ncart, npoints))

    idx = 0
    for i in range(L + 1):
        l = L - i + 2
        for j in range(i + 1):
            m = i - j + 2
            n = j + 2

            a_lp = l - 2
            a_mp = m - 2
            a_np = n - 2
            
        
            A = xc_pow[l] * yc_pow[m] * zc_pow[n];
            AX = a_lp * xc_pow[l - 1] * yc_pow[m] * zc_pow[n];
            AY = a_mp * xc_pow[l] * yc_pow[m - 1] * zc_pow[n];
            AZ = a_np * xc_pow[l] * yc_pow[m] * zc_pow[n - 1];
            AXY = a_lp * a_mp * xc_pow[l - 1] * yc_pow[m - 1] * zc_pow[n];
            AXZ = a_lp * a_np * xc_pow[l - 1] * yc_pow[m] * zc_pow[n - 1];
            AYZ = a_mp * a_np * xc_pow[l] * yc_pow[m - 1] * zc_pow[n - 1];
            AXX = a_lp * (a_lp - 1) * xc_pow[l - 2] * yc_pow[m] * zc_pow[n];
            AYY = a_mp * (a_mp - 1) * xc_pow[l] * yc_pow[m - 2] * zc_pow[n];
            AZZ = a_np * (a_np - 1) * xc_pow[l] * yc_pow[m] * zc_pow[n - 2];

            output["PHI"][idx] = S * A;
            output["PHI_X"][idx] = S * AX + SX * A;
            output["PHI_Y"][idx] = S * AY + SY * A;
            output["PHI_Z"][idx] = S * AZ + SZ * A;
            output["PHI_XX"][idx] = SXX * A + SX * AX + SX * AX + S * AXX;
            output["PHI_YY"][idx] = SYY * A + SY * AY + SY * AY + S * AYY;
            output["PHI_ZZ"][idx] = SZZ * A + SZ * AZ + SZ * AZ + S * AZZ;
            output["PHI_XY"][idx] = SXY * A + SX * AY + SY * AX + S * AXY;
            output["PHI_XZ"][idx] = SXZ * A + SX * AZ + SZ * AX + S * AXZ;
            output["PHI_YZ"][idx] = SYZ * A + SY * AZ + SZ * AY + S * AYZ;
            idx += 1

    return output

# Sum up g2g points
g2g_results = {k : [] for k in list(points)}
for shell in py_basis:
    for k, v in compute_shell(cart_x, cart_y, cart_z, shell).items():
        g2g_results[k].append(v)

# Transform both results
g2g_results = {k : np.vstack(v) for k, v in g2g_results.items()}
psi_results = {k : np.array(v).T for k, v in points.items()}

# Test each points
for k in list(g2g_results):
    close = np.allclose(g2g_results[k], psi_results[k])

    print("Psi4 and Numpy %6s same: %s" % (k, close))

# ... 
