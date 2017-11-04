"""
NumPy based example project.
"""

import numpy as np
import psi4
import time
np.set_printoptions(linewidth=120)

mol = psi4.geometry("""
#He 0 0 0
He 0 0 2
no_com
no_reorient
""")

# Tweakers
npoints = 1000000
basis = "cc-pV5Z"
#basis = "cc-pVQZ"


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
t = time.time()
p4_points.compute_functions(block)
print("Psi4 : Collocation build %5.3f" % (time.time() - t))
points = p4_points.basis_values()

sttime = 0
order = []
def compute_shell(x, y, z, shell, grad=0):
    global sttime

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
        tmp = - alpha[K] * R2
        np.exp(tmp, out=tmp)
        tmp *= norm[K]
        S0 += tmp 
        #S0 += norm[K] * np.exp(- alpha[K] * R2)
    
    # SX, SY, SZ, SXX, SXZ, SXZ, SYY, SYZ, SZZ
    
    xc_pow = np.zeros((L + 1, npoints))
    yc_pow = np.zeros((L + 1, npoints))
    zc_pow = np.zeros((L + 1, npoints))
    
    xc_pow[:] = 1.0
    yc_pow[:] = 1.0
    zc_pow[:] = 1.0
  
    for LL in range(1, L + 1):
        xc_pow[LL:] *= xc
        yc_pow[LL:] *= yc
        zc_pow[LL:] *= zc
        #xc_pow[LL] = xc_pow[LL - 1] * xc
        #yc_pow[LL] = yc_pow[LL - 1] * yc
        #zc_pow[LL] = zc_pow[LL - 1] * zc
   
    ncart = int((L + 1) * (L + 2) / 2) 
    output = {}
    output["PHI"] = np.zeros((ncart, npoints))
    phi_out = output["PHI"]

    t = time.time()
    idx = 0
    if L == 0:
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[0], zc_pow[0], out=phi_out[0])

    elif L == 1:
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[0], zc_pow[0], out=phi_out[0])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[1], zc_pow[0], out=phi_out[1])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[0], zc_pow[1], out=phi_out[2])
    elif L == 2:
        np.einsum('p,p,p,p->p', S0, xc_pow[2], yc_pow[0], zc_pow[0], out=phi_out[0])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[1], zc_pow[0], out=phi_out[1])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[0], zc_pow[1], out=phi_out[2])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[2], zc_pow[0], out=phi_out[3])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[1], zc_pow[1], out=phi_out[4])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[0], zc_pow[2], out=phi_out[5])
    elif L == 3:
        np.einsum('p,p,p,p->p', S0, xc_pow[3], yc_pow[0], zc_pow[0], out=phi_out[0])
        np.einsum('p,p,p,p->p', S0, xc_pow[2], yc_pow[1], zc_pow[0], out=phi_out[1])
        np.einsum('p,p,p,p->p', S0, xc_pow[2], yc_pow[0], zc_pow[1], out=phi_out[2])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[2], zc_pow[0], out=phi_out[3])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[1], zc_pow[1], out=phi_out[4])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[0], zc_pow[2], out=phi_out[5])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[3], zc_pow[0], out=phi_out[6])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[2], zc_pow[1], out=phi_out[7])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[1], zc_pow[2], out=phi_out[8])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[0], zc_pow[3], out=phi_out[9])

    elif L == 4:
        np.einsum('p,p,p,p->p', S0, xc_pow[4], yc_pow[0], zc_pow[0], out=phi_out[0])
        np.einsum('p,p,p,p->p', S0, xc_pow[3], yc_pow[1], zc_pow[0], out=phi_out[1])
        np.einsum('p,p,p,p->p', S0, xc_pow[3], yc_pow[0], zc_pow[1], out=phi_out[2])
        np.einsum('p,p,p,p->p', S0, xc_pow[2], yc_pow[2], zc_pow[0], out=phi_out[3])
        np.einsum('p,p,p,p->p', S0, xc_pow[2], yc_pow[1], zc_pow[1], out=phi_out[4])
        np.einsum('p,p,p,p->p', S0, xc_pow[2], yc_pow[0], zc_pow[2], out=phi_out[5])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[3], zc_pow[0], out=phi_out[6])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[2], zc_pow[1], out=phi_out[7])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[1], zc_pow[2], out=phi_out[8])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[0], zc_pow[3], out=phi_out[9])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[4], zc_pow[0], out=phi_out[10])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[3], zc_pow[1], out=phi_out[11])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[2], zc_pow[2], out=phi_out[12])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[1], zc_pow[3], out=phi_out[13])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[0], zc_pow[4], out=phi_out[14])

    elif L == 5:
        np.einsum('p,p,p,p->p', S0, xc_pow[5], yc_pow[0], zc_pow[0], out=phi_out[0])
        np.einsum('p,p,p,p->p', S0, xc_pow[4], yc_pow[1], zc_pow[0], out=phi_out[1])
        np.einsum('p,p,p,p->p', S0, xc_pow[4], yc_pow[0], zc_pow[1], out=phi_out[2])
        np.einsum('p,p,p,p->p', S0, xc_pow[3], yc_pow[2], zc_pow[0], out=phi_out[3])
        np.einsum('p,p,p,p->p', S0, xc_pow[3], yc_pow[1], zc_pow[1], out=phi_out[4])
        np.einsum('p,p,p,p->p', S0, xc_pow[3], yc_pow[0], zc_pow[2], out=phi_out[5])
        np.einsum('p,p,p,p->p', S0, xc_pow[2], yc_pow[3], zc_pow[0], out=phi_out[6])
        np.einsum('p,p,p,p->p', S0, xc_pow[2], yc_pow[2], zc_pow[1], out=phi_out[7])
        np.einsum('p,p,p,p->p', S0, xc_pow[2], yc_pow[1], zc_pow[2], out=phi_out[8])
        np.einsum('p,p,p,p->p', S0, xc_pow[2], yc_pow[0], zc_pow[3], out=phi_out[9])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[4], zc_pow[0], out=phi_out[10])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[3], zc_pow[1], out=phi_out[11])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[2], zc_pow[2], out=phi_out[12])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[1], zc_pow[3], out=phi_out[13])
        np.einsum('p,p,p,p->p', S0, xc_pow[1], yc_pow[0], zc_pow[4], out=phi_out[14])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[5], zc_pow[0], out=phi_out[15])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[4], zc_pow[1], out=phi_out[16])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[3], zc_pow[2], out=phi_out[17])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[2], zc_pow[3], out=phi_out[18])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[1], zc_pow[4], out=phi_out[19])
        np.einsum('p,p,p,p->p', S0, xc_pow[0], yc_pow[0], zc_pow[5], out=phi_out[20])
    sttime += time.time() - t
        

    return output

t = time.time()
PHI_tmp = []
for shell in py_basis:
    tmp = compute_shell(cart_x, cart_y, cart_z, shell)
    PHI_tmp.append(tmp["PHI"])
print("NumPy: Collocation build %5.3f" % (time.time() - t))

PHI = np.vstack(PHI_tmp)
psi_PHI = points["PHI"].np.T
print("Psi4 vs NumPy Matches? %s" % np.allclose(PHI, psi_PHI))
#print(np.where(np.any(np.abs(PHI - psi_PHI) > 1.e-14, axis=1)))
#wrong = np.any(np.abs(PHI - psi_PHI) > 1.e-14, axis=1)
#print(np.sum(np.array(order)[wrong], axis=1))
##print(PHI[wrong])
##print(psi_PHI[wrong])
#
#print(np.allclose(PHI, psi_PHI))

# ... 
