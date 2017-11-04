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
#npoints = 4
npoints = 100000
#basis = "sto-3g"
basis = "cc-pV5Z"


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
p4_points.set_deriv(1)
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
    
    S0 = np.zeros((npoints), dtype=np.float)
    V1 = np.zeros((npoints), dtype=np.float)
    V2 = np.zeros((npoints), dtype=np.float)
    for K in range(nprim):
        #tmp = - alpha[K] * R2
        #np.exp(tmp, out=tmp)
        #tmp *= norm[K]
        T1 = norm[K] * np.exp(-alpha[K] * R2)
        T2 = -2.0 * alpha[K] * T1

        V1 += T1 
        V2 += T2
        #S0 += norm[K] * np.exp(- alpha[K] * R2)
    
    S0 = V1
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
        #xc_pow[LL:] *= xc
        #yc_pow[LL:] *= yc
        #zc_pow[LL:] *= zc
        xc_pow[LL] = xc_pow[LL - 1] * xc
        yc_pow[LL] = yc_pow[LL - 1] * yc
        zc_pow[LL] = zc_pow[LL - 1] * zc
   
    ncart = int((L + 1) * (L + 2) / 2) 
    output = {}
    output["PHI"] = np.zeros((ncart, npoints))
    output["PHI_X"] = np.zeros((ncart, npoints))
    output["PHI_Y"] = np.zeros((ncart, npoints))
    output["PHI_Z"] = np.zeros((ncart, npoints))

    phi_out = output["PHI"]
    phi_out_x = output["PHI_X"]
    phi_out_y = output["PHI_Y"]
    phi_out_z = output["PHI_Z"]

    t = time.time()
    idx = 0
    if L == 0:
        # Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[1], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[0])
        # Derivatives
        phi_out_x[0] += SX * xyz
        phi_out_y[0] += SY * xyz
        phi_out_z[0] += SZ * xyz
    if L == 1:
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[1], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[0])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[1], zc_pow[1], out=phi_out_x[0])
        phi_out_x[0] += SX * xyz
        phi_out_y[0] += SY * xyz
        phi_out_z[0] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[2], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[1])
# Derivatives
        phi_out_x[1] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[1], out=phi_out_y[1])
        phi_out_y[1] += SY * xyz
        phi_out_z[1] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[1], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[2])
# Derivatives
        phi_out_x[2] += SX * xyz
        phi_out_y[2] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[1], zc_pow[2], out=phi_out_z[2])
        phi_out_z[2] += SZ * xyz
    elif L == 2:
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[3], yc_pow[1], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[0])
# Derivatives
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[1], zc_pow[1], out=phi_out_x[0])
        phi_out_x[0] += SX * xyz
        phi_out_y[0] += SY * xyz
        phi_out_z[0] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[2], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[1])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[1], out=phi_out_x[1])
        phi_out_x[1] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[1], out=phi_out_y[1])
        phi_out_y[1] += SY * xyz
        phi_out_z[1] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[1], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[2])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[1], zc_pow[2], out=phi_out_x[2])
        phi_out_x[2] += SX * xyz
        phi_out_y[2] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[1], zc_pow[2], out=phi_out_z[2])
        phi_out_z[2] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[3], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[3])
# Derivatives
        phi_out_x[3] += SX * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[1], out=phi_out_y[3])
        phi_out_y[3] += SY * xyz
        phi_out_z[3] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[2], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[4])
# Derivatives
        phi_out_x[4] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[2], out=phi_out_y[4])
        phi_out_y[4] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[2], out=phi_out_z[4])
        phi_out_z[4] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[1], zc_pow[3])
        np.einsum('p,p->p', S0, xyz, out=phi_out[5])
# Derivatives
        phi_out_x[5] += SX * xyz
        phi_out_y[5] += SY * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[1], zc_pow[3], out=phi_out_z[5])
        phi_out_z[5] += SZ * xyz
    elif L == 3:
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[4], yc_pow[1], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[0])
# Derivatives
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[1], zc_pow[1], out=phi_out_x[0])
        phi_out_x[0] += SX * xyz
        phi_out_y[0] += SY * xyz
        phi_out_z[0] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[3], yc_pow[2], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[1])
# Derivatives
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[2], zc_pow[1], out=phi_out_x[1])
        phi_out_x[1] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[1], out=phi_out_y[1])
        phi_out_y[1] += SY * xyz
        phi_out_z[1] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[3], yc_pow[1], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[2])
# Derivatives
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[1], zc_pow[2], out=phi_out_x[2])
        phi_out_x[2] += SX * xyz
        phi_out_y[2] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[1], zc_pow[2], out=phi_out_z[2])
        phi_out_z[2] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[3], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[3])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[3], zc_pow[1], out=phi_out_x[3])
        phi_out_x[3] += SX * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[1], out=phi_out_y[3])
        phi_out_y[3] += SY * xyz
        phi_out_z[3] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[2], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[4])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[2], out=phi_out_x[4])
        phi_out_x[4] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[2], out=phi_out_y[4])
        phi_out_y[4] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[2], out=phi_out_z[4])
        phi_out_z[4] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[1], zc_pow[3])
        np.einsum('p,p->p', S0, xyz, out=phi_out[5])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[1], zc_pow[3], out=phi_out_x[5])
        phi_out_x[5] += SX * xyz
        phi_out_y[5] += SY * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[1], zc_pow[3], out=phi_out_z[5])
        phi_out_z[5] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[4], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[6])
# Derivatives
        phi_out_x[6] += SX * xyz
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[4], zc_pow[1], out=phi_out_y[6])
        phi_out_y[6] += SY * xyz
        phi_out_z[6] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[3], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[7])
# Derivatives
        phi_out_x[7] += SX * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[2], out=phi_out_y[7])
        phi_out_y[7] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[3], zc_pow[2], out=phi_out_z[7])
        phi_out_z[7] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[2], zc_pow[3])
        np.einsum('p,p->p', S0, xyz, out=phi_out[8])
# Derivatives
        phi_out_x[8] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[3], out=phi_out_y[8])
        phi_out_y[8] += SY * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[2], zc_pow[3], out=phi_out_z[8])
        phi_out_z[8] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[1], zc_pow[4])
        np.einsum('p,p->p', S0, xyz, out=phi_out[9])
# Derivatives
        phi_out_x[9] += SX * xyz
        phi_out_y[9] += SY * xyz
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[1], zc_pow[4], out=phi_out_z[9])
        phi_out_z[9] += SZ * xyz
    elif L == 4:
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[5], yc_pow[1], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[0])
# Derivatives
        np.einsum(',p,p,p,p->p', 4, S0, xc_pow[4], yc_pow[1], zc_pow[1], out=phi_out_x[0])
        phi_out_x[0] += SX * xyz
        phi_out_y[0] += SY * xyz
        phi_out_z[0] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[4], yc_pow[2], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[1])
# Derivatives
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[2], zc_pow[1], out=phi_out_x[1])
        phi_out_x[1] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[1], out=phi_out_y[1])
        phi_out_y[1] += SY * xyz
        phi_out_z[1] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[4], yc_pow[1], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[2])
# Derivatives
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[1], zc_pow[2], out=phi_out_x[2])
        phi_out_x[2] += SX * xyz
        phi_out_y[2] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[1], zc_pow[2], out=phi_out_z[2])
        phi_out_z[2] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[3], yc_pow[3], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[3])
# Derivatives
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[1], out=phi_out_x[3])
        phi_out_x[3] += SX * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[1], out=phi_out_y[3])
        phi_out_y[3] += SY * xyz
        phi_out_z[3] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[3], yc_pow[2], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[4])
# Derivatives
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[2], zc_pow[2], out=phi_out_x[4])
        phi_out_x[4] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[2], out=phi_out_y[4])
        phi_out_y[4] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[2], out=phi_out_z[4])
        phi_out_z[4] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[3], yc_pow[1], zc_pow[3])
        np.einsum('p,p->p', S0, xyz, out=phi_out[5])
# Derivatives
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[1], zc_pow[3], out=phi_out_x[5])
        phi_out_x[5] += SX * xyz
        phi_out_y[5] += SY * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[1], zc_pow[3], out=phi_out_z[5])
        phi_out_z[5] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[4], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[6])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[4], zc_pow[1], out=phi_out_x[6])
        phi_out_x[6] += SX * xyz
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[4], zc_pow[1], out=phi_out_y[6])
        phi_out_y[6] += SY * xyz
        phi_out_z[6] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[3], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[7])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[3], zc_pow[2], out=phi_out_x[7])
        phi_out_x[7] += SX * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[2], out=phi_out_y[7])
        phi_out_y[7] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[3], zc_pow[2], out=phi_out_z[7])
        phi_out_z[7] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[2], zc_pow[3])
        np.einsum('p,p->p', S0, xyz, out=phi_out[8])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[3], out=phi_out_x[8])
        phi_out_x[8] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[3], out=phi_out_y[8])
        phi_out_y[8] += SY * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[2], zc_pow[3], out=phi_out_z[8])
        phi_out_z[8] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[1], zc_pow[4])
        np.einsum('p,p->p', S0, xyz, out=phi_out[9])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[1], zc_pow[4], out=phi_out_x[9])
        phi_out_x[9] += SX * xyz
        phi_out_y[9] += SY * xyz
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[1], zc_pow[4], out=phi_out_z[9])
        phi_out_z[9] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[5], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[10])
# Derivatives
        phi_out_x[10] += SX * xyz
        np.einsum(',p,p,p,p->p', 4, S0, xc_pow[4], yc_pow[5], zc_pow[1], out=phi_out_y[10])
        phi_out_y[10] += SY * xyz
        phi_out_z[10] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[4], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[11])
# Derivatives
        phi_out_x[11] += SX * xyz
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[4], zc_pow[2], out=phi_out_y[11])
        phi_out_y[11] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[4], zc_pow[2], out=phi_out_z[11])
        phi_out_z[11] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[3], zc_pow[3])
        np.einsum('p,p->p', S0, xyz, out=phi_out[12])
# Derivatives
        phi_out_x[12] += SX * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[3], out=phi_out_y[12])
        phi_out_y[12] += SY * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[3], out=phi_out_z[12])
        phi_out_z[12] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[2], zc_pow[4])
        np.einsum('p,p->p', S0, xyz, out=phi_out[13])
# Derivatives
        phi_out_x[13] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[4], out=phi_out_y[13])
        phi_out_y[13] += SY * xyz
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[2], zc_pow[4], out=phi_out_z[13])
        phi_out_z[13] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[1], zc_pow[5])
        np.einsum('p,p->p', S0, xyz, out=phi_out[14])
# Derivatives
        phi_out_x[14] += SX * xyz
        phi_out_y[14] += SY * xyz
        np.einsum(',p,p,p,p->p', 4, S0, xc_pow[4], yc_pow[1], zc_pow[5], out=phi_out_z[14])
        phi_out_z[14] += SZ * xyz
    elif L == 5:
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[6], yc_pow[1], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[0])
# Derivatives
        np.einsum(',p,p,p,p->p', 5, S0, xc_pow[5], yc_pow[1], zc_pow[1], out=phi_out_x[0])
        phi_out_x[0] += SX * xyz
        phi_out_y[0] += SY * xyz
        phi_out_z[0] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[5], yc_pow[2], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[1])
# Derivatives
        np.einsum(',p,p,p,p->p', 4, S0, xc_pow[4], yc_pow[2], zc_pow[1], out=phi_out_x[1])
        phi_out_x[1] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[1], out=phi_out_y[1])
        phi_out_y[1] += SY * xyz
        phi_out_z[1] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[5], yc_pow[1], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[2])
# Derivatives
        np.einsum(',p,p,p,p->p', 4, S0, xc_pow[4], yc_pow[1], zc_pow[2], out=phi_out_x[2])
        phi_out_x[2] += SX * xyz
        phi_out_y[2] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[1], zc_pow[2], out=phi_out_z[2])
        phi_out_z[2] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[4], yc_pow[3], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[3])
# Derivatives
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[3], zc_pow[1], out=phi_out_x[3])
        phi_out_x[3] += SX * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[1], out=phi_out_y[3])
        phi_out_y[3] += SY * xyz
        phi_out_z[3] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[4], yc_pow[2], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[4])
# Derivatives
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[2], zc_pow[2], out=phi_out_x[4])
        phi_out_x[4] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[2], out=phi_out_y[4])
        phi_out_y[4] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[2], out=phi_out_z[4])
        phi_out_z[4] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[4], yc_pow[1], zc_pow[3])
        np.einsum('p,p->p', S0, xyz, out=phi_out[5])
# Derivatives
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[1], zc_pow[3], out=phi_out_x[5])
        phi_out_x[5] += SX * xyz
        phi_out_y[5] += SY * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[1], zc_pow[3], out=phi_out_z[5])
        phi_out_z[5] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[3], yc_pow[4], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[6])
# Derivatives
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[4], zc_pow[1], out=phi_out_x[6])
        phi_out_x[6] += SX * xyz
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[4], zc_pow[1], out=phi_out_y[6])
        phi_out_y[6] += SY * xyz
        phi_out_z[6] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[3], yc_pow[3], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[7])
# Derivatives
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[2], out=phi_out_x[7])
        phi_out_x[7] += SX * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[2], out=phi_out_y[7])
        phi_out_y[7] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[3], zc_pow[2], out=phi_out_z[7])
        phi_out_z[7] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[3], yc_pow[2], zc_pow[3])
        np.einsum('p,p->p', S0, xyz, out=phi_out[8])
# Derivatives
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[2], zc_pow[3], out=phi_out_x[8])
        phi_out_x[8] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[3], out=phi_out_y[8])
        phi_out_y[8] += SY * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[2], zc_pow[3], out=phi_out_z[8])
        phi_out_z[8] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[3], yc_pow[1], zc_pow[4])
        np.einsum('p,p->p', S0, xyz, out=phi_out[9])
# Derivatives
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[1], zc_pow[4], out=phi_out_x[9])
        phi_out_x[9] += SX * xyz
        phi_out_y[9] += SY * xyz
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[1], zc_pow[4], out=phi_out_z[9])
        phi_out_z[9] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[5], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[10])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[5], zc_pow[1], out=phi_out_x[10])
        phi_out_x[10] += SX * xyz
        np.einsum(',p,p,p,p->p', 4, S0, xc_pow[4], yc_pow[5], zc_pow[1], out=phi_out_y[10])
        phi_out_y[10] += SY * xyz
        phi_out_z[10] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[4], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[11])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[4], zc_pow[2], out=phi_out_x[11])
        phi_out_x[11] += SX * xyz
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[4], zc_pow[2], out=phi_out_y[11])
        phi_out_y[11] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[4], zc_pow[2], out=phi_out_z[11])
        phi_out_z[11] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[3], zc_pow[3])
        np.einsum('p,p->p', S0, xyz, out=phi_out[12])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[3], zc_pow[3], out=phi_out_x[12])
        phi_out_x[12] += SX * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[3], out=phi_out_y[12])
        phi_out_y[12] += SY * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[3], out=phi_out_z[12])
        phi_out_z[12] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[2], zc_pow[4])
        np.einsum('p,p->p', S0, xyz, out=phi_out[13])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[4], out=phi_out_x[13])
        phi_out_x[13] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[4], out=phi_out_y[13])
        phi_out_y[13] += SY * xyz
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[2], zc_pow[4], out=phi_out_z[13])
        phi_out_z[13] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[2], yc_pow[1], zc_pow[5])
        np.einsum('p,p->p', S0, xyz, out=phi_out[14])
# Derivatives
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[1], zc_pow[5], out=phi_out_x[14])
        phi_out_x[14] += SX * xyz
        phi_out_y[14] += SY * xyz
        np.einsum(',p,p,p,p->p', 4, S0, xc_pow[4], yc_pow[1], zc_pow[5], out=phi_out_z[14])
        phi_out_z[14] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[6], zc_pow[1])
        np.einsum('p,p->p', S0, xyz, out=phi_out[15])
# Derivatives
        phi_out_x[15] += SX * xyz
        np.einsum(',p,p,p,p->p', 5, S0, xc_pow[5], yc_pow[6], zc_pow[1], out=phi_out_y[15])
        phi_out_y[15] += SY * xyz
        phi_out_z[15] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[5], zc_pow[2])
        np.einsum('p,p->p', S0, xyz, out=phi_out[16])
# Derivatives
        phi_out_x[16] += SX * xyz
        np.einsum(',p,p,p,p->p', 4, S0, xc_pow[4], yc_pow[5], zc_pow[2], out=phi_out_y[16])
        phi_out_y[16] += SY * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[5], zc_pow[2], out=phi_out_z[16])
        phi_out_z[16] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[4], zc_pow[3])
        np.einsum('p,p->p', S0, xyz, out=phi_out[17])
# Derivatives
        phi_out_x[17] += SX * xyz
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[4], zc_pow[3], out=phi_out_y[17])
        phi_out_y[17] += SY * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[4], zc_pow[3], out=phi_out_z[17])
        phi_out_z[17] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[3], zc_pow[4])
        np.einsum('p,p->p', S0, xyz, out=phi_out[18])
# Derivatives
        phi_out_x[18] += SX * xyz
        np.einsum(',p,p,p,p->p', 2, S0, xc_pow[2], yc_pow[3], zc_pow[4], out=phi_out_y[18])
        phi_out_y[18] += SY * xyz
        np.einsum(',p,p,p,p->p', 3, S0, xc_pow[3], yc_pow[3], zc_pow[4], out=phi_out_z[18])
        phi_out_z[18] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[2], zc_pow[5])
        np.einsum('p,p->p', S0, xyz, out=phi_out[19])
# Derivatives
        phi_out_x[19] += SX * xyz
        np.einsum(',p,p,p,p->p', 1, S0, xc_pow[1], yc_pow[2], zc_pow[5], out=phi_out_y[19])
        phi_out_y[19] += SY * xyz
        np.einsum(',p,p,p,p->p', 4, S0, xc_pow[4], yc_pow[2], zc_pow[5], out=phi_out_z[19])
        phi_out_z[19] += SZ * xyz
# Density
        xyz = np.einsum('p,p,p->p', xc_pow[1], yc_pow[1], zc_pow[6])
        np.einsum('p,p->p', S0, xyz, out=phi_out[20])
# Derivatives
        phi_out_x[20] += SX * xyz
        phi_out_y[20] += SY * xyz
        np.einsum(',p,p,p,p->p', 5, S0, xc_pow[5], yc_pow[1], zc_pow[6], out=phi_out_z[20])
        phi_out_z[20] += SZ * xyz
    sttime += time.time() - t
        

    return output

t = time.time()
PHI_tmp = []
PHI_tmp_x = []
for shell in py_basis:
    tmp = compute_shell(cart_x, cart_y, cart_z, shell)
    PHI_tmp.append(tmp["PHI"])
    PHI_tmp_x.append(tmp["PHI_X"])
print("NumPy: Collocation build %5.3f" % (time.time() - t))

PHI = np.vstack(PHI_tmp)
psi_PHI = points["PHI"].np.T
PHI_x = np.vstack(PHI_tmp)
psi_PHI_x = points["PHI_X"].np.T
#print(PHI_x)
#print(psi_PHI_x)
print("Psi4 vs NumPy Matches? %s" % np.allclose(PHI, psi_PHI))
#print(np.allclose(PHI_x, psi_PHI_x))
#print(np.where(np.any(np.abs(PHI - psi_PHI) > 1.e-14, axis=1)))
#wrong = np.any(np.abs(PHI - psi_PHI) > 1.e-14, axis=1)
#print(np.sum(np.array(order)[wrong], axis=1))
##print(PHI[wrong])
##print(psi_PHI[wrong])
#
#print(np.allclose(PHI, psi_PHI))

# ... 
