import numpy as np

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
        yc_pow[LL] = yc_pow[LL - 1] * yc
        zc_pow[LL] = zc_pow[LL - 1] * zc

    ncart = int((L + 1) * (L + 2) / 2)
    output = {}
    output["PHI"] = np.zeros((ncart, npoints))

    idx = 0
    for i in range(L + 1):
        l = L - i
        for j in range(i + 1):
            m = i - j
            n = j
            output["PHI"][idx] = S0 * xc_pow[l] * yc_pow[m] * zc_pow[n]
            idx += 1

    return output
