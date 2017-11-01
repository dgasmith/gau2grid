"""
NumPy based example project.
"""

import numpy as np

# Generator
print("Order generator")
L = 3
for i in range(L + 1):
    l = L - i
    for j in range(i + 1):
        m = i - j
        n = j
        print(l, m, n)




npoints = 10
nprim = 3
alpha = [5, 6, 7]
norm = [0.25, 0.5, 0.25]
L = 2


alpha *= -1

x = np.random.rand(npoints)
y = np.random.rand(npoints)
z = np.random.rand(npoints)

R2 = x * x + y * y + z * z

S0 = np.zeros(npoints)
for K in range(nprim):
    S0 += norm[K] * np.exp(- alpha[K] * R2)

# SX, SY, SZ, SXX, SXZ, SXZ, SYY, SYZ, SZZ

xc_pow = np.zeros(L, npoints)
yc_pow = np.zeros(L, npoints)
zc_pow = np.zeros(L, npoints)

xc_pow[0] = 1.0
yc_pow[0] = 1.0
zc_pow[0] = 1.0

for LL in range(1, L + 1):
    xc_pow[LL] = xc_pow[LL - 1] * xc
    yc_pow[LL] = xc_pow[LL - 1] * yc
    zc_pow[LL] = xc_pow[LL - 1] * zc

output = np.zeros((LL * 2 + 1, npoints))

# S 
output[0] = S0 * xc_pow[0] * yc_pow[0] * zc_pow[0]

# P
output[0] = S0 * xc_pow[1] * yc_pow[0] * zc_pow[0]
output[1] = S0 * xc_pow[0] * yc_pow[1] * zc_pow[0]
output[2] = S0 * xc_pow[0] * yc_pow[0] * zc_pow[1]

# D
output[0] = S0 * xc_pow[2] * xc_pow[0] * xc_pow[0]
output[1] = S0 * xc_pow[1] * xc_pow[1] * xc_pow[0]
output[2] = S0 * xc_pow[1] * xc_pow[0] * xc_pow[1]
output[3] = S0 * xc_pow[0] * xc_pow[2] * xc_pow[0]
output[4] = S0 * xc_pow[0] * xc_pow[1] * xc_pow[1]
output[5] = S0 * xc_pow[0] * xc_pow[0] * xc_pow[2]


# ... 
