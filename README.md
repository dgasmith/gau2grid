<p align="center">
<a href="https://travis-ci.org/dgasmith/gau2grid">
  <img src="https://travis-ci.org/dgasmith/gau2grid.svg?branch=master" alt="Travis CI"/>
</a>
<a href="https://codecov.io/gh/dgasmith/gau2grid">
  <img src="https://codecov.io/gh/dgasmith/gau2grid/branch/master/graph/badge.svg" alt="Codecov" />
</a>
</p>

# gau2grid
A collocation code for computing gaussians on a grid of the form:
```
    ret_Lp = x^l y^m z^n \sum_i coeff_i e^(exponent_i * (|center - p|)^2)
```
Where the return matrix in the angular momentum by number of requested points.

```python
import gau2grid
import numpy as np

# Build coordinates along the Z axis
>>> xyz = np.zeros((3, 5))
>>> xyz[2] = np.arange(5)

# Compute a 's' gaussian with a scaling and exponent of one at the origin
>>> ret = gau2grid.collocation(xyz, 0, [1], [1], [0, 0, 0])
>>> print(ret["PHI"])
[[  1.00000e+00   3.67879e-01   1.83156e-02   1.23409e-04   1.12535e-07]]

# Compute a 'p' gaussian with a scaling and exponent of one at the origin
>>> ret = gau2grid.collocation(xyz, 1, [1], [1], [0, 0, 0], spherical=False)
>>> print(ret["PHI"])
[[  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
 [  0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00   0.00000e+00]
 [  0.00000e+00   3.67879e-01   3.66312e-02   3.70229e-04   4.50140e-07]]

# Note that the X and Y components are zero as they are orthogonal to our Z vector. 
```

The returned matrix can be in either cartesian or regular solid harmonics. There currently
exists three algorithms in which to compute these collocation matrices:
 - NumPy Reference: A 
 - Optimized/Generated NumPy: A exploratory tool to 
 - Optimize C: A simple dy

## Building Gau2Grid
The C library is currently built with CMake and currently has C no required dependancies
other than the standard library. A simple CMake and build example can found below:

```python
cmake -H. -Bobjdir \
    -DPYTHON_EXECUTABLE='python_path'/bin/python \
    -DMAX_AM=6 \
    -DCMAKE_INSTALL_PREFIX="/usr/local/gau2grid"

cd objdir; make -j2
```

## Python installation
The gau2grid program (without the optimized C library) can be installed using
the canonical `setup.py` script,
```
python setup.py install
```

# Authors
This code was inspired by a number of folks and quite a few provided excellent advice.

 - Daniel G. A. Smith - Code author
 - Rob M. Parrish - Author of the Psi4 section which contains the original equations
 - Lori A. Burns - CMake, building, and library linking
 - Andy C. Simmonett - RSH coefficients
 - Ben Pritchard - Generator and vectorization recommendations

