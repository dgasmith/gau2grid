Collocation Example
===================

Single Basis Functions
----------------------

A collocation grid between a single basis and a Cartesian grid can be computed
with the :c:func:`~gg_collocation` function. For example, we will use a grid
starting at the origin along the ``z`` axis and a ``S`` shell at the origin:

.. code-block:: C

  #include <stdio.h>
  #include "gau2grid.h"

  int main() {
      // Generate grid
      long int npoints = 5;
      double x[5] = {0, 0, 0, 0, 0};
      double y[5] = {0, 0, 0, 0, 0};
      double z[5] = {0, 1, 2, 3, 4};

      // Gaussian data
      int nprim = 1;
      double coef[1] = {1};
      double exp[1] = {1};
      double center[3] = {0, 0, 0};
      int spherical = 0; // Use cartesian components

      double s_output[5] = {0};
      gg_collocation(0,                                   // The angular momentum
                     npoints, x, y, z,                    // Grid data
                     nprim, coef, exp, center, spherical, // Gaussian data
                     s_output);                           // Output

      // Print output to stdout
      for (int i = 0; i < npoints; i += 1) {
          printf("%lf  ", s_output[i]);
      }
      printf("\n");
  }

The resulting output should be:

.. code-block:: bash

  1.000000  0.367879  0.018316  0.000123  0.000000

For higher angular momentum functions that output size should ``ncomponents x
npoints`` in size. Where each component is on a unique row or the ``X``
component starts at position ``0``, the ``Y`` component starts at position
``5``, and the ``Z`` component starts at position ``10`` as out grid is of
length ``5``. See `Gaussian Component Orders`_ for more details or order output.

Multiple Basis Functions
------------------------

Often collocation matrices are computed for multiple basis functions at once.
The below is an example of usage:

.. code-block:: C

  #include <stdio.h>
  #include "gau2grid.h"

  int main() {
      // Generate grid
      long int npoints = 5;
      double x[5] = {0, 0, 0, 0, 0};
      double y[5] = {0, 0, 0, 0, 0};
      double z[5] = {0, 1, 2, 3, 4};

      // Gaussian data
      int nprim = 1;
      double coef[1] = {1};
      double exp[1] = {1};
      double center[3] = {0, 0, 0};
      int spherical = 1; // Use spherical components

      // Size ncomponents * npoints, (1 + 3 + 5) * 5
      double output[45] = {0};
      int row = 0;
      for (int L = 0; L < 3; L++) {
          gg_collocation(L,                                    // The angular momentum
                         npoints, x, y, z,                     // Grid data
                         nprim, coef, exp, center, spherical,  // Gaussian data
                         output + (row * npoints));            // Output, shift pointer

          row += gg_ncomponents(L, spherical); // Increment rows skipped
      }

      // Print out by row
      for (int i = 0; i < row; i += 1) {
          for (int j = 0; j < npoints; j += 1) {
              printf("%lf  ", output[i * npoints + j]);
          }
          printf("\n");
      }
  }

The resulting output should be:

.. code-block:: bash

  1.000000  0.367879  0.018316  0.000123  0.000000 // S
  0.000000  0.367879  0.036631  0.000370  0.000000 // P_0
  0.000000  0.000000  0.000000  0.000000  0.000000 // P^+_0
  0.000000  0.000000  0.000000  0.000000  0.000000 // P^-_0
  0.000000  0.367879  0.073263  0.001111  0.000002 // D_0
  0.000000  0.000000  0.000000  0.000000  0.000000 // D^+_1
  0.000000  0.000000  0.000000  0.000000  0.000000 // D^-_1
  0.000000  0.000000  0.000000  0.000000  0.000000 // D^+_2
  0.000000  0.000000  0.000000  0.000000  0.000000 // D^-_2
