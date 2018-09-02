Collocation Example
===================

Single Collocation
------------------

A collocation grid between a single basis and a Cartesian grid can be computed
with the :func:`~gau2grid.collocation` function. For example, we will use a grid
starting at the origin along the ``z`` axis:

.. code-block:: C

    // Generate grid
    long int npoints = 5; 
    double x[5] = {0, 0, 0, 0, 0}; 
    double y[5] = {0, 0, 0, 0, 0}; 
    double z[5] = {0, 1, 2, 3, 4}; 

    
    // Gaussian data
    int nprim = 0
    double coef[1] = {1};
    double exp[1] = {1};
    double center[3] = {0, 0, 0};
    int spherical = 0;

To compute the collocation matrix for the above gaussian and grid we can call ``gg_collocation`` as follows:

.. code-block:: C

    double s_output[5] = {0}; 
    gg_collocation(0,

                   // Grid data
                   npoints, x, y, z,
                    
                   // Gaussian data
                   nprim, coef, exp, center, spherical,
                    
                   // Output data
                   s_output)

For higher angular momentum functions that output size should ``ncomponents x
npoints`` in size. Where each component is on a unique row or the ``X``
component starts at position ``0``, the ``Y`` component starts at position
``5``, and the ``Z`` component starts at position ``10`` as out grid is of
length ``5``. 

.. code-block:: C

    double p_output[15] = {0}; 
    gg_collocation(1,

                   // Grid data
                   npoints, x, y, z,
                    
                   // Gaussian data
                   nprim, coef, exp, center, spherical,
                    
                   // Output data
                   p_output)
