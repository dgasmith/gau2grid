API Reference
=============


.. code-block:: C

    // Return the maximum compiled angular momentum
    int max_L();
   
    // Returns the cartesian order ("row", "molden") 
    const char* cartesian_order();
    
    // Returns the spherical order ("cca", "gaussian") 
    const char* spherical_order();

.. code-block:: C

    void gg_naive_transpose(unsigned long n, unsigned long m,
                       const double* PRAGMA_RESTRICT input,
                       double* PRAGMA_RESTRICT output);

    void gg_fast_transpose(unsigned long n, unsigned long m,
                           const double* PRAGMA_RESTRICT input,
                           double* PRAGMA_RESTRICT output);

.. code-block:: C

    void gg_orbitals(int L,

        // Orbital data (norbs x nlength)
        const double* PRAGMA_RESTRICT orb,
        const unsigned long norbs,

        // XYZ grid points
        const unsigned long npoints,
        const double* PRAGMA_RESTRICT x,
        const double* PRAGMA_RESTRICT y,
        const double* PRAGMA_RESTRICT z,

        // Gaussian information
        const int nprim,
        const double* PRAGMA_RESTRICT coeffs,
        const double* PRAGMA_RESTRICT exponents,
        const double* PRAGMA_RESTRICT center,

        // Spherical transform true (1)/ false (0)
        const int spherical,

        // Ouput
        double* PRAGMA_RESTRICT phi_out);

.. code-block:: C

    void gg_collocation(int L,

        // XYZ grid points
        const unsigned long npoints,
        const double* PRAGMA_RESTRICT x,
        const double* PRAGMA_RESTRICT y,
        const double* PRAGMA_RESTRICT z,

        // Gaussian information
        const int nprim,
        const double* PRAGMA_RESTRICT coeffs,
        const double* PRAGMA_RESTRICT exponents,
        const double* PRAGMA_RESTRICT center,

        // Spherical transform true (1)/ false (0)
        const int spherical,

        // Ouput
        double* PRAGMA_RESTRICT phi_out);

.. code-block:: C

    void gg_collocation_deriv1(int L,

        // XYZ grid points
        const unsigned long npoints,
        const double* PRAGMA_RESTRICT x,
        const double* PRAGMA_RESTRICT y,
        const double* PRAGMA_RESTRICT z,

        // Gaussian information
        const int nprim,
        const double* PRAGMA_RESTRICT coeffs,
        const double* PRAGMA_RESTRICT exponents,
        const double* PRAGMA_RESTRICT center,

        // Spherical transform true (1)/ false (0)
        const int spherical,

        // Ouput
        double* PRAGMA_RESTRICT phi_out,
        double* PRAGMA_RESTRICT phi_x_out,
        double* PRAGMA_RESTRICT phi_y_out,
        double* PRAGMA_RESTRICT phi_z_out);

.. code-block:: C

    void gg_collocation_deriv2(int L,

        // XYZ grid points
        const unsigned long npoints,
        const double* PRAGMA_RESTRICT x,
        const double* PRAGMA_RESTRICT y,
        const double* PRAGMA_RESTRICT z,

        // Gaussian information
        const int nprim,
        const double* PRAGMA_RESTRICT coeffs,
        const double* PRAGMA_RESTRICT exponents,
        const double* PRAGMA_RESTRICT center,

        // Spherical transform true (1)/ false (0)
        const int spherical,

        // Ouput
        double* PRAGMA_RESTRICT phi_out,
        double* PRAGMA_RESTRICT phi_x_out,
        double* PRAGMA_RESTRICT phi_y_out,
        double* PRAGMA_RESTRICT phi_z_out,
        double* PRAGMA_RESTRICT phi_xx_out,
        double* PRAGMA_RESTRICT phi_xy_out,
        double* PRAGMA_RESTRICT phi_xz_out,
        double* PRAGMA_RESTRICT phi_yy_out,
        double* PRAGMA_RESTRICT phi_yz_out,
        double* PRAGMA_RESTRICT phi_zz_out);
