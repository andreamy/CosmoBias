1) May need to remove the function **int recompute_cosmo3D_CLASS(cosmopara C)** in the file CosmoCov/cosmolike_core/theory/recompute.c (in order to avoid errors when running CosmoCov with CLASS).

2) The file CosmoCov/cosmolike_core/theory/run_covariances_real_bin_fft.c  has to be substituted by the one in here in order to include the modified function **void run_cov_shear_shear_real_binned(char *SIMULATION, char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm1, int pm2, int start)**. This allows to print the xipm output. 

3) The file cosmo2D_fullsky_xipm.c computes the xipm from C_ells. It needs to be added to the folder theory/

4) The file CosmoCov/cosmolike_core/theory/cosmo3D.c has to be substituted by the one in here in order to use CLASS for the matter power spectra computation.
