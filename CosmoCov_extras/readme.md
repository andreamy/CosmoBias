**Files to be added/replaced**:

- The file CosmoCov/cosmolike_core/theory/run_covariances_real_bin_fft.c  needs to be replaced by the one in here in order to include the modified function **void run_cov_shear_shear_real_binned(char *SIMULATION, char *OUTFILE, char *PATH, double *t, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm1, int pm2, int start)**. This allows to print the xipm output. 

- The file cosmo2D_fullsky_xipm.c computes the xipm from C_ells. It needs to be added to the folder theory/

- The file CosmoCov/cosmolike_core/theory/cosmo3D.c needs to be replaced by the one in here in order to use CLASS for the matter power spectra computation.

- The Makefile in Cosmolike/covs/ needs to be replaced by the one in here in order to call CLASS from CosmoCov. 

- The file Cosmolike/covs/init.c needs to be replaced by the one here in order to compute the angular binning as in KiDS-1000, i.e. logarithmic spaced bins and centered angles instead of using lower bound thetas for each bin. It also adds the code snippet necessary to pass the variable SIMULATION from the ini_file. 

**Minor modifications:**

-May need to remove the function **int recompute_cosmo3D_CLASS(cosmopara C)** in the file CosmoCov/cosmolike_core/theory/recompute.c (in order to avoid errors when running CosmoCov with CLASS).
 
- The following line must be added in the typedef struc covpar of the file CosmoCov/cosmolike_core/theory/structs.c

    char simulation[200]; /* simulation name containing info about changed parameters*/
    
- In the file CosmoCov/covs/compute_covariances_real_flat_fft.c, the following modifications have to be made:
  - Add to the include section: #include "../cosmolike_core/theory/cosmo2D_fullsky_xipm.c" 
  - To follow the new notation of the modified function **run_cov_shear_shear_real_binned()**, replace the line below by the new one 
    (and repeat for any time that you call the function):  
    
      OLD: run_cov_shear_shear_real_bin(OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,1,1,k);
      
      NEW: run_cov_shear_shear_real_bin(covparams.simulation, OUTFILE,covparams.outdir,thetamin,dtheta,Ntheta,l,m,1,1,k);
        
 
