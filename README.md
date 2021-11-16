# Biasing_covmat

**Prerequisites**: 
- Both CosmoCov and CLASS v=>3.0 need to be installed in order to use CosmoBias. Some files (CLASS_wrapper.c, cosmo2D_fullsky_xipm.c, run_covariances_real_bin_fft.c) have to be added/modified in CosmoCov/cosmolike_core/theory/ in order to get the xipm data vector and to run HMCode via CLASS. These can be found in the folder 'CosmoCov_extras'.

- Cosmocov needs to be placed in the main directory CosmoBias/ (e.g. at the same level as notebooks), and CLASS should be installed in CosmoCov/cosmolike_core/.

**Quick step guide to get all the data needed to run the task list scripts:**

** Note: all Python scripts contain information regarding the data that the user wishes to use for an specific simulation (i.e. the data vector and data covariance to be loaded, the CosmoCov inifile, etc.) . There are standard values already inserted but they will have to be changed by the user if they wish to run something different.

- xipm: run get_xipms.py. This will run a bash script with the details specified in the inifile that will compute all the necessary blocks in CosmoCov in order to get all the possible zbin pairs. 

- Covmat: run in the terminal script.sh, which that calls ng_cov_kids.ini(to include non-Gaussian contributions) for all blocks (1..465 for KiDS-1000). This script contains a variable N that specifies the number of CPUs available in your computer and may have to be changed (default is 12)

- X matrix : run create_X_matrix.py

- shifted_covariance (to be used in point 4, 5, 6): run shift_covmat.py

