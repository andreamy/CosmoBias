import sys
import importlib
sys.path.append('../pycharmProject1/')

import utils.function_library
importlib.reload(utils.function_library)
from utils.function_library import *

#DATA COVARIANCE MATRIX
sigma_kids = np.loadtxt("../KiDS1000_cosmis_shear_data_release/data_fits/kids_covariance_matrix")

#DATA VECTOR
kids_file = np.loadtxt("../KiDS1000_cosmis_shear_data_release/data_fits/data_xipm")
data_xipm = np.reshape(kids_file[:,1], (-1,1))

#change output and simulation variables in inifile, then setup_simulation,
#then run script from terminal, script_for_xipm.sh or script.sh
simulation = "K1000BestFit"
xipm_directory = "../CosmoCov/covs/xipm/" + simulation +"/"
script = './script_for_xipm.sh' #the script to run by Cosmolike
#setup_simulation(simulation, non_gaussian=False, euklids=False) #creates directories and saves inifile to simulation_settings as well

# Inifile: this is just for running purposes, as this is the file set to run in the bash script.
# We copy over it the inifile of interest.
bash_inifile = "../CosmoCov/covs/ini_files/cov_kids_xipm_runs.ini"
inifile_to_copy = "../CosmoCov/covs/simulation_settings/inifile_Class3.0_closer_derivative.txt"

#This step allows to use different inifiles (inifile_to_copy),  without having to change
#  the bash script whenever we switch.
copy_inifile = subprocess.run(['cp', inifile_to_copy, bash_inifile],
                                 stdout=subprocess.PIPE,
                                 universal_newlines=True)

setup_and_run_CosmoLike(simulation, bash_inifile, script) #sets up & runs cosmolike

#Creates a python array from the output files containing xipm's, that can be later loaded by any file
create_theory_vector(xipm_directory, simulation)
mu_cosmolike = np.reshape(np.load("../CosmoCov/covs/xipm/{}/"
                                  "vector_xipm_kids_{}.npy".format(simulation,simulation)), (-1,1))

chi_square_test(sigma_kids, mu_cosmolike, data_xipm)
