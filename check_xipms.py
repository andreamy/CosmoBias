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
setup_simulation(simulation, non_gaussian=False, euklids=False) #creates directories and saves inifile to simulation_settings as well

setup_and_run_CosmoLike(simulation,euklids=False,non_gaussian=False) #sets up & runs cosmolike

create_theory_vector(xipm_directory, simulation)
mu_cosmolike = np.reshape(np.load("../CosmoCov/covs/xipm/{}/"
                                  "vector_xipm_kids_{}.npy".format(simulation,simulation)), (-1,1))

chi_square_test(sigma_kids, mu_cosmolike, data_xipm)