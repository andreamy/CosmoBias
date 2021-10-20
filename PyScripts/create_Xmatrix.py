# !/usr/bin/env python
# coding: utf-8
# Created by Andrea Manrique Yus -  June 2021 - copyrights

from utils.function_library import *

# DATA LOADING FOR TESTS

# number of cosmological parameters
N_parameters = 10

# THEORY VECTOR
mu_cosmolike = np.reshape(np.load("../CosmoCov/covs/xipm/Class3.0_closer_derivative/vector_xipm_kids_Class3.0_closer_derivative.npy"), (-1, 1))

# THEORY COVARIANCE MATRIX
sigma = np.load("../notebooks/CovMatrix_class3.0_closer_derivative.npy")

# DATA VECTOR
kids_file = np.loadtxt("../KiDS1000_cosmis_shear_data_release/data_fits/data_xipm")
data_xipm = np.reshape(kids_file[:, 1], (-1, 1))
Ndim = len(data_xipm)

# DATA COVARIANCE MATRIX
sigma_kids = np.loadtxt("../KiDS1000_cosmis_shear_data_release/data_fits/kids_covariance_matrix")

# Inifile: does not matter that we use a different one later, this is just for running purporses
# as this is the file set to run in the bash script. We copy over it the inifile of interest.
cosmolike_inifile = "../CosmoCov/covs/ini_files/cov_kids_xipm_runs.ini"

default_inifile = "../CosmoCov/covs/simulation_settings/inifile_Class3.0_closer_derivative.txt"
#default_inifile = '../CosmoCov/covs/ini_files/cov_kids_xipm_runs_copy.ini'
copy_default_inifile = subprocess.run(['cp', default_inifile, cosmolike_inifile],
                                 stdout=subprocess.PIPE,
                                 universal_newlines=True)
# The step above ensures that the ini file is always in the format needed for X matrix

simulation = "class3.0_closer_derivative"
xmatrix_filename = 'X_' + simulation + '.npy'
script = './script_for_xipm.sh' #the script to run by Cosmolike

if (os.path.isfile("../notebooks/Xmatrices/"+xmatrix_filename) == False):
    print("Xmatrix does not exists, computing...")
    Xmatrix = create_X_matrix(simulation, mu_cosmolike, cosmolike_inifile,
                              script, N_parameters)
    print("Program executed. Saving Xmatrix to file.")

else:
    print("Xmatrix already exists, loading from file.")
    Xmatrix = np.load("../notebooks/Xmatrices/"+xmatrix_filename)


