# !/usr/bin/env python
# coding: utf-8
# Created by Andrea Manrique Yus -  June 2021 - copyrights

import sys
import importlib
sys.path.append('.')

import utils.function_library 
importlib.reload(utils.function_library)
from utils.function_library import *

# DATA LOADING FOR TESTS

# number of cosmological parameters
N_parameters = 10

# THEORY VECTOR
mu_cosmolike = np.reshape(np.load("CosmoCov/covs/xipm/Class3.0_closer_derivative/vector_xipm_kids_Class3.0_closer_derivative.npy"), (-1, 1))

# THEORY COVARIANCE MATRIX
sigma = np.load("Output/covmats/CovMatrix_class3.0_closer_derivative.npy")

# DATA VECTOR
kids_file = np.loadtxt("KiDS1000_cosmis_shear_data_release/data_fits/data_xipm")
data_xipm = np.reshape(kids_file[:, 1], (-1, 1))
Ndim = len(data_xipm)

# DATA COVARIANCE MATRIX
sigma_kids = np.loadtxt("KiDS1000_cosmis_shear_data_release/data_fits/kids_covariance_matrix")

# Inifile: this is just for running purposes, as this is the file set to run in the bash script.
# We copy over it the inifile of interest.
bash_inifile = "CosmoCov/covs/ini_files/cov_kids_xipm_runs.ini"
inifile_to_copy = "CosmoCov/covs/simulation_settings/inifile_Class3.0_closer_derivative.txt"

#This step allows to use different inifiles (inifile_to_copy),  without having to change
#  the bash script whenever we switch.
copy_inifile = subprocess.run(['cp', inifile_to_copy, bash_inifile],
                                 stdout=subprocess.PIPE,
                                 universal_newlines=True)


simulation = "class3.0_closer_derivative"
xmatrix_filename = 'X_' + simulation + '.npy'
script = './script_for_xipm.sh' #the script to run by Cosmolike

if (os.path.isfile("Output/Xmatrices/"+xmatrix_filename) == False):
    print("Xmatrix does not exists, computing...")
    Xmatrix = create_X_matrix(simulation, mu_cosmolike, bash_inifile,
                              script, N_parameters)
    print("Program executed. Saving Xmatrix to file.")

else:
    print("Xmatrix already exists, loading from file.")
    Xmatrix = np.load("Output/Xmatrices/"+xmatrix_filename)


