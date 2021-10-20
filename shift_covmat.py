#!/usr/bin/env python
# coding: utf-8
# Created by Andrea Manrique Yus -  June 2021 - copyrights

from utils.function_library import *

# DATA LOADING FOR TESTS

# number of cosmological parameters
N_parameters = 10

# THEORY VECTOR
mu_cosmolike = np.reshape(np.load("../CosmoCov/covs/xipm/Class3.0_closer_derivative/vector_xipm_kids_Class3.0_closer_derivative.npy"), (-1, 1))

# DATA VECTOR
kids_file = np.loadtxt("../KiDS1000_cosmis_shear_data_release/data_fits/data_xipm")
data_xipm = np.reshape(kids_file[:, 1], (-1, 1))
Ndim = len(data_xipm)

# DATA COVARIANCE MATRIX
sigma_kids = np.loadtxt("../KiDS1000_cosmis_shear_data_release/data_fits/kids_covariance_matrix")

# Inifile
cosmolike_inifile = "../CosmoCov/covs/ini_files/cov_kids_xipm_runs.ini"

# cosmological parameters
cosmo_params = get_cosmological_parameters(cosmolike_inifile, N_parameters)
params_without_bestfit = dict.fromkeys(['Omega_v', 'w0', 'wa', 'IA'])
list(map(cosmo_params.pop, params_without_bestfit))

#Xmatrix - computed previously with create_Xmatrix.py
Xmatrix = np.load("../notebooks/X_class3.0_closer_derivative.npy")

shift_matrix_up = shift_covmat_elements('class3.0_closer_derivative', Xmatrix.T,
                                        sigma_kids, mu_cosmolike, data_xipm, cosmo_params, 'up')
shift_matrix_down = shift_covmat_elements('class3.0_closer_derivative', Xmatrix.T,
                                          sigma_kids, mu_cosmolike, data_xipm, cosmo_params, 'down')