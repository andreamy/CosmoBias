#!/usr/bin/env python
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

# Inifile
cosmolike_inifile = "../CosmoCov/covs/ini_files/cov_kids_xipm_runs.ini"

# cosmological parameters
cosmo_params = get_cosmological_parameters(cosmolike_inifile, N_parameters)
params_without_bestfit = dict.fromkeys(['Omega_v', 'w0', 'wa', 'IA'])
list(map(cosmo_params.pop, params_without_bestfit))

#Xmatrix - computed previously with create_Xmatrix.py
Xmatrix = np.load("../notebooks/X_class3.0_closer_derivative.npy")

#1 sigma errors
omega_m_symerr = 0.5*(0.065+0.033)
sigma8_symerr = 0.5*(0.084+0.107)
n_s_symerr = 0.5*(0.093+0.049)
omega_b_symerr = 0.5*(0.005+0.000)#Omega_b*h**2 = 0.046*h**2
h_symerr = 0.5*(0.110+0.001)
A_ia_symerr = 0.5 *(0.321 + 0.374)
denom_list = [omega_m_symerr, sigma8_symerr, n_s_symerr, omega_b_symerr, h_symerr, A_ia_symerr]

print("linear theta with cosmolike matrix")
lin_theta = linear_theta(Xmatrix.T, sigma, mu_cosmolike, data_xipm, cosmo_params)
print(lin_theta)
#print(round(lin_theta[0,0],3), 'is', round(sigmas_away(mean_value,lin_theta,denom_list[4])[0,0],3),
#      'sigma away from mean value', mean_value)

print("\n with sigma kids:")
kids_lin_theta = linear_theta(Xmatrix.T, sigma_kids, mu_cosmolike, data_xipm,cosmo_params)
print(kids_lin_theta)

