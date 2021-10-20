import os
import sys
import natsort
import fileinput
import subprocess
import shutil
import matplotlib
import numpy as np

import matplotlib.pyplot as plt
from numpy.linalg import inv
from scipy.stats import chi2
from scipy.stats import multivariate_normal
from matplotlib import colors


def structure_matrix(A, i, j):
    """
    See matrix cookbook section on structure matrix.
    Case for symmetric matrices, otherwise the structure matrix is the
    single entry matrix: all elements are zero except for ij, equal to 1
    """
    Jij = np.zeros(A.shape)
    Jji = np.zeros(A.shape)
    Jij[i][j] = 1
    Jji[j][i] = 1
    Sij = Jij + Jji - np.matmul(Jij, Jij)
    return Sij


def partial_derivative_old(mu, mu_shifted, delta_theta, param, i, j):
    # load data vector and covariance matrix (sigma)
    data = np.loadtxt("../KiDS1000_cosmis_shear_data_release/data_fits/data_xipm")
    d = np.reshape(data[:, 1], (-1, 1))
    sigma = np.load("../notebooks/CovMatrix_kids.npy")

    # Compute the equation by parts
    aux1 = forward_derivative_old(mu_shifted, mu)
    X = np.reshape(aux1, (-1, 1))
    aux2 = np.matmul(X.T, inv(sigma))
    bracket = 1 / (np.matmul(aux2, X))
    dsigma = structure_matrix(sigma, i, j)
    result = bracket * X.T @ inv(sigma) @ dsigma @ (inv(sigma) @ X * bracket * X.T @ inv(sigma) -
                                                    inv(sigma)) @ (d - np.reshape(mu, (-1, 1)) + X * param)
    return result


def linear_theta_old(mu, mu_shifted, sigma, d, param, inverted_sigma):
    # sigma must be the shifted one here
    # mu and mu_shifted must have been reshaped to 2d column arrays (N, 1)
    X = np.reshape(forward_derivative_old(mu_shifted, mu), (-1, 1))

    if inverted_sigma:
        theta_lin = 1 / (X.T @ sigma @ X) @ X.T @ sigma @ (d - mu + X * param)

    else:
        theta_lin = 1 / (X.T @ inv(sigma) @ X) @ X.T @ inv(sigma) @ (d - mu + X * param)

    return theta_lin


def forward_derivative_old(mu_shifted, mu):
    N = len(mu)
    X = np.zeros((N,))
    h = 0.01  # had 0.0001 but cosmolike (or other theory codes) do not have such resolution

    for i in range(N): X[i] = (mu_shifted[i] - mu[i]) / h
    return X


def pos_neg_colormap(matrix):
    # Mask an array where smaller/greater than a given value.
    pos_matrix = np.ma.masked_less(matrix, 0)
    neg_matrix = np.ma.masked_greater(matrix, 0)

    fig, (ax1, ax2, ax3) = plt.subplots(figsize=(26, 6), ncols=3)

    # plot just the positive data and save the
    # color "mappable" object returned by ax1.imshow
    pos = ax1.imshow(pos_matrix, cmap='Blues', interpolation='none')

    # add the colorbar using the figure's method,
    # telling which mappable we're talking about and
    # which axes object it should be near
    fig.colorbar(pos, ax=ax1)

    # repeat everything above for the negative data
    neg = ax2.imshow(neg_matrix, cmap='Reds_r', interpolation='none')
    fig.colorbar(neg, ax=ax2)

    # Flatten matrix to find vmin and vmax
    flat_mat = matrix.flatten()

    # Plot both positive and negative values between +/- 1.2
    # pos_neg_clipped = ax3.imshow(Z, cmap='RdBu', vmin=-1.2, vmax=1.2,interpolation='none')
    pos_neg = ax3.imshow(matrix, vmin=np.min(flat_mat), vmax=np.max(flat_mat), cmap='RdBu', interpolation='none')
    # Add minorticks on the colorbar to make it easy to read the
    # values off the colorbar.
    # cbar = fig.colorbar(pos_neg_clipped, ax=ax3, extend='both')
    cbar = fig.colorbar(pos_neg, ax=ax3, extend='both')
    cbar.minorticks_on()
    plt.show()


def largest_values(matrix, number_of_values):
    # Make array starting from 0 and ending at the dimensions of the matrix, multiplied (270x270)
    index_matrix = np.arange(0, len(matrix.flatten()), 1, dtype=np.int)

    # Find the location of the 10 highest elements of the matrix (flattened out)
    top_indices = matrix.flatten().argsort()[-number_of_values:]

    # Reshape the index matrix to the dimensions of the input matrix (270x270). Then, put True where the index matrix
    # matches the indices of the top 10 ids, and false everywhere else (np.isin looks whether top10ids are in
    # index_matrix.reshape(omega_matrix.shape)). Then find the location where this matrix is True.
    # This will be the location of the top 10 elements of the input matrix.
    top_indices_matrix = np.where(np.isin(index_matrix.reshape(matrix.shape), top_indices) == True)
    # print(top10matrix_ids)

    return top_indices_matrix


def reduce_dimension(covmat, theory_vector, data_vector):
    """
    Loop: for every xi_minus vector per z-bin combination, keep first three indices with
    np.s_[0:3] (to later discard first 3 theta bins). Finally reshape data and theory vectors
    to (225,1) just like above with the covariance matrix (covmat).
    """
    xim_indices = np.arange(135, 270)
    xim_indices_per_zbin = np.split(xim_indices, 15)
    indices_to_discard = np.array(())

    for xim_zbin in xim_indices_per_zbin:
        indices_to_discard = np.concatenate((indices_to_discard, xim_zbin[np.s_[0:3]]))

    indices_to_discard = indices_to_discard.astype(int)

    reshaped_covmat = covmat
    reshaped_covmat = np.delete(reshaped_covmat, indices_to_discard, axis=0)
    reshaped_covmat = np.delete(reshaped_covmat, indices_to_discard, axis=1)
    # print(reshaped_covmat.shape)

    reshaped_data = np.delete(data_vector, indices_to_discard, axis=0)
    reshaped_theory = np.delete(theory_vector, indices_to_discard, axis=0)

    return reshaped_covmat, reshaped_theory, reshaped_data


def reduce_vector_dimension(vector):
    """
    Loop: for every xi_minus vector per z-bin combination, keep first three indices with
    np.s_[0:3] (to later discard first 3 theta bins). Finally reshape vector (data or theory)
    to (225,1) just like above with the covariance matrix (covmat).
    """
    xim_indices = np.arange(135, 270)
    xim_indices_per_zbin = np.split(xim_indices, 15)
    indices_to_discard = np.array(())

    for xim_zbin in xim_indices_per_zbin:
        indices_to_discard = np.concatenate((indices_to_discard, xim_zbin[np.s_[0:3]]))

    indices_to_discard = indices_to_discard.astype(int)
    reshaped_vector = np.delete(vector, indices_to_discard, axis=0)

    return reshaped_vector


def reduce_covmat_dimension(covmat):
    """
    Loop: for every xi_minus vector per z-bin combination, keep first three indices with
    np.s_[0:3] (to later discard first 3 theta bins). Finally reshape covmat to (225,225).
    """
    xim_indices = np.arange(135, 270)
    xim_indices_per_zbin = np.split(xim_indices, 15)
    indices_to_discard = np.array(())

    for xim_zbin in xim_indices_per_zbin:
        indices_to_discard = np.concatenate((indices_to_discard, xim_zbin[np.s_[0:3]]))

    indices_to_discard = indices_to_discard.astype(int)

    reshaped_covmat = covmat
    reshaped_covmat = np.delete(reshaped_covmat, indices_to_discard, axis=0)
    reshaped_covmat = np.delete(reshaped_covmat, indices_to_discard, axis=1)
    # print(reshaped_covmat.shape)

    return reshaped_covmat


def get_cosmological_parameters(ini_filename, N_parameters):
    k = 0
    cosmo_parameters = {}
    f = open(ini_filename, "r")

    while (k != N_parameters):
        line = f.readline()
        if not line.startswith('#'):
            parameter_name = line.split(' : ')[0]
            value = line.split(' : ')[1][:-1]
            cosmo_parameters[parameter_name] = float(value)
            k = k + 1

    f.close()

    return cosmo_parameters


# make a function of editing the file by replacing, useful for running cosmolike anytime.
# can be called for every thing we want to change
def edit_ini_file(ini_filename, string_to_replace, new_string):
    # Read in the file, replace the target string and write the file out again
    with open(ini_filename, 'r') as file:
        filedata = file.read()

    filedata = filedata.replace(string_to_replace, new_string)

    with open(ini_filename, 'w') as file:
        file.write(filedata)


def create_theory_vector(xipm_path, simulation_name):
    xip_vector = np.array(())
    xim_vector = np.array(())
    listdir = []

    for file in os.listdir(xipm_path):
        if file.startswith("xipm_"):
            listdir.append(file)

    listdir = natsort.natsorted(listdir)

    for file in listdir:
        data = np.loadtxt(xipm_path + '/' + file)
        xip = data[:, 1]
        xim = data[:, 2]
        xip_vector = np.append(xip_vector, xip)
        xim_vector = np.append(xim_vector, xim)

    xipm_vector = np.reshape(np.hstack((xip_vector, xim_vector)), (-1, 1))
    np.save(xipm_path + '/vector_xipm_kids_' + simulation_name, xipm_vector)

    return xipm_vector


def create_X_matrix(simulation_name, fiducial_mu, ini_filename, script, N_parameters=None,
                    Ndim_parameters=None, kids=True, reduced_dimension=False):
    if Ndim_parameters:
        cosmo_parameters = Ndim_parameters #must be a dictionary with keys and values

    else:
        cosmo_parameters = get_cosmological_parameters(ini_filename, N_parameters)
        if kids:
            params_without_bestfit = dict.fromkeys(['Omega_v', 'w0', 'wa', 'IA'])
            #'Omega_m', 'sigma_8', 'n_spec', 'omb', 'h0'
            list(map(cosmo_parameters.pop, params_without_bestfit))
            # choose only the ones that have a best-fit with 1sigma error

    k = 0
    N_dimension = len(fiducial_mu)
    N_parameters = len(cosmo_parameters)

    if reduced_dimension == True:
        N_dimension = 225
        X_matrix = np.zeros((N_parameters, N_dimension))

    else:
        X_matrix = np.zeros((N_parameters, N_dimension))

    for i in cosmo_parameters:
        print('Running simulation for parameter :', i)
        h = 2/100*cosmo_parameters[i]
        #h = K1000_standard_dev_symmetrized()[k]  # for taylor test, we want taylor_theta=1sigma
        if i == 'wa': h = 0.005  # otherwise 2% of zero is zero

        new_value = round((cosmo_parameters[i] + h), 4)

        edit_ini_file(ini_filename, i + ' : ' + str(cosmo_parameters[i]), i + ' : ' + str(new_value))
        edit_ini_file(ini_filename, 'output/default', 'output/' + simulation_name + '_derivative_shift_' + i)
        edit_ini_file(ini_filename, 'simulation : default',
                      'simulation : ' + simulation_name + '_derivative_shift_' + i)

        X_matrix[k] = forward_derivative(simulation_name, i, fiducial_mu, h, N_dimension, ini_filename, script)

        edit_ini_file(ini_filename, i + ' : ' + str(new_value), i + ' : ' + str(cosmo_parameters[i]))
        edit_ini_file(ini_filename, 'output/' + simulation_name + '_derivative_shift_' + i, 'output/default')
        edit_ini_file(ini_filename, 'simulation : ' + simulation_name + '_derivative_shift_' + i,
                      'simulation : default')
        k = k + 1
        np.save("../notebooks/Xmatrices/X_" + simulation_name, X_matrix)

    return X_matrix.T


def forward_derivative(simulation_name, parameter, fiducial_mu, h, N_dimension, ini_file, script, reduced_dimension=False):
    # Parameter must be a string
    # h = 0.01
    derivative = np.zeros((N_dimension,))

    if type(parameter) != str:
        sys.exit(
            'Error: <parameter> is not a string. Please insert the name of a cosmological parameter in string format.')

    mu_shifted = setup_and_run_CosmoLike(simulation_name + '_derivative_shift_' + parameter, ini_file, script )

    if reduced_dimension == True:
        fiducial_mu = reduce_vector_dimension(fiducial_mu)
        mu_shifted = reduce_vector_dimension(mu_shifted)

    for i in range(N_dimension):
        derivative[i] = (mu_shifted[i] - fiducial_mu[i]) / h

    return derivative


# do the for loop out, in case that I need to use the derivative for something else
# do another function called X that returns the whole thing together. Much better
def setup_and_run_CosmoLike(simulation_name, ini_file, script):

    # save inifile to know later which settings were used.
    target = '../CosmoCov/covs/simulation_settings/inifile_' + simulation_name + '.txt'
    shutil.copyfile(ini_file, target)

    xipm_path = os.path.join("../CosmoCov/covs/xipm/", simulation_name)
    covmat_path = os.path.join("../CosmoCov/covs/output/", simulation_name)

    if (os.path.exists(xipm_path) == False):
        os.mkdir(xipm_path)
        print("Directory '%s' created " % (xipm_path))
    else:
        print("Directory '%s' already exists " % (xipm_path))

    if (os.path.exists(covmat_path) == False):
        os.mkdir(covmat_path)
        print("Directory '%s' created " % (covmat_path))
    else:
        print("Directory '%s'  already exists " % (covmat_path))

    # run_covs = subprocess.run([script, 'chmod', '+x', 'script_for_xipm.sh'], ...)
    run_covs = subprocess.run([script],
                              cwd='../CosmoCov/covs/',
                              stdout=subprocess.PIPE,
                              universal_newlines=True)

    # If the program has run correctly, it will print the last part of the output, saying PROGRAM EXECUTED.
    print(run_covs.stdout[-54:])

    return create_theory_vector(xipm_path, simulation_name)


def setup_and_run_CosmoLike_old(simulation_name, euklids, non_gaussian=False):
    if euklids == True:
        ini_file = 'cov_euklids.ini'
        script = './script_for_xipm_euklids.sh'
    else:
        if non_gaussian == True:
            ini_file = 'ng_cov_kids.ini'
            script = './script.sh'
        else:
            ini_file = 'cov_kids_xipm_runs.ini'
            script = './script_for_xipm.sh'

    # save inifile to know later which settings were used.
    original = '../CosmoCov/covs/ini_files/' + ini_file
    target = '../CosmoCov/covs/simulation_settings/inifile_' + simulation_name + '.txt'
    shutil.copyfile(original, target)

    xipm_path = os.path.join("../CosmoCov/covs/xipm/", simulation_name)
    covmat_path = os.path.join("../CosmoCov/covs/output/", simulation_name)

    if (os.path.exists(xipm_path) == False):
        os.mkdir(xipm_path)
        print("Directory '%s' created " % (xipm_path))
    else:
        print("Directory '%s' already exists " % (xipm_path))

    if (os.path.exists(covmat_path) == False):
        os.mkdir(covmat_path)
        print("Directory '%s' created " % (covmat_path))
    else:
        print("Directory '%s'  already exists " % (covmat_path))

    # run_covs = subprocess.run([script, 'chmod', '+x', 'script_for_xipm.sh'], ...)
    run_covs = subprocess.run([script],
                              cwd='../CosmoCov/covs/',
                              stdout=subprocess.PIPE,
                              universal_newlines=True)

    # If the program has run correctly, it will print the last part of the output, saying PROGRAM EXECUTED.
    print(run_covs.stdout[-54:])

    return create_theory_vector(xipm_path, simulation_name)


def setup_simulation_old(simulation_name, non_gaussian=False, euklids=False):
    # this can be used to create all directories and setup the file without
    # runnning Cosmolike. That can be done separately

    if euklids == True:
        ini_file = 'cov_euklids.ini'
        script = './euklids_script.sh'
    else:
        if non_gaussian == True:
            ini_file = 'ng_cov_kids.ini'
            script = './script.sh'
        else:
            ini_file = 'cov_kids_xipm_runs.ini'
            script = './script_for_xipm.sh'

    # save inifile to know later which settings were used.
    original = '../CosmoCov/covs/ini_files/' + ini_file
    target = '../CosmoCov/covs/simulation_settings/inifile_' + simulation_name + '.txt'
    shutil.copyfile(original, target)

    xipm_path = os.path.join("../CosmoCov/covs/xipm/", simulation_name)
    covmat_path = os.path.join("../CosmoCov/covs/output/", simulation_name)

    if (os.path.exists(xipm_path) == False):
        os.mkdir(xipm_path)
        print("Directory '%s' created " % (xipm_path))
    else:
        print("Directory '%s' already exists " % (xipm_path))

    if (os.path.exists(covmat_path) == False):
        os.mkdir(covmat_path)
        print("Directory '%s' created " % (covmat_path))
    else:
        print("Directory '%s'  already exists " % (covmat_path))

    return script

def setup_simulation(simulation_name, ini_file):
    # save inifile to know later which settings were used.
    target = '../CosmoCov/covs/simulation_settings/inifile_' + simulation_name + '.txt'
    shutil.copyfile(ini_file, target)

    #create xipm directory
    xipm_path = os.path.join("../CosmoCov/covs/xipm/", simulation_name)
    covmat_path = os.path.join("../CosmoCov/covs/output/", simulation_name)

    if (os.path.exists(xipm_path) == False):
        os.mkdir(xipm_path)
        print("Directory '%s' created " % (xipm_path))
    else:
        print("Directory '%s' already exists " % (xipm_path))

    if (os.path.exists(covmat_path) == False):
        os.mkdir(covmat_path)
        print("Directory '%s' created " % (covmat_path))
    else:
        print("Directory '%s'  already exists " % (covmat_path))

def run_CosmoLike(script):
    # run_covs = subprocess.run([script, 'chmod', '+x', 'script_for_xipm.sh'], ...)
    run_covs = subprocess.run([script],
                              cwd='../CosmoCov/covs/',
                              stdout=subprocess.PIPE,
                              universal_newlines=True)

    # If the program has run correctly, it will print the last part of the output, saying PROGRAM EXECUTED.
    print(run_covs.stdout[-54:])


def chi_square_test(data_covariance_matrix, theory_vector, data_vector):
    """ REDUCE COVMAT AND VECTORS DIMENSION FOR CHI-SQUARED TEST
    Following the KiDS-1000 cosmic shear paper (asgari 2020) where they discard first three values of
    xi_m for the chi-square test"""

    reshaped_covmat, mu_reshaped, data_reshaped = reduce_dimension(data_covariance_matrix,
                                                                   theory_vector, data_vector)

    chi = (mu_reshaped - data_reshaped).T @ inv(reshaped_covmat) @ (mu_reshaped - data_reshaped)

    # dof = degrees of freedom (see Asgari 2020 page 9 and table 3)
    Ndim = 270 - (3 * 15)  # first 3 values of each xipm are discarded due
                            # to their sensitivity to small physical scales
    Nparams = 4.5
    dof = Ndim - Nparams

    print('Minimum chi-square to reject null hypothesis: %.3f\n' %(chi2.ppf(0.95, dof)))
    print('chi-squared = %.3f, p-value = %r' % (chi, 1 - chi2.cdf(chi, dof)))
    print('chi-squared reduced from kids is 260, and p value is 0.034\n')

    return chi[0,0]


def linear_theta(X, sigma, mu, data, cosmo_parameters, inverted_sigma=False):
    parameter_list = np.reshape(list(cosmo_parameters.values()), (-1, 1))

    if inverted_sigma == True:
        theta_lin = 1 / (X.T @ sigma @ X) @ X.T @ sigma @ (data - mu + X @ parameter_list)

    else:
        theta_lin = 1/(X.T @ inv(sigma) @ X) @ X.T @ inv(sigma) @ (data - mu + X @ parameter_list)

    return theta_lin


def sigmas_away(mean_value, value, sigma):
    return abs((value - mean_value)/sigma)


def shift_covmat_elements(simulation, X, sigma, mu, data, cosmo_parameters, type_of_shift='up', inverted_sigma=False):
    k = 0
    N_dim = sigma.shape[0]
    #N_dim=10
    N_parameters = len(cosmo_parameters)
    parameter_shifts = np.zeros(shape=(N_dim, N_dim, N_parameters))

    if inverted_sigma == True:
        sigma = inv(sigma)

    if type_of_shift == 'up':
        shift = 1.05
    elif type_of_shift == 'down':
        shift = 0.95
    else:
        sys.exit('Error: type_of_shift must be \'up\'(default) or \'down\'')

    initial_linear_parameters = linear_theta(X, sigma, mu, data, cosmo_parameters, inverted_sigma)

    # ATTENTION, 36315 ELEMENTS TO COMPUTE, TAKES LONG TIME AND ALL CPU's
    for i in range(N_dim):
        for j in range(i, N_dim):
            element = sigma[i, j]
            sigma[i, j] = element * shift
            parameter_shifts[i, j] = linear_theta(X, sigma, mu, data, cosmo_parameters, inverted_sigma)[:, 0] \
                                     - initial_linear_parameters[:, 0]
            sigma[i, j] = element
            k = k + 1
            if k % 1000 == 0: print(k)

    parameter_shifts_transposed = np.transpose(parameter_shifts.copy(), (1, 0, 2))

    for i in range(N_parameters):
        np.fill_diagonal(parameter_shifts_transposed[:, :, i], 0)

    parameter_shifts_sym = parameter_shifts_transposed + parameter_shifts

    np.save("../notebooks/parameter_shifts_" + str(shift)+"_"+simulation, parameter_shifts_sym)
    return parameter_shifts_sym


# NOT NEEDED IN THE END
def shift_covmat_free_parameters(free_parameter, X, sigma_shifted, mu, data, cosmo_parameters):
    k = 0
    N_dim = sigma_shifted.shape[0]
    N_parameters = len(cosmo_parameters)
    parameter_shifts = np.zeros(shape=(N_dim, N_dim, N_parameters))

    initial_linear_parameters = linear_theta(X, sigma_shifted, mu, data, cosmo_parameters)

    # ATTENTION, 36315 ELEMENTS TO COMPUTE, TAKES LONG TIME AND ALL CPU's
    for i in range(N_dim):
        for j in range(i, N_dim):
            parameter_shifts[i, j] = linear_theta(X, sigma_shifted, mu, data, cosmo_parameters)[:, 0] \
                                     - initial_linear_parameters[:, 0]
            k = k + 1
            if k % 1000 == 0: print(k)

    parameter_shifts_transposed = np.transpose(parameter_shifts.copy(), (1, 0, 2))

    for i in range(N_parameters):
        np.fill_diagonal(parameter_shifts_transposed[:, :, i], 0)

    parameter_shifts_sym = parameter_shifts_transposed + parameter_shifts

    np.save("parameter_shifts_modified_" + free_parameter, parameter_shifts_sym)
    return parameter_shifts_sym


def feature_scaling_norm(vector, vector_max, vector_min):
    a = -1
    b = +1
    return a + (((vector - vector_min) * (b - a)) / (vector_max - vector_min))


def precision_recommendation(X, sigma, mu, data, cosmo_parameters, survey='K1000'):
    """
        Returns a matrix with percentage errors per element to avoid biasing Euclid precision,
        by showing the maximum percent deviation a covariance element could have while keeping
        the parameters at 1sigma.

        The symmetrized errors come from Asgari(2020), cosmic shear, table A.2, 2PCF, MAP values
        The symmetrized error is 1 sigma
    """
    N_dim = sigma.shape[0]
    p = 0.1
    # one sigma errors KiDS1000
    omega_m_symerr = 0.5 * (0.065 + 0.033) * p
    sigma8_symerr = 0.5 * (0.084 + 0.107) * p
    n_s_symerr = 0.5 * (0.093 + 0.049) * p
    omega_b_symerr = (0.5 * (0.005 + 0.000)) / (cosmo_parameters['h0']) ** 2 * p
    # A_ia_symerr = 0.5*(0.321+0.374)*p
    h_symerr = 0.5 * (0.110 + 0.001) * p

    if survey == 'K1000':
        one_sigma_errors = {"Omega_m": omega_m_symerr, "sigma_8": sigma8_symerr, "n_spec": n_s_symerr,
                            "omb": omega_b_symerr, "h0": h_symerr}

    # one sigma errors Euclid precision, only weak lensing LambdaCDM, non-flat (A&A 642, A191 (2020))
    if survey == 'EUCLID':
        one_sigma_errors = {"Omega_m": 0.012, "Omega_v": 0.032, "h0": 0.2,
                            "n_spec": 0.032, "omb": 0.26, "sigma_8": 0.0074}

    N_parameters = len(one_sigma_errors)
    parameters = {"Omega_m": 0, "sigma_8": 1, "n_spec": 2, "omb": 3, "h0": 4}  # mimic one_sigma_errors keys
    recommended_precision_matrix = np.zeros(shape=(N_parameters, N_dim, N_dim))
    initial_linear_parameters = linear_theta(X, sigma, mu, data, cosmo_parameters)

    for i in range(N_dim):
        for j in range(i, N_dim):
            shift = 0.05  # 1.0 equals 100% shift of a matrix element
            epsilon = shift / 10
            reduce_shift = False

            while (True):
                element = sigma[i, j]
                sigma[i, j] = element * (1 + shift)
                parameter_shifts = abs(linear_theta(X, sigma, mu, data, cosmo_parameters)[:, 0] \
                                       - initial_linear_parameters[:, 0])
                sigma[i, j] = element

                # to ensure that cosmo_parameters has same parameters as one_sigma_errors
                shifts_dict = dict(zip(cosmo_parameters.keys(), parameter_shifts))

                if (one_sigma_errors.keys() == shifts_dict.keys()) == False:
                    different_keys = shifts_dict.keys() - one_sigma_errors.keys()  # returns the different keys
                    list(map(shifts_dict.pop, different_keys))

                # print("shifts = ",shifts_dict.values(), "\n 1sigmas = ", one_sigma_errors.values())

                for param, s, o in zip(parameters.items(), list(shifts_dict.values()), list(one_sigma_errors.values())):
                    if s < o and recommended_precision_matrix[param[1], i, j] == 0:
                        recommended_precision_matrix[param[1], i, j] = shift * 100
                        print("{}% error at ({},{}) does not bias {} at {}'s "
                              "precision".format(round(shift * 100, 3), i, j, param[0], survey))
                        # instead of param, parameter.keys()[param]
                    elif s < o:
                        pass

                    else:
                        if round(shift, 3) <= epsilon:
                            print("\n{}% error at ({},{}) still "
                                  "biases {} at {}'s precision.\n".format(epsilon * 100, i, j, param[0], survey))
                        else:
                            reduce_shift = True

                if (reduce_shift == True) and (round(shift, 3) > epsilon):
                    shift = shift - epsilon
                    # print("Reducing shift to ", round(shift, 3))
                    reduce_shift = False
                else:
                    break

    np.save("recommended_covariance_{}_precision.npy".format(survey), recommended_precision_matrix)
    return recommended_precision_matrix


def redshift_Gauss_distribution(z, z_0, std_dev, A):
    return np.round_(A * np.exp(-0.5 * (z - z_0) ** 2 / std_dev ** 2) / np.sqrt(2 * np.pi * std_dev ** 2), 12)


def normalize_area(f, x):
    area = np.trapz(f, x)
    return f / area


def fisher_matrix(X, sigma):
    return  X.T @ inv(sigma) @ X


def precision_matrix(sigma, X):
    # top-hat priors from K1000 (methodology paper)
    h = [0.64, 0.82]
    omega_ch2 = [0.051, 0.255]
    omega_bh2 = [0.019, 0.026]
    S8 = [0.1, 1.3]
    n_s = [0.84, 1.1]
    A_ia = [-6, 6]

    # Sample uniform prior distributions with previous info and convert some parameters to cosmolike parameters
    h_uniform = np.random.uniform(h[0], h[1], 100000)

    omch2_uniform = np.random.uniform(omega_ch2[0], omega_ch2[1], 100000)
    ombh2_uniform = np.random.uniform(omega_bh2[0], omega_bh2[1], 100000)
    Omegab_uniform = ombh2_uniform / h_uniform ** 2

    Omegam_uniform = (omch2_uniform + ombh2_uniform) / (h_uniform ** 2)

    S8_uniform = np.random.uniform(S8[0], S8[1], 100000)
    sigma8_uniform = S8_uniform / np.sqrt(Omegam_uniform / 0.3)

    Aia_uniform = np.random.uniform(A_ia[0], A_ia[1], 100000)

    Ns_uniform = np.random.uniform(n_s[0], n_s[1], 100000)

    # stack prior distributions together
    parameter_priors = np.vstack((Omegam_uniform, sigma8_uniform, Ns_uniform, Omegab_uniform, h_uniform, Aia_uniform))

    # Compute covariance matrix of priors
    prior_covariances = np.cov(parameter_priors)

    # Add prior information to Fisher matrix (see notes from elena - to be written on LateX)
    precision_matrix = inv(inv(prior_covariances) + fisher_matrix(X.T, sigma))

    # To get just the sigmas
    # fisher_sigmas = np.sqrt(np.diag(precision_matrix))

    return precision_matrix


def standardDeviation_Fisher(X, sigma):
    fisher_matrix = X.T @ inv(sigma) @ X
    return np.sqrt(np.diag(fisher_matrix))


def sigmas_away(mean_value, value, sigma):
    return abs((value - mean_value)/sigma)


def gaussian_2D(x, y, x_mean, y_mean, sigma_x, sigma_y, norm=True):
    z = np.exp(-( (x-x_mean)**2/(2*sigma_x**2) + (y-y_mean)**2/(2*sigma_y**2) ))
    if norm:
        # Normalize
        x_grid = np.linspace(-5*sigma_x, 5*sigma_x, 3000)
        y_grid = np.linspace(-5*sigma_y, 5*sigma_y, 3000)
        x_grid, y_grid = np.meshgrid(x_grid, y_grid)
        z_grid  = np.exp(- ( (x_grid-x_mean)**2/(2*sigma_x**2)+(y_grid-y_mean)**2/(2*sigma_y**2) ) )
        Integral = np.trapz(np.trapz(z_grid, x_grid[0,:]), y_grid[:,0])
        return z/Integral
    else:
        return z


def multivariate_gaussian(pos, mu, Sigma):
    """Return the multivariate Gaussian distribution on array pos.

    pos is an array constructed by packing the meshed arrays of variables
    x_1, x_2, x_3, ..., x_k into its _last_ dimension.

    """
    n = mu.shape[0]
    Sigma_det = np.linalg.det(Sigma)
    Sigma_inv = np.linalg.inv(Sigma)
    N = np.sqrt((2*np.pi)**n * Sigma_det)
    # This einsum call calculates (x-mu)T.Sigma-1.(x-mu) in a vectorized
    # way across all the input variables.
    fac = np.einsum('...k,kl,...l->...', pos-mu, Sigma_inv, pos-mu)

    return np.exp(-fac / 2) / N


def plot_matrix_colormap(matrix, colormap=None, min_value=None, max_value=None, title="Default title",
                         ylabel=None, pngfilename="changefilename"):
    fig, (ax1) = plt.subplots(figsize=(10, 7), ncols=1)
    # fig.suptitle(title, fontsize=20)
    ax1.set_title(title, fontdict={'fontsize': 15, 'fontweight': 'medium'})
    plt.ylabel(ylabel)
    axes = plt.gca()
    axes.yaxis.label.set_size(20)
    # define your scale, with white at zero
    discrete_cmap = matplotlib.cm.get_cmap("Reds_r", 10)
    figplot = ax1.imshow(matrix, vmax=max_value, vmin=min_value, cmap=discrete_cmap,
                         interpolation='none')
    fig.colorbar(figplot, ax=ax1)
    plt.savefig(pngfilename + ".pdf", bbox_inches='tight')
    plt.margins(0, 0)
    plt.show()
    return ax1


def plot_matrix_colormap_old(matrix, min_value=None, max_value=None,
                             title="Default title"):
    fig, (ax1) = plt.subplots(figsize=(10, 7), ncols=1)
    plt.title(title)
    norm = colors.TwoSlopeNorm(vcenter=0)
    discrete_cmap = matplotlib.cm.get_cmap("RdBu", 10)
    figplot = ax1.imshow(matrix, vmax=max_value, vmin=min_value, norm=norm,
                         cmap="RdBu", interpolation='none')
    fig.colorbar(figplot, ax=ax1)
    plt.show()


def plot_derivatives(simulation):
    xmatrix_filename = 'X_' + simulation + '.npy'
    #xmatrix_filename = 'X_h_{}percent_ofsigma.npy'.format(h_percent)
    Xmatrix = np.load(xmatrix_filename)
    fig, ax = plt.subplots(5, 1, figsize=(6, 18), sharex=True)
    label_list = [r"$\Omega_m$", r"$\sigma_8$", r"$n_s$", r"$\Omega_b$", r"$h$"]

    for i in range(len(Xmatrix)):
        ax[i].plot(Xmatrix[i])
        ax[i].set_xlabel(label_list[i], fontsize=20)
        ax[i].set_ylabel(r'$d\mu/d$' + label_list[i], fontsize=15)

    plt.subplots_adjust(top=0.98, bottom=0.02, right=1, left=0, hspace=0.07, wspace=0.05)
    fig.tight_layout()
    #fig.savefig('derivatives_h_{}percent.pdf'.format(h_percent))
    fig.savefig('derivatives_{}.pdf'.format(simulation))


def K1000_standard_dev_symmetrized():
    # 1 sigma errors from KiDS-1000 best fit of 2PCF, symmetrized
    omega_m_symerr = 0.5 * (0.065 + 0.033)
    sigma8_symerr = 0.5 * (0.084 + 0.107)
    n_s_symerr = 0.5 * (0.093 + 0.049)
    omega_b_symerr = 0.5 * (0.005 + 0.000)  # Omega_b*h**2 = 0.046*h**2
    h_symerr = 0.5 * (0.110 + 0.001)
    A_ia_symerr = 0.5 * (0.321 + 0.374)
    return np.array([omega_m_symerr, sigma8_symerr, n_s_symerr, omega_b_symerr, h_symerr, A_ia_symerr])


def insert_new_values(parameters_dict, filename):
    """Note: this function is case sensitive, be mindful of how variables are named in CosmoCov"""

    # read file once per parameter and check which line starts with i (e.g omega_m, h0, etc)
    for i in parameters_dict:
        for line in fileinput.input([filename], inplace=True):
            if line.strip().startswith(i):
                line = i + ' : ' + str(parameters_dict[i]) + '\n'
            sys.stdout.write(line)



def RunCosmolike_Nparameters(mini_dict):
    inifile = '../CosmoCov/covs/ini_files/Ndim_kids.ini'

    # modify value(s) in file.The mini_dict can get any number of parameters
    insert_new_values(mini_dict, inifile)

    # to get simulation name
    dict_keys = list(mini_dict.keys())
    N = len(dict_keys)

    # create directories and save inifile to covs/simulation_settings/
    simulation = str(N) + 'D_' + ('-'.join(dict_keys))
    setup_simulation(simulation, inifile)

    # set directory and simulation name in ini_file
    edit_ini_file(inifile, 'output/default', 'output/' + simulation)
    edit_ini_file(inifile, 'simulation : default', 'simulation : ' + simulation)

    # Run CosmoLike
    run_CosmoLike('./Ndim_script.sh')

    # Create theory vector
    xipm_directory = "../CosmoCov/covs/xipm/" + simulation + "/"
    create_theory_vector(xipm_directory, simulation)
    mu_cosmolike = np.reshape(np.load("../CosmoCov/covs/xipm/{}/"
                                      "vector_xipm_kids_{}.npy".format(simulation, simulation)), (-1, 1))

    # Compute Xmatrix
    edit_ini_file(inifile, 'output/' + simulation, 'output/default')
    edit_ini_file(inifile, 'simulation : ' + simulation, 'simulation : default')
    Xmatrix = create_X_matrix(simulation, mu_cosmolike, inifile,
                              './Ndim_script.sh', Ndim_parameters=mini_dict)

    # Load data
    kids_covmat = np.loadtxt("../KiDS1000_cosmis_shear_data_release/data_fits/kids_covariance_matrix")
    theory_covmat = np.load("CovMatrix_class3.0_closer_derivative.npy")
    fake_data_vector = np.load("fake_datavectors.npy")[225]  # for example, Number 225
    data_xipm = np.reshape(fake_data_vector, (-1, 1))

    # Compute linear approximation of parameters
    lin_theta = linear_theta(Xmatrix, theory_covmat, mu_cosmolike,
                             data_xipm, mini_dict)

    kids_lin_theta = linear_theta(Xmatrix, kids_covmat, mu_cosmolike,
                                  data_xipm, mini_dict)

    print('parameter \t original \t theory_cov  \t kids_cov')
    for j in range(len(lin_theta)):
        print('{} \t {} \t\t {} \t\t {}\n'.format(dict_keys[j],
                                                  list(mini_dict.values())[j], round(lin_theta[j][0], 3),
                                                  round(kids_lin_theta[j][0], 3)))

    # Restore original values
    close_to_bestfit_cosmo = {'Omega_m': 0.24, 'sigma_8': 0.85, 'n_spec': 0.9,
                              'omb': 0.041, 'h0': 0.68, 'A_ia': 0.41}
    for i in mini_dict:
        insert_new_values({i: close_to_bestfit_cosmo[i]}, inifile)
