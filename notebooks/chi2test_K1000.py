import sys
import importlib
sys.path.append('../PyScripts/')

import utils.function_library 
importlib.reload(utils.function_library)
from utils.function_library import *


def chi_square_test_K1000(data_covariance_matrix, theory_vector, data_vector):
    """ 
    Reduced chi-square test for the KiDS-1000 cosmic shear paper (asgari 2020) 
    where they discard first three values of xi_m for the chi-square test.
    """

    reshaped_covmat, mu_reshaped, data_reshaped = reduce_dimension(data_covariance_matrix,
                                                                   theory_vector, data_vector)

    chi = (mu_reshaped - data_reshaped).T @ inv(reshaped_covmat) @ (mu_reshaped - data_reshaped)

    #dof = degrees of freedom (see Asgari 2020 page 9 and table 3)
    Ndim = 270 - (3 * 15)  # first 3 values of each xipm are discarded due
                            # to their sensitivity to small physical scales
    Nparams = 4.5
    dof = Ndim - Nparams

    print('Minimum chi-square to reject null hypothesis: %.3f\n' %(chi2.ppf(0.95, dof)))
    print('chi-squared = %.3f, p-value = %r' % (chi, 1 - chi2.cdf(chi, dof)))
    print('chi-squared reduced from kids is 260, and p value is 0.034\n')

    return chi[0,0]
  
  
  def reduce_dimension(covmat, theory_vector, data_vector):
    """
    Reduce the dimension of the covmat and theory/data vectors, following the KiDS-1000 cosmic 
    shear paper (asgari 2020) where they discard first three values of xi_m for the chi-square test.
    
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
