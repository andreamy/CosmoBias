def plot_covmat_differences(standard_covmat, shifted_covmat)
    """
    Plots the percentual deviation of each covariance matrix element, if the free 
    parameters of the covariance are varied within their error bars. (e.g. n_eff 
    or A_sky, or the fiducial cosmology).
    
    Parameters
    ----------
    standard_covmat : 2darray
                     covariance matrix computed with standard (fiducial) 
                     free parameters, without changes.
                     
    shifted_covmat : 2darray
                     covariance matrix computed with one (or more) free 
                     parameters varied within their error bars.
    """
    
    Ndim = standard_covmat.shape[0]
    percent_deviation = np.zeros(shape=(Ndim, Ndim))
    
    for i in range(Ndim):
        for j in range(Ndim):
            percent_deviation[i,j] = abs((shifted_covmat[i,j] - 
                                          standard_covmat[i,j])/standard_covmat[i,j])*100

    fig, (ax1) = plt.subplots(figsize=(10, 7), ncols=1)
    plt.title(r'% deviation from fiducial covariance')
    figplot = ax1.imshow(percent_deviation, vmax = 20, cmap='viridis', 
                         interpolation='none')
    fig.colorbar(figplot, ax=ax1)
    plt.savefig("covmat_deviation.pdf")
    plt.show()
