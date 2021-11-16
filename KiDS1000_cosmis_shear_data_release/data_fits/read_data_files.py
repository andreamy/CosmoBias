import astropy.io.fits as fits


cat_version = 'V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid'


#############################################################################
# Band powers

stats_name ='bp'

filename = stats_name+'_KIDS1000_BlindC_with_m_bias_'+cat_version+'.fits'
F=fits.open(filename)

# Gives the names of all the extension
F.info()

# Covariance matrix
ext=F["COVMAT"]
covariance= ext.data

# cosmic shear data vector
ext=F['PeeE']
data_vector=ext.data['VALUE']

# GGL data vector
ext=F['PneE']
data_vector=ext.data['VALUE']

# nofz sources (shear)
ext=F["nz_source"]

# nofz lenses (positions)
ext=F["nz_lens"]

F.close()

#############################################################################
# COSEBIs

stats_name ='cosebis'

filename = stats_name+'_KIDS1000_BlindC_with_m_bias_'+cat_version+'.fits'
F=fits.open(filename)

# Covariance matrix
ext=F["COVMAT"]
covariance= ext.data

# cosmic shear data vector
ext=F['En']
data_vector=ext.data['VALUE']

F.close()

#############################################################################
# 2PCFs

stats_name ='xipm'

filename = stats_name+'_KIDS1000_BlindC_with_m_bias_'+cat_version+'.fits'
F=fits.open(filename)

# Covariance matrix includes covariances for xip, xim and their cross covariance
ext=F["COVMAT"]
covariance= ext.data

# cosmic shear data vector: xi_plus
ext=F['xip']
xip=ext.data['VALUE']

# cosmic shear data vector: xi_minus
ext=F['xim']
xim=ext.data['VALUE']


F.close()