import astropy.io.fits as fits
import numpy as np

cat_version = 'V1.0.0A_ugriZYJHKs_photoz_SG_mask_LF_svn_309c_2Dbins_v2_goldclasses_Flag_SOM_Fid'


#############################################################################
# Band powers

stats_name ='bp'

filename = stats_name+'_KIDS1000_BlindC_with_m_bias_'+cat_version+'.fits'
F=fits.open(filename)

# Gives the names of all the extension
F.info()


# Covariance matrix
#ext=F["COVMAT"]
#covariance= ext.data

# cosmic shear data vector
ext=F['PeeE']
data_vector=ext.data
#data_vector=ext.data['VALUE']


# nofz sources (shear)
ext=F["nz_source"]
noz = ext.data

F.close()

#############################################################################

print(data_vector.dtype)
print(noz.dtype)
np.savetxt("shear_data_vector.txt", data_vector)
#np.savetxt("shear_noz.txt", noz)
