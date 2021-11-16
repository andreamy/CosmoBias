import numpy as np
import os

# reads in multinest chains from a file and makes a chain with the given columns
# If nSample is set to True only takes the last
# nSamples read from the last line of the file
# if cols set to all returns all columns
def makeChainMultinest(filename,cols,nSampleOn):
	file=open(filename)
	chain_in=np.loadtxt(file,comments='#')
	# print(cols)
	# print(chain_in.shape)
	if(nSampleOn):
		file=open(filename)
		lines=file.readlines()
		file.close()
		nSample_str=lines[-3]
		nSample=int(nSample_str[9:-1])
	else:
		nSample=len(chain_in)
	weight=chain_in[-nSample:,-1]
	if(cols=='all'):
		return chain_in[-nSample:,:], weight, nSample
	else:
		return chain_in[-nSample:,cols], weight, nSample

def mkdir_mine(dirName):
    try:
        # Create target Directory
        os.mkdir(dirName)
        # print("Directory " , dirName ,  " Created ") 
    except FileExistsError:
    	x=1
        # print("Directory " , dirName ,  " already exists")


def find_cols_for_params(parameter_list,input_names,parameter_names,my_names):
	cols=[]
	param_names=[]
	for param in my_names:
		cols.append(parameter_list.index(input_names[param]))
		param_names.append(parameter_names[param])
	return cols,param_names

def find_best_fit_maxpost(filename,input_names,parameter_names,multinest=False):
	file=open(filename)
	line=file.readline()
	parameter_list = (line.replace('#','')).split()
	if(multinest):
		chain,weights, n= makeChainMultinest(filename,'all',True)
	else:
		chain  = np.loadtxt(file,comments='#')
	cols, param_names=find_cols_for_params(parameter_list,input_names,parameter_names,['post','like'])
	if(chain.ndim==1):
		best_fit = chain
		post     = chain[cols[0]]
		like     = chain[cols[1]]
	else:
		minimum_chi = np.argmax(chain[:,cols[0]])
		best_fit = chain[minimum_chi,:]
		post     = chain[minimum_chi,cols[0]]
		like     = chain[minimum_chi,cols[1]]
	return best_fit,post,like


def find_best_fit_maxlike(filename,input_names,parameter_names,multinest=False):
	file=open(filename)
	line=file.readline()
	parameter_list = (line.replace('#','')).split()
	if(multinest):
		chain,weights, n= makeChainMultinest(filename,'all',True)
	else:
		chain  = np.loadtxt(file,comments='#')
	cols, param_names=find_cols_for_params(parameter_list,input_names,parameter_names,['post','like'])
	minimum_chi = np.argmax(chain[:,cols[1]])
	best_fit = chain[minimum_chi,:]
	post     = chain[minimum_chi,cols[0]]
	like     = chain[minimum_chi,cols[1]]
	return best_fit,post,like


# numerical Hellinger distance based on two samples (histogram)
def hellinger_sample(sample1,weight1,sample2,weight2):
	import scipy.stats
	from scipy.spatial.distance import euclidean
	bin_low = np.floor(np.minimum(min(sample1),min(sample2)))
	bin_upp = np.ceil(np.maximum(max(sample1),max(sample2)))
	bin_n = int(np.ceil(np.sqrt(np.maximum(len(sample1),len(sample2)))))
	bins = np.linspace(bin_low,bin_upp,bin_n)

	histogram1, bins = np.histogram(sample1, bins=bins, weights=weight1, density=True)
	histogram2, bins = np.histogram(sample2, bins=bins, weights=weight2,density=True)
	histogram1 = histogram1 / np.sum(histogram1)  #normalise
	histogram2 = histogram2 / np.sum(histogram2)  #normalise
	return(euclidean(np.sqrt(histogram1), np.sqrt(histogram2)) / np.sqrt(2.))


def hellinger_sample_kde(sample1,weight1,sample2,weight2):
	import scipy.stats
	from scipy.spatial.distance import euclidean
	bin_low = np.floor(np.minimum(min(sample1),min(sample2)))
	bin_upp = np.ceil(np.maximum(max(sample1),max(sample2)))
	bin_n = int(np.ceil(np.sqrt(np.maximum(len(sample1),len(sample2)))))
	bins = np.linspace(bin_low,bin_upp,bin_n)

	kde1 = scipy.stats.gaussian_kde(sample1, weights=weight1)
	kde2 = scipy.stats.gaussian_kde(sample2, weights=weight2)
	pdf1 = scipy.stats.gaussian_kde.pdf(kde1,bins)
	pdf2 = scipy.stats.gaussian_kde.pdf(kde2,bins)
	pdf1 = pdf1 / np.sum(pdf1)  # normalise
	pdf2 = pdf2 / np.sum(pdf2)  # normalise
	return(euclidean(np.sqrt(pdf1), np.sqrt(pdf2)) / np.sqrt(2.))

# inverse function for given std deviations
def inverse_hellinger_gauss(H,std1,std2):
    s1sq = np.square(std1)
    s2sq = np.square(std2)
    root = np.sqrt(2*std1*std2 / (s1sq+s2sq))
    return( np.sqrt( (-4.)*(s1sq+s2sq) * np.log((1-np.square(H))/root) ) )

# Hellinger distance-based tension estimate for two samples, each with associated weights
def hellinger_tension(sample1,weight1,sample2,weight2,htype=3):
	import scipy.stats
	mean1 = np.average(sample1,weights=weight1)
	mean2 = np.average(sample2,weights=weight2)
	var1 = len(sample1)/(len(sample1)-1.) * np.average(np.square(sample1-mean1),weights=weight1)
	var2 = len(sample2)/(len(sample2)-1.) * np.average(np.square(sample2-mean2),weights=weight2)

	if (htype==1):
		d_hell = hellinger_sample(sample1,weight1,sample2,weight2)
	elif (htype==2):
		d_hell = hellinger_sample_kde(sample1,weight1,sample2,weight2)
	elif (htype==3):
		d_hell = 0.5*(hellinger_sample(sample1,weight1,sample2,weight2)+hellinger_sample_kde(sample1,weight1,sample2,weight2))
	else:
		print("Incorrect value for htype.")
		exit(-1)

	# find mean_diff for which hellinger_gauss(mean_diff,std1,std2)=d_hell
	mean_diff = inverse_hellinger_gauss(d_hell,np.sqrt(var1),np.sqrt(var2))
	return(mean_diff/np.sqrt(var1+var2))


input_names={ 'omch2' :'cosmological_parameters--omch2', 
			  'ombh2':'cosmological_parameters--ombh2', 
			  'h'    :'cosmological_parameters--h0', 
			  'h_out'     :'COSMOLOGICAL_PARAMETERS--H0',
			  'n_s'   :'cosmological_parameters--n_s', 
			  's8_in':'cosmological_parameters--s_8_input', 
			  'a_bar'      :'halo_model_parameters--a', 
			  'a_ia':'intrinsic_alignment_parameters--a', 
			  'alpha_ia': 'intrinsic_alignment_parameters--alpha',
			  'deltaz_uncorr_1' :'nofz_shifts--uncorr_bias_1', 
			  'deltaz_uncorr_2' :'nofz_shifts--uncorr_bias_2', 
			  'deltaz_uncorr_3' :'nofz_shifts--uncorr_bias_3', 
			  'deltaz_uncorr_4' :'nofz_shifts--uncorr_bias_4', 
			  'deltaz_uncorr_5' :'nofz_shifts--uncorr_bias_5', 
			  'uncorr_m1' :'shear_calibration_parameters--uncorr_m1', 
			  'uncorr_m2' :'shear_calibration_parameters--uncorr_m2', 
			  'uncorr_m3' :'shear_calibration_parameters--uncorr_m3', 
			  'uncorr_m4' :'shear_calibration_parameters--uncorr_m4', 
			  'uncorr_m5' :'shear_calibration_parameters--uncorr_m5', 
			  'delta_c' : 'shear_c_bias--delta_c',
			  's8':    'COSMOLOGICAL_PARAMETERS--S_8', 
			  'sigma_8':  'COSMOLOGICAL_PARAMETERS--SIGMA_8', 
			  'a_s': 'COSMOLOGICAL_PARAMETERS--A_S', 
			  'omega_m':    'COSMOLOGICAL_PARAMETERS--OMEGA_M', 
			  'omega_nu':    'COSMOLOGICAL_PARAMETERS--OMEGA_NU', 
			  'omega_lambda':       'COSMOLOGICAL_PARAMETERS--OMEGA_LAMBDA', 
			  'theta_mc':      'COSMOLOGICAL_PARAMETERS--COSMOMC_THETA', 
			  'deltaz_1':'NOFZ_SHIFTS--BIAS_1', 
			  'deltaz_2':'NOFZ_SHIFTS--BIAS_2', 
			  'deltaz_3':'NOFZ_SHIFTS--BIAS_3', 
			  'deltaz_4':'NOFZ_SHIFTS--BIAS_4', 
			  'deltaz_5':'NOFZ_SHIFTS--BIAS_5', 
			  'deltaz_out_1':'DELTA_Z_OUT--BIN_1', 
			  'deltaz_out_2':'DELTA_Z_OUT--BIN_2', 
			  'deltaz_out_3':'DELTA_Z_OUT--BIN_3', 
			  'deltaz_out_4':'DELTA_Z_OUT--BIN_4', 
			  'deltaz_out_5':'DELTA_Z_OUT--BIN_5',
			  'm1':'SHEAR_CALIBRATION_PARAMETERS--M1',
			  'm2':'SHEAR_CALIBRATION_PARAMETERS--M2',
			  'm3':'SHEAR_CALIBRATION_PARAMETERS--M3',
			  'm4':'SHEAR_CALIBRATION_PARAMETERS--M4',
			  'm5':'SHEAR_CALIBRATION_PARAMETERS--M5',
			  'prior':'prior', 
			  'like':'like', 
			  'post':'post', 
			  'weight':'weight'}

parameter_names={ 'omch2' : r'$\Omega_{\rm c}h^2$', 
			      'ombh2':r'$\Omega_{\rm b}h^2$', 
			      'h'    :r'$h$', 
			      'h_out':r'$h$', 
				  'n_s'   :r'$n_{\rm s}$', 
				  's8_in':r'$S_8^{\rm in}$', 
				  'a_bar' :r'$A_{\rm bary}$', 
				  'a_ia'  :r'$A_{\rm IA}$', 
				  'alpha_ia': r'$\alpha_{\rm IA}$',
				  'deltaz_uncorr_1' :r'$\Delta z_1^{\rm uncorr}$', 
				  'deltaz_uncorr_2' :r'$\Delta z_2^{\rm uncorr}$', 
				  'deltaz_uncorr_3' :r'$\Delta z_3^{\rm uncorr}$', 
				  'deltaz_uncorr_4' :r'$\Delta z_4^{\rm uncorr}$', 
				  'deltaz_uncorr_5' :r'$\Delta z_5^{\rm uncorr}$', 
				  'uncorr_m1' :r'$m_1^{\rm uncorr}$', 
				  'uncorr_m2' :r'$m_2^{\rm uncorr}$', 
				  'uncorr_m3' :r'$m_3^{\rm uncorr}$', 
				  'uncorr_m4' :r'$m_4^{\rm uncorr}$', 
				  'uncorr_m5' :r'$m_5^{\rm uncorr}$', 
				  'delta_c' : r'$\delta_{\rm c}$',
				  's8':    r'$S_8$', 
				  'sigma_8':  r'$\sigma_8$', 
				  'a_s': r'$A_{\rm s}$', 
				  'omega_m':  r'$\Omega_{\rm m}$'  , 
				  'omega_nu':    r'$\Omega_{\nu}$' , 
				  'omega_lambda': r'$\Omega_{\Lambda}$' , 
				  'theta_mc':     r'$\theta_{\rm MC}$' , 
				  'deltaz_1':r'$\Delta z_1^{\rm in}$', 
				  'deltaz_2':r'$\Delta z_2^{\rm in}$', 
				  'deltaz_3':r'$\Delta z_3^{\rm in}$', 
				  'deltaz_4':r'$\Delta z_4^{\rm in}$', 
				  'deltaz_5':r'$\Delta z_5^{\rm in}$', 
				  # 'deltaz_out_1':r'$\Delta z_1^{\rm out}$', 
				  # 'deltaz_out_2':r'$\Delta z_2^{\rm out}$', 
				  # 'deltaz_out_3':r'$\Delta z_3^{\rm out}$', 
				  # 'deltaz_out_4':r'$\Delta z_4^{\rm out}$', 
				  # 'deltaz_out_5':r'$\Delta z_5^{\rm out}$', 
				  'deltaz_out_1':r'$\delta_{z,1}$', 
				  'deltaz_out_2':r'$\delta_{z,2}$', 
				  'deltaz_out_3':r'$\delta_{z,3}$', 
				  'deltaz_out_4':r'$\delta_{z,4}$', 
				  'deltaz_out_5':r'$\delta_{z,5}$', 
				  'm1':r'$m_1$',
				  'm2':r'$m_2$',
				  'm3':r'$m_3$',
				  'm4':r'$m_4$',
				  'm5':r'$m_5$',
				  'prior':'prior', 
				  'like':'like', 
				  'post':'post', 
				  'weight':'weight'}