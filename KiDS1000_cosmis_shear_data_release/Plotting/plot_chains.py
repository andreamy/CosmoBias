import matplotlib.pylab as plt
import numpy as np
from chainconsumer import ChainConsumer
from matplotlib import colors as mcolors
from chain_functions import makeChainMultinest
from chain_functions import input_names, parameter_names, find_cols_for_params, mkdir_mine

colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
# Sort colors by hue, saturation, value and name.
by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                for name, color in colors.items())
sorted_names = [name for hsv, name in by_hsv]


#####################################################

output="plots"
# where the plots will be saved
mkdir_mine(output)
outputfolder = output

name_tag = "KiDS-1000"


chain_name_cosebis = "COSEBIs"
chain_name_xipm    = "2PCFs"
chain_name_bp      = "Band Powers"

color_cosebis = colors['darkorange']
color_bp      = colors['hotpink']
color_xipm    = colors['darkturquoise']

FolderName = '../chains_and_config_files/main_chains_iterative_covariance/'

##########################################################################################
# Read in the multinest chain for COSEBIs
filename=FolderName+'/cosebis/chain/output_multinest_C.txt' 
file=open(filename)
line=file.readline()
parameter_list_cosebis=(line.replace('#','')).split()
file.close()
chain_cosebis,weights_cosebis, n= makeChainMultinest(filename,'all',True)


##########################################################################################
# read in the multinest chain bp
filename=FolderName+'/bp/chain/output_multinest_C.txt' 
file=open(filename)
line=file.readline()
parameter_list_bp=(line.replace('#','')).split()
file.close()
chain_bp,weights_bp, n= makeChainMultinest(filename,'all',True)

##########################################################################################
# read in the multinest chain xipm
filename=FolderName+'/xipm/chain/output_multinest_C.txt' 
file=open(filename)
line=file.readline()
parameter_list_xipm=(line.replace('#','')).split()
file.close()
chain_xipm,weights_xipm, n= makeChainMultinest(filename,'all',True)

###############################
# 
# now plotting, which parameter do you want to plot?
# for list of parameter names see parameter_names in chain_functions.py
# 
###############################
# plot S_8
param = ["s8"]
name = ''
for p in param:
	name+=p
	name+='_'

name+=name_tag
# 
c = ChainConsumer()
# Find the relevant columns
cols, param_names=find_cols_for_params(parameter_list_cosebis,input_names,parameter_names,param)
c.add_chain(chain_cosebis[:,cols], weights=weights_cosebis,linestyle='-',color=color_cosebis,parameters=param_names, name=chain_name_cosebis)
cols, param_names=find_cols_for_params(parameter_list_bp,input_names,parameter_names,param)
c.add_chain(chain_bp[:,cols], weights=weights_bp,linestyle='-',color=color_bp,parameters=param_names, name=chain_name_bp,shade_alpha=0.8)
cols, param_names=find_cols_for_params(parameter_list_xipm,input_names,parameter_names,param)
c.add_chain(chain_xipm[:,cols], weights=weights_xipm,linestyle='-',color=color_xipm,parameters=param_names, name=chain_name_xipm,shade_alpha=0.8)
c.configure(kde=1.0,statistics="max",shade_gradient=1.0,shade_alpha=0.8,bar_shade=True)
c.plotter.plot(filename=outputfolder+"/"+name+".pdf", figsize="column")
plt.close()
#
#################################
# plot omega_m S_8
param = ["omega_m","s8"]
name = ''
for p in param:
	name+=p
	name+='_'

name+=name_tag
#
cols, param_names=find_cols_for_params(parameter_list_cosebis,input_names,parameter_names,param)
c = ChainConsumer()
c.add_chain(chain_cosebis[:,cols], weights=weights_cosebis,linestyle='-',color=color_cosebis,parameters=param_names, name=chain_name_cosebis,shade_alpha=0.8)
cols, param_names=find_cols_for_params(parameter_list_bp,input_names,parameter_names,param)
c.add_chain(chain_bp[:,cols], weights=weights_bp,linestyle='-',color=color_bp,parameters=param_names, name=chain_name_bp,shade_alpha=0.7)
cols, param_names=find_cols_for_params(parameter_list_xipm,input_names,parameter_names,param)
c.add_chain(chain_xipm[:,cols], weights=weights_xipm,linestyle='-',color=color_xipm,parameters=param_names, name=chain_name_xipm,shade_alpha=0.6)
c.configure(kde=1.5,plot_hists=False,shade_gradient=1.0,diagonal_tick_labels=False,label_font_size=14,tick_font_size=13,serif=True,legend_color_text=True,linewidths=1.5,statistics="max",shade=True)
c.plotter.plot(filename=outputfolder+"/"+name+".pdf", figsize=1.5)
plt.close()
# 
#################################
# plot omega_m sigma_8
param = ["omega_m","sigma_8"]
name = ''
for p in param:
	name+=p
	name+='_'

name+=name_tag

cols, param_names=find_cols_for_params(parameter_list_cosebis,input_names,parameter_names,param)

c = ChainConsumer()
c.add_chain(chain_cosebis[:,cols], weights=weights_cosebis,linestyle='-',color=color_cosebis,parameters=param_names, name=chain_name_cosebis,shade_alpha=1.0)
cols, param_names=find_cols_for_params(parameter_list_bp,input_names,parameter_names,param)
c.add_chain(chain_bp[:,cols], weights=weights_bp,linestyle='-',color=color_bp,parameters=param_names, name=chain_name_bp,shade_alpha=0.7)
cols, param_names=find_cols_for_params(parameter_list_xipm,input_names,parameter_names,param)
c.add_chain(chain_xipm[:,cols], weights=weights_xipm,linestyle='-',color=color_xipm,parameters=param_names, name=chain_name_xipm,shade_alpha=0.6)
c.configure(kde=1.5,plot_hists=False,shade_gradient=1.0,diagonal_tick_labels=False,label_font_size=14,tick_font_size=13,serif=True,legend_color_text=True,linewidths=1.5,statistics="max",shade=True)
c.plotter.plot(filename=outputfolder+"/"+name+".pdf", figsize=1.5)
plt.close()
# 
