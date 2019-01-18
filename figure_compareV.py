# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 14:09:05 2019

@author: dirge
"""

from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import multivariate_normal
from numpy import inf
import matplotlib.ticker as mtick

import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab


# load processed matlab data for figure
# -----------------------------------------------------------------------------
pickle_file_name = 'fig_compareVdata-01.pickle'
pickle_file = open(pickle_file_name,'rb') 
[N,s,U,v,parameters,grand_means,sarry,Uarry,v1_data,v2_data] = pickle.load(pickle_file)
pickle_file.close()

v_comp = np.log10(v1_data/v2_data)
arry_dim = len(sarry)
my_xlabel = [str(np.round(np.log10(sarry[i,0]),1)) for i in range(len(sarry))]
my_ylabel = [str(np.round(np.log10(Uarry[i,0]),1)) for i in range(len(Uarry))]

# create heatmap of v reduction
# -----------------------------------------------------------------------------
fig1, ax1 = plt.subplots(1,1,figsize=[10,8])
fit_distr_2d = ax1.pcolormesh(v_comp.transpose(),cmap=plt.cm.gray_r)
cbar = plt.colorbar(fit_distr_2d)
ax1.axis('tight')        
ax1.set_xticks(np.arange(arry_dim)+0.5)
ax1.set_yticks(np.arange(arry_dim)+0.5)        
ax1.set_xticklabels(my_xlabel)
ax1.set_yticklabels(my_ylabel)        
ax1.set_xlabel('Trait 1 (s Log10 Scale)',fontsize=18,labelpad=20)
ax1.set_ylabel('Trait 2 (U Log10 Scale)',fontsize=18,labelpad=10)
ax1.tick_params(axis='both',labelsize=14)        
cbar.ax.text(2.5,0.65,'Log$_{10}$ of $v_1/v_2$',rotation=90,fontsize=18)
fig1.savefig('fig_compareVdata.pdf')