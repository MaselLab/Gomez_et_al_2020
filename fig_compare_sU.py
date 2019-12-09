# -*- coding: utf-8 -*-
"""
Created on Tue Jan 15 14:09:05 2019
Masel Lab
Project: Mutation-driven Adaptation
@author: Kevin Gomez

Description:
Script for creating a figure that compares the rates of adaptation in two 
traits of an asexual population, evolving together in the concurrent mutations
regime.
"""

#libraries
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import multivariate_normal
from numpy import inf
import matplotlib.ticker as mtick

import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab
    
# set file name of data
# -----------------------------------------------------------------------------
pickle_file_name = 'fig_compare_sU_data-03.pickle'    # data for comparison with fixed v

# load processed matlab data for figure
# -----------------------------------------------------------------------------
pickle_file = open("data/" + pickle_file_name,'rb') 
[N,s,U,v,parameters,grand_means,sarry,Uarry,v1_data,v2_data] = pickle.load(pickle_file)
pickle_file.close()

Uarry = np.flip(Uarry,0)
rate_comp = v1_data/v

# get rate_comp array dimensions 
[m,n] = rate_comp.shape
[m,n] = [int(m),int(n)]
arry_dim = len(sarry)

# set labels for axes
my_slabel = ['$10^{'+str(np.round(np.log10(sarry[i,0]),1))+'}$' for i in range(len(sarry))]
my_Ulabel = ['$10^{'+str(np.round(np.log10(Uarry[i,0]),1))+'}$' for i in range(len(Uarry))]

# set some labels blank to have them fit on graph
for i in range(len(my_slabel)):
    if ((i%3==1) | (i%3==2)):
        my_Ulabel[i]=''

for i in range(len(my_slabel)):
    if ((i%4==1) | (i%4==2) | (i%4==3)):
        my_slabel[i]=''
        
x_border = [0.0+m*i/1000.0 for i in range(1001)]
y_border = [min(np.floor(1.0+m*i/1000.0),m) for i in range(1001)]

# create heatmap of v reduction
# -----------------------------------------------------------------------------
fig1, ax1 = plt.subplots(1,1,figsize=[11,10])
fit_distr_2d = ax1.pcolormesh(rate_comp,cmap=plt.cm.gray_r)

cbar = plt.colorbar(fit_distr_2d,pad = 0.03)
cbar.ax.tick_params(labelsize=16)

# *****************************************************************************

#ax1.plot(x_border,y_border,color="black")
ax1.axis('tight')        
ax1.set_xticks(np.arange(arry_dim)+0.5)
ax1.set_yticks(np.arange(arry_dim)+0.5)        
ax1.set_xticklabels(my_slabel)
ax1.set_yticklabels(my_Ulabel[::-1])        
ax1.set_xlabel('Selection coefficient trait 2',multialignment='center',fontsize=18,labelpad=10)
ax1.set_ylabel('Mutation rate trait 2',multialignment='center',fontsize=18,labelpad=10)
ax1.tick_params(axis='both',labelsize=16)
# *****************************************************************************
cbar.ax.text(3,0.55,'Ratio $v_1/v$',rotation=270,fontsize=22)    # use this label of comparing v
# *****************************************************************************
plt.text(18,31.5,r'$s_1 = 10^{-2}$,',fontsize=18)
plt.text(24.5,31.5,r'$U_1 = 10^{-5}$',fontsize=18)     # use this label of comparing v
plt.text(-6.0,30,'(a)',fontsize=20)
# *****************************************************************************
fig1.subplots_adjust(bottom=0.2,left=0.2)
fig1.savefig('figures/fig_two_trait_compare_sU_data_v1.pdf',bbox_inches='tight')


# *****************************************************************************
#                Figure identifying tradeoff curve
# *****************************************************************************

min_v2_data = 1

for i in range(v2_data.shape[0]):
    for j in range(v2_data.shape[1]):
        if(v2_data[i,j] >0):
            min_v2_data = min([v2_data[i,j],min_v2_data])

for i in range(v2_data.shape[0]):
    for j in range(v2_data.shape[1]):
        if(v2_data[i,j] <= 0):
            v2_data[i,j] = min_v2_data

for i in range(v2_data.shape[0]):
    for j in range(v2_data.shape[1]):
        rate_comp[i,j] = min([v1_data[i,j],v2_data[i,j]])/v
        
# create heatmap of v reduction
# -----------------------------------------------------------------------------
fig2, ax2 = plt.subplots(1,1,figsize=[11,10])
fit_distr_2d = ax2.pcolormesh(rate_comp,cmap=plt.cm.gray_r)

cbar = plt.colorbar(fit_distr_2d,pad = 0.03)
cbar.ax.tick_params(labelsize=16)

# *****************************************************************************
#ax1.plot(x_border,y_border,color="black")
ax2.axis('tight')        
ax2.set_xticks(np.arange(arry_dim)+0.5)
ax2.set_yticks(np.arange(arry_dim)+0.5)        
ax2.set_xticklabels(my_slabel)
ax2.set_yticklabels(my_Ulabel[::-1])        
ax2.set_xlabel('Selection coefficient trait 2',multialignment='center',fontsize=18,labelpad=10)
ax2.set_ylabel('Mutation rate trait 2',multialignment='center',fontsize=18,labelpad=10)
ax2.tick_params(axis='both',labelsize=16)
# *****************************************************************************
cbar.ax.text(3,0.7,'Ratio $\min(v_1,v_2)/v$',rotation=270,fontsize=22)    # use this label of comparing v
# *****************************************************************************
plt.text(17,31.5,r'$N = 10^9$,',fontsize=18)
plt.text(22.5,31.5,r'$v = 5.3\times 10^{-5}$',fontsize=18)     # use this label of comparing v
plt.text(-6,30,'(b)',fontsize=20)
# *****************************************************************************

fig2.subplots_adjust(bottom=0.2,left=0.2)            
fig2.savefig('figures/fig_two_trait_compare_sU_data_tradeoff.pdf',bbox_inches='tight')
