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
    
# load processed matlab data for figure
# -----------------------------------------------------------------------------
pickle_file_name = 'fig_compareVdata-01.pickle'
pickle_file = open("data\\" + pickle_file_name,'rb') 
[N,s,U,v,parameters,grand_means,sarry,Uarry,v1_data,v2_data] = pickle.load(pickle_file)
pickle_file.close()

v_comp = np.log10(v1_data/v2_data)
[m,n] = v_comp.shape
[m,n] = [int(m),int(n)]

for i in range(n-1):
    for j in range(m-1-i):
        v_comp[i,m-1-j]= 0

arry_dim = len(sarry)
my_slabel = ['$'+str(np.round(np.log10(sarry[i,0]),1))+'$' for i in range(len(sarry))]
my_Ulabel = ['$'+str(np.round(np.log10(Uarry[i,0]),1))+'$' for i in range(len(Uarry))]

x_border = [0.0+m*i/1000.0 for i in range(1001)]
y_border = [min(np.floor(1.0+m*i/1000.0),m) for i in range(1001)]

# create heatmap of v reduction
# -----------------------------------------------------------------------------
fig1, ax1 = plt.subplots(1,1,figsize=[15,10])
fit_distr_2d = ax1.pcolormesh(v_comp.transpose(),cmap=plt.cm.gray_r)
cbar = plt.colorbar(fit_distr_2d,pad = 0.15)
#ax1.plot(x_border,y_border,color="black")
#ax1.axis('tight')        
ax1.set_xticks(np.arange(arry_dim)+0.5)
ax1.set_yticks(np.arange(arry_dim)+0.5)        
ax1.set_xticklabels(my_slabel)
ax1.set_yticklabels(my_slabel)        
ax1.set_xlabel('Selection coefficient trait 1 ($\log_{10}$)',multialignment='center',fontsize=18,labelpad=10)
ax1.set_ylabel('Selection coefficient trait 2 ($\log_{10}$)',multialignment='center',fontsize=18,labelpad=10)
ax1.tick_params(axis='both',labelsize=16)
cbar.ax.text(1.8,0.6,'$\log_{10}$ of $v_1/v_2$',rotation=90,fontsize=22)

plt.text(1,m-1,r'$N = 10^9$',fontsize=18)
plt.text(1,m-1.5,r'$v = 5.3\times 10^{-5}$',fontsize=18)

ax2 = ax1.twinx()
ax2.yaxis.set_ticks(np.arange(0+0.5/(len(my_Ulabel)),1-0.5/(len(my_Ulabel)),1.0/(len(my_Ulabel))))
ax2.set_yticklabels(my_Ulabel)
ax2.tick_params(labelsize=16)
ax2.set_ylabel('Mutation rate trait 2 ($\log_{10}$)',fontsize=18,labelpad=20)

ax3 = ax1.twiny()
ax3.xaxis.set_ticks(np.arange(0+0.5/(len(my_Ulabel)),1-0.5/(len(my_Ulabel)),1.0/(len(my_Ulabel))))
ax3.set_xticklabels(my_Ulabel)
ax3.tick_params(labelsize=16)
ax3.set_xlabel('Mutation rate trait 1 ($\log_{10}$)',fontsize=18,labelpad=20)
#plt.tight_layout()
fig1.subplots_adjust(bottom=0.2,left=0.2)

plt.text(.21,-1.9,'Trait 1 favored by selection',fontsize=22)
plt.text(-0.2,9.2,'Trait 2 favored by mutation',rotation=90,fontsize=22)

fig1.savefig('figures\\fig_compareVdata.pdf')

# need to creat a figure with v_g rather than v_w used above