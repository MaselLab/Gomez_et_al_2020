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

#pickle_file_name = 'fig_compareVdata-16.pickle'    # data for comparison with fixed v
#pickle_file_name = 'fig_compareVdata-19.pickle'    # data for comparison with fixed R
#pickle_file_name = 'fig_compareVdata-20.pickle'    # data for comparison with fixed R
#pickle_file_name = 'fig_compareVdata-100.pickle'    # data for comparison with fixed v
#pickle_file_name = 'fig_compareVdata-101.pickle'    # data for comparison with fixed v
pickle_file_name = 'fig_compareVdata-102.pickle'    # data for comparison with fixed v

# load processed matlab data for figure
# -----------------------------------------------------------------------------

pickle_file = open("data/" + pickle_file_name,'rb') 
[N,s,U,v,parameters,grand_means,sarry,Uarry,v1_data,v2_data] = pickle.load(pickle_file)
pickle_file.close()

# comment out lines depending on R or v (these sections indicated by ****)

rate_comp = v1_data/v2_data
rate_comp = np.transpose(rate_comp[::-1])[::-1]

[m,n] = rate_comp.shape
[m,n] = [int(m),int(n)]

for i in range(n-1):                # remove lower portion of grid
    for j in range(m-1-i):
        rate_comp[i,m-1-j]= 1

for i in range(n):                  #set R1/R2 or v1/v2 < 1 to 1 for cbar map correction  
    for j in range(m-i):
# *****************************************************************************        
#        rate_comp[i+j,i]=rate_comp[i+j,i]*(sarry[i]/sarry[i+j])
# *****************************************************************************
        if (rate_comp[i+j,i]< 1):
            rate_comp[i+j,i]=2-rate_comp[i+j,i]

# flip the entries to get right matching of values with ordinal axis U vs. s
for i in range(n):
    rate_comp[i] = rate_comp[i][::-1]

arry_dim = len(sarry)
my_slabel = ['$10^{'+str(np.round(np.log10(sarry[i,0]),1))+'}$' for i in range(len(sarry))]
my_Ulabel = ['$10^{'+str(np.round(np.log10(Uarry[i,0]),1))+'}$' for i in range(len(Uarry))]

#my_slabel = ['$'+str(np.round(np.log10(sarry[i,0]),2))+'$' for i in range(len(sarry))]
#my_Ulabel = ['$'+str(np.round(np.log10(Uarry[i,0]),2))+'$' for i in range(len(Uarry))]

for i in range(len(my_slabel)):
    if (i%2==1):
        my_slabel[i]=''
        my_Ulabel[i]=''

x_border = [0.0+m*i/1000.0 for i in range(1001)]
y_border = [min(np.floor(1.0+m*i/1000.0),m) for i in range(1001)]

# create heatmap of v reduction
# -----------------------------------------------------------------------------
fig1, ax1 = plt.subplots(1,1,figsize=[15,10])

# *****************************************************************************
#fit_distr_2d = ax1.pcolormesh(rate_comp.transpose(),cmap=plt.cm.gray_r)
#cbar = plt.colorbar(fit_distr_2d,pad = 0.15,ticks=[1.0+.1*i for i in range(13)])

fit_distr_2d = ax1.pcolormesh(rate_comp,cmap=plt.cm.gray_r)
cbar = plt.colorbar(fit_distr_2d,pad = 0.15)
# *****************************************************************************

#ax1.plot(x_border,y_border,color="black")
ax1.axis('tight')        
ax1.set_xticks(np.arange(arry_dim)+0.5)
ax1.set_yticks(np.arange(arry_dim)+0.5)        
ax1.set_xticklabels(my_slabel)
ax1.set_yticklabels(my_Ulabel[::-1])        
ax1.set_xlabel('Selection coefficient trait 1',multialignment='center',fontsize=18,labelpad=10)
ax1.set_ylabel('Mutation rate trait 2',multialignment='center',fontsize=18,labelpad=10)
ax1.tick_params(axis='both',labelsize=18)
# *****************************************************************************
#cbar.ax.text(2.8,0.6,'Ratio of $R_1/R_2$',rotation=270,fontsize=22)     # use this label of comparing R
cbar.ax.text(2.8,0.6,'Ratio of $v_1/v_2$',rotation=270,fontsize=22)    # use this label of comparing v
# *****************************************************************************

plt.text(1,2,r'$N = 10^9$',fontsize=18)
# *****************************************************************************
#plt.text(1,m-1.7,r'$R = 5.3\times 10^{-3}$',fontsize=18)     # use this label of comparing R
plt.text(1,1,r'$v = 5.3\times 10^{-5}$',fontsize=18)     # use this label of comparing v
# *****************************************************************************

#ax1.add_patch(Rectangle((someX - 0.1, someY - 0.1), 0.2, 0.2, alpha=1, facecolor='none'))

ax2 = ax1.twinx()
ax2.yaxis.set_ticks(np.arange(0+0.5/(len(my_Ulabel)),1+0.5/(len(my_Ulabel)),1.0/(len(my_Ulabel))))
ax2.set_yticklabels(my_slabel[::-1])
ax2.tick_params(labelsize=18)
ax2.set_ylabel('Selection coefficient trait 2',fontsize=18,labelpad=25,rotation=270)

ax3 = ax1.twiny()
ax3.xaxis.set_ticks(np.arange(0+0.5/(len(my_Ulabel)),1+0.5/(len(my_Ulabel)),1.0/(len(my_Ulabel))))
ax3.set_xticklabels(my_Ulabel)
ax3.tick_params(labelsize=18)
ax3.set_xlabel('Mutation rate trait 1',fontsize=18,labelpad=20)
#plt.tight_layout()

fig1.subplots_adjust(bottom=0.2,left=0.2)

plt.text(.2,-2.1,'Trait 1 favored by selection',fontsize=22)
plt.text(-0.213,11.5,'Trait 2 favored by mutation',rotation=90,fontsize=22)
        
# color x axis
ax1.annotate("",
            xy=(0,-.70), xycoords='data',
            xytext=(5,-0.70), textcoords='data',
            arrowprops=dict(arrowstyle="-",connectionstyle="arc3",color='lime',lw=5),
            annotation_clip=False)

ax1.annotate("",
            xy=(5,-.70), xycoords='data',
            xytext=(10,-.70), textcoords='data',
            arrowprops=dict(arrowstyle="-",connectionstyle="arc3",color='yellow',lw=5),
            annotation_clip=False)

ax1.annotate("",
            xy=(10,-.70), xycoords='data',
            xytext=(15,-.70), textcoords='data',
            arrowprops=dict(arrowstyle="-",connectionstyle="arc3",color='cyan',lw=5),
            annotation_clip=False)
            
# color y axis
ax1.annotate("",
            xy=(-1.8,0), xycoords='data',
            xytext=(-1.8,5), textcoords='data',
            arrowprops=dict(arrowstyle="-",connectionstyle="arc3",color='cyan',lw=5),
            annotation_clip=False)

ax1.annotate("",
            xy=(-1.8,5), xycoords='data',
            xytext=(-1.8,10), textcoords='data',
            arrowprops=dict(arrowstyle="-",connectionstyle="arc3",color='yellow',lw=5),
            annotation_clip=False)

ax1.annotate("",
            xy=(-1.8,10), xycoords='data',
            xytext=(-1.8,15), textcoords='data',
            arrowprops=dict(arrowstyle="-",connectionstyle="arc3",color='lime',lw=5),
            annotation_clip=False)

# add lines to separate regimes
ax1.annotate("",
            xy=(5,0), xycoords='data',
            xytext=(5,15), textcoords='data',
            arrowprops=dict(arrowstyle="-",connectionstyle="arc3",color='black',lw=1.5),
            annotation_clip=False)

ax1.annotate("",
            xy=(10,0), xycoords='data',
            xytext=(10,15), textcoords='data',
            arrowprops=dict(arrowstyle="-",connectionstyle="arc3",color='black',lw=1.5),
            annotation_clip=False)

ax1.annotate("",
            xy=(0,5), xycoords='data',
            xytext=(15,5), textcoords='data',
            arrowprops=dict(arrowstyle="-",connectionstyle="arc3",color='black',lw=1.5),
            annotation_clip=False)

ax1.annotate("",
            xy=(0,10), xycoords='data',
            xytext=(15,10), textcoords='data',
            arrowprops=dict(arrowstyle="-",connectionstyle="arc3",color='black',lw=1.5),
            annotation_clip=False)
            
fig1.savefig('figures/fig_two_trait_compare_v.pdf',bbox_inches='tight')
