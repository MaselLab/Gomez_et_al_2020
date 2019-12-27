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
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab

import fig_functions as myfun            # my functions in a seperate file
    
def get_heatmap_data(pickle_file_name):
    # load processed matlab data for figure panels
    # returns the relavant data for figures (see return statement)

    pickle_file = open("data/" + pickle_file_name,'rb') 
    [N,s,U,v,parameters,grand_means,sarry,Uarry,v1_data,v2_data] = pickle.load(pickle_file)
    pickle_file.close()    
    
    v1_min = 1
    v2_min = 1
    
    for i in range(v1_data.shape[0]):
        for j in range(v1_data.shape[1]):
            if (v1_data[i,j]>0):
                v1_min = np.min([v1_min,v1_data[i,j]])
            if (v2_data[i,j]>0):
                v2_min = np.min([v2_min,v2_data[i,j]])

    for i in range(v1_data.shape[0]):
        for j in range(v1_data.shape[1]):
            if (v1_data[i,j]<=0):
                v1_data[i,j] = v1_min
            if (v2_data[i,j]<=0):
                v2_data[i,j] = v2_min
                
    Uarry = np.flipud(Uarry)
    v_max = np.amax(v1_data)
    rate_comp = np.log10(v1_data/v_max)
#    v_tot = v1_data+v2_data
#    rate_comp = v1_data/v_tot

#    for i in range(rate_comp.shape[0]):
#        for j in range(rate_comp.shape[1]):
#            if (rate_comp[i,j]<=-4):
#                rate_comp[i,j] = -4
            
    log_s_min = np.log10(np.amin(sarry))
    log_s_max = np.log10(np.amax(sarry))
    log_U_min = np.log10(np.amin(Uarry))
    log_U_max = np.log10(np.amax(Uarry))
    log_s_crd = np.log10(parameters[0][1])    
    log_U_crd = np.log10(parameters[0][2])
    
    s1_coord = 31*(log_s_crd-log_s_min)/(log_s_max-log_s_min)
    U1_coord = 61*(log_U_crd-log_U_min)/(log_U_max-log_U_min)
        
    return [sarry, Uarry, rate_comp,s1_coord,U1_coord]

def get_vContourThry(N,s,v,log_s_lbd1,log_s_lbd2,log_s_lbd3,\
                    hlsmin,hlsmax,hlumin,hlumax,gridsize_s,gridsize_u):
    # computes tradeoff curves given parameters. Also, the coordinates have to
    # converted into those for the bounds given for the heatmap.
    #
    # inputs:
    # s = selection coefficient
    # N = Population size
    # v = rate of adaptation
    # log_s_lbd1 = lower bound s
    # log_s_lbd2 = transition to diff regime
    # log_s_lbd3 = upper bound s
    # hlsmin = lower bound s of heatmap log10
    # hlsmax = upper bound s of heatmap log10
    # hlumin = lower bound U of heatmap log10
    # hlumax = upper bound U of heatmap log10
    #
    # outputs:
    # curves for plots
    
    no_div,no_div1,no_div2 = [100,50,75]        # spacing between s points
    
    # setting bounds for the window and computing their log10 values for the log-plot
    [s_min,s_max,U_min,U_max,sc_max,sc_trans] = myfun.sU_bounds(N,v)          
    log10_s_min = np.log10(s_min)
    log10_s_max = np.log10(s_max)
    log10_U_min = np.log10(U_min)
    log10_U_max = np.log10(10*U_max)
    log10_sc_max = np.log10(sc_max)
    log10_sc_trans = np.log10(sc_trans)
    
    # Define range for s and U
    s1 = np.logspace(np.log10(s_min), np.log10(s_max), no_div)
    u1 = np.logspace(np.log10(U_min), np.log10(U_max), no_div)
    [log10_s1,log10_u1] = [np.log10(s1),np.log10(u1)]
    
    # special set of s values to help shade the drift barrier in sU space
    # s-thresholds and curves in phenotype space
    log10_sd = np.log10(np.logspace(np.log10(s_min), np.log10(20*s_min), no_div/10))
    
    # set s values between thresholds
    s_reg1 = np.logspace(log_s_lbd1,log10_sc_trans,no_div1)  # solid OF thry curve
    s_reg2 = np.logspace(log10_sc_trans,log10_s_max,no_div1) # dashed OF thry curve

    s_reg3 = np.logspace(log_s_lbd2,log10_sc_trans,no_div1)  # solid MM thry curve
    s_reg4 = np.logspace(log_s_lbd1,log10_s_max,no_div1)     # dashed MM thry curve
    
    s_reg5 = np.logspace(log_s_lbd2,log_s_lbd3,no_div1)     # solid OFMM thry curve
    s_reg6 = np.logspace(log_s_lbd1,log_s_lbd3,no_div1)     # dashed OFMM thry curve
    
    s_reg7 = np.logspace(log_s_lbd1,1.15*log_s_lbd2,no_div1)    # solid DM thry curve
    s_reg8 = np.logspace(log_s_lbd1,0.96*log_s_lbd2,no_div1)    # dashed DM thry curve
    
    # caluculate v-isoquant for successional regime
    vCont_OF1 = np.log10(np.asarray([[s_reg1[i],myfun.vContour_OF(s_reg1[i],N,v)] for i in range(no_div1)]))    
    vCont_OF2 = np.log10(np.asarray([[s_reg2[i],myfun.vContour_OF(s_reg2[i],N,v)] for i in range(no_div1)]))    
    
    # caluculate v-isoquant for concurrent regime
    vCont_MM1 = np.log10(np.asarray([[s_reg3[i],myfun.vContour_MM(s_reg3[i],N,v)] for i in range(no_div1)]))    
    vCont_MM2 = np.log10(np.asarray([[s_reg4[i],myfun.vContour_MM(s_reg4[i],N,v)] for i in range(no_div1)]))    

    # caluculate piecewise v-isoquant for combined regimes
    vCont_OFMM1 = np.log10(np.asarray([[s_reg5[i],myfun.vContour_OFMM(s_reg5[i],N,v)] for i in range(no_div1)]))    
    vCont_OFMM2 = np.log10(np.asarray([[s_reg6[i],myfun.vContour_OFMM(s_reg6[i],N,v)] for i in range(no_div1)]))    

    # caluculate v-isoquant for using halletschek approximations (Hallatschek 20011) 
    vCont_DM1 = np.log10(np.asarray([[s_reg7[i],myfun.vContour_DM(s_reg7[i],N,v)] for i in range(no_div1)]))    
    vCont_DM2 = np.log10(np.asarray([[s_reg8[i],myfun.vContour_DM(s_reg8[i],N,v)] for i in range(no_div1)]))    
    
    # convert all values to coordinates heatmap 
    
    log_s_min = np.log10(np.amin(sarry))
    log_s_max = np.log10(np.amax(sarry))
    log_U_min = np.log10(np.amin(Uarry))
    log_U_max = np.log10(np.amax(Uarry))
    
    vCont_MM1[:,0] = gridsize_s*(vCont_MM1[:,0]-hlsmin)/(hlsmax-hlsmin)
    vCont_MM1[:,1] = gridsize_u*(vCont_MM1[:,1]-hlumin)/(hlumax-hlumin)
    
    vCont_MM2[:,0] = gridsize_s*(vCont_MM2[:,0]-hlsmin)/(hlsmax-hlsmin)
    vCont_MM2[:,1] = gridsize_u*(vCont_MM2[:,1]-hlumin)/(hlumax-hlumin)
    
    vCont_OF1[:,0] = gridsize_s*(vCont_OF1[:,0]-hlsmin)/(hlsmax-hlsmin)
    vCont_OF1[:,1] = gridsize_u*(vCont_OF1[:,1]-hlumin)/(hlumax-hlumin)
    
    vCont_OF2[:,0] = gridsize_s*(vCont_OF2[:,0]-hlsmin)/(hlsmax-hlsmin)
    vCont_OF2[:,1] = gridsize_u*(vCont_OF2[:,1]-hlumin)/(hlumax-hlumin)
    
    vCont_OFMM1[:,0] = gridsize_s*(vCont_OFMM1[:,0]-hlsmin)/(hlsmax-hlsmin)
    vCont_OFMM1[:,1] = gridsize_u*(vCont_OFMM1[:,1]-hlumin)/(hlumax-hlumin)
    
    vCont_OFMM2[:,0] = gridsize_s*(vCont_OFMM2[:,0]-hlsmin)/(hlsmax-hlsmin)
    vCont_OFMM2[:,1] = gridsize_u*(vCont_OFMM2[:,1]-hlumin)/(hlumax-hlumin)
    
    vCont_DM1[:,0] = gridsize_s*(vCont_DM1[:,0]-hlsmin)/(hlsmax-hlsmin)
    vCont_DM1[:,1] = gridsize_u*(vCont_DM1[:,1]-hlumin)/(hlumax-hlumin)
    
    vCont_DM2[:,0] = gridsize_s*(vCont_DM2[:,0]-hlsmin)/(hlsmax-hlsmin)
    vCont_DM2[:,1] = gridsize_u*(vCont_DM2[:,1]-hlumin)/(hlumax-hlumin)
    
    return [vCont_MM1, vCont_MM2,vCont_OF1, \
            vCont_OF2,vCont_OFMM1,vCont_OFMM2,vCont_DM1,vCont_DM2]    

# set file name of data and load it into script
# note sarry and Uarry (values of s and U) are the same in all cases
pickle_file_name = 'fig_compare_sU_data-01a-DF.pickle'    # data for comparison with fixed v
[sarry,Uarry,rate_comp_MM,s1_coord_MM,U1_coord_MM]=get_heatmap_data(pickle_file_name)

pickle_file_name = 'fig_compare_sU_data-01a-OF.pickle'    # data for comparison with fixed v
[sarry,Uarry,rate_comp_OF,s1_coord_OF,U1_coord_OF]=get_heatmap_data(pickle_file_name)

pickle_file_name = 'fig_compare_sU_data-01a-HR.pickle'    # data for comparison with fixed v
[sarry,Uarry,rate_comp_DM,s1_coord_DM,U1_coord_DM]=get_heatmap_data(pickle_file_name)

# get rate_comp array dimensions 
[m,n] = rate_comp_MM.shape
[m,n] = [int(m),int(n)]
arry_dim_s = len(sarry)
arry_dim_u = len(Uarry)

hlsmin = np.log10(np.amin(sarry))
hlsmax = np.log10(np.amax(sarry))
hlumin = np.log10(np.amin(Uarry))
hlumax = np.log10(np.amax(Uarry))

# set labels for axes
my_slabel = ['$10^{'+str(np.round(np.log10(sarry[i,0]),2))+'}$' for i in range(len(sarry))]
my_Ulabel = ['$10^{'+str(np.round(np.log10(Uarry[i,0]),1))+'}$' for i in range(len(Uarry))]

# set some labels blank to have them fit on graph
for i in range(len(my_Ulabel)):
    if (i%10!=0):
        my_Ulabel[i]=''

for i in range(len(my_slabel)):
    if ((i%6!=3)):
        my_slabel[i]=''
        
x_border = [0.0+m*i/1000.0 for i in range(1001)]
y_border = [min(np.floor(1.0+m*i/1000.0),m) for i in range(1001)]

# -----------------------------------------------------------------------------
# create v-contours to add

[Nm,s,U] = [1e9,1e-2,1e-5]
gridsize_s = 31
gridsize_u = 61

vm = myfun.get_vDF(Nm,s,U)

[log_s_lbd1,log_s_lbd2,log_s_lbd3] = [-3.5,-2.3*0.9,-0.5]
         
[vCont_MM1, vCont_MM2,vCont_OF1,vCont_OF2, \
        vCont_OFMM1,vCont_OFMM2,vCont_DM1,vCont_DM2] \
        = get_vContourThry(Nm,s,vm,log_s_lbd1,log_s_lbd2,log_s_lbd3, \
                    hlsmin,hlsmax,hlumin,hlumax,gridsize_s,gridsize_u)

# -----------------------------------------------------------------------------
# create heatmaps of v reduction
fig = plt.figure(figsize=[7.5,17])

ax=plt.subplot(311)    # first panel-trait 1 in origin-fixation regime OF

#fit_distr_2d = ax.pcolormesh(rate_comp_OF,cmap=plt.cm.bwr)
fit_distr_2d = ax.pcolormesh(rate_comp_OF)
#cbar = plt.colorbar(fit_distr_2d,pad = 0.03,
#                    ticks=[0+i/10.0 for i in range(12)],
#                           norm=mpl.colors.Normalize(vmin=0.0, vmax=1.0))
cbar = plt.colorbar(fit_distr_2d,pad = 0.03,ticks=[0.0-i for i in range(9)])
                           
#cbar.set_clim(0.0, 1.0)
cbar.ax.set_yticklabels(['$1.00$','$0.10$','$0.01$'] + ['$10^{'+str(-3-i)+'}$' for i in range(6)])
cbar.ax.tick_params(labelsize=18)
ax.axis('tight')        
ax.set_xticks(np.arange(arry_dim_s)+0.5)
ax.set_yticks(np.arange(arry_dim_u)+0.5)        
ax.set_xticklabels([])
ax.set_yticklabels(my_Ulabel[::-1])        
ax.set_ylabel('Mutation rate trait 2',multialignment='center',fontsize=18,labelpad=10)
ax.tick_params(axis='both',labelsize=18)
#cbar.ax.text(3.8,0.80,'Ratio $\log_{10}(v_1/v)$',rotation=270,fontsize=22)    # use this label of comparing v
cbar.ax.text(4.4,0.65,'Ratio $v_1/v$',rotation=270,fontsize=22)    # use this label of comparing v
#cbar.ax.text(3.5,0.75,'Ratio $v_1/(v_1+v_2)$',rotation=270,fontsize=22)    # use this label of comparing v
plt.text(-8.0,60,'(a)',fontsize=20)

# plot isoquant calculated from theory on figures and location of s1,U1
ax.plot(vCont_OFMM1[:,0],vCont_OFMM1[:,1],color="black",linewidth=2,linestyle="-",label=r'$v=5.31\times 10^{-5}$')
ax.plot(vCont_OFMM2[:,0],vCont_OFMM2[:,1],color="black",linewidth=2,linestyle=":")
ax.plot(vCont_DM1[:,0],vCont_DM1[:,1],color="black",linewidth=2,linestyle="-")
ax.plot(vCont_DM2[:,0],vCont_DM2[:,1],color="black",linewidth=2,linestyle=":")
ax.scatter(s1_coord_OF,U1_coord_OF, facecolors='none',linewidth=2, edgecolors='k',s=80)


# *****************************************************************************
ax=plt.subplot(312)    # second panel-trait 1 in multiple mutations (U<s) MM

#fit_distr_2d = ax.pcolormesh(rate_comp_MM,cmap=plt.cm.bwr)
fit_distr_2d = ax.pcolormesh(rate_comp_MM)
#cbar = plt.colorbar(fit_distr_2d,pad = 0.03,
#                    ticks=[0+i/10.0 for i in range(12)],
#                           norm=mpl.colors.Normalize(vmin=0.0, vmax=1.0))
cbar = plt.colorbar(fit_distr_2d,pad = 0.03,ticks=[0.0-i for i in range(9)])
                           
#cbar.set_clim(0.0, 1.0)
cbar.ax.set_yticklabels(['$1.00$','$0.10$','$0.01$'] + ['$10^{'+str(-3-i)+'}$' for i in range(6)])
cbar.ax.tick_params(labelsize=18)
ax.axis('tight')        
ax.set_xticks(np.arange(arry_dim_s)+0.5)
ax.set_yticks(np.arange(arry_dim_u)+0.5)        
ax.set_xticklabels([])
ax.set_yticklabels(my_Ulabel[::-1])        
ax.set_ylabel('Mutation rate trait 2',multialignment='center',fontsize=18,labelpad=10)
ax.tick_params(axis='both',labelsize=18)
#cbar.ax.text(3.8,0.80,'Ratio $\log_{10}(v_1/v)$',rotation=270,fontsize=22)    # use this label of comparing v
cbar.ax.text(4.4,0.65,'Ratio $v_1/v$',rotation=270,fontsize=22)    # use this label of comparing v
#cbar.ax.text(3.5,0.75,'Ratio $v_1/(v_1+v_2)$',rotation=270,fontsize=22)    # use this label of comparing v
plt.text(-8.0,60,'(b)',fontsize=20)

# plot isoquant calculated from theory on figures and location of s1,U1 
ax.plot(vCont_OFMM1[:,0],vCont_OFMM1[:,1],color="black",linewidth=2,linestyle="-",label=r'$v=5.31\times 10^{-5}$')
ax.plot(vCont_OFMM2[:,0],vCont_OFMM2[:,1],color="black",linewidth=2,linestyle=":")
ax.plot(vCont_DM1[:,0],vCont_DM1[:,1],color="black",linewidth=2,linestyle="-")
ax.plot(vCont_DM2[:,0],vCont_DM2[:,1],color="black",linewidth=2,linestyle=":")
ax.scatter(s1_coord_MM,U1_coord_MM, facecolors='none',linewidth=2, edgecolors='k',s=80)


# *****************************************************************************
ax=plt.subplot(313)    # third panel-trait 1 in diffusive mutations regime (DM)

#fit_distr_2d = ax.pcolormesh(rate_comp_DM,cmap=plt.cm.bwr)
fit_distr_2d = ax.pcolormesh(rate_comp_DM)
#cbar = plt.colorbar(fit_distr_2d,pad = 0.03,
#                    ticks=[0+i/10.0 for i in range(12)],
#                           norm=mpl.colors.Normalize(vmin=0.0, vmax=1.0))
cbar = plt.colorbar(fit_distr_2d,pad = 0.03,ticks=[0.0-i for i in range(9)])
                           
#cbar.set_clim(0.0, 1.0)
cbar.ax.set_yticklabels(['$1.00$','$0.10$','$0.01$'] + ['$10^{'+str(-3-i)+'}$' for i in range(6)])
cbar.ax.tick_params(labelsize=18)
ax.axis('tight')        
ax.set_xticks(np.arange(arry_dim_s)+0.5)
ax.set_yticks(np.arange(arry_dim_u)+0.5)        
ax.set_xticklabels(my_slabel)
ax.set_yticklabels(my_Ulabel[::-1])        
ax.set_xlabel('Selection coefficient trait 2',multialignment='center',fontsize=18,labelpad=10)
ax.set_ylabel('Mutation rate trait 2',multialignment='center',fontsize=18,labelpad=10)
ax.tick_params(axis='both',labelsize=18)
#cbar.ax.text(3.8,0.80,'Ratio $\log_{10}(v_1/v)$',rotation=270,fontsize=22)    # use this label of comparing v
cbar.ax.text(4.4,0.65,'Ratio $v_1/v$',rotation=270,fontsize=22)    # use this label of comparing v
#cbar.ax.text(3.5,0.75,'Ratio $v_1/(v_1+v_2)$',rotation=270,fontsize=22)    # use this label of comparing v
plt.text(-8.0,60,'(c)',fontsize=20)

# plot isoquant calculated from theory on figures and location of s1,U1         
ax.plot(vCont_OFMM1[:,0],vCont_OFMM1[:,1],color="black",linewidth=2,linestyle="-",label=r'$v=5.31\times 10^{-5}$')
ax.plot(vCont_OFMM2[:,0],vCont_OFMM2[:,1],color="black",linewidth=2,linestyle=":")
ax.plot(vCont_DM1[:,0],vCont_DM1[:,1],color="black",linewidth=2,linestyle="-")
ax.plot(vCont_DM2[:,0],vCont_DM2[:,1],color="black",linewidth=2,linestyle=":")
ax.scatter(s1_coord_DM,U1_coord_DM, facecolors='none',linewidth=2, edgecolors='k',s=80)

plt.tight_layout()
fig.savefig('figures/fig_two_trait_compare_sU.pdf',bbox_inches='tight')

# -----------------------------------------------------------------------------
# OLD CODE THAT MIGHT BE USED
#plt.text(17.5,28.5,r'$N = 10^9$',fontsize=20)
#plt.text(17.5,26.5,r'$v = 5.3\times 10^{-5}$',fontsize=20)     # use this label of comparing v
#plt.text(23,28.5,r'$s_1 = 10^{-2}$',fontsize=18)
#plt.text(22.7,26.5,r'$U_1 = 10^{-5}$',fontsize=18)     # use this label of comparing v
#ax.plot(x_border,y_border,color="black")
#fig2.subplots_adjust(bottom=0.2,left=0.2)            