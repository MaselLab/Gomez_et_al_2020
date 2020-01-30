#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 14:15:28 2019
Masel Lab
Project: Mutation-driven Adaptation
@author: Kevin Gomez

Description:
Script for creating phase plot, and sample trajectories, for the relationship
between U and s that preserves the total rate of adaptation given a population
size N.
"""

#libraries
import matplotlib.pyplot as plt
import numpy as np
import fig_functions as myfun
import pickle

import fig_functions as myfun            # my functions in a seperate file

## uncomment if you make changes to fig_functions and run these lines to reload
#import fig_functions
#reload(fig_functions)

# functions
# -----------------------------------------------------------------------------

def get_graph_data(N,s,v,log_s_lbd1,log_s_lbd2,log_s_lbd3):
    # computes tradeoff curves given parameters, function is only used here.
    #
    # inputs:
    # s = selection coefficient
    # N = Population size
    # v = rate of adaptation
    # log_s_lbd1 = lower bound s
    # log_s_lbd2 = transition to diff regime
    # log_s_lbd3 = upper bound s
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
    s_reg1 = np.logspace(log_s_lbd1,log10_sc_trans,no_div1)
    s_reg2 = np.logspace(log10_sc_trans,log10_s_max,no_div1)
    s_reg3 = np.logspace(log_s_lbd2,log10_sc_trans,no_div1)
    s_reg4 = np.logspace(log_s_lbd1,log10_s_max,no_div1)
    s_reg5 = np.logspace(log_s_lbd2,log_s_lbd3,no_div1)
    s_reg6 = np.logspace(log_s_lbd1,log_s_lbd3,no_div1)
    
    s_reg7 = np.logspace(log_s_lbd1,1.15*log_s_lbd2,no_div1)
    s_reg8 = np.logspace(log_s_lbd1,1.1*log_s_lbd2,no_div1)
    
    
    # successional-concurrent barrier (same in pheno and genotype sU space)
    succ_shade = np.log10(np.asarray([[s1[i],U_min,10*U_max] for i in range(no_div)]))
    conc_shade = np.log10(np.asarray([[s1[i],myfun.succ_conc_barrier(s1[i],N,v),10*U_max] for i in range(no_div)]))
    disc_shade = np.log10(np.asarray([[s1[i],min(0.1*s1[i],10*U_max),10*U_max] for i in range(no_div)]))
    
    # caluculate v-isoquant for concurrent regime
    vCont_MM1 = np.log10(np.asarray([[s_reg3[i],myfun.vContour_MM(s_reg3[i],N,v)] for i in range(no_div1)]))    
    vCont_MM2 = np.log10(np.asarray([[s_reg4[i],myfun.vContour_MM(s_reg4[i],N,v)] for i in range(no_div1)]))    
    
    # caluculate v-isoquant for successional regime
    vCont_OF1 = np.log10(np.asarray([[s_reg1[i],myfun.vContour_OF(s_reg1[i],N,v)] for i in range(no_div1)]))    
    vCont_OF2 = np.log10(np.asarray([[s_reg2[i],myfun.vContour_OF(s_reg2[i],N,v)] for i in range(no_div1)]))    

    # caluculate piecewise v-isoquant for combined regimes
    vCont_OFMM1 = np.log10(np.asarray([[s_reg5[i],myfun.vContour_OFMM(s_reg5[i],N,v)] for i in range(no_div1)]))    
    vCont_OFMM2 = np.log10(np.asarray([[s_reg6[i],myfun.vContour_OFMM(s_reg6[i],N,v)] for i in range(no_div1)]))    

    # caluculate v-isoquant for using halletschek approximations (Hallatschek 20011) 
    vCont_DM1 = np.log10(np.asarray([[s_reg7[i],myfun.vContour_DM(s_reg7[i],N,v)] for i in range(no_div1)]))    
    vCont_DM2 = np.log10(np.asarray([[s_reg8[i],myfun.vContour_DM(s_reg8[i],N,v)] for i in range(no_div1)]))    
    
    return [succ_shade,conc_shade,disc_shade,vCont_MM1, vCont_MM2,vCont_OF1, \
            vCont_OF2,vCont_OFMM1,vCont_OFMM2,vCont_DM1,vCont_DM2]

# set file names with saved data
# -----------------------------------------------------------------------------

#my_saved_data = 'data/fig_sUtradeoff_simdata-01.pickle'
#my_saved_data = 'data/fig_sUtradeoff_simdata-06.pickle'
#my_saved_data = 'data/fig_sUtradeoff_simdata-22.pickle'
#my_saved_data = 'data/fig_sUtradeoff_simdata-100.pickle'
my_saved_data = 'data/fig_sUtradeoff_simdata-101.pickle'

# load estimates for v isoquants from simulations
# -----------------------------------------------------------------------------
pickle_file = open(my_saved_data,'rb') 
[Nv_param,sU_data] = pickle.load(pickle_file)
pickle_file.close()

for i in range(6):
    sU_data[i]=np.log10(sU_data[i])

# -----------------------------------------------------------------------------
#               v isoquants from stochastic approximations
# -----------------------------------------------------------------------------

# set basic parameters of the figure for varying v
# -----------------------------------------------------------------------------
[Nl,Nm,Nh,s,U] = [1e6,1e9,1e18,1e-2,1e-5]

vm = myfun.get_vDF(Nm,s,U)
vl = vm*10**(-1)
vh = vm*10**(1)

[log_s_lbd1,log_s_lbd2,log_s_lbd3] = [-3.5,-2.3,-0.5]

[succ_sh_vl,conc_sh_vl,disc_sh_vl,c_curve1_vl,c_curve2_vl,s_curve1_vl, \
     s_curve2_vl,sU_curve1_vl,sU_curve2_vl,sU_curve1h_vl,sU_curve2h_vl] \
     = get_graph_data(Nm,s,vl,log_s_lbd1,1.15*log_s_lbd2,log_s_lbd3)
         
[succ_sh_vm,conc_sh_vm,disc_sh_vm,c_curve1_vm,c_curve2_vm,s_curve1_vm, \
     s_curve2_vm,sU_curve1_vm,sU_curve2_vm,sU_curve1h_vm,sU_curve2h_vm] \
     = get_graph_data(Nm,s,vm,log_s_lbd1,1.0*log_s_lbd2,log_s_lbd3)

[succ_sh_vh,conc_sh_vh,disc_sh_vh,c_curve1_vh,c_curve2_vh,s_curve1_vh, \
     s_curve2_vh,sU_curve1_vh,sU_curve2_vh,sU_curve1h_vh,sU_curve2h_vh] \
     = get_graph_data(Nm,s,vh,log_s_lbd1,0.85*log_s_lbd2,log_s_lbd3)

[s_min,s_max,U_min,U_max,sc_max,sc_trans] = myfun.sU_bounds(Nm,vm)
[log10_s_min,log10_s_max,log10_U_min,log10_U_max,log10_sc_max,log10_sc_trans] = \
    [np.log10(s_min),np.log10(s_max),np.log10(U_min),np.log10(10*U_max),np.log10(sc_max),np.log10(sc_trans)]

# create figure 1    
# -----------------------------------------------------------------------------
fig = plt.figure(figsize=[7,11])

#fig1, ax1 = plt.subplots(2,1,figsize=[16,8])

ax1=plt.subplot(211)

# set colors identifying regions
# -----------------------------------------------------------------------------
ax1.fill_between(succ_sh_vm[:,0],succ_sh_vm[:,1],succ_sh_vm[:,2],facecolor="cyan")
ax1.fill_between(conc_sh_vm[:,0],conc_sh_vm[:,1],conc_sh_vm[:,2],facecolor="yellow")
ax1.fill_between(disc_sh_vm[:,0],disc_sh_vm[:,1],disc_sh_vm[:,2],facecolor="lime")

# plot three isoquants calculated from theory
# -----------------------------------------------------------------------------
ax1.plot(sU_curve1_vl[:,0],sU_curve1_vl[:,1],color="blue",linewidth=2,linestyle="-",label=r'$v=5.31\times 10^{-6}$')
ax1.plot(sU_curve2_vl[:,0],sU_curve2_vl[:,1],color="blue",linewidth=2,linestyle=":")
ax1.plot(sU_curve1_vm[:,0],sU_curve1_vm[:,1],color="purple",linewidth=2,linestyle="-",label=r'$v=5.31\times 10^{-5}$')
ax1.plot(sU_curve2_vm[:,0],sU_curve2_vm[:,1],color="purple",linewidth=2,linestyle=":")
ax1.plot(sU_curve1_vh[:,0],sU_curve1_vh[:,1],color="red",linewidth=2,linestyle="-",label=r'$v=5.31\times 10^{-4}$')
ax1.plot(sU_curve2_vh[:,0],sU_curve2_vh[:,1],color="red",linewidth=2,linestyle=":")

# plot three isoquants calculated from theory with Hallatschek
# -----------------------------------------------------------------------------
ax1.plot(sU_curve1h_vl[:,0],sU_curve1h_vl[:,1],color="blue",linewidth=2,linestyle="-")
ax1.plot(sU_curve2h_vl[:,0],sU_curve2h_vl[:,1],color="blue",linewidth=2,linestyle=":")
ax1.plot(sU_curve1h_vm[:,0],sU_curve1h_vm[:,1],color="purple",linewidth=2,linestyle="-")
ax1.plot(sU_curve2h_vm[:,0],sU_curve2h_vm[:,1],color="purple",linewidth=2,linestyle=":")
ax1.plot(sU_curve1h_vh[:,0],sU_curve1h_vh[:,1],color="red",linewidth=2,linestyle="-")
ax1.plot(sU_curve2h_vh[:,0],sU_curve2h_vh[:,1],color="red",linewidth=2,linestyle=":")

# plot simulated data points of sU tradeoff concurrent/successional curves
# -----------------------------------------------------------------------------
ax1.scatter(sU_data[0][:,0],sU_data[0][:,1],color="blue",linewidth=2)
ax1.scatter(sU_data[1][:,0],sU_data[1][:,1],color="purple",linewidth=2)
ax1.scatter(sU_data[2][3:,0],sU_data[2][3:,1],color="red",linewidth=2)

# set figure dimensions and labels
# -----------------------------------------------------------------------------
new_xtick_labels=['',r'$10^{-3.5}$',r'$10^{-3.0}$',r'$10^{-2.5}$',r'$10^{-2.0}$',r'$10^{-1.5}$',r'$10^{-1.0}$']
new_ytick_labels=['',r'$10^{-12}$',r'$10^{-10}$',r'$10^{-8}$',r'$10^{-6}$',r'$10^{-4}$',r'$10^{-2}$']

ax1.set_xlim([1.2*log10_sc_max,log10_s_max])
ax1.set_ylim([log10_U_min,log10_U_max])
ax1.set_xticklabels([])
ax1.set_yticklabels(new_ytick_labels)

#ax1.set_xlabel(r'Selection coefficient',fontsize=18,labelpad=20)
ax1.set_ylabel(r'Beneficial mutation rate',fontsize=18,labelpad=8)

#locs, labels = xticks()            # Get locations and labels
#xticks(ticks, [labels], **kwargs)  # Set locations and labels

ax1.tick_params(labelsize=18)
ax1.legend(loc=3,fontsize=18)

# annotations to graphs
# -----------------------------------------------------------------------------
xh_loc = (log10_s_max-1.2*log10_sc_max)
yh_loc = (log10_U_max-log10_U_min)
plt.text(1.1*log10_sc_max-0.082*xh_loc,log10_U_min+0.30*yh_loc,r'$N = 10^9$',fontsize=22)
plt.text(1.15*log10_sc_max,0.55*log10_U_min,"Multiple Mutations\n Regime (U << s)",fontsize=16)
plt.text(0.72*log10_sc_max,0.93*log10_U_min,"Origin-fixation\n     Regime",fontsize=16)
plt.text(0.67*log10_sc_max,0.13*log10_U_min,"Diffusive Mutations\n        Regime",color="black",fontsize=16)
plt.text(1.2*log10_sc_max-.16*xh_loc,log10_U_min+.96*yh_loc,'(a)',fontsize=20)

#plt.close()

# save figure
#fig1.savefig('figures/fig_v_isoquants_vary_v.pdf')

#--------------------------------------------------------------------------------
#           Figure 2 - v-isoquants with varying N 
#--------------------------------------------------------------------------------

# set 
[Nl,Nm,Nh,s,U] = [1e7,1e9,1e11,1e-2,1e-5]
[succ_sh_Nl,conc_sh_Nl,disc_sh_Nl,c_curve1_Nl,c_curve2_Nl,s_curve1_Nl,s_curve2_Nl, \
                             sU_curve1_Nl,sU_curve2_Nl,sU_curve1h_Nl,sU_curve2h_Nl] \
                             = get_graph_data(Nl,s,vm,log_s_lbd1,0.9*log_s_lbd2,log_s_lbd3)
                             
[succ_sh_Nm,conc_sh_Nm,disc_sh_Nm,c_curve1_Nm,c_curve2_Nm,s_curve1_Nm,s_curve2_Nm, \
                             sU_curve1_Nm,sU_curve2_Nm,sU_curve1h_Nm,sU_curve2h_Nm] \
                             = get_graph_data(Nm,s,vm,log_s_lbd1,0.9*log_s_lbd2,log_s_lbd3)
                             
[succ_sh_Nh,conc_sh_Nh,disc_sh_Nh,c_curve1_Nh,c_curve2_Nh,s_curve1_Nh,s_curve2_Nh, \
                             sU_curve1_Nh,sU_curve2_Nh,sU_curve1h_Nh,sU_curve2h_Nh] \
                             = get_graph_data(Nh,s,vm,log_s_lbd1,0.9*log_s_lbd2,log_s_lbd3)

[s_min,s_max,U_min,U_max,sc_max,sc_trans] = myfun.sU_bounds(Nm,vm)
[log10_s_min,log10_s_max,log10_U_min,log10_U_max,log10_sc_max,log10_sc_trans] = \
    [np.log10(s_min),np.log10(s_max),np.log10(U_min),np.log10(10*U_max),np.log10(sc_max),np.log10(sc_trans)]
    
# create figure 2    
# -----------------------------------------------------------------------------
#fig2, ax2 = plt.subplots(1,1,figsize=[8,8])
ax2=plt.subplot(212)

# plot three isoquants calculated from theory
# -----------------------------------------------------------------------------
ax2.plot(sU_curve1_Nl[:,0],sU_curve1_Nl[:,1],color="mediumblue",linewidth=2,linestyle="-",label=r'$N=10^7$')
ax2.plot(sU_curve2_Nl[:,0],sU_curve2_Nl[:,1],color="mediumblue",linewidth=2,linestyle=":")
ax2.plot(sU_curve1_Nm[:,0],sU_curve1_Nm[:,1],color="purple",linewidth=2,linestyle="-",label=r'$N=10^9$')
ax2.plot(sU_curve2_Nm[:,0],sU_curve2_Nm[:,1],color="purple",linewidth=2,linestyle=":")
ax2.plot(sU_curve1_Nh[:,0],sU_curve1_Nh[:,1],color="red",linewidth=2,linestyle="-",label=r'$N=10^{11}$')
ax2.plot(sU_curve2_Nh[:,0],sU_curve2_Nh[:,1],color="red",linewidth=2,linestyle=":")

## plot three isoquants calculated from theory with Hallatschek
## -----------------------------------------------------------------------------
#ax2.plot(sU_curve1h_Nl[:,0],sU_curve1h_Nl[:,1],color="blue",linewidth=2,linestyle="-")
#ax2.plot(sU_curve2h_Nl[:,0],sU_curve2h_Nl[:,1],color="blue",linewidth=2,linestyle=":")
#ax2.plot(sU_curve1h_Nm[:,0],sU_curve1h_Nm[:,1],color="purple",linewidth=2,linestyle="-")
#ax2.plot(sU_curve2h_Nm[:,0],sU_curve2h_Nm[:,1],color="purple",linewidth=2,linestyle=":")
#ax2.plot(sU_curve1h_Nh[:,0],sU_curve1h_Nh[:,1],color="red",linewidth=2,linestyle="-")
#ax2.plot(sU_curve2h_Nh[:,0],sU_curve2h_Nh[:,1],color="red",linewidth=2,linestyle=":")

# plot estimates of isoquants from stochastic approximation
# -----------------------------------------------------------------------------
ax2.scatter(sU_data[3][:,0],sU_data[3][:,1],color="mediumblue",linewidth=2)
ax2.scatter(sU_data[4][:,0],sU_data[4][:,1],color="purple",linewidth=2)
ax2.scatter(sU_data[5][:,0],sU_data[5][:,1],color="red",linewidth=2)

# set figure dimensions and add labels
# -----------------------------------------------------------------------------

new_xtick_labels=['',r'$10^{-3.5}$',r'$10^{-3.0}$',r'$10^{-2.5}$',r'$10^{-2.0}$',r'$10^{-1.5}$',r'$10^{-1.0}$']
new_ytick_labels=['',r'$10^{-12}$',r'$10^{-10}$',r'$10^{-8}$',r'$10^{-6}$',r'$10^{-4}$',r'$10^{-2}$']

ax2.set_xlim([1.2*log10_sc_max,log10_s_max])
ax2.set_ylim([log10_U_min,log10_U_max])
ax2.set_xlabel(r'Selection coefficient',fontsize=18,labelpad=20)
ax2.set_ylabel(r'Beneficial mutation rate',fontsize=18,labelpad=8)
ax2.set_xticklabels(new_xtick_labels)
ax2.set_yticklabels(new_ytick_labels)

ax2.tick_params(labelsize=18)
ax2.legend(loc=3,fontsize=18)

# add annotations
# -----------------------------------------------------------------------------
xh_loc = (log10_s_max-1.2*log10_sc_max)
yh_loc = (log10_U_max-log10_U_min)
plt.text(1.2*log10_sc_max+0.025*xh_loc,log10_U_min+0.30*yh_loc,r'$v = 5.3\times 10^{-5}$',fontsize=22)
plt.text(1.2*log10_sc_max-.16*xh_loc,log10_U_min+.96*yh_loc,'(b)',fontsize=20)
plt.tight_layout()
#plt.close()

fig.savefig('figures/fig_v_isoquants.pdf')

# -----------------------------------------------------------------------------
#                           OLD CODE
# -----------------------------------------------------------------------------
# piecewise concurrent/successional curves
#ax1.plot(c_curve1_vl[:,0],c_curve1_vl[:,1],color="mediumblue",linewidth=2,linestyle="-",label="Concurrent")
#ax1.plot(c_curve2_vl[:,0],c_curve2_vl[:,1],color="mediumblue",linewidth=2,linestyle=":")
#ax1.plot(s_curve1_vl[:,0],s_curve1_vl[:,1],color="red",linewidth=2,linestyle=":")
#ax1.plot(s_curve2_vl[:,0],s_curve2_vl[:,1],color="red",linewidth=2,linestyle="-",label="Origin-fixation")
#ax1.plot(c_curve1_vm[:,0],c_curve1_vm[:,1],color="mediumblue",linewidth=2,linestyle="-",label="Concurrent")
#ax1.plot(c_curve2_vm[:,0],c_curve2_vm[:,1],color="mediumblue",linewidth=2,linestyle=":")
#ax1.plot(s_curve1_vm[:,0],s_curve1_vm[:,1],color="red",linewidth=2,linestyle=":")
#ax1.plot(s_curve2_vm[:,0],s_curve2_vm[:,1],color="red",linewidth=2,linestyle="-",label="Origin-fixation")
#ax1.plot(c_curve1_vh[:,0],c_curve1_vh[:,1],color="mediumblue",linewidth=2,linestyle="-",label="Concurrent")
#ax1.plot(c_curve2_vh[:,0],c_curve2_vh[:,1],color="mediumblue",linewidth=2,linestyle=":")
#ax1.plot(s_curve1_vh[:,0],s_curve1_vh[:,1],color="red",linewidth=2,linestyle=":")
#ax1.plot(s_curve2_vh[:,0],s_curve2_vh[:,1],color="red",linewidth=2,linestyle="-",label="Origin-fixation")
