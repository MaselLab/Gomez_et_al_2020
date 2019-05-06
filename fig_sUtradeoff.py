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

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# functions
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

def sU_bounds(N,v):
    # computes appropriate bounds of sU space for looking at tradeoff with 
    # phenotype adaptation
    #
    # inputs:
    # s = selection coefficient
    # N = Population size
    # v = rate of adaptation
    #
    # outputs:
    # list of bounds for s and U

    sm = 1.1*myfun.s_max_Uc(N,v)                      # approx s for max Uc
    st = 1.2*myfun.s_transition2succ(N,v)             # approx s for transition conc to succ
    
    s_min = 1/N
    s_max = 5*st     
    U_min = 0.01*myfun.sU_tradeoff_succ(s_max,N,v)
    U_max = 100*myfun.sU_tradeoff_conc(sm,N,v)
    
    return [s_min,s_max,U_min,U_max,sm,st]

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
    
def get_graph_data(N,s,v,log_s_lbd1,log_s_lbd2,log_s_lbd3):
    # computes tradeoff curves given parameters
    #
    # inputs:
    # s = selection coefficient
    # N = Population size
    # v = rate of adaptation
    #
    # outputs:
    # curves for plots
    
    no_div,no_div1,no_div2 = [100,50,75]
    
    # setting bounds for the window and computing their log10 values for the log-plot
    [s_min,s_max,U_min,U_max,sc_max,sc_trans] = sU_bounds(N,v)
    [log10_s_min,log10_s_max,log10_U_min,log10_U_max,log10_sc_max,log10_sc_trans] = \
        [np.log10(s_min),np.log10(s_max),np.log10(U_min),np.log10(10*U_max),np.log10(sc_max),np.log10(sc_trans)]
    
    # Define range for s and U
    s1 = np.logspace(np.log10(s_min), np.log10(s_max), no_div)
    u1 = np.logspace(np.log10(U_min), np.log10(U_max), no_div)
    [log10_s1,log10_u1] = [np.log10(s1),np.log10(u1)]
    
    # special set of s values to help shade the drift barrier in sU space
    # s-thresholds and curves in phenotype space
    log10_sd = np.log10(np.logspace(np.log10(s_min), np.log10(20*s_min), no_div/10))
    
    # pick s values between thresholds
    s_reg1 = np.logspace(log_s_lbd1,log10_sc_trans,no_div1)
    s_reg2 = np.logspace(log10_sc_trans,log10_s_max,no_div1)
    s_reg3 = np.logspace(log_s_lbd2,log10_sc_trans,no_div1)
    s_reg4 = np.logspace(log_s_lbd1,log10_s_max,no_div1)
    s_reg5 = np.logspace(log_s_lbd2,log_s_lbd3,no_div1)
    s_reg6 = np.logspace(log_s_lbd1,log_s_lbd3,no_div1)
    
    # successional-concurrent barrier (same in pheno and genotype sU space)
    succ_shade = np.log10(np.asarray([[s1[i],U_min,10*U_max] for i in range(no_div)]))
    conc_shade = np.log10(np.asarray([[s1[i],myfun.succ_conc_barrier(s1[i],N,v),10*U_max] for i in range(no_div)]))
    disc_shade = np.log10(np.asarray([[s1[i],min(0.1*s1[i],10*U_max),10*U_max] for i in range(no_div)]))
    
    # caluculate U along sU tradeoff for concurrent regime
    sU_tradeoff_conc_curve1 = np.log10(np.asarray([[s_reg3[i],myfun.sU_tradeoff_conc(s_reg3[i],N,v)] for i in range(no_div1)]))    
    sU_tradeoff_conc_curve2 = np.log10(np.asarray([[s_reg4[i],myfun.sU_tradeoff_conc(s_reg4[i],N,v)] for i in range(no_div1)]))    
    
    # caluculate U along sU tradeoff for successional regime
    sU_tradeoff_succ_curve1 = np.log10(np.asarray([[s_reg1[i],myfun.sU_tradeoff_succ(s_reg1[i],N,v)] for i in range(no_div1)]))    
    sU_tradeoff_succ_curve2 = np.log10(np.asarray([[s_reg2[i],myfun.sU_tradeoff_succ(s_reg2[i],N,v)] for i in range(no_div1)]))    

    # caluculate U along sU tradeoff for combined regimes
    sU_curve1 = np.log10(np.asarray([[s_reg5[i],myfun.sU_tradeoff(s_reg5[i],N,v)] for i in range(no_div1)]))    
    sU_curve2 = np.log10(np.asarray([[s_reg6[i],myfun.sU_tradeoff(s_reg6[i],N,v)] for i in range(no_div1)]))    
    
    return [succ_shade,conc_shade,disc_shade,sU_tradeoff_conc_curve1, \
            sU_tradeoff_conc_curve2,sU_tradeoff_succ_curve1, \
            sU_tradeoff_succ_curve2,sU_curve1,sU_curve2]


# -----------------------------------------------------------------------------    
# load data for simulation estimates
# -----------------------------------------------------------------------------

# simulation estimates for N-medium, v-medium
pickle_file = open('data/fig_sUtradeoff_simdata-01.pickle','rb') 
[Nv_param,sU_data] = pickle.load(pickle_file)
pickle_file.close()

for i in range(6):
    sU_data[i]=np.log10(sU_data[i])

#indx = [4,6,8,10,12,14,16,18,20,22,24,26,29,31,33,36,38,40]
#sU_pair_log = np.log10(sU_pair)
#sU_comp = []
#
#for i in range(len(indx)):
#    sU_comp += [[np.log10(sU_pair[indx[i],0]),np.log10(sU_pair[indx[i],1])]]
#
#sU_comp = np.asarray(sU_comp)
#
## manual estimates
#sU_comp = np.asarray([[-3.2,-0.975],[-3.1,-1.175],[-3,-1.375],
#                      [-2.9,-1.58],[-2.8,-1.81],[-2.7,-2.04],[-2.6,-2.32],[-2.5,-2.65],
#                      [-2.4,-2.98],[-2.3,-3.375],[-2.2,-3.765],[-2.1,-4.275],[-2,-4.93],
#                      [-1.9,-5.54],[-1.8,-6.429],[-1.7,-7.34],[-1.6,-8.32],[-1.5,-9.4]])

# -----------------------------------------------------------------------------
# set basic parameters of the figure for varying v
# -----------------------------------------------------------------------------

[Nl,Nm,Nh,s,U] = [1e6,1e9,1e18,1e-2,1e-5]
vl = myfun.get_vDF(Nl,s,U)
vm = myfun.get_vDF(Nm,s,U)
vh = myfun.get_vDF(Nh,s,U)

[log_s_lbd1,log_s_lbd2,log_s_lbd3] = [-3.5,-2.4,-0.5]
[succ_sh_vl,conc_sh_vl,disc_sh_vl,c_curve1_vl,c_curve2_vl,s_curve1_vl,s_curve2_vl,sU_curve1_vl,sU_curve2_vl] = get_graph_data(Nm,s,vl,log_s_lbd1,log_s_lbd2,log_s_lbd3)
[succ_sh_vm,conc_sh_vm,disc_sh_vm,c_curve1_vm,c_curve2_vm,s_curve1_vm,s_curve2_vm,sU_curve1_vm,sU_curve2_vm] = get_graph_data(Nm,s,vm,log_s_lbd1,0.95*log_s_lbd2,log_s_lbd3)
[succ_sh_vh,conc_sh_vh,disc_sh_vh,c_curve1_vh,c_curve2_vh,s_curve1_vh,s_curve2_vh,sU_curve1_vh,sU_curve2_vh] = get_graph_data(Nm,s,vh,log_s_lbd1,0.85*log_s_lbd2,log_s_lbd3)

[s_min,s_max,U_min,U_max,sc_max,sc_trans] = sU_bounds(Nm,vm)
[log10_s_min,log10_s_max,log10_U_min,log10_U_max,log10_sc_max,log10_sc_trans] = \
    [np.log10(s_min),np.log10(s_max),np.log10(U_min),np.log10(10*U_max),np.log10(sc_max),np.log10(sc_trans)]

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# sU-tradeoff for phenotpe adaptation
# plot that includes thresholds and sU-tradeoff curve
# -----------------------------------------------------------------------------
    
fig1, ax1 = plt.subplots(1,1,figsize=[8,8])
ax1.fill_between(succ_sh_vm[:,0],succ_sh_vm[:,1],succ_sh_vm[:,2],facecolor="deepskyblue")
ax1.fill_between(conc_sh_vm[:,0],conc_sh_vm[:,1],conc_sh_vm[:,2],facecolor="gold")
ax1.fill_between(disc_sh_vm[:,0],disc_sh_vm[:,1],disc_sh_vm[:,2],facecolor="limegreen")

ax1.plot(sU_curve1_vl[:,0],sU_curve1_vl[:,1],color="mediumblue",linewidth=2,linestyle="-",label='v='+'%e' % vl)
ax1.plot(sU_curve2_vl[:,0],sU_curve2_vl[:,1],color="mediumblue",linewidth=2,linestyle=":")
ax1.plot(sU_curve1_vm[:,0],sU_curve1_vm[:,1],color="purple",linewidth=2,linestyle="-",label='v='+'%e' % vm)
ax1.plot(sU_curve2_vm[:,0],sU_curve2_vm[:,1],color="purple",linewidth=2,linestyle=":")
ax1.plot(sU_curve1_vh[:,0],sU_curve1_vh[:,1],color="red",linewidth=2,linestyle="-",label='v='+'%e' % vh)
ax1.plot(sU_curve2_vh[:,0],sU_curve2_vh[:,1],color="red",linewidth=2,linestyle=":")

# plot of piecewise concurrent/successional curves
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

# plot of piecewise concurrent/successional curves
ax1.scatter(sU_data[0][:,0],sU_data[0][:,1],color="mediumblue",marker='.')
ax1.scatter(sU_data[1][:,0],sU_data[1][:,1],color="purple",marker='.')
ax1.scatter(sU_data[2][:,0],sU_data[2][:,1],color="red",marker='.')
#ax1.scatter(sU_comp[:,0],sU_comp[:,1],color="black",linewidth=4)

# set figure dimensions
ax1.set_xlim([1.2*log10_sc_max,log10_s_max])
ax1.set_ylim([log10_U_min,log10_U_max])
ax1.set_xlabel(r'Selection coefficient ($\log_{10}s$)',fontsize=18,labelpad=20)
ax1.set_ylabel(r'Mutation rate ($\log_{10}U$)',fontsize=18,labelpad=10)
ax1.tick_params(labelsize=16)
ax1.legend(loc=3)

# annotations to graphs
xh_loc = (log10_s_max-1.2*log10_sc_max)
yh_loc = (log10_U_max-log10_U_min)

plt.text(1.1*log10_sc_max+0.70*xh_loc,log10_U_min+0.65*yh_loc,r'$N = 10^9$',fontsize=16)
#plt.text(1.2*log10_sc_max+0.70*xh_loc,log10_U_min+0.60*yh_loc,r'$v = 5.3\times 10^{-5}$',fontsize=16)
plt.text(1.0*log10_sc_max,0.55*log10_U_min,"Traveling Wave\n    Regime",fontsize=16)
plt.text(0.75*log10_sc_max,0.93*log10_U_min,"Origin-fixation\n     Regime",fig1ontsize=16)
plt.text(0.85*log10_sc_max,0.13*log10_U_min,"Discontinuous\n     Regime",color="black",fontsize=16)

fig1.savefig('figures/fig_sUtradeoff_pheno_adapt_v.pdf')

#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# selecting parameters for sU tradeoff curves with varying N

[Nl,Nm,Nh,s,U] = [1e7,1e9,1e11,1e-2,1e-5]
[succ_sh_Nl,conc_sh_Nl,disc_sh_Nl,c_curve1_Nl,c_curve2_Nl,s_curve1_Nl,s_curve2_Nl,sU_curve1_Nl,sU_curve2_Nl] = get_graph_data(Nl,s,vm,log_s_lbd1,log_s_lbd2)
[succ_sh_Nm,conc_sh_Nm,disc_sh_Nm,c_curve1_Nm,c_curve2_Nm,s_curve1_Nm,s_curve2_Nm,sU_curve1_Nm,sU_curve2_Nm] = get_graph_data(Nm,s,vm,log_s_lbd1,log_s_lbd2)
[succ_sh_Nh,conc_sh_Nh,disc_sh_Nh,c_curve1_Nh,c_curve2_Nh,s_curve1_Nh,s_curve2_Nh,sU_curve1_Nh,sU_curve2_Nh] = get_graph_data(Nh,s,vm,log_s_lbd1,log_s_lbd2)

[s_min,s_max,U_min,U_max,sc_max,sc_trans] = sU_bounds(Nm,vm)
[log10_s_min,log10_s_max,log10_U_min,log10_U_max,log10_sc_max,log10_sc_trans] = \
    [np.log10(s_min),np.log10(s_max),np.log10(U_min),np.log10(10*U_max),np.log10(sc_max),np.log10(sc_trans)]
    
# -----------------------------------------------------------------------------
# sU-tradeoff for phenotpe adaptation
# plot that includes thresholds and sU-tradeoff curve
# -----------------------------------------------------------------------------
    
fig2, ax2 = plt.subplots(1,1,figsize=[8,8])

ax2.plot(sU_curve1_Nl[:,0],sU_curve1_Nl[:,1],color="mediumblue",linewidth=2,linestyle="-",label='N='+'%e' % Nl)
ax2.plot(sU_curve2_Nl[:,0],sU_curve2_Nl[:,1],color="mediumblue",linewidth=2,linestyle=":")
ax2.plot(sU_curve1_Nm[:,0],sU_curve1_Nm[:,1],color="purple",linewidth=2,linestyle="-",label='N='+'%e' % Nm)
ax2.plot(sU_curve2_Nm[:,0],sU_curve2_Nm[:,1],color="purple",linewidth=2,linestyle=":")
ax2.plot(sU_curve1_Nh[:,0],sU_curve1_Nh[:,1],color="red",linewidth=2,linestyle="-",label='N='+'%e' % Nh)
ax2.plot(sU_curve2_Nh[:,0],sU_curve2_Nh[:,1],color="red",linewidth=2,linestyle=":")

#ax1.plot(sU_pair_log[0:32,0],sU_pair_log[0:32,1],color="black",linewidth=2,linestyle="-")
#ax2.scatter(sU_pair_log[:,0],sU_pair_log[:,1],color="black",marker='.')
#for i in range(len(indx)):
#ax2.scatter(sU_comp[:,0],sU_comp[:,1],color="black",linewidth=4)
#ax1.errorbar(sU_pair_log[:,0],sU_pair_log[:,1],2*new_verr[:],color="black",linewidth=2,linestyle="-")

#ax1.plot(sU_pair_log[:,0],logU_check+1,color="magenta",linewidth=2,linestyle="--")
#ax1.plot(sU_pair_log[:,0],logU_check-1,color="magenta",linewidth=2,linestyle="--")

ax2.set_xlim([1.2*log10_sc_max,log10_s_max])
ax2.set_ylim([log10_U_min,log10_U_max])
ax2.set_xlabel(r'Selection coefficient ($\log_{10}s$)',fontsize=18,labelpad=20)
ax2.set_ylabel(r'Mutation rate ($\log_{10}U$)',fontsize=18,labelpad=10)
ax2.tick_params(labelsize=16)
ax2.legend(loc=3)

xh_loc = (log10_s_max-1.2*log10_sc_max)
yh_loc = (log10_U_max-log10_U_min)

plt.text(1.2*log10_sc_max+0.70*xh_loc,log10_U_min+0.65*yh_loc,r'$N = 10^9$',fontsize=16)
plt.text(1.2*log10_sc_max+0.70*xh_loc,log10_U_min+0.60*yh_loc,r'$v = 5.3\times 10^{-5}$',fontsize=16)

#plt.text(0.9*log10_sc_max,0.55*log10_U_min,"Traveling Wave\n    Regime",fontsize=16)
#plt.text(0.85*log10_sc_max,0.93*log10_U_min,"Origin-fixation\n     Regime",fontsize=16)
#plt.text(0.85*log10_sc_max,0.13*log10_U_min,"Discontinuous\n     Regime",color="black",fontsize=16)
#plt.arrow(1.09*log10_sc_max,0.33*log10_U_min,-.1,1.5,linewidth=2,head_width=.07,color="black")

fig2.savefig('figures/fig_sUtradeoff_pheno_adapt_N.pdf')



#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
# region with smaller v




#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------




#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------
#--------------------------------------------------------------------------------