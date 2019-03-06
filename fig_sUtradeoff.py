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

# functions

def sU_bounds(N,v):
    # computes mutation rate U that maintains fixed v given s 
    #
    # inputs:
    # s = selection coefficient
    # N = Population size
    # v = rate of adaptation
    #
    # outputs:
    # list of bounds for s and U

    sm = 0.5*np.sqrt(v)/np.sqrt(np.log(0.5*N*np.sqrt(v)))
    st = np.sqrt(v*np.log(N*np.sqrt(v)))
    
    s_min = 1/N
    s_max = 5*np.sqrt(v*np.log(N*np.sqrt(v)))     # accuracy depends on Nv choice 
    U_min = min(0.1/(N*np.log(N*np.sqrt(v))),v/(N*s_max**2))
    U_max = 10*sm*np.exp(-(0.5*sm**2/v)*(np.sqrt(8*myfun.theta(sm,N,v)+1)-1))
    
    return [s_min,s_max,U_min,U_max,sm,st]

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
    
# set basic parameters of the figure
[N0,s0,U0] = [1e9,1e-2,1e-5]
v0 = myfun.get_vDF(N0,s0,U0)    # rate of adpatation (concurrent mutation regime)

sp = np.sqrt(v0)/(2*np.sqrt(np.log(N0*np.sqrt(v0)/2)))
Up = myfun.sU_tradeoff(sp,N0,v0)

# setting bounds for the window and computing their log10 values for the log-plot
[s_min,s_max,U_min,U_max,sc_max,sc_trans] = sU_bounds(N0,v0)

log10_s_min = np.log10(s_min)
log10_s_max = np.log10(s_max)
log10_U_min = np.log10(U_min)
log10_U_max = np.log10(U_max)
log10_sc_max = np.log10(sc_max)
log10_sc_trans = np.log10(sc_trans)

# Define range for s and U
no_div = 100
no_div1 = no_div/2
no_div2 = 3*no_div/4

s1 = np.logspace(np.log10(s_min), np.log10(s_max), no_div)
u1 = np.logspace(np.log10(U_min), np.log10(U_max), no_div)

log10_s1 = np.log10(s1)
log10_u1 = np.log10(u1)

# special set of s values to help shade the drift barrier in sU space
sd = np.logspace(np.log10(s_min), np.log10(20*s_min), no_div/10) 
log10_sd = np.log10(sd)

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# drift barrier 
drift_shade = np.log10(np.asarray([[sd[i], U_min, U_max] for i in range(no_div/10)]))
drift_barrier = np.log10(np.asarray([[10*s_min, u1[i]] for i in range(no_div)]))

# successional-concurrent barrier
succ_shade = np.log10(np.asarray([[s1[i],U_min,U_max] for i in range(no_div)]))
conc_shade = np.log10(np.asarray([[s1[i],myfun.succ_conc_barrier(s1[i],N0,v0),U_max] for i in range(no_div)]))
conc_barrier = np.log10(np.asarray([[s1[i],myfun.succ_conc_barrier(s1[i],N0,v0)] for i in range(no_div)]))

# discontinuous-concurrent barrier
disc_shade = np.log10(np.asarray([[s1[i],min(s1[i],U_max),U_max] for i in range(no_div)]))
disc_barrier = np.log10(np.asarray([[s1[i],min(0.1*s1[i],U_max)] for i in range(no_div)]))

# sU-tradeoff curves in concurrent and successional regimes
log10_sc_max = np.log10(sc_max)
log10_sc_trans = np.log10(sc_trans)

s_reg = np.logspace(log10_sc_max,log10_s_max,no_div2)
s_reg1 = np.logspace(log10_sc_max,log10_sc_trans,no_div1)
s_reg2 = np.logspace(log10_sc_trans,log10_s_max,no_div1)

sU_tradeoff_curve = np.log10(np.asarray([[s_reg[i],myfun.sU_tradeoff(s_reg[i],N0,v0)] for i in range(no_div2)]))    

sU_tradeoff_conc_curve1 = np.log10(np.asarray([[s_reg1[i],myfun.sU_tradeoff_conc(s_reg1[i],N0,v0)] for i in range(no_div1)]))    
sU_tradeoff_conc_curve2 = np.log10(np.asarray([[s_reg2[i],myfun.sU_tradeoff_conc(s_reg2[i],N0,v0)] for i in range(no_div1)]))    

sU_tradeoff_succ_curve1 = np.log10(np.asarray([[s_reg1[i],myfun.sU_tradeoff_succ(s_reg1[i],N0,v0)] for i in range(no_div1)]))    
sU_tradeoff_succ_curve2 = np.log10(np.asarray([[s_reg2[i],myfun.sU_tradeoff_succ(s_reg2[i],N0,v0)] for i in range(no_div1)]))    

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# plot that includes thresholds and sU-tradeoff curve
fig1, ax1 = plt.subplots(1,1,figsize=[8,8])
ax1.fill_between(succ_shade[:,0],succ_shade[:,1],succ_shade[:,2],facecolor="deepskyblue")
ax1.fill_between(conc_shade[:,0],conc_shade[:,1],conc_shade[:,2],facecolor="gold")
ax1.fill_between(disc_shade[:,0],disc_shade[:,1],disc_shade[:,2],facecolor="limegreen")

ax1.plot(sU_tradeoff_conc_curve1[:,0],sU_tradeoff_conc_curve1[:,1],color="mediumblue",linewidth=2,linestyle="-")
ax1.plot(sU_tradeoff_conc_curve2[:,0],sU_tradeoff_conc_curve2[:,1],color="mediumblue",linewidth=2,linestyle="--")
ax1.plot(sU_tradeoff_succ_curve1[:,0],sU_tradeoff_succ_curve1[:,1],color="red",linewidth=2,linestyle="--")
ax1.plot(sU_tradeoff_succ_curve2[:,0],sU_tradeoff_succ_curve2[:,1],color="red",linewidth=2,linestyle="-")

ax1.set_xlim([1.2*log10_sc_max,log10_s_max])
ax1.set_ylim([log10_U_min,log10_U_max])

ax1.set_xlabel(r'Selection coefficient ($\log_{10}s$)',fontsize=18,labelpad=20)
ax1.set_ylabel(r'Mutation rate ($\log_{10}U$)',fontsize=18,labelpad=10)

plt.text(0.83*log10_sc_max,0.6*log10_U_min,"Concurrent\n   Regime",fontsize=16)
plt.text(0.85*log10_sc_max,0.96*log10_U_min,"Origin-Fixation\n     Regime",fontsize=16)
plt.text(1.15*log10_sc_max,0.4*log10_U_min,"Discontinuous\n   Regime",fontsize=16)
plt.arrow(1.09*log10_sc_max,0.33*log10_U_min,-.1,.65,linewidth=2,head_width=.07,color="black")

fig1.savefig('fig_sUtradeoff_pheno_adapt.pdf')

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# sU-tradeoff for genotype adaptation

fig2, ax2 = plt.subplots(1,1,figsize=[8,8])
ax2.plot(theta_curve[:,0],theta_curve[:,1],label="Succ Regime")
#ax2.plot(np.log10(np.asarray(s)/s0),NsUreg2,label="Succ Regime")
#ax2.legend()


