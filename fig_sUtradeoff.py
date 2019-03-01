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
    
    s_min = 1/N 
    s_max = np.sqrt(v*np.log(N*np.sqrt(v)))     # accuracy depends on Nv choice 
    U_min = 0.1/(N*np.log(N*np.sqrt(v)))
    U_max = 10/N0
    
    return [s_min,s_max,U_min,U_max]

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
    
# set basic parameters of the figure
[N0,s0,U0] = [1e9,1e-2,1e-5]
v0 = myfun.get_vDF(N0,s0,U0)    # rate of adpatation (concurrent mutation regime)

# setting bounds for the window and computing their log10 values for the log-plot
[s_min,s_max,U_min,U_max] = sU_bounds(N0,v0)

log10_s_min = np.log10(s_min)
log10_s_max = np.log10(s_max)
log10_U_min = np.log10(U_min)
log10_U_max = np.log10(U_max)

# Define range for s and U
no_div = 100
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
disc_barrier = np.log10(np.asarray([[s1[i],min(s1[i],U_max)] for i in range(no_div)]))
    
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

# plot that includes thresholds and sU-tradeoff curve
fig1, ax1 = plt.subplots(1,1,figsize=[8,8])
ax1.fill_between(succ_shade[:,0],succ_shade[:,1],succ_shade[:,2],facecolor="blue")
ax1.fill_between(conc_shade[:,0],conc_shade[:,1],conc_shade[:,2],facecolor="orange")
#ax1.plot(conc_barrier[:,0],conc_barrier[:,1],c="black",label="conc/succ threshold")
ax1.fill_between(disc_shade[:,0],disc_shade[:,1],disc_shade[:,2],facecolor="green")
#ax1.plot(disc_barrier[:,0],disc_barrier[:,1],c="black",label="U(s)")

#ax1.set_xticklabels(my_xlabel)
#ax1.set_yticklabels(my_ylabel)        
ax1.set_xlabel('Selection coefficient (log10)',fontsize=18,labelpad=20)
ax1.set_ylabel('Mutation rate (log10)',fontsize=18,labelpad=10)
ax1.tick_params(axis='both',labelsize=14)        
ax1.legend()
ax1.axis('tight')        
ax2 = plt.twinx(ax1)
ax2.plot(np.log10(s),np.log10(wr),c="blue",label="U/s")
ax2.set_ylabel('U/s (log10)',fontsize=18,labelpad=10)
fig1.savefig('fig_sUtradeoff.pdf')

#fig2, ax2 = plt.subplots(1,1,figsize=[8,8])
#ax2.plot(np.log10(np.asarray(s)/s0),NsUreg1,label="Succ Regime")
#ax2.plot(np.log10(np.asarray(s)/s0),NsUreg2,label="Succ Regime")
#ax2.legend()


# PARAMETERS NEED TO BE SPACED ON A LOG SCALE TO HIGHLIGHT TRANSITIONS

#s_arry = [(s_min+i*(s_max-s_min)/no_div)/s  for i in range(no_div+1)]
#U_arry_succ = [(v/(N*s_arry[i]**2))/U for i in range(no_div+1)]
#U_arry_conc = [(conc_sU_tradeoff(s_arry[i],N,v))/U for i in range(no_div+1)]
#NsU_regime_succ = [np.log(N*U_arry_succ[i]*np.log(N*s_arry[i])) for i in range(no_div+1)]
#NsU_regime_conc = [np.log(N*U_arry_conc[i]*np.log(N*s_arry[i])) for i in range(no_div+1)]


# testing phase portait plotting


x

plt.figure(figsize=(18,6))
plt.quiver(apv, avv, dapv, davv, color='b', alpha=.75)
plt.box('off')
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.xlabel('Radians', fontsize=14)
plt.ylabel('Radians/Second', fontsize=14)
plt.title('Phase portrait for a simple pendulum', fontsize=16);



def f(sU_vect, N, v):
    S, U = sU_vect
    Theta = v*np.log(N*S)/S**2
    if():
        dU = -2*U/S
        dS = 1
    else:
        dU = ( (S + 8*S*Theta)/v + (2 - 4*np.log(N*S))/(S*np.sqrt(1+8*Theta)) )*U
        dS = 1
    return [dS, dU]

S1, U1 = np.meshgrid(s1, u1)
dS1, dU1 = np.zeros(S1.shape), np.zeros(U1.shape)

NI, NJ = S1.shape

for i in range(NI):
    for j in range(NJ):
        x, y = S1[i, j], U1[i, j]
        dS1[i,j], dU1[i,j] = f([x, y], N0, v0)
        
Q = plt.quiver(S1, U1, dS1, dU1, color='r')
plt.xlabel('$S1$')
plt.ylabel('$U1$')