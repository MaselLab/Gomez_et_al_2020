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
import fig_functions.py as myfun

def get_theta(N,s,v):
    theta = v*np.log(N*s)/s**2
    return theta
    
# basic parameters
[N0,s0,U0] = [1e9,1e-2,1e-5]
v0=myfun.get_v(N0,s0,U0)

[s_min,s_max] = [1e-6, 5e-1]
[U_min,U_max] = [1e-13, 5e-3]
no_div = 100

def f(Y, t):
    y1, y2 = Y
    return [y2, -np.sin(y1)]

y1 = np.linspace(-2.0, 8.0, 20)
y2 = np.linspace(-2.0, 2.0, 20)
Y1, Y2 = np.meshgrid(y1, y2)
t = 0
u, v = np.zeros(Y1.shape), np.zeros(Y2.shape)
NI, NJ = Y1.shape
for i in range(NI):
    for j in range(NJ):
        x = Y1[i, j]
        y = Y2[i, j]
        yprime = f([x, y], t)
        u[i,j] = yprime[0]
        v[i,j] = yprime[1]
Q = plt.quiver(Y1, Y2, u, v, color='r')
plt.xlabel('$y_1$')
plt.ylabel('$y_2$')
plt.xlim([-2, 8])
plt.ylim([-4, 4])


s = [(s_min*10**(i*np.log10(s_max/s_min)/no_div)) for i in range(no_div+1)]

Us = [v0/(N0*s[i]**2) for i in range(no_div+1)]
Uc = [conc_sU_tradeoff(s[i],N0,v0) for i in range(no_div+1)]
U = get_sUtransition(s,Us,Uc,N0)

vs = [N0*Us[i]*s[i]**2 for i in range(no_div+1)]
vc = [s[i]**2*(2*np.log(N0*(s[i]))-np.log(s[i]/Uc[i]))/(np.log(s[i]/Uc[i])**2) for i in range(no_div+1)]

wr = [U[i]/s[i] for i in range(no_div+1)]

fig1, ax1 = plt.subplots(1,1,figsize=[8,8])
#ax1.plot(np.log10(s),np.log10(Us),c="red",label="Succ Regime")
#ax1.plot(np.log10(s),np.log10(Uc),c="blue",label="Conc Regime")
ax1.plot(np.log10(s),np.log10(U),c="black",label="U(s)")
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