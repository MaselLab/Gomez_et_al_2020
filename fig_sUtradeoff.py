#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 14:15:28 2019

@author: kgomez81
"""
from mpl_toolkits.mplot3d import Axes3D
from numpy import inf
import matplotlib.ticker as mtick

import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import matplotlib.mlab as mlab

def conc_sU_tradeoff(s,N,v):
    U = s*np.exp(-(0.5*s**2/v)*(np.sqrt(8*v*np.log(N*s)/s**2+1)-1))
    return U

# basic parameters
[N0,s0,U0] = [1e9,1e-2,1e-5]
v0=s0**2*(2*np.log(N0*s0)-np.log(s0/U0))/(np.log(s0/U0)**2)

[s_min,s_max,no_div] = [1e-5, 3e-1, 100]

s = [(s_min*10**(i*np.log10(s_max/s_min)/no_div)) for i in range(no_div+1)]
Us = [v0/(N0*s[i]**2) for i in range(no_div+1)]
Uc = [conc_sU_tradeoff(s[i],N0,v0) for i in range(no_div+1)]

vs = [N0*Us[i]*s[i]**2 for i in range(no_div)]
vc = [s[i]**2*(2*np.log(N*(s[i]))-np.log(s[i]/Uc[i]))/(np.log(s[i]/Uc[i])**2) for i in range(no_div)]

# checking regime
NsUreg1 = [np.log(N0*Us[i]*np.log(N0*s[i])) for i in range(no_div+1)]
NsUreg2 = [np.log(N0*Uc[i]*np.log(N0*s[i])) for i in range(no_div+1)]

fig, ax = plt.subplots(1,1,figsize=[8,8])
ax.plot(np.log10(s_larry),np.log10(U_larry_succ),label="succ")
ax.plot(np.log10(s_larry),np.log10(U_larry_conc),label="conc")
ax.legend()
ax2 = plt.twinx(ax)
ax2.plot(np.log10(s_larry),NsU_lregime_succ,label="succ")
ax2.plot(np.log10(s_larry),NsU_lregime_conc,label="conc")
ax2.legend()


# PARAMETERS NEED TO BE SPACED ON A LOG SCALE TO HIGHLIGHT TRANSITIONS

#s_arry = [(s_min+i*(s_max-s_min)/no_div)/s  for i in range(no_div+1)]
#U_arry_succ = [(v/(N*s_arry[i]**2))/U for i in range(no_div+1)]
#U_arry_conc = [(conc_sU_tradeoff(s_arry[i],N,v))/U for i in range(no_div+1)]
#NsU_regime_succ = [np.log(N*U_arry_succ[i]*np.log(N*s_arry[i])) for i in range(no_div+1)]
#NsU_regime_conc = [np.log(N*U_arry_conc[i]*np.log(N*s_arry[i])) for i in range(no_div+1)]