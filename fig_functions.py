# -*- coding: utf-8 -*-
"""
Created on Sat Feb 02 11:25:45 2019
Masel Lab
Project: Mutation-driven Adaptation
@author: Kevin Gomez

Description:
Defines the basics functions used in all scripts that process matlab
data and create figures in the mutation-driven adaptation manuscript.
"""
# libraries
import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
from numpy import inf

# functions
def get_vDF(N,s,U):
    # Calculates the rate of adaptation v, derived in Desai and Fisher 2007
    # for a given population size (N), selection coefficient (s) and beneficial
    # mutation rate (U)
    #
    # Inputs:
    # N - population size
    # s - selection coefficient
    # U - beneficial mutation rate
    #
    # Output: 
    # v - rate of adaptation
        
    v = s**2*(2*np.log(N*s)-np.log(s/U))/(np.log(s/U)**2)
    
    return v

def sU_tradeoff(s,N,v):
    # Computes the U-s trade-off function that preserves rate of adaptation (v)
    # for a chosen population size (N)
    #    
    # Inputs:
    # s - selection coefficient
    # N - populations size
    # v - fixed rate of adaptation
    #    
    # Outputs:
    # U - beneficial mutation rate yielding v, given N and s

    if (N*Uc[i]*np.log(N*s[i]) < 1):    # successional regime 
        U = v/N*s**2
    else:                               # continuous concurrent mutations regime                              
        U = s*np.exp(-(0.5*s**2/v)*(np.sqrt(8*v*np.log(N*s)/s**2+1)-1))
        
    return U
    
def get_min_s_disc_reg(N,v):
    # Computes the U-s trade-off function that preserves rate of adaptation (v)
    # for a chosen population size (N)
    #    
    # Inputs:
    # N - populations size
    # v - fixed rate of adaptation
    #    
    # Outputs:
    # min_s - smallest s yielding v, given N, that remains in conc-regime

    sp.minimize()
    return min_s
