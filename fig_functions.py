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

def theta(s,N,v):
    # computes the characteristic width of the fitness distribution
    #
    # inputs:
    # s = selection coefficient
    # N = Population size
    # v = rate of adaptation
    #
    # outputs:
    # theta = characteristic width of the fitness distribution
    
    theta = v*np.log(N*s)/s**2
    
    return theta

def s_transition2succ(N,v):
    # approximates the transition point from concurrent to successional 
    #
    # inputs:
    # s = selection coefficient
    # N = Population size
    # v = rate of adaptation
    #
    # outputs:
    # s_t = estimated s marking transition to successional from concurrent
    
    s_t = np.sqrt(v*np.log(N*np.sqrt(v)))
    
    return s_t

def succ_conc_barrier(s,N,v):
    # computes mutation rate cutoff between succ and conc regime, given a 
    # selection coefficient s.
    #
    # inputs:
    # s = selection coefficient
    # N = Population size
    # v = rate of adaptation
    #
    # outputs:
    # U = mutation rate cutoff    
    
    U = 1/(N*np.log(N*s))
    
    return U

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

    if (theta(s,N,v) < 1):              # successional regime 
        U = v/(N*s**2)
    else:                               # concurrent mutations regime                              
        U = s*np.exp(-(0.5*s**2/v)*(np.sqrt(8*theta(s,N,v)+1)-1))
        
    return U

def sU_tradeoff_succ(s,N,v):
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
    
    U = v/(N*s**2)
    
    return U
    
def sU_tradeoff_conc(s,N,v):
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

    U = s*np.exp(-(0.5*s**2/v)*(np.sqrt(8*theta(s,N,v)+1)-1))
        
    return U