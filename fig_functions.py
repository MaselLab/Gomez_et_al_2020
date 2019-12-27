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
import scipy.optimize as op
import numpy as np
from numpy import inf


# *****************************************************************************
# FUNCTIONS TO GET QUANTITIES FROM DESAI AND FISHER 2007
# *****************************************************************************


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

def get_qDF(N,s,U):
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
        
    q = 2*np.log(N*s)/np.log(s/U)
    
    return q

# *****************************************************************************
# FUNCTIONS TO APPROXIMATE TRANSITIONS BETWEEN REGIMES IN s,U SPACE
# *****************************************************************************

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

def s_max_Uc(N,v):
    # approximates the s value that maximizes U in concurrent regime
    #
    # inputs:
    # s = selection coefficient
    # N = Population size
    # v = rate of adaptation
    #
    # outputs:
    # s_m = estimated s marking transition to successional from concurrent
    
    s_m = 0.5*np.sqrt(v)/np.sqrt(np.log(0.5*N*np.sqrt(v)))
    
    return s_m

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
    
def sU_bounds(N,v):
    # computes appropriate bounds of sU space for looking at tradeoff with 
    # phenotype adaptation
    # inputs:
    # s = selection coefficient
    # N = Population size
    # v = rate of adaptation
    # outputs:
    # list of bounds for s and U

    sm = 1.1*s_max_Uc(N,v)                      # approx s for max Uc
    st = 1.2*s_transition2succ(N,v)             # approx s for transition conc to succ
    s_min = 1/N
    s_max = 5*st     
    U_min = 0.01*vContour_OF(s_max,N,v)
    U_max = 100*vContour_MM(sm,N,v)
    
    return [s_min,s_max,U_min,U_max,sm,st]

# *****************************************************************************
# FUNCTIONS FOR U(s) PARAMETERIZATIONS
# *****************************************************************************

#def vContour_full(sarry,N,v):
#    # Computes the U(s) parameteriziation of a v-contour
#    # for a chosen population size (N). The curve is returns as Uarry in two
#    # pieces because the theory gives a discontinuous transition between the 
#    # the diffusive mutations regime (DM) and the multiple mutation regime (MM).
#    # Origin fixation (OF) and mutltiple mutations intersect so those can be 
#    # returned as one piece.
#    #    
#    # Inputs:
#    # sarry - array of selection coefficients, should be column
#    # N - populations size
#    # v - fixed rate of adaptation
#    #    
#    # Outputs:
#    # Uarry - beneficial mutation rate yielding v, given N and s OF and MM
#
#
#    Uarry1 = np.zeros(sarry.shape)   # stores the calculated
#     = np.zeros(sarry.shape)   # stores the calculated 
#    
#    n = int(sarry.shape[0])     # assuming sarry stored as column
#
#    for i in range (n):    
#        U1 = vContour_OF(s,N,v)
#        U2 = vContour_OF(s,N,v)
#        U3 = vContour_OF(s,N,v)
#    
#    U = max(U1,U2)
#    
#    return U

def vContour_OFMM(s,N,v):
    # Computes the U-s trade-off function that preserves rate of adaptation (v)
    # for a chosen population size (N) for OF and MM as one curve
    #    
    # Inputs:
    # s - selection coefficient
    # N - populations size
    # v - fixed rate of adaptation
    #    
    # Outputs:
    # U - beneficial mutation rate yielding v, given N and s
    theta  = v*np.log(N*s)/s**2
    
    U1 = v*(1+s)**2/(2*N*s**2)                 # successional regime 
    U2 = s*np.exp(-(0.5*s**2/v)*(np.sqrt(8*theta+1)-1))  # concurrent mutations regime      
    
    U = max(U1,U2)
    
    return U
    
def vContour_OF(s,N,v):
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
    
    U = v*(1+s)**2/(2*N*s**2)
    
    return U

def vContour_MM(s,N,v):   
    # Computes the U-s trade-off function that preserves rate of adaptation (v)
    # for a chosen population size (N) for the multiple mutations regime (U<<s)
    #    
    # Inputs:
    # s - selection coefficient
    # N - populations size
    # v - fixed rate of adaptation
    #    
    # Outputs:
    # U - beneficial mutation rate yielding v, given N and s     
    
    theta = v*np.log(N*s)/s**2    # special ratio marking transition to MM from OF
    
    U = s*np.exp(-(0.5*s**2/v)*(np.sqrt( 8*theta + 1 ) - 1 ))  # concurrent mutations regime      

    return U
    
def vContour_DM(s,N,v):
    # Computes the U-s trade-off function that preserves rate of adaptation (v)
    # for a chosen population size (N) for the diffusive mutations regime
    #    
    # Inputs:
    # s - selection coefficient
    # N - populations size
    # v - fixed rate of adaptation
    #    
    # Outputs:
    # U - beneficial mutation rate yielding v, given N and s    
    #
    # hallatschek approximation has to be adjusted by a constant (~0.3)

    U = (0.3*v**(1.5)/s**2)/np.log(6**(1/6.0)*N*np.sqrt(v))**(1/2.0)  

    return U 



