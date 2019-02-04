# -*- coding: utf-8 -*-
"""
Created on Tue Jan 01 11:38:07 2019
Masel Lab
Project: Mutation-driven Adaptation
@author: Kevin Gomez

Description:
Code for processing matlab simulation data and converting it to pickle format
for use in python scripts for figures.
"""

cd ~/Documents/mutBiasCI/data/ # must be evaluated separately

#libraries
import pickle
import scipy as sp
import numpy as np
import copy as cpy

# basic functions needed for processing data
# -----------------------------------------------------------------------------


# read matlab outputs and create python data for figure using grand means
# -----------------------------------------------------------------------------

data_file=open('mutBiasCI_data_all_simulation_grand_means_ml-01-1.dat')
grand_means = data_file.read().splitlines()
data_file.close()
data_file=open('mutBiasCI_data_all_simulation_parameters_ml-01-0.dat')
parameters = data_file.read().splitlines()
data_file.close()
del data_file

num_of_sims = len(grand_means)     # number of simulations
for i in range(num_of_sims):
    grand_means[i]='grand_means[i]=np.array(['+grand_means[i].replace('\t',',')+'])'
    exec(grand_means[i])
    parameters[i]='parameters[i]=np.array(['+parameters[i].replace('\t',',')+'])'
    exec(parameters[i])
    
grand_means = np.asarray(grand_means)
parameters = np.asarray(parameters)

# basic parameters (these should be exported with paremeters)
[N,s,U] = [1e9, 1e-2, 1e-5]
v = s**2*(2*np.log(N*s)-np.log(s/U))/(np.log(s/U)**2)

# reconstruct array of parameters 
arry_dim = int(0.5*(np.sqrt(8*num_of_sims+1)-1))
sarry = np.zeros([arry_dim, 1])
Uarry = np.zeros([arry_dim, 1])
v1_data = np.zeros([arry_dim, arry_dim])
v2_data = np.zeros([arry_dim, arry_dim])

for i in range(arry_dim):
    for j in range(arry_dim-i):
        if (i==0):
            sarry[j,0] = parameters[j,3]
            Uarry[j,0] = parameters[j,4]    
        indx = int(i*arry_dim+j-i*(i-1)/2)
        v1_data[i,j+i] = grand_means[indx,1]
        v2_data[i,j+i] = grand_means[indx,2]
        if (i!=j+i):
            v1_data[j+i,i] = grand_means[indx,2]
            v2_data[j+i,i] = grand_means[indx,1]
            
pickle_file_name = 'fig_compareVdata-01.pickle'
pickle_file = open(pickle_file_name,'wb') 
pickle.dump([N,s,U,v,parameters,grand_means,sarry,Uarry,v1_data,v2_data],pickle_file,pickle.HIGHEST_PROTOCOL)
pickle_file.close()

    
    