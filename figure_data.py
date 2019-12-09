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

#libraries
import pickle
import scipy as sp
import numpy as np
import copy as cpy

# -----------------------------------------------------------------------------
#            Process matlab simulation data for two traits
# -----------------------------------------------------------------------------

# set file names to read and store data
# -----------------------------------------------------------------------------

#sim_data_grandmeans = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_grand_means_ml-01-1.dat'
#sim_data_parameters = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_grand_means_ml-01-0.dat'
#pickle_file_name = 'data/fig_compareVdata-01.pickle'

#sim_data_grandmeans = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_grand_means_ml-08-1-v2.dat'
#sim_data_parameters = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_parameters_ml-08-0-v2.dat'
#pickle_file_name = 'data/fig_compareVdata-08.pickle'

#sim_data_grandmeans = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_grand_means_ml-15-1.dat')
#sim_data_parameters = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_parameters_ml-15-0.dat')
#pickle_file_name = 'data/fig_compareVdata-15.pickle'

#sim_data_grandmeans = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_grand_means_ml-16-1.dat'
#sim_data_parameters = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_parameters_ml-16-0.dat'
#pickle_file_name = 'data/fig_compareVdata-16.pickle'

#sim_data_grandmeans = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_grand_means_ml-19-1.dat'
#sim_data_parameters = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_parameters_ml-19-0.dat'
#pickle_file_name = 'data/fig_compareVdata-19.pickle'

#sim_data_grandmeans = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_grand_means_ml-20-1.dat'
#sim_data_parameters = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_parameters_ml-20-0.dat'
#pickle_file_name = 'data/fig_compareVdata-20.pickle'

#sim_data_grandmeans = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_grand_means_ml-100-1.dat'
#sim_data_parameters = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_parameters_ml-100-0.dat'
#pickle_file_name = 'data/fig_compareVdata-100.pickle'

#sim_data_grandmeans = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_grand_means_ml-101-1.dat'
#sim_data_parameters = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_parameters_ml-101-0.dat'
#pickle_file_name = 'data/fig_compareVdata-101.pickle'

sim_data_grandmeans = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_grand_means_ml-102-1.dat'
sim_data_parameters = 'data/TwoTraitSim/mutBiasCI_data_all_simulation_parameters_ml-102-0.dat'
pickle_file_name = 'data/fig_compareVdata-102.pickle'

# read matlab outputs and create python data for figure using grand means
# -----------------------------------------------------------------------------
data_file=open(sim_data_grandmeans)
grand_means = data_file.read().splitlines()
data_file.close()
data_file=open(sim_data_parameters)
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
# -----------------------------------------------------------------------------
[N,s,U] = [1e9, 1e-2, 1e-5]     # these are set in the simulation
v = s**2*(2*np.log(N*s)-np.log(s/U))/(np.log(s/U)**2)

# reconstruct array of parameters 
# -----------------------------------------------------------------------------
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
            
pickle_file = open(pickle_file_name,'wb') 
pickle.dump([N,s,U,v,parameters,grand_means,sarry,Uarry,v1_data,v2_data],pickle_file,pickle.HIGHEST_PROTOCOL)
pickle_file.close()

# -----------------------------------------------------------------------------    
# process data for U estimates to construct figure 1 (v isoquants)
# -----------------------------------------------------------------------------    

#libraries
import pickle
import scipy as sp
import numpy as np
import copy as cpy

# functions
# -----------------------------------------------------------------------------
def read_my_data(mydata):
    n = len(mydata)     
    for i in range(n):
        mydata[i]='mydata[i]=['+mydata[i].replace('\t',',')+']'
        exec(mydata[i])
    return np.asarray(mydata)

# set file names to import data from Stochastic approximation algorithm
# -----------------------------------------------------------------------------

#sim_data_parameters = 'data/SAapprox/mutBiasCI_estimate_U_ml-2-'
#sim_data_sU_pairs = 'data/SAapprox/mutBiasCI_estimate_U_ml-2-'
#pickle_file_name = 'data/fig_sUtradeoff_simdata-01.pickle'

#sim_data_parameters = 'data/SAapprox/mutBiasCI_estimate_U_ml-6-'
#sim_data_sU_pairs = 'data/SAapprox/mutBiasCI_estimate_U_ml-6-'
#pickle_file_name = 'data/fig_sUtradeoff_simdata-06.pickle'

#sim_data_parameters = 'data/SAapprox/mutBiasCI_estimate_U_ml-19-'
#sim_data_sU_pairs = 'data/SAapprox/mutBiasCI_estimate_U_ml-19-'
#pickle_file_name = 'data/fig_sUtradeoff_simdata-19.pickle'

#sim_data_parameters = 'data/SAapprox/mutBiasCI_estimate_U_ml-20-'
#sim_data_sU_pairs = 'data/SAapprox/mutBiasCI_estimate_U_ml-20-'
#pickle_file_name = 'data/fig_sUtradeoff_simdata-20.pickle'

#sim_data_parameters = 'data/SAapprox/mutBiasCI_estimate_U_ml-21-'
#sim_data_sU_pairs = 'data/SAapprox/mutBiasCI_estimate_U_ml-21-'
#pickle_file_name = 'data/fig_sUtradeoff_simdata-21.pickle'

#sim_data_parameters = 'data/SAapprox/mutBiasCI_estimate_U_ml-22-'
#sim_data_sU_pairs = 'data/SAapprox/mutBiasCI_estimate_U_ml-22-'
#pickle_file_name = 'data/fig_sUtradeoff_simdata-22.pickle'

sim_data_parameters = 'data/SAapprox/mutBiasCI_estimate_U_ml-100-'
sim_data_sU_pairs = 'data/SAapprox/mutBiasCI_estimate_U_ml-100-'
pickle_file_name = 'data/fig_sUtradeoff_simdata-100.pickle'

# read, process and save data from files
# -----------------------------------------------------------------------------
Nv_param=[]
sU_data=[]
   
for i in range(6):
    data_file=open(sim_data_parameters+str(i+1)+'-0.dat')
    temp_param = data_file.read().splitlines()
    Nv_param += [read_my_data(temp_param)] 
    data_file.close()
    
    data_file=open(sim_data_sU_pairs+str(i+1)+'-1.dat')
    temp_data = data_file.read().splitlines()
    sU_data += [read_my_data(temp_data)]
    data_file.close()

del data_file

pickle_file = open(pickle_file_name,'wb') 
pickle.dump([Nv_param,sU_data],pickle_file,pickle.HIGHEST_PROTOCOL)
pickle_file.close()

# -----------------------------------------------------------------------------
#            Process matlab simulation data for two traits
#             to compare v with one having fixed s and U
# -----------------------------------------------------------------------------

# set file names to read and store data
# -----------------------------------------------------------------------------
#libraries
import pickle
import scipy as sp
import numpy as np
import copy as cpy

#sim_data_grandmeans = 'data/FixedsU/mutBiasCI_data_all_simulation_grand_means_fixed_sU_ml-01-1.dat'
#sim_data_parameters = 'data/FixedsU/mutBiasCI_data_all_simulation_parameters_fixed_sU_ml-01-0.dat'
#pickle_file_name = 'data/fig_compare_sU_data-01.pickle'

#sim_data_grandmeans = 'data/FixedsU/mutBiasCI_data_all_simulation_grand_means_fixed_sU_ml-02-1.dat'
#sim_data_parameters = 'data/FixedsU/mutBiasCI_data_all_simulation_parameters_fixed_sU_ml-02-0.dat'
#pickle_file_name = 'data/fig_compare_sU_data-02.pickle'

sim_data_grandmeans = 'data/FixedsU/mutBiasCI_data_all_simulation_grand_means_fixed_sU_ml-03-1.dat'
sim_data_parameters = 'data/FixedsU/mutBiasCI_data_all_simulation_parameters_fixed_sU_ml-03-0.dat'
pickle_file_name = 'data/fig_compare_sU_data-03.pickle'

# read matlab outputs and create python data for figure using grand means
# -----------------------------------------------------------------------------
data_file=open(sim_data_grandmeans)
grand_means = data_file.read().splitlines()
data_file.close()
data_file=open(sim_data_parameters)
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
# -----------------------------------------------------------------------------
[N,s,U] = [1e9, 1e-2, 1e-5]     # these are set in the simulation
v = s**2*(2*np.log(N*s)-np.log(s/U))/(np.log(s/U)**2)

# reconstruct array of parameters 
# -----------------------------------------------------------------------------
arry_dim = int(np.sqrt(num_of_sims))
sarry = np.zeros([arry_dim, 1])
Uarry = np.zeros([arry_dim, 1])
v1_data = np.zeros([arry_dim, arry_dim])
v2_data = np.zeros([arry_dim, arry_dim])

for i in range(arry_dim):
    for j in range(arry_dim):
        indx = int(arry_dim*i+j)
        if (i==0):
            Uarry[j,0] = parameters[j,4]    
        
        sarry[i,0] = parameters[indx,3]
        
        v1_data[i,j] = grand_means[indx,1]
        v2_data[i,j] = grand_means[indx,2]
            
v1_data = np.transpose(v1_data)
v2_data = np.transpose(v2_data)

pickle_file = open(pickle_file_name,'wb') 
pickle.dump([N,s,U,v,parameters,grand_means,sarry,Uarry,v1_data,v2_data],pickle_file,pickle.HIGHEST_PROTOCOL)
pickle_file.close()