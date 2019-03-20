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

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#libraries
import pickle
import scipy as sp
import numpy as np
import copy as cpy

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

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------

#libraries
import pickle
import scipy as sp
import numpy as np
import copy as cpy

# read matlab outputs and create python data for figure using grand means
# -----------------------------------------------------------------------------

data_file=open('data/mutBiasCI_data_all_simulation_grand_means_ml-08-1.dat')
grand_means = data_file.read().splitlines()
data_file.close()
data_file=open('data/mutBiasCI_data_all_simulation_parameters_ml-08-0.dat')
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
            
pickle_file_name = 'data/fig_compareVdata-08.pickle'
pickle_file = open(pickle_file_name,'wb') 
pickle.dump([N,s,U,v,parameters,grand_means,sarry,Uarry,v1_data,v2_data],pickle_file,pickle.HIGHEST_PROTOCOL)
pickle_file.close()

# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------    
# checking the discontinuous regime

#libraries
import pickle
import scipy as sp
import numpy as np
import copy as cpy

# basic parameters (these should be exported with paremeters)
[N,s,U] = [1e9, 1e-2, 1e-5]
v = s**2*(2*np.log(N*s)-np.log(s/U))/(np.log(s/U)**2)

# read matlab outputs and create python data for figure using grand means
data_file=open('data/mutBiasCI_data_all_simulation_grand_means_ml-05-1.dat')
grand_means1 = data_file.read().splitlines()
data_file.close()
data_file=open('data/mutBiasCI_data_all_simulation_parameters_ml-05-0.dat')
parameters1 = data_file.read().splitlines()
data_file.close()
del data_file

num_of_sims = len(grand_means1)     # number of simulations
for i in range(num_of_sims):
    grand_means1[i]=grand_means1[i].replace('NaN','0')
    grand_means1[i]='grand_means1[i]=np.array(['+grand_means1[i].replace('\t',',')+'])'
    exec(grand_means1[i])
    parameters1[i]='parameters1[i]=np.array(['+parameters1[i].replace('\t',',')+'])'
    exec(parameters1[i])

data_file=open('data/mutBiasCI_data_all_simulation_grand_means_ml-06-1.dat')
grand_means2 = data_file.read().splitlines()
data_file.close()
data_file=open('data/mutBiasCI_data_all_simulation_parameters_ml-06-0.dat')
parameters2 = data_file.read().splitlines()
data_file.close()
del data_file

num_of_sims = len(grand_means2)     # number of simulations
for i in range(num_of_sims):    
    grand_means2[i]=grand_means2[i].replace('NaN','0')
    grand_means2[i]='grand_means2[i]=np.array(['+grand_means2[i].replace('\t',',')+'])'
    exec(grand_means2[i])
    parameters2[i]='parameters2[i]=np.array(['+parameters2[i].replace('\t',',')+'])'
    exec(parameters2[i])

grand_means = grand_means1 + grand_means2    # extra empty data for some reason
parameters = parameters1 + parameters2
    
grand_means = np.asarray(grand_means)
parameters = np.asarray(parameters)

[dim_s1,dim_s2,dim_U] = [41,43,41]       # these are the number of samples I took

# reconstruct array of parameters 
sarry = np.zeros([dim_s1+dim_s2, 1])
Uarry = np.zeros([dim_s1+dim_s2, dim_U])
sU_pair = np.zeros([dim_s1+dim_s2,2])
v_err = np.zeros([dim_s1+dim_s2,4])
v1_data = np.zeros([dim_s1+dim_s2, dim_U])

for i in range(dim_s1):
    sarry[i,0] = parameters[i*dim_U,1]
    for j in range(dim_U):
        v1_data[i,j] = grand_means[i*dim_U+j,1]
        Uarry[i,j] = parameters[i*dim_U+j,2]

for i in range(dim_s2):
    sarry[i+dim_s1,0] = parameters[(i+dim_s1)*dim_U,1]
    for j in range(dim_U):
        v1_data[i+dim_s1,j] = grand_means[(i+dim_s1)*dim_U+j,1]
        Uarry[i+dim_s1,j] = parameters[(i+dim_s1)*dim_U+j,2]
        
# construct U(s) tradeoff in disc regime
for i in range(dim_s1):
    sU_pair[i,:] = [sarry[i],Uarry[i,0]]
    v_err[i,:] = np.asarray([i,v1_data[i,0],v,np.abs(v1_data[i,0]-v)/v])
    for j in range(dim_U-1):
        if (np.abs(v1_data[i,j+1]-v) < np.abs(v_err[i,1]-v)):
             sU_pair[i,1] = Uarry[i,j+1]
             v_err[i,:] = np.asarray([i,v1_data[i,j+1],v,np.abs(v1_data[i,j+1]-v)/v])

for i in range(dim_s2):
    sU_pair[i+dim_s1,:] = [sarry[i+dim_s1],Uarry[i+dim_s1,0]]
    v_err[i+dim_s1,:] = np.asarray([i+dim_s1,v1_data[i+dim_s1,0],v,np.abs(v1_data[i+dim_s1,0]-v)/v])
    for j in range(dim_U-1):
        if (np.abs(v1_data[i+dim_s1,j+1]-v) < np.abs(v_err[i+dim_s1,1]-v)):
             sU_pair[i+dim_s1,1] = Uarry[i+dim_s1,j+1]
             v_err[i+dim_s1,:] = np.asarray([i+dim_s1,v1_data[i+dim_s1,j+1],v,np.abs(v1_data[i+dim_s1,j+1]-v)/v])

# --------- get best estimates for sU tradeoff --------------------------------
sU_pair = np.concatenate((sU_pair[41:52,:],sU_pair[1:13,:],sU_pair[52:67,:],sU_pair[26:30,:],sU_pair[67:,:]),axis=0)
v_err = np.concatenate((v_err[41:52,:],v_err[1:13,:],v_err[52:67,:],v_err[26:30,:],v_err[67:,:]),axis=0)            

# -------------------------------------------------------------------------------
# saving data [s,U,v_data,parameters,grand_means]
pickle_file_name = 'data/fig_discVdata-06.pickle'
pickle_file = open(pickle_file_name,'wb') 
pickle.dump([sU_pair,v_err,parameters,grand_means],pickle_file,pickle.HIGHEST_PROTOCOL)
pickle_file.close()
    
# function defining regions scanned in to get target U
    
#logU_check = np.zeros([dim_s,1])
#
#for i in range(dim_s):
#    if(np.log10(sarry[i])<-2.5):
#        logU_check[i,0] = -2.5*(np.log10(sarry[i])+2.95)-1.5
#    elif((np.log10(sarry[i])>=-2.5) and  (np.log10(sarry[i])<-1.5)):
#        logU_check[i,0] = -6.3*(np.log10(sarry[i])+1.5)-8.2
#    else:
#        logU_check[i,0] = -2*(np.log10(sarry[i])+1.5)-10.3    

#logU_check = np.zeros([dim_s,1])
#for i in range(dim_s):
#    if(np.log10(sarry[i])<-1.5):
#        logU_check[i,0] = -6.3*(np.log10(sarry[i])+1.5)-9.6
#    else:
#        logU_check[i,0] = -2*(np.log10(sarry[i])+1.5)-10