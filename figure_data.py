# -*- coding: utf-8 -*-
"""
Created on Tue Jan 01 11:38:07 2019
@author: Kevin Gomez (Masel Lab)
Library of functions used in plots.py and plotdata.py
"""

#--------FUNCTIONS REQUIRE PACKAGES LISTED:------------------------------------
from scipy.stats import multivariate_normal
import pickle
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
def get_sample_window(times,start_time,end_time):
# returns: indeces of times that correspond to start_time and end_time
 
    [num_pts,start_indx,end_indx] = [len(times),0,0]
    
    for i in range(num_pts):
        if times[start_indx] <= start_time:
            start_indx = start_indx + 1
        if times[end_indx] <= end_time:
            end_indx = end_indx + 1
    
    return [start_indx,end_indx]

# -----------------------------------------------------------------------------
def get_trait_mean_var(genotypes,abundances,traitno):
    
    if ((traitno == 1) | (traitno == 2)):
        mean = (genotypes[:,traitno-1].dot(abundances[0]))/sum(abundances[0])
        var = (((genotypes[:,traitno-1]-mean*np.ones(np.shape(genotypes[:,traitno-1])))**2).dot(abundances[0]))/sum(abundances[0])
    if (traitno == 0):
        mean1 = (genotypes[:,0].dot(abundances[0]))/sum(abundances[0])
        mean2 = (genotypes[:,1].dot(abundances[0]))/sum(abundances[0])
        means_arry = np.asarray([[mean1,mean2] for i in range(len(genotypes[:,0]))])
        mean = mean1 + mean2
        var = (abundances.dot((((genotypes - means_arry)**2).dot(np.ones([2,1]))))[0][0])/sum(abundances[0])
    return [mean, var]
       
# -----------------------------------------------------------------------------
def get_1D_proj(genotypes,abundances,traitno):

    if ((traitno == 1) | (traitno == 2)):
        trait_min = np.min(genotypes[:,traitno-1])
        trait_max = np.max(genotypes[:,traitno-1])
    
        trait_classes = [trait_min+i-3 for i in range(trait_max-trait_min+1+6)]
        trait_totals = [0 for i in range(trait_max-trait_min+1+6)]
    
        for i in range(len(genotypes[:,traitno-1])):
            indx = genotypes[i,traitno-1]-trait_min+3
            trait_totals[indx] = trait_totals[indx]+abundances[0][i] 
    
    if (traitno == 0):
        genotype_fitnesses = genotypes[:,0]+genotypes[:,1]
        trait_min = np.min(genotype_fitnesses)
        trait_max = np.max(genotype_fitnesses)
    
        trait_classes = [trait_min+i-3 for i in range(trait_max-trait_min+1+6)]
        trait_totals = [0 for i in range(trait_max-trait_min+1+6)]
    
        for i in range(len(genotype_fitnesses)):
            indx = genotype_fitnesses[i]-trait_min+3
            trait_totals[indx] = trait_totals[indx]+abundances[0][i] 
    
    return [trait_classes, trait_totals]

# -----------------------------------------------------------------------------
def get_2D_distr(genotypes,abundances,box_dim):
# box_dim = gives array data to bound distr correponding to genotypes & abund.
#           [[width1,margin1],[width2,margin2]]
# returns: an array whose elements are the abundances of the fit classes
    
    hhfgenotypes = np.asarray(get_hifit_front_genos(genotypes))
    hhf_points = []
    tot_pop_size = sum(abundances[0])  # be careful with your sums of arrays!!!    
    dim1_data = [np.min(genotypes[:,0]),np.max(genotypes[:,0])]
    dim2_data = [np.min(genotypes[:,1]),np.max(genotypes[:,1])]
    
    if((box_dim[0][0] < dim1_data[1]-dim1_data[0]) | (box_dim[1][0] < dim2_data[1]-dim2_data[0])):
        print "Error with box dimensions"
        end()
    
    my_distr = np.zeros([box_dim[0][0],box_dim[1][0]])
    
    for i in range(len(genotypes)):
        indx1 = genotypes[i][0] - dim1_data[0] + box_dim[0][1]
        indx2 = genotypes[i][1] - dim2_data[0] + box_dim[1][1]        
        my_distr[indx1,indx2] = max([abundances[0][i]/tot_pop_size,1/tot_pop_size])
        
    xlabels = [(dim1_data[0] - box_dim[0][1]-1 + i) for i in range(box_dim[0][0]+1)]
    ylabels = [(dim2_data[0] - box_dim[1][1]-1 + i) for i in range(box_dim[1][0]+1)]

    for i in range(len(hhfgenotypes)):
        indx1 = hhfgenotypes[i][0] - dim1_data[0] + box_dim[0][1]
        indx2 = hhfgenotypes[i][1] - dim2_data[0] + box_dim[1][1]        
        hhf_points.append([indx1+.5,indx2+0.5])
    
    hhf_points = np.asarray(hhf_points)
    return [my_distr,xlabels,ylabels,hhf_points]

# -----------------------------------------------------------------------------
def get_cov_by_fitness_line(genotypes,abundances,s):
    
    mean_fit = get_trait_mean_var(genotypes,abundances,0)[0]
    num_genotypes = len(abundances[0])
    
    fit1D = [genotypes[i,0]+genotypes[i,1] for i in range(num_genotypes)]
    fit1Dshrt = list(set(fit1D))
    fit1Dshrt.sort()
    fit1Dcovs = []    
    tempcov = 0
    tempfreq = 0
    tempmean1 = 0
    tempmean2 = 0
    popsize = sum(abundances[0])
    
    for i in range(len(fit1Dshrt)):
        tempcov = 0
        tempfreq = 0
        tempmean1 = 0
        tempmean2 = 0
        for j in range(num_genotypes):
            if(fit1D[j]==fit1Dshrt[i]):
                tempmean1 += s*genotypes[j,0]*abundances[0][j]/popsize
                tempmean2 += s*genotypes[j,1]*abundances[0][j]/popsize                
                tempfreq += abundances[0][j]/popsize
                tempcov += s**2*genotypes[j,0]*genotypes[j,1]*(abundances[0][j]/popsize)
        tempcov = tempcov/tempfreq - (tempmean1/tempfreq)*(tempmean2/tempfreq)
        fit1Dcovs = fit1Dcovs+[[fit1Dshrt[i]-mean_fit,tempfreq,tempcov]]
    
    # should return [[fit_i,p_i,cov_i,] for i = min_fit,...,max_fit]
    return fit1Dcovs

# -----------------------------------------------------------------------------    
def get_vNsU_perChg(N,s,U,n):
    vrate = ((2*np.log(N*s)-np.log(s/n/U))/np.log(s/n/U)**2/n) / ((2*np.log(N*s)-np.log(s/U))/np.log(s/U)**2) 
    return vrate

# -----------------------------------------------------------------------------
def get_vNsU(N,s,U):
    vrate = s**2*(2*np.log(N*s)-np.log(s/U))/(np.log(s/U)**2) 
    return vrate

# -----------------------------------------------------------------------------    
def get_cov_cov(times,nose_cov,fit_cov,N,s,U):
    tau_q = (np.log(s/U))**2/(s*(2*np.log(N*s)-np.log(s/U)))
    q = (2*np.log(N*s))/(np.log(s/U))
    time_d = int(np.floor(q*tau_q))
    
    new_times = []
    new_covs = []
    new_ncovs = []
    
    for i in range(len(times)):
        if(np.mod(times[i],1)<0.00000001):
            new_times = new_times+[times[i]]
            new_covs = new_covs + [fit_cov[i]]
            new_ncovs = new_ncovs + [nose_cov[i]]
    
    new_covs = np.asarray(new_covs)
    new_ncovs = np.asarray(new_ncovs)
    t_off = [i+1 for i in range(2*time_d)]
    t_cov = [0 for i in range(2*time_d)]
    
    for i in range(2*time_d):
        t_cov[i] = (np.cov(np.vstack((new_covs[i+1:],new_ncovs[:-(1+i)])))[0,1])/np.std(new_covs[i+1:])/np.std(new_ncovs[:-(1+i)])
    
    return [t_off,t_cov,new_times,new_covs,new_ncovs]
    
# -----------------------------------------------------------------------------
def get_subset_times(N,s,U,times,scaling):
    tau_q = scaling*((np.log(s/U))**2)/(s*(2*np.log(N*s)-np.log(s/U)))
    indx_list = [0]
    indx = 0
    
    while(times[indx]+tau_q < times[-1]):
        indx = get_sample_window(times,times[indx],times[indx]+tau_q)[1]
        indx_list = indx_list + [indx]
    return indx_list
# -----------------------------------------------------------------------------    

def get_hifit_front_line(genotypes,num_points,box_dim):
    num_geno = len(genotypes)
    min_x = min(genotypes[:,0])
    min_y = min(genotypes[:,1])
    
    hffrt = []    
    L1 = 1+max([genotypes[i][0]-min_x+genotypes[i][1]-min_y+2*box_dim[0][1] for i in range(num_geno)])
    
    x_start = L1-box_dim[0][0]
    x_end = box_dim[0][0]
    
    xl = np.asarray([1.0*(x_end-x_start)*i/num_points + x_start for i in range(num_points+1)])
    yl = np.asarray([L1-1.0*xl[i] for i in range(num_points+1)])
    return [xl,yl]
    
# -----------------------------------------------------------------------------  
    
def get_hifit_front_genos(genotypes):
    num_geno = len(genotypes)
    hhfgenotypes = []
    L = 1+np.max([genotypes[i][0]+genotypes[i][1] for i in range(num_geno)])
    
    for i in range(num_geno):
        if(genotypes[i][0]+genotypes[i][1]+1 == L):
            if [genotypes[i][0]+1,genotypes[i][1]] not in hhfgenotypes:
                hhfgenotypes.append([genotypes[i][0]+1,genotypes[i][1]])
            if [genotypes[i][0],genotypes[i][1]+1] not in hhfgenotypes:
                hhfgenotypes.append([genotypes[i][0],genotypes[i][1]+1])
    
    return hhfgenotypes

# -----------------------------------------------------------------------------  

def get_stoch_genotypes(genotypes,abundances,cutoff):
# returns the set of classes that are smaller than the given cutoff
# and the those classes that are at the high fitness front.
    nosefitness = np.max(np.matmul(genotypes,np.ones([2,1])))
    hhf_points = []
    stoch_points = []
    
    for i in range(len(genotypes)):
        if (abundances[0][i]<cutoff):
            stoch_points = stoch_points+[genotypes[i]]
            if(genotypes[i][0]+genotypes[i][1]==nosefitness):
                hhf_points = hhf_points + [genotypes[i]]
    
    hhf_points = np.asarray(hhf_points)
    stoch_points = np.asarray(stoch_points)
    return [hhf_points,stoch_points]

def get_2D_distr2(genotypes,abundances,box_dim,cutoff):
# box_dim = gives array data to bound distr correponding to genotypes & abund.
#           [[width1,margin1],[width2,margin2]]
# returns: an array whose elements are the abundances of the fit classes
    
    [hhfgenotypes,stochgenotypes] = get_stoch_genotypes(genotypes,abundances,cutoff)
    hhf_points = []
    stoch_points = []
    
    tot_pop_size = sum(abundances[0])  # be careful with your sums of arrays!!!    
    dim1_data = [np.min(genotypes[:,0]),np.max(genotypes[:,0])]
    dim2_data = [np.min(genotypes[:,1]),np.max(genotypes[:,1])]
    
    if((box_dim[0][0] < dim1_data[1]-dim1_data[0]) | (box_dim[1][0] < dim2_data[1]-dim2_data[0])):
        print "Error with box dimensions"
        end()
    
    my_distr = np.zeros([box_dim[0][0],box_dim[1][0]])
    
    for i in range(len(genotypes)):
        indx1 = genotypes[i][0] - dim1_data[0] + box_dim[0][1]
        indx2 = genotypes[i][1] - dim2_data[0] + box_dim[1][1]        
        my_distr[indx1,indx2] = max([abundances[0][i]/tot_pop_size,1/tot_pop_size])
        
    xlabels = [(dim1_data[0] - box_dim[0][1]-1 + i) for i in range(box_dim[0][0]+1)]
    ylabels = [(dim2_data[0] - box_dim[1][1]-1 + i) for i in range(box_dim[1][0]+1)]

    for i in range(len(hhfgenotypes)):
        indx1 = hhfgenotypes[i][0] - dim1_data[0] + box_dim[0][1]
        indx2 = hhfgenotypes[i][1] - dim2_data[0] + box_dim[1][1]        
        hhf_points.append([indx1+.5,indx2+0.5])

    for i in range(len(stochgenotypes)):
        indx1 = stochgenotypes[i][0] - dim1_data[0] + box_dim[0][1]
        indx2 = stochgenotypes[i][1] - dim2_data[0] + box_dim[1][1]        
        stoch_points.append([indx1+0.5,indx2+0.5])
        
    hhf_points = np.asarray(hhf_points)
    stoch_points = np.asarray(stoch_points)
    
    return [my_distr,xlabels,ylabels,hhf_points,stoch_points]

def get_normlzd_thry_indv_var(N,s,U):
    sigma1sqrd = 0.25 * (get_vNsU(N,s,2*U)/get_vNsU(N,s,U)) * ( 1 + np.log(s/(2*U)) + (s / np.sqrt(np.pi*get_vNsU(N,s,2*U))))
    return sigma1sqrd

def get_normlzd_thry_cov(N,s,U):
    sigma12 = 0.25 * (get_vNsU(N,s,2*U)/get_vNsU(N,s,U)) * ( 1 - np.log(s/(2*U)) - (s / np.sqrt(np.pi*get_vNsU(N,s,2*U))))
    return sigma12

def get_q(N,s,U):
    q_est = 2*np.log(N*s)/np.log(s/U) 
    return q_est
