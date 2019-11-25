%% setup array for parameters N,s1,s2,u1,u2
% Simulate the s-U relationship in the U>s regime with
% deleterious mutations included

N = 1e9;
s = 1e-2;
U = 1e-5;
digits(16)

% rng(7);                                                     % set seed for random number generator

% simulations with two traits either run for a fixed R or a fixed v, which
% is determined by the combinations of s and U being fed in. 

sU = dlmread('data/SAapprox/mutBiasCI_estimate_U_ml-6-2-1.dat',','); % this is the original data set for figure 3
% sU = dlmread('data/SAapprox/mutBiasCI_estimate_U_ml-22-1-1.dat',','); % this is the new data for getting figure 3 along the bottom contour

% the selected s and U pairs below give either v=5.308e-5 or R=5.308e-3 with N=1e9
% depending on the flag used (fixed R for rate_flag==0, and fixed v for R==1)
sarry = sU(1:13,1);        
Uarry = sU(1:13,2);

steps = 5.0e5;
start_time = 5e3;                     % collect data on distribution at start time
end_time = 5e3;                       % collect data on distribution at end time
outputfile = '~/Documents/mutBiasCI/data/mutBiasCI_data_for_2d_distribution_ml-22'; 

data_pts_s = length(sarry);
number_of_sims = floor(0.5*data_pts_s*(data_pts_s+1));

collect_distribution_data = zeros(number_of_sims,1);
indx_of_collected_data = [];

NsU = zeros(number_of_sims,7);          % array that stores the parameters [N,s1,u1,s2,u2]
sim_data = zeros(number_of_sims,6);     % data collected [v,v1,v2,varx,vary,cov]

indx = 0;
tic
for i=1:data_pts_s
    for j=i:data_pts_s
        indx = indx + 1
        NsU(indx,:)=[N,sarry(i),Uarry(i),sarry(j),Uarry(j),Uarry(i),Uarry(j)];
        [sim_data(indx,1),sim_data(indx,2),sim_data(indx,3),sim_data(indx,4),sim_data(indx,5),sim_data(indx,6)] ...
            = stochastic_simulation_two_traits(N,sarry(i),Uarry(i),sarry(j),Uarry(j),Uarry(i),Uarry(j),steps, ...
                collect_distribution_data(indx),start_time,end_time,[outputfile '-' num2str(indx)]);
    end
    
end
toc

dlmwrite('~/Documents/mutBiasCI/data/mutBiasCI_data_all_simulation_parameters_ml-22-0.dat',NsU,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/mutBiasCI_data_all_simulation_grand_means_ml-22-1.dat',sim_data,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/mutBiasCI_data_all_simulation_indx_of_collected_data_ml-22-2.dat',indx_of_collected_data,'delimiter',',','precision',16);
