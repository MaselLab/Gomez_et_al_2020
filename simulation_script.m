%% setup array for parameters N,s1,s2,u1,u2
% Simulate the s-U relationship in the discontinuous regime with
% deleterious mutations included

N = 1e9;
s = 1e-2;
U = 1e-5;
v = s^2*(2*log(N*s)-log(s/U))/(log(s/U)^2);                 % extend this to a range of v
digits(16)

rng(7);                                                     % set seed for random number generator

% s and U pairs that give v=5.308e-5 with N=1e9
sU = dlmread('data/mutBiasCI_estimate_U_ml-6-2-1.dat',',');
sarry = sU(:,1);
Uarry = sU(:,2);

steps = 1.0e6;
start_time = 5e3;                     % collect data on distribution at start time
end_time = 5e3;                       % collect data on distribution at end time
outputfile = '~/Documents/mutBiasCI/data/mutBiasCI_data_for_2d_distribution_ml-16'; 

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
            = stochastic_simulation_two_traits2(N,sarry(i),Uarry(i),sarry(j),Uarry(j),Uarry(i),Uarry(j),steps, ...
                collect_distribution_data(indx),start_time,end_time,[outputfile '-' num2str(indx)]);
    end
    
end
toc

dlmwrite('~/Documents/mutBiasCI/data/mutBiasCI_data_all_simulation_parameters_ml-16-0.dat',NsU,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/mutBiasCI_data_all_simulation_grand_means_ml-16-1.dat',sim_data,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/mutBiasCI_data_all_simulation_indx_of_collected_data_ml-16-2.dat',indx_of_collected_data,'delimiter',',','precision',16);
