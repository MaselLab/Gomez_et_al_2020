%% setup array for parameters N,s1,s2,u1,u2
% Simulate the s-U relationship in the U>s regime with
% deleterious mutations included

N = 1e9;
s = 1e-2;
U = 1e-5;

s_min = 1e-3; s_max = 1e-1;
U_min = 1e-2; U_max = 1e-14;

% the variable start time can be changed to sample a trajectory in detail,
% but currently it is set ot the last time point to sample the distribution
% so that the simulations can be continued from that point.
steps = 5e5;
start_time = steps;                         % collect data on distribution at start time
end_time = steps;    % collect data on distribution at end time

digits(16)
% rng(7);    % set seed for random number generator

n = 31; % must be an odd number to ensure that s that U of trait 1 are used with 2nd trait

sarry = logspace(log10(s_min),log10(s_max),n);
Uarry = logspace(log10(U_min),log10(U_max),n+30);

data_pts_s = length(sarry);    
data_pts_U = length(Uarry);     

number_of_sims = floor(data_pts_s*data_pts_U);  % define trait 1 vs trait 2 grid
collect_distribution_data = ones(number_of_sims,1);
indx_of_collected_data = [];

%% trait one origin fixation regime

% name the output files and location to store them
outputfile = '~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_fixed_sU_OF_ml-200'; 

s1 = 0.07943282347;
U1 = 5.22e-12;

NsU = zeros(number_of_sims,7);          % array that stores the parameters [N,s1,u1,s2,u2]
sim_data = zeros(number_of_sims,6);     % data collected [v,v1,v2,varx,vary,cov]
new_sim_flag = true;  %set to true if new simulations, otherwise set false if continuing simulations
init_flag = false;  % set to true if some continue simulations, otherwise set to false
init_time = 0;
init_pop = N;
init_fit = 0;
init_fitx = 0;
init_fity = 0;
init_means = 0;
init_summr = 0;

% indicate which simulations should be initialized from prior simulations
init_filename = '';
init_indx = [];      % set indx values here

indx = 0;
tic
for i=1:data_pts_s
    for j=1:data_pts_U
        indx = indx + 1
        
        if (new_sim_flag)
            NsU(indx,:)=[N,s1,U1,sarry(i),Uarry(j),U1,Uarry(j)];
            [sim_data(indx,1),sim_data(indx,2),sim_data(indx,3),sim_data(indx,4),sim_data(indx,5),sim_data(indx,6)] ...
                = stochastic_simulation_two_traits(N,s,U,sarry(i),Uarry(j),U,Uarry(j),steps, ...
                collect_distribution_data(indx),start_time,end_time,[outputfile '-' num2str(indx)], ...
                init_flag,init_time,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr);
            
        elseif ( ~new_sim_flag & ( sum(init_indx == indx) == 1 ) )
            input_file_param = [init_filename '-' num2str(indx) '-0.txt'];
            input_file_class = [init_filename '-' num2str(indx) '-2.txt'];
            input_file_abund = [init_filename '-' num2str(indx) '-3.txt'];
            input_file_summr = [init_filename '-' num2str(indx) '-1.txt'];
            input_file_means = [init_filename '-' num2str(indx) '-4.txt'];
            [init_time,init_param,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr] = ...
                get_population_data(input_file_param,input_file_class,input_file_abund,input_file_summr,input_file_means);
            
            NsU(indx,:)=[init_param init_param(1,3) init_param(1,5)];
            [sim_data(indx,1),sim_data(indx,2),sim_data(indx,3),sim_data(indx,4),sim_data(indx,5),sim_data(indx,6)] ...
                = stochastic_simulation_two_traits(NsU(indx,1),NsU(indx,2),NsU(indx,3),NsU(indx,4),NsU(indx,5),NsU(indx,6),NsU(indx,7),steps, ...
                collect_distribution_data(indx),start_time,end_time,[outputfile '-' num2str(indx)], ...
                init_flag,init_time,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr);
        else
            input_file_param = [init_filename '-' num2str(indx) '-0.txt'];
            input_file_class = [init_filename '-' num2str(indx) '-2.txt'];
            input_file_abund = [init_filename '-' num2str(indx) '-3.txt'];
            input_file_summr = [init_filename '-' num2str(indx) '-1.txt'];
            input_file_means = [init_filename '-' num2str(indx) '-4.txt'];
            [init_time,init_param,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr] = ...
                get_population_data(input_file_param,input_file_class,input_file_abund,input_file_summr,input_file_means);
            
            NsU(indx,:)=[init_param init_param(1,3) init_param(1,5)];
            sim_data(indx,:) = init_means;
            
            copyfile([input_filename '-' num2str(indx) '-0.txt'], [outputfile '-' num2str(indx) '-0.txt']);
            copyfile([input_filename '-' num2str(indx) '-1.txt'], [outputfile '-' num2str(indx) '-1.txt']);
            copyfile([input_filename '-' num2str(indx) '-2.txt'], [outputfile '-' num2str(indx) '-2.txt']);
            copyfile([input_filename '-' num2str(indx) '-3.txt'], [outputfile '-' num2str(indx) '-3.txt']);
            copyfile([input_filename '-' num2str(indx) '-4.txt'], [outputfile '-' num2str(indx) '-4.txt']);
                
        end
        
    end
    
end
toc

dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_parameters_fixed_sU_OF_ml-200-0.dat',NsU,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_grand_means_fixed_sU_OF_ml-200-1.dat',sim_data,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_indx_of_collected_data_fixed_sU_OF_ml-200-2.dat',indx_of_collected_data,'delimiter',',','precision',16);

% clean up outputs
save_file = '~/Documents/mutBiasCI/data/FixedsU/compare_sU_sim_data_exp200_OF.mat';
clean_up_output_files(outputfile,number_of_sims,save_file)

%% trait one multiple mutations regime U<<s

% name the output files and location to store them
outputfile = '~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_fixed_sU_DF_ml-200'; 

s1 = 1e-2;
U1 = 1e-5;

NsU = zeros(number_of_sims,7);          % array that stores the parameters [N,s1,u1,s2,u2]
sim_data = zeros(number_of_sims,6);     % data collected [v,v1,v2,varx,vary,cov]
new_sim_flag = true;  %set to true if new simulations, otherwise set false if continuing simulations
init_flag = false;  % set to true if some continue simulations, otherwise set to false
init_time = 0;
init_pop = N;
init_fit = 0;
init_fitx = 0;
init_fity = 0;
init_means = 0;
init_summr = 0;

% indicate which simulations should be initialized from prior simulations
init_filename = '';
init_indx = [];      % set indx values here

indx = 0;
tic
for i=1:data_pts_s
    for j=1:data_pts_U
        indx = indx + 1
        
        if (new_sim_flag)
            NsU(indx,:)=[N,s1,U1,sarry(i),Uarry(j),U1,Uarry(j)];
            [sim_data(indx,1),sim_data(indx,2),sim_data(indx,3),sim_data(indx,4),sim_data(indx,5),sim_data(indx,6)] ...
                = stochastic_simulation_two_traits(N,s,U,sarry(i),Uarry(j),U,Uarry(j),steps, ...
                collect_distribution_data(indx),start_time,end_time,[outputfile '-' num2str(indx)], ...
                init_flag,init_time,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr);
            
        elseif ( ~new_sim_flag & ( sum(init_indx == indx) == 1 ) )
            input_file_param = [init_filename '-' num2str(indx) '-0.txt'];
            input_file_class = [init_filename '-' num2str(indx) '-2.txt'];
            input_file_abund = [init_filename '-' num2str(indx) '-3.txt'];
            input_file_summr = [init_filename '-' num2str(indx) '-1.txt'];
            input_file_means = [init_filename '-' num2str(indx) '-4.txt'];
            [init_time,init_param,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr] = ...
                get_population_data(input_file_param,input_file_class,input_file_abund,input_file_summr,input_file_means);
            
            NsU(indx,:)=[init_param init_param(1,3) init_param(1,5)];
            [sim_data(indx,1),sim_data(indx,2),sim_data(indx,3),sim_data(indx,4),sim_data(indx,5),sim_data(indx,6)] ...
                = stochastic_simulation_two_traits(NsU(indx,1),NsU(indx,2),NsU(indx,3),NsU(indx,4),NsU(indx,5),NsU(indx,6),NsU(indx,7),steps, ...
                collect_distribution_data(indx),start_time,end_time,[outputfile '-' num2str(indx)], ...
                init_flag,init_time,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr);
        else
            input_file_param = [init_filename '-' num2str(indx) '-0.txt'];
            input_file_class = [init_filename '-' num2str(indx) '-2.txt'];
            input_file_abund = [init_filename '-' num2str(indx) '-3.txt'];
            input_file_summr = [init_filename '-' num2str(indx) '-1.txt'];
            input_file_means = [init_filename '-' num2str(indx) '-4.txt'];
            [init_time,init_param,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr] = ...
                get_population_data(input_file_param,input_file_class,input_file_abund,input_file_summr,input_file_means);
            
            NsU(indx,:)=[init_param init_param(1,3) init_param(1,5)];
            sim_data(indx,:) = init_means;
            
            copyfile([input_filename '-' num2str(indx) '-0.txt'], [outputfile '-' num2str(indx) '-0.txt']);
            copyfile([input_filename '-' num2str(indx) '-1.txt'], [outputfile '-' num2str(indx) '-1.txt']);
            copyfile([input_filename '-' num2str(indx) '-2.txt'], [outputfile '-' num2str(indx) '-2.txt']);
            copyfile([input_filename '-' num2str(indx) '-3.txt'], [outputfile '-' num2str(indx) '-3.txt']);
            copyfile([input_filename '-' num2str(indx) '-4.txt'], [outputfile '-' num2str(indx) '-4.txt']);
                
        end
        
    end
    
end
toc

dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_parameters_fixed_sU_DF_ml-200-0.dat',NsU,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_grand_means_fixed_sU_DF_ml-200-1.dat',sim_data,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_indx_of_collected_data_fixed_sU_DF_ml-200-2.dat',indx_of_collected_data,'delimiter',',','precision',16);

save_file = '~/Documents/mutBiasCI/data/FixedsU/compare_sU_sim_data_exp200_DF.mat'
clean_up_output_files(outputfile,number_of_sims,save_file)

%% trait one diffusive mutations regime

% name the output files and location to store them
outputfile = '~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_fixed_sU_HR_ml-200'; 

s = 1e-2;
U = 1e-5;

NsU = zeros(number_of_sims,7);          % array that stores the parameters [N,s1,u1,s2,u2]
sim_data = zeros(number_of_sims,6);     % data collected [v,v1,v2,varx,vary,cov]
new_sim_flag = true;  %set to true if new simulations, otherwise set false if continuing simulations
init_flag = false;  % set to true if some continue simulations, otherwise set to false
init_time = 0;
init_pop = N;
init_fit = 0;
init_fitx = 0;
init_fity = 0;
init_means = 0;
init_summr = 0;

% indicate which simulations should be initialized from prior simulations
init_filename = '';
init_indx = [];      % set indx values here

indx = 0;
tic
for i=1:data_pts_s
    for j=1:data_pts_U
        indx = indx + 1
        
        if (new_sim_flag)
            NsU(indx,:)=[N,s1,U1,sarry(i),Uarry(j),U1,Uarry(j)];
            [sim_data(indx,1),sim_data(indx,2),sim_data(indx,3),sim_data(indx,4),sim_data(indx,5),sim_data(indx,6)] ...
                = stochastic_simulation_two_traits(N,s,U,sarry(i),Uarry(j),U,Uarry(j),steps, ...
                collect_distribution_data(indx),start_time,end_time,[outputfile '-' num2str(indx)], ...
                init_flag,init_time,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr);
            
        elseif ( ~new_sim_flag & ( sum(init_indx == indx) == 1 ) )
            input_file_param = [init_filename '-' num2str(indx) '-0.txt'];
            input_file_class = [init_filename '-' num2str(indx) '-2.txt'];
            input_file_abund = [init_filename '-' num2str(indx) '-3.txt'];
            input_file_summr = [init_filename '-' num2str(indx) '-1.txt'];
            input_file_means = [init_filename '-' num2str(indx) '-4.txt'];
            [init_time,init_param,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr] = ...
                get_population_data(input_file_param,input_file_class,input_file_abund,input_file_summr,input_file_means);
            
            NsU(indx,:)=[init_param init_param(1,3) init_param(1,5)];
            [sim_data(indx,1),sim_data(indx,2),sim_data(indx,3),sim_data(indx,4),sim_data(indx,5),sim_data(indx,6)] ...
                = stochastic_simulation_two_traits(NsU(indx,1),NsU(indx,2),NsU(indx,3),NsU(indx,4),NsU(indx,5),NsU(indx,6),NsU(indx,7),steps, ...
                collect_distribution_data(indx),start_time,end_time,[outputfile '-' num2str(indx)], ...
                init_flag,init_time,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr);
        else
            input_file_param = [init_filename '-' num2str(indx) '-0.txt'];
            input_file_class = [init_filename '-' num2str(indx) '-2.txt'];
            input_file_abund = [init_filename '-' num2str(indx) '-3.txt'];
            input_file_summr = [init_filename '-' num2str(indx) '-1.txt'];
            input_file_means = [init_filename '-' num2str(indx) '-4.txt'];
            [init_time,init_param,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr] = ...
                get_population_data(input_file_param,input_file_class,input_file_abund,input_file_summr,input_file_means);
            
            NsU(indx,:)=[init_param init_param(1,3) init_param(1,5)];
            sim_data(indx,:) = init_means;
            
            copyfile([input_filename '-' num2str(indx) '-0.txt'], [outputfile '-' num2str(indx) '-0.txt']);
            copyfile([input_filename '-' num2str(indx) '-1.txt'], [outputfile '-' num2str(indx) '-1.txt']);
            copyfile([input_filename '-' num2str(indx) '-2.txt'], [outputfile '-' num2str(indx) '-2.txt']);
            copyfile([input_filename '-' num2str(indx) '-3.txt'], [outputfile '-' num2str(indx) '-3.txt']);
            copyfile([input_filename '-' num2str(indx) '-4.txt'], [outputfile '-' num2str(indx) '-4.txt']);
                
        end
        
    end
    
end
toc

dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_parameters_fixed_sU_HR_ml-200-0.dat',NsU,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_grand_means_fixed_sU_HR_ml-200-1.dat',sim_data,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_indx_of_collected_data_fixed_sU_HR_ml-200-2.dat',indx_of_collected_data,'delimiter',',','precision',16);

save_file = '~/Documents/mutBiasCI/data/FixedsU/compare_sU_sim_data_exp200_HR.mat'
clean_up_output_files(outputfile,number_of_sims,save_file)