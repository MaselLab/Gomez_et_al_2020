%% setup array for parameters N,s1,s2,u1,u2
% Simulate the s-U relationship in the U>s regime with
% deleterious mutations included

N = 1e9;

%trait one parameters
s = 1e-2;
U = 1e-5;

s1 = 0.07943282347;
U1 = 5.22e-12;

s_min = s/10; s_max = s*10;
U_min = U/1000; U_max = U*1000;

digits(16)
% rng(7);    % set seed for random number generator

n = 31; % must be an odd number to ensure that s that U of trait 1 are used with 2nd trait
grid = -(n-1)/2:1:(n-1)/2;

sarry = logspace(log10(s_min),log10(s_max),n);
% Uarry = logspace(log10(U_min),log10(U_max),n);
Uarry = logspace(-14.0,log10(U_min)-.2,30);


% name the output files and location to store them
outputfile = '~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_fixed_sU_OF_ml-02'; 

data_pts_s = length(sarry);    
data_pts_U = length(Uarry);     

number_of_sims = floor(data_pts_s*data_pts_U);  % define trait 1 vs trait 2 grid

collect_distribution_data = ones(number_of_sims,1);
indx_of_collected_data = [];

% the variable start time can be changed to sample a trajectory in detail,
% but currently it is set ot the last time point to sample the distribution
% so that the simulations can be continued from that point.
steps = 1.5e5;
start_time = steps;                         % collect data on distribution at start time
end_time = steps;    % collect data on distribution at end time

NsU = zeros(number_of_sims,7);          % array that stores the parameters [N,s1,u1,s2,u2]
sim_data = zeros(number_of_sims,6);     % data collected [v,v1,v2,varx,vary,cov]

% ARRANGMENT OF SIMULATIONS BY indx (comment out if not used)
% indx = 0;
% sim_table = zeros(data_pts_s,data_pts_U);
% for i=1:data_pts_s
%     for j=1:data_pts_U
%         indx = indx + 1
%         sim_table(i,j)=indx;
%     end
% end
% heatmap(sim_table,'XData',round(log10(sarry),2),'YData',round(log10(Uarry),1),'XLabel','s-values','YLabel','U-values');

% The original simulations were setup to run trait 1 with index i, and
% trait 2 with index j. However, now trait 2 is index i and trait 1 is
% index j because figure 3 was switched from lower diagonal to upper
% diagonal.

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
init_indx = [1:data_pts_s*data_pts_U];      % set indx values here

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

dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_parameters_fixed_sU_OF_ml-02-0.dat',NsU,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_grand_means_fixed_sU_OF_ml-02-1.dat',sim_data,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_indx_of_collected_data_fixed_sU_OF_ml-02-2.dat',indx_of_collected_data,'delimiter',',','precision',16);



%% 
N = 1e9;

%trait one parameters
s = 1e-2;
U = 1e-5;

s1 = 1e-2;
U1 = 1e-5;

s_min = s/10; s_max = s*10;
U_min = U/1000; U_max = U*1000;

digits(16)
% rng(7);    % set seed for random number generator

n = 31; % must be an odd number to ensure that s that U of trait 1 are used with 2nd trait
grid = -(n-1)/2:1:(n-1)/2;

sarry = logspace(log10(s_min),log10(s_max),n);
% Uarry = logspace(log10(U_min),log10(U_max),n);
Uarry = logspace(-14.0,log10(U_min)-.2,30);

% name the output files and location to store them
outputfile = '~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_fixed_sU_OF_ml-01'; 

data_pts_s = length(sarry);    
data_pts_U = length(Uarry);     

number_of_sims = floor(data_pts_s*data_pts_U);  % define trait 1 vs trait 2 grid

collect_distribution_data = ones(number_of_sims,1);
indx_of_collected_data = [];

% the variable start time can be changed to sample a trajectory in detail,
% but currently it is set ot the last time point to sample the distribution
% so that the simulations can be continued from that point.
steps = 1.5e5;
start_time = steps;                         % collect data on distribution at start time
end_time = steps;    % collect data on distribution at end time

NsU = zeros(number_of_sims,7);          % array that stores the parameters [N,s1,u1,s2,u2]
sim_data = zeros(number_of_sims,6);     % data collected [v,v1,v2,varx,vary,cov]

% ARRANGMENT OF SIMULATIONS BY indx (comment out if not used)
% indx = 0;
% sim_table = zeros(data_pts_s,data_pts_U);
% for i=1:data_pts_s
%     for j=1:data_pts_U
%         indx = indx + 1
%         sim_table(i,j)=indx;
%     end
% end
% heatmap(sim_table,'XData',round(log10(sarry),2),'YData',round(log10(Uarry),1),'XLabel','s-values','YLabel','U-values');

% The original simulations were setup to run trait 1 with index i, and
% trait 2 with index j. However, now trait 2 is index i and trait 1 is
% index j because figure 3 was switched from lower diagonal to upper
% diagonal.

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
init_indx = [1:data_pts_s*data_pts_U];      % set indx values here

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

dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_parameters_fixed_sU_DF_ml-02-0.dat',NsU,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_grand_means_fixed_sU_DF_ml-02-1.dat',sim_data,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_indx_of_collected_data_fixed_sU_DF_ml-02-2.dat',indx_of_collected_data,'delimiter',',','precision',16);


%%

N = 1e9;

%trait one parameters
s = 1e-2;
U = 1e-5;

s1 = 0.001995262315;
U1 = 0.004466209113;

s_min = s/10; s_max = s*10;
U_min = U/1000; U_max = U*1000;

digits(16)
% rng(7);    % set seed for random number generator

n = 31; % must be an odd number to ensure that s that U of trait 1 are used with 2nd trait
grid = -(n-1)/2:1:(n-1)/2;

sarry = logspace(log10(s_min),log10(s_max),n);
% Uarry = logspace(log10(U_min),log10(U_max),n);
Uarry = logspace(-14.0,log10(U_min)-.2,30);

% name the output files and location to store them
outputfile = '~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_fixed_sU_HR_ml-02'; 

data_pts_s = length(sarry);    
data_pts_U = length(Uarry);     

number_of_sims = floor(data_pts_s*data_pts_U);  % define trait 1 vs trait 2 grid

collect_distribution_data = ones(number_of_sims,1);
indx_of_collected_data = [];

% the variable start time can be changed to sample a trajectory in detail,
% but currently it is set ot the last time point to sample the distribution
% so that the simulations can be continued from that point.
steps = 1.5e5;
start_time = steps;                         % collect data on distribution at start time
end_time = steps;    % collect data on distribution at end time

NsU = zeros(number_of_sims,7);          % array that stores the parameters [N,s1,u1,s2,u2]
sim_data = zeros(number_of_sims,6);     % data collected [v,v1,v2,varx,vary,cov]

% ARRANGMENT OF SIMULATIONS BY indx (comment out if not used)
% indx = 0;
% sim_table = zeros(data_pts_s,data_pts_U);
% for i=1:data_pts_s
%     for j=1:data_pts_U
%         indx = indx + 1
%         sim_table(i,j)=indx;
%     end
% end
% heatmap(sim_table,'XData',round(log10(sarry),2),'YData',round(log10(Uarry),1),'XLabel','s-values','YLabel','U-values');

% The original simulations were setup to run trait 1 with index i, and
% trait 2 with index j. However, now trait 2 is index i and trait 1 is
% index j because figure 3 was switched from lower diagonal to upper
% diagonal.

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
init_indx = [1:data_pts_s*data_pts_U];      % set indx values here

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

dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_parameters_fixed_sU_HR_ml-02-0.dat',NsU,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_grand_means_fixed_sU_HR_ml-02-1.dat',sim_data,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/FixedsU/mutBiasCI_data_all_simulation_indx_of_collected_data_fixed_sU_HR_ml-02-2.dat',indx_of_collected_data,'delimiter',',','precision',16);