%% setup array for parameters N,s1,s2,u1,u2
% Simulate the s-U relationship in the U>s regime with
% deleterious mutations included

N = 1e9;
s = 1e-2;
U = 1e-5;
digits(16)

% rng(7);    % set seed for random number generator

% simulations with two traits either run for a fixed R or a fixed v, which
% is determined by the combinations of s and U being fed in. 

% Output file holds the s-U pairs of v contour (stoch. Approx).
% sU = dlmread('data/SAapprox/mutBiasCI_estimate_U_ml-6-2-1.dat',','); % data from figure 3
% sU = dlmread('data/SAapprox/mutBiasCI_estimate_U_ml-22-1-1.dat',','); % data from figure 3 (bottom contour)
sU = dlmread('data/SAapprox/mutBiasCI_estimate_U_ml-100-2-1.dat',','); % data from figure 3

% the selected s and U pairs below give either v=5.308e-5 or R=5.308e-3 with N=1e9
% depending on the flag used (fixed R for rate_flag==0, and fixed v for R==1)
sarry = sU(:,1);        
Uarry = sU(:,2);

% name the output files and location to store them
outputfile = '~/Documents/mutBiasCI/data/TimeSeries/mutBiasCI_data_for_2d_distribution_ml-102'; 

data_pts_s = length(sarry);     % define the grid size, given by number of s values sampled.
number_of_sims = floor(0.5*data_pts_s*(data_pts_s+1));  % define trait 1 vs trait 2 grid

collect_distribution_data = ones(number_of_sims,1);
indx_of_collected_data = [];

% the variable start time can be changed to sample a trajectory in detail,
% but currently it is set ot the last time point to sample the distribution
% so that the simulations can be continued from that point.
steps = 2.0e6;
start_time = steps;                         % collect data on distribution at start time
end_time = steps*ones(number_of_sims,1);    % collect data on distribution at end time

NsU = zeros(number_of_sims,7);          % array that stores the parameters [N,s1,u1,s2,u2]
sim_data = zeros(number_of_sims,6);     % data collected [v,v1,v2,varx,vary,cov]

% % ARRANGMENT OF SIMULATIONS BY indx (comment out if not used)
% indx = 0;
% sim_table = zeros(data_pts_s,data_pts_s);
% for i=1:data_pts_s
%     for j=i:data_pts_s
%         indx = indx + 1
%         sim_table(i,j)=indx;
%     end
% end
% heatmap(sim_table,'XData',log10(sarry),'YData',round(log10(Uarry),1),'XLabel','s-values','YLabel','U-values');

% The original simulations were setup to run trait 1 with index i, and
% trait 2 with index j. However, now trait 2 is index i and trait 1 is
% index j because figure 3 was switched from lower diagonal to upper
% diagonal.

init_flag = true;
init_time = 0;
init_pop = N;
init_fit = 0;
init_fitx = 0;
init_fity = 0;

% indicate which simulations should be initialized from prior simulations
init_filename = '~/Documents/mutBiasCI/data/TimeSeries/mutBiasCI_data_for_2d_distribution_ml-101';
init_indx = [10:15 24:29 37:42 49:54 60:65 70:75 79:84 87:92 94:120];      % set indx values here

indx = 0;
tic
for i=1:data_pts_s
    for j=i:data_pts_s
        indx = indx + 1
        
        % initial population data from a prior simulation
        if ( sum(init_indx == indx) == 1 )
            input_file_param = [init_filename '-' num2str(indx) '-0.txt'];
            input_file_class = [init_filename '-' num2str(indx) '-2.txt'];
            input_file_abund = [init_filename '-' num2str(indx) '-3.txt'];
            input_file_summr = [init_filename '-' num2str(indx) '-1.txt'];
            input_file_means = [init_filename '-' num2str(indx) '-4.txt'];
            [init_time,init_param,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr] = ...
                get_population_data(input_file_param,input_file_class,input_file_abund,input_file_summr,input_file_means);
            
            NsU(indx,:)=[N,sarry(i),Uarry(i),sarry(j),Uarry(j),Uarry(i),Uarry(j)];
            [sim_data(indx,1),sim_data(indx,2),sim_data(indx,3),sim_data(indx,4),sim_data(indx,5),sim_data(indx,6)] ...
                = stochastic_simulation_two_traits(N,sarry(i),Uarry(i),sarry(j),Uarry(j),Uarry(i),Uarry(j),steps, ...
                collect_distribution_data(indx),start_time,end_time(indx),[outputfile '-' num2str(indx)], ...
                init_flag,init_time,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr);
        else
            input_file_param = [init_filename '-' num2str(indx) '-0.txt'];
            input_file_class = [init_filename '-' num2str(indx) '-2.txt'];
            input_file_abund = [init_filename '-' num2str(indx) '-3.txt'];
            input_file_summr = [init_filename '-' num2str(indx) '-1.txt'];
            input_file_means = [init_filename '-' num2str(indx) '-4.txt'];
            [init_time,init_param,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr] = ...
                get_population_data(input_file_param,input_file_class,input_file_abund,input_file_summr,input_file_means);
            
            NsU(indx,:)=[N,sarry(i),Uarry(i),sarry(j),Uarry(j),Uarry(i),Uarry(j)];
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

dlmwrite('~/Documents/mutBiasCI/data/TwoTraitSim/mutBiasCI_data_all_simulation_parameters_ml-102-0.dat',NsU,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/TwoTraitSim/mutBiasCI_data_all_simulation_grand_means_ml-102-1.dat',sim_data,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/TwoTraitSim/mutBiasCI_data_all_simulation_indx_of_collected_data_ml-102-2.dat',indx_of_collected_data,'delimiter',',','precision',16);

%% script to copy old simulations to continue simulations

% input_filename = '~/Documents/mutBiasCI/data/TimeSeries/mutBiasCI_data_for_2d_distribution_ml-100';
% output_filename = '~/Documents/mutBiasCI/data/TimeSeries/mutBiasCI_data_for_2d_distribution_ml-101';
% init_indx = [10:15 24:29 37:42 49:54 60:65 70:75 79:84 87:92 94:120];      % set indx values here
% 
% indx = 0;
% for i=1:data_pts_s
%     for j=i:data_pts_s
%         indx = indx + 1
%         
%         % copy and rename files for simulations that were not continued to
%         % complete the set of simulations.
%         if ( sum(init_indx == indx) == 0 )
%             copyfile([input_filename '-' num2str(indx) '-0.txt'], [output_filename '-' num2str(indx) '-0.txt']);
%             copyfile([input_filename '-' num2str(indx) '-1.txt'], [output_filename '-' num2str(indx) '-1.txt']);
%             copyfile([input_filename '-' num2str(indx) '-2.txt'], [output_filename '-' num2str(indx) '-2.txt']);
%             copyfile([input_filename '-' num2str(indx) '-3.txt'], [output_filename '-' num2str(indx) '-3.txt']);
%             copyfile([input_filename '-' num2str(indx) '-4.txt'], [output_filename '-' num2str(indx) '-4.txt']);
%         end
%         
%     end
%     
% end
