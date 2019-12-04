function [adj_time,param,pop,fit,fitx,fity,means_data,summr_data] = get_population_data(input_file_param,input_file_class,input_file_abund,input_file_summr,input_file_means)
% get_population_data take ending distribution data from a previous simlation and
% recreates the arrays to use in a new simulation that will continue the
% prior one.

% input files:
% -input_file_param = set of parameters used in prior simulation
% -input_file_class = set of classes at end of prior simulation
% -input_file_abund = set of abundances at end of prior simulation
% output files:
% -param = array of parameters from prior simulation
% -pop = array of abundances at end of prior simulation
% -fit = array of abundances at end of prior simulation
% -fitx = array of abundances at end of prior simulation
% -fity = array of abundances at end of prior simulation

    pop = 0;
    fit = 0;
    fitx = 0;
    fity = 0;

    % open files holding population data from a prior simulation
    fileID1 = fopen(input_file_param,'r');
    fileID2 = fopen(input_file_class,'r');
    fileID3 = fopen(input_file_abund,'r');
    fileID4 = fopen(input_file_summr,'r');
    fileID5 = fopen(input_file_means,'r');
    
    % reading data from prior simulation
    param_data = fgetl(fileID1);
    class_data = fgetl(fileID2); 
        class_data = class_data(1:end-1);
        class_data = replace(class_data,"],[",";");
    abund_data = fgetl(fileID3); 
        abund_data = abund_data(1:end-1);
    summr_data = fgetl(fileID4);
    means_data = fgetl(fileID5);
        
    fclose(fileID1);
    fclose(fileID2);
    fclose(fileID3);
    fclose(fileID4);
    fclose(fileID5);

    eval(['param = [' param_data '];'])
    eval(['temp_class = ' class_data ';'])
    eval(['temp_abund = [' abund_data '];'])
    eval(['temp_summr = [' summr_data '];'])
    eval(['means_data = [' means_data '];'])

    % set time adjustment
    adj_time = temp_summr(1,1);
    summr_data = temp_summr(1,2:end);
    
    % build fitness arrays
    s1 = param(1,2);
    s2 = param(1,4);
    fitx = min(temp_class(:,1)):1:max(temp_class(:,1));
    fity = min(temp_class(:,2)):1:max(temp_class(:,2));
    dim = [length(fitx) length(fity)];
    fitx_arry = s1*fitx'*ones(1,dim(2));           
    fity_arry = s2*ones(dim(1),1)*fity;            
    fit = fitx_arry + fity_arry;                   
    
    % build array with abundances
    pop = zeros(dim(1),dim(2));
    
    dim=size(temp_class);
    
    for i=1:dim(1)
        i1 = temp_class(i,1)-min(fitx)+1;
        i2 = temp_class(i,2)-min(fity)+1;
        pop(i1,i2) = temp_abund(i);
    end
    
end