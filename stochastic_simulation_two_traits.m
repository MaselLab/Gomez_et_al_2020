function [v,v1,v2,varx,vary,cov] = stochastic_simulation_two_traits(N,s1,u1,s2,u2,ud1,ud2,steps,...
                                    collect_data,start_time,end_time,outputfile,init_flag,init_time,init_pop,init_fit,init_fitx,init_fity,init_means,init_summr)
% The code below has been modified from the source code made
% availabe by Pearce MT and Fisher DS, obtained from:
% 
% https://datadryad.org/resource/doi:10.5061/dryad.f36v6 
%
% Stochastic simulation for two chromosome model. 
% Each step is one generation. Rates (such as s,u,r) are per generation. 
% The expected size of each subpop after selection, mutation, and mating is computed.
% If expected size < :cutoff: then the size is drawn from Poisson
% distribution with mean equal to the expected size. Otherwise the
% subpopulation size = expected size. 

% output :v, v1, v2:        total rate of adaptation, and rate of adaptation in traits 1 and 2. 
% output :varx, vary, cov:  time averaged variances and covariance in traits 1 and 2.

% input :N: population size
% input :s1: effect size of beneficial mutation in trait 1. 
% input :u1: mutation rate per locus of trait 1.
% input :s2: effect size of beneficial mutation in trait 2. 
% input :u2: mutation rate per locus of trait 2.
% input :ud1: deleterious mutation rate per locus of trait 1.
% input :ud2: deleterious mutation rate per locus of trait 2.
% input :steps: number of steps for simulation.
% input :collect_data: true/false - collect detailed data on 2d distr. per generation
% input :start_time: start time for collecting detailed data on 2d distribution
% input :end_time: end time for collecting detailed data on 2d distribution
% input :init_flag: inicates if simulations should be initialized prior one
% input :init_time: end time from prior simulation
% input :init_pop: abundances from prior simulation
% input :init_fit: total fitness from prior simulation
% input :init_fitx: trait one fitness from prior simulation
% input :init_fity: trait two fitness from prior simulation
% input :outputfile: string with filename where detailed data will be stored

digits(16);

% initialize variables
pop=N;                      % abundances of a classes
Na = 0;                     % actual population size
fit=0;                      % total fitness of a class
fitx=0;                     % fitness in trait 1 of a class
fity=0;                     % fitness in trait 2 of a class
nosefitness = 0;            % total fitness of the front
meanfitness = 0;            % mean fitness of the population
meanfitx = 0;               % mean fitness in x 
meanfity = 0;               % mean fitness in y
meanfit_s = 0;              % mean fitness of the population
meanfitx_s = 0;             % mean fitness in x 
meanfity_s = 0;             % mean fitness in y
varx = 0;                   % variance in trait 1
vary = 0;                   % variance in trait 2
cov = 0;                    % covariance between trait 1 and 2
adj_time = 0;               % adjustment to timestamp (cont. simulations)
burn_time = 1;

if (s2 == 0)    % simulation with one trait should be done with setting s2 = 0, U2=0, U2d=0
    cutoff = 10/s1;     % cutoff for stochastic dynamics should be set by first s
else
    cutoff = max([10/s1 10/s2]);       % cutoff for stochastic dynamics will use highest of the two cuttoffs
end

if (init_flag)  % initialize population details to those from a prior simulation
    pop = init_pop;
    adj_time = init_time;
    fit = init_fit;
    fitx = init_fitx;
    fity = init_fity;
    varx = init_means(1,4);
    vary = init_means(1,5);
    cov = init_means(1,6);
    meanfit_s = init_summr(1,4) - init_means(1,1)*(adj_time-burn_time);
    meanfitx_s = init_summr(1,5) - init_means(1,2)*(adj_time-burn_time);
    meanfity_s = init_summr(1,6) - init_means(1,3)*(adj_time-burn_time);
end

if (collect_data)           % store parameters used in simulation
    fileID = fopen([outputfile '-0.txt'],'w');
    fprintf(fileID,'%.16f,%.16f,%.16f,%.16f,%.16f',N, s1, s2, u1, u2);
    fclose(fileID);
    fileID1 = fopen([outputfile '-1.txt'],'w'); %file for all other 2d wave data per generation
    fileID2 = fopen([outputfile '-2.txt'],'w'); %file for data on classes per generation
    fileID3 = fopen([outputfile '-3.txt'],'w'); %file for data on abundances per generation
    fileID4 = fopen([outputfile '-4.txt'],'w'); %file for grand means
end

% Main loop for simulation of each generation
for timestep=1:steps   

    %%%%%%%%%%%%%%%%%%%
    % Remove columns of zeros of decreasing fitness classes
    while any(pop(:,1))==0  % remove y-fitness classes
        pop(:,1)=[];
        fit(:,1)=[];
        fity(1)=[];
    end
    
    while any(pop(1,:))==0  % remove x-fitness classes
        pop(1,:)=[];
        fit(1,:)=[];
        fitx(1)=[];
    end
    
    % Add columns for padding, i.e. for new class produced by mutations
    dim=size(pop);
    if any(pop(:,dim(2)))==1    % check for expansion of front in direction of trait two (y)
        pop(:,dim(2)+1)=zeros(dim(1),1);
        fit(:,dim(2)+1)=fit(:,dim(2))+ones(dim(1),1);
        fity(dim(2)+1)=fity(dim(2))+1;
    end
    
    dim=size(pop);
    if any(pop(dim(1),:))==1    % check for expansion of front in direction of trait one (x)
        pop(dim(1)+1,:)=zeros(1,dim(2));
        fit(dim(1)+1,:)=fit(dim(1),:)+ones(1,dim(2));
        fitx(dim(1)+1)=fitx(dim(1))+1;
    end
    
    %%%%%%%%%%%% 
    % Find expected frequencies after selection and mutation
    dim=size(pop);  
    freq=pop/N;                                 % array with approximate frequencies of each class
    
    % calculate growth due to selection
    fitx_arry = s1*fitx'*ones(1,dim(2));           % array with # of mutations in trait 1
    fity_arry = s2*ones(dim(1),1)*fity;            % array with # of mutations in trait 2
    fit = fitx_arry + fity_arry;                   % total fitness

    % not actual mean fitness since freq are divided by N. When actual
    % population size exceeds N, then value results in decrese of popsize.
    % This allows selection to control population size.
    meanfitness = sum(sum(times(freq,fit)));        
    
    newpop=times(exp(fit-meanfitness),pop);  % after selection
    newpop=N*newpop/sum(sum(newpop));        
    
    % calculate changes in abundances due to mutations
    
    % beneficial mutations
    z1=zeros(1,dim(2));
    mutatex=[z1; newpop];
    mutatex(dim(1)+1,:)=[];                     % newfreq already has padding, get rid of extra padding from shift
    
    z2=zeros(dim(1),1);
    mutatey=[z2 newpop];
    mutatey(:,dim(2)+1)=[];                     % newfreq already has padding, get rid of extra padding from shift
    
    % deletirious mutations (note: del mutations in lowest classes ignored)
    z1=zeros(1,dim(2));
    mutatexdel=[newpop; z1];
    mutatexdel(1,:)=[];                     % newfreq already has padding, get rid of extra padding from shift
    
    z2=zeros(dim(1),1);
    mutateydel=[newpop z2];
    mutateydel(:,1)=[];                     % newfreq already has padding, get rid of extra padding from shift
    
    % Sum of the mutation rates cannot exceed 
    nomutate=(1-u1-u2-ud1-ud2)*newpop;
    newpop=nomutate+u1*mutatex+u2*mutatey+ud1*mutatexdel+ud2*mutateydel;    
    
    % For subpopulations with size less than the stoch_cutoff, draw size
    % from poisson distribution. Otherwise, size = N*(expected frequency).
    
    stoch=newpop<cutoff;
    stochpop=poissrnd(newpop(stoch));           % sample poisson for abundances below stochasticity cutoff
    newpop(stoch)=stochpop;
    
    newpop=round(newpop);    
    Na = sum(sum(newpop));
    Nas = sum(sum(newpop(stoch)));
    
    meanfitness = sum(sum(times(newpop,fit)))/Na;
    meanfitx = sum(sum(times(newpop,fitx_arry)))/Na;
    meanfity = sum(sum(times(newpop,fity_arry)))/Na;
    
    pop = newpop;
    
    % recompute time-average of variances and covariances
    if (timestep+adj_time > burn_time)
        varx = ( 1/(timestep+adj_time) )*((timestep+adj_time-1)*varx + sum(sum(times(newpop,(fitx_arry-meanfitx).^2)))/Na);
        vary = ( 1/(timestep+adj_time) )*((timestep+adj_time-1)*vary + sum(sum(times(newpop,(fity_arry-meanfity).^2)))/Na);
        cov = ( 1/(timestep+adj_time) )*((timestep+adj_time-1)*cov + sum(sum(times(newpop,(fitx_arry-meanfitx).*(fity_arry-meanfity))))/Na);
    else
        if (timestep+adj_time)==burn_time
            varx = sum(sum(times(newpop,(fitx_arry-meanfitx).^2)))/Na;
            vary = sum(sum(times(newpop,(fity_arry-meanfity).^2)))/Na;
            cov = sum(sum(times(newpop,(fitx_arry-meanfitx).*(fity_arry-meanfity))))/Na;
            meanfit_s = meanfitness;
            meanfitx_s = meanfitx;
            meanfity_s = meanfity;
        end
    end
    
    if( collect_data && (timestep >= start_time) && (timestep <= end_time) )
        % compute variances, covarainces and population load
        sigmax2 = sum(sum(times(newpop,(fitx_arry-meanfitx).^2)))/Na;
        sigmay2 = sum(sum(times(newpop,(fity_arry-meanfity).^2)))/Na;
        sigmaxy = sum(sum(times(newpop,(fitx_arry-meanfitx).*(fity_arry-meanfity))))/Na;

        % print data to output files, need: times,mean_fit,fit_var,fit cov
        fprintf(fileID1,'%i,%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n',timestep+adj_time,sigmax2,sigmay2,sigmaxy,meanfitness,meanfitx,meanfity);
        
        for i=1:size(pop,1)
            for j=1:size(pop,2)
                if(pop(i,j)>0)
                    fprintf(fileID2,'[%i,%i],',fitx(i),fity(j));
                    fprintf(fileID3,'%f,',pop(i,j));
                end
            end
        end
        
        fprintf(fileID2,'\n');
        fprintf(fileID3,'\n');
    end
    
end

v = (meanfitness-meanfit_s)/(steps+adj_time-burn_time);
v1 = (meanfitx-meanfitx_s)/(steps+adj_time-burn_time);
v2 = (meanfity-meanfity_s)/(steps+adj_time-burn_time);

fprintf(fileID4,'%.16e,%.16e,%.16e,%.16e,%.16e,%.16e\n',v,v1,v2,varx,vary,cov);

% close output files
if(collect_data)
    fclose(fileID1);
    fclose(fileID2);
    fclose(fileID3);
    fclose(fileID4);
end

end