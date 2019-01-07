function [v,v1,v2,varx,vary,cov] = stochastic_simulation_two_traits(N,s1,u1,s2,u2,steps,...
                                    collect_data,start_time,end_time,outputfile)
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
% input :steps: number of steps for simulation.
% input :collect_data: true/false - collect detailed data on 2d distr. per generation
% input :start_time: start time for collecting detailed data on 2d distribution
% input :end_time: end time for collecting detailed data on 2d distribution
% input :outputfile: string with filename where detailed data will be stored

% initialize variables
pop=N;                      %abundances of a classes
Na = 0;                     %actual population size
fit=0;                      %total fitness of a class
fitx=0;                     %fitness in trait 1 of a class
fity=0;                     %fitness in trait 2 of a class
nosefitness = 0;            %total fitness of the front
meanfitness = 0;            %mean fitness of the population
meanfitx = 0;               %mean fitness in x 
meanfity = 0;               %mean fitness in y
varx = 0;                   %variance in trait 1
vary = 0;                   %variance in trait 2
cov = 0;                    %covariance between trait 1 and 2
cutoff=10/min(s1,s2);       %population cutoff for stochasticity

if (collect_data)   %store parameters used in simulation
    fileID = fopen([outputfile '-0.txt'],'w');
    fprintf(fileID,'%f,%f,%f,%f,%f',N, s1, s2, u1, u2);
    fclose(fileID);
end

fileID1 = fopen([outputfile '-1.txt'],'w'); %file for all other 2d wave data per generation
fileID2 = fopen([outputfile '-2.txt'],'w'); %file for data on classes per generation
fileID3 = fopen([outputfile '-3.txt'],'w'); %file for data on abundances per generation

% Main loop for simulation of each generation
for timestep=1:steps   
        
    %%%%%%%%%%%%%%%%%%%
    % Remove columns of zeros of decreasing fitness classes
    while any(pop(:,1))==0
        pop(:,1)=[];
        fit(:,1)=[];
        fity(1)=[];
    end
    
    while any(pop(1,:))==0 
        pop(1,:)=[];
        fit(1,:)=[];
        fitx(1)=[];
    end
    
    % Add columns for padding, i.e. for new class produced by mutations
    dim=size(pop); 
    if any(pop(:,dim(2)))==1    % check for expansion of front in direction of trait 1
        pop(:,dim(2)+1)=zeros(dim(1),1);
        fit(:,dim(2)+1)=fit(:,dim(2))+ones(dim(1),1);
        fity(dim(2)+1)=fity(dim(2))+1;
    end
    
    dim=size(pop);
    if any(pop(dim(1),:))==1    % check for expansion of front in direction of trait 2
        pop(dim(1)+1,:)=zeros(1,dim(2));
        fit(dim(1)+1,:)=fit(dim(1),:)+ones(1,dim(2));
        fitx(dim(1)+1)=fitx(dim(1))+1;
    end
    
    %%%%%%%%%%%% 
    % Find expected frequencies after selection and mutation
    dim=size(pop);  
    freq=pop/N;                                 % array with frequencies of the each class
    
    % calculate growth due to selection
    fitx_arry = fitx'*ones(1,dim(2));           % array with # of mutations in trait 1
    fity_arry = ones(dim(1),1)*fity;            % array with # of mutations in trait 2
    fit = s1*fitx_arry + s2*fity_arry;          % total fitness
    meanfitness = sum(sum(times(freq,fit)));
    meanfit_arry = meanfitness*ones(size(fit));
    
    newfreq=times(exp(fit-meanfit_arry),freq);  %after selection
    newfreq=newfreq/sum(sum(newfreq));          %make sure frequencies still add to one.
    
    % calculate changes in abundances due to mutations
    z1=zeros(dim(1),1);
    mutatex=[z1 newfreq];
    mutatex(:,dim(2)+1)=[];                     % newfreq already has padding, get rid of extra padding from shift
    
    z2=zeros(1,dim(2));
    mutatey=[z2; newfreq];
    mutatey(dim(1)+1,:)=[];                     % newfreq already has padding, get rid of extra padding from shift
    
    nomutate=(1-u1-u2)*newfreq;
    postmutate=nomutate+u1*mutatex+u2*mutatey;    
    newfreq=nomutate+postmutate;
    
    % For subpopulations with size less than the stoch_cutoff, draw size
    % from poisson distribution. Otherwise, size = N*(expected frequency).
    
    newpop=N*newfreq;
    stoch=newpop<cutoff;
    stochpop=poissrnd(newpop(stoch));           % sample poisson for abundances below stochasticity cutoff
    newpop(stoch)=stochpop;
    
    newpop=round(newpop);    
    Na = sum(sum(newpop));
    
    meanfitness = sum(sum(times(newpop,fit)))/Na;
    meanfitx = s1*sum(sum(newpop,2).*fitx')/Na;
    meanfity = s2*sum(sum(newpop,1).*fity)/Na;
    
    nosefitness = max(max(times(fit,sign(newpop))));    % calculate most fitness of most fit class
    
    % recompute time-average of variances and covariances
    if timestep == 1
        varx = sum(sum(times(newpop,(fitx_arry-meanfitx).^2)))/Na;
        vary = sum(sum(times(newpop,(fity_arry-meanfity).^2)))/Na;
        cov = sum(sum(times(newpop,(fitx_arry-meanfitx).*(fity_arry-meanfity(timestep)))))/Na;
    else
        varx = (1/timestep)*((timestep-1)*varx + sum(sum(times(newpop,(fitx_arry-meanfitx).^2)))/Na);
        vary = (1/timestep)*((timestep-1)*vary + sum(sum(times(newpop,(fity_arry-meanfity).^2)))/Na);
        cov = (1/timestep)*((timestep-1)*cov + sum(sum(times(newpop,(fitx_arry-meanfitx).*(fity_arry-meanfity))))/Na);
    end
    
    pop=newpop;
    
    if( collect_data && (timestep >= start_time) && (timestep <= end_time) )
        
        % compute the covariance using only classes at the front
        indx_front = (fit==nosefitness);
        meanfitx_front = sum(pop(indx_front).*fitx_arry(indx_front))/sum(pop(indx_front));
        meanfity_front = sum(pop(indx_front).*fity_arry(indx_front))/sum(pop(indx_front));
        front_cov = sum(pop(indx_front).*(fitx_arry(indx_front)-meanfitx_front).*(fity_arry(indx_front)-meanfity_front))/sum(pop(indx_front));
        
        % compute variances, covarainces and population load
        sigmax2 = s1^2*sum(sum(times(newpop,(fitx_arry-meanfitx).^2)))/Na;
        sigmay2 = s2^2*sum(sum(times(newpop,(fity_arry-meanfity).^2)))/Na;
        sigmaxy = s1*s2*sum(sum(times(newpop,(fitx_arry-meanfitx).*(fity_arry-meanfity))))/Na;
        pop_load = nosefitness - meanfitness;
        
        % compute v1, v2, sigma12, G eigenvalues and orientation 
        [D,L] = eig([sigmax2 sigmaxy; sigmaxy sigmay2]);
        evec1 = [cosd(45) sind(45); -sind(45) cosd(45)]*(sign(D(1,1))*D(:,1)); 
        Gang = atan2d(evec1(2),evec1(1)); 
        
        % print data to output files, need: times,mean_fit,fit_var,fit_cov,pop_load,dcov_dt,vU_thry,v2U_thry
        fprintf(fileID1,'%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n',timestep,sigmax2,sigmay2,sigmaxy,front_cov,pop_load,L(2,2),L(1,1),Gang,meanfitness,meanfitx,meanfity);
        
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

v = meanfitness/steps;
v1 = meanfitx/steps;
v2 = meanfity/steps;

% close output files
fclose(fileID1);
fclose(fileID2);
fclose(fileID3);

end

