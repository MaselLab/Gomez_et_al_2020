function v = stochastic_simulation_one_trait(N,s1,u1,ud1,steps)
% The code below has been modified from the source code made
% availabe by Pearce MT and Fisher DS, obtained from:
% 
% https://datadryad.org/resource/doi:10.5061/dryad.f36v6 
%
% Stochastic simulation for two chromosome model. 
% Each step is one generation. Rates (such as s,u) are per generation. 
% The expected size of each subpop after selection, mutation, and mating is computed.
% If expected size < :cutoff: then the size is drawn from Poisson
% distribution with mean equal to the expected size. Otherwise the
% subpopulation size = expected size. 

% output :v     total rate of adaptation, and rate of adaptation in traits 1 and 2. 

% input :N: population size
% input :s1: effect size of beneficial mutation in trait 1. 
% input :u1: mutation rate per locus of trait 1.
% input :ud1: deleterious mutation rate per locus of trait 1.
% input :steps: number of steps for simulation.

% fitness distribution stored as a column vector in pop

digits(16);

% initialize variables
pop=N;                      % abundances of a classes
Na = 0;                     % actual population size
fit=0;                      % total fitness of a class
fitx=0;                     % fitness in trait 1 of a class

meanfitness = 0;            % mean fitness of the population
meanfit_s = 0;              % mean fitness of the population
cutoff=10/s1;       % population cutoff for stochasticity

% Main loop for simulation of each generation
for timestep=1:steps   

    %%%%%%%%%%%%%%%%%%%
    % Remove columns of zeros of decreasing fitness classes
    while any(pop(1,:))==0 
        pop(1,:)=[];
        fit(1,:)=[];
        fitx(1)=[];
    end

    % Add columns for padding, i.e. for new class produced by mutations
    dim=size(pop); 
    if any(pop(dim(1),:))==1    % check for expansion of front in direction of trait 2
        pop(dim(1)+1,:)=0;
        fit(dim(1)+1,:)=fit(dim(1),:)+1;
        fitx(dim(1)+1)=fitx(dim(1))+1;
    end

    %%%%%%%%%%%% 
    % Find expected frequencies after selection and mutation
    dim=size(pop);  
    freq=pop/sum(sum(pop));                                 % array with frequencies of the each class

    % calculate growth due to selection
    fit = s1*fitx';                                % total fitness
    meanfitness = sum(sum(times(freq,fit)));

    newfreq=times(exp(fit-meanfitness),freq);  % after selection
    newfreq=newfreq/sum(sum(newfreq));          % make sure frequencies still add to one.

    % calculate changes in abundances due to mutations

    % beneficial mutations
    z1=zeros(1,dim(2));
    mutatex=[z1; newfreq];
    mutatex(dim(1)+1,:)=[];                     % newfreq already has padding, get rid of extra padding from shift

    % deletirious mutations (note: del mutations in lowest classes ignored)
    z1=zeros(1,dim(2));
    mutatexdel=[newfreq; z1];
    mutatexdel(1,:)=[];                     % newfreq already has padding, get rid of extra padding from shift

    nomutate=(1-u1-ud1)*newfreq;
    postmutate=nomutate+u1*mutatex+ud1*mutatexdel;    
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
    pop = newpop*(N/Na);

    % recompute time-average of variances and covariances

    if timestep==3000       % setting burn in period at 3000 generations
        meanfit_s = meanfitness;
    end

end

v = (meanfitness-meanfit_s)/(steps-3000);

end
