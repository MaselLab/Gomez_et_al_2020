function [Uarry,Uthry,rate] = get_U_estimates(N,trgt_rate,sarry,steps,n)
% main function for estimating U's to achieve target v
% inputs:
% N - pop size
% trgt_rate - target rate of evolution (R or v depending on rate_flag)
% sarry - array with selection coefficients (column vector)
% steps - number of generations to use in simulation
% n - number of iterations to use in solving for U
% rate_flag - flag to search for R (0) or v (1) isoclines
% 
% outputs:
% Uarry - corresponding U values (column vector)

digits(16)
rng(7);     % set seed for random number generator

Uarry = ones(size(sarry));
Uthry = ones(size(sarry));
rate = ones(size(sarry));
data_pts = length(sarry);

for i=1:data_pts
    Uthry(i) = initial_U(N,trgt_rate,sarry(i));
    
    % need to bound U, otherwise we don't get a probability of mutation
    if (i>=9) % origin-fixation estimates need to run longer to converge
        [Uarry(i),rate(i)] = approximate_U(N,trgt_rate,sarry(i),Uthry(i),16*steps,n); % more steps
    else
        [Uarry(i),rate(i)] = approximate_U(N,trgt_rate,sarry(i),Uthry(i),steps,n);
    end
end

%%    
    function [Ut,est_rate] = approximate_U(Ni,trgt_rate,si,U0,steps,n)
    % Script attempts to find Ut such that trgt_rate = v given using the
    % N, s, and a target rate of evolution (trgt_rate), and starting at U0
    % robbins-monro algorithm.
    % 
    % Inputs:
    % N - population size
    % s - selection coefficient
    % trgt_rate - target rate of evolution/adaptation
    % U0 - starting mutation rate
    % n - number of iterations to take
    % 
    % Outputs:
    % Ut - estimate of U that provides target rate of adaptation vt
    % est_rate - return the last rate v associated with Ut 

        Ut = U0;
        log10Ut = log10(U0);
        init_flag = false;
        
        pop = Ni;
        fitx = 0;

        for k=1:n
            %set initial U value and start simulation with no initial pop
            Ut = 10^log10Ut;
            [vn,pop,fitx] = stochastic_simulation_one_trait(Ni,si,Ut,Ut,steps,init_flag,pop,fitx);
            
            % adjust power of U based on error. Here I don't allow it to
            % exceed mutation rates to go higher that 10^-3.1 because then
            % the mutation rate would not be a probability (2*10^-3.1>1).
            log10Ut = min([log10(0.49) log10Ut+(1/k^0.6)*log10(trgt_rate/vn)]);
            init_flag = true;
        end
        
        est_rate = vn;
    end
%%
    function Ui = initial_U(Ni,trgt_rate_v,si)
    % function to get obtain initial estimate to use in the algorithm   
      
        % find initial U for fixed v
        U1 = trgt_rate_v/(2*Ni*si*(sella_hirsh_pfix(Ni,si)));  % origin-fixation v
        U2 = si*exp(-si^2/(2*trgt_rate_v)*(sqrt(1+8*trgt_rate_v*log(Ni*si)/si^2)-1));       % Desai and Fisher v
        U3 = 2*trgt_rate_v^(1.5)/(si^2*sqrt(log(Ni*sqrt(trgt_rate_v)/sqrt(log(Ni*sqrt(trgt_rate_v)*6^(1/6))))));    % Hallatschek v
        Ui = max(U1,U2);
        
        if (Ui>si/10) % check is s is in diffusive mutations regime
            Ui = U3;
        end
        
        Ui = min(Ui,0.49);
    end
%%
    function pfix = sella_hirsh_pfix(Ni,si)
        pfix = ( 1-(1/(1+si))^2 )/( 1-(1/(1+si))^Ni );
    end
end
