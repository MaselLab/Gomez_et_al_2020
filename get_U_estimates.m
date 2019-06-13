function [Uarry,Uthry] = get_U_estimates(N,trgt_rate,sarry,steps,n,rate_flag)
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
data_pts = length(sarry);

for i=1:data_pts
    U0 = initial_U(N,trgt_rate,sarry(i),rate_flag);
    Uthry(i) = U0;
    Uarry(i) = approximate_U(N,trgt_rate,sarry(i),U0,steps,n,rate_flag);
end

%%    
    function Ut = approximate_U(N,trgt_rate,si,U0,steps,n,rate_flag)
    % Script attempts to find Ut such that trgt_rate = R or v, given 
    % N, s, and a target rate of evolution (trgt_rate), and starting at U0
    % algorithm based on robbins-monro approach
    % 
    % Inputs:
    % N - population size
    % s - selection coefficient
    % trgt_rate - target rate of evolution/adaptation
    % U0 - starting mutation rate
    % n - number of iterations to take
    % rate_flag - flag to search for R (0) or v (1) isoclines
    % 
    % Outputs:
    % Ut - estimate of U that provides target rate of adaptation vt

        Ut = U0;
        log10Ut = log10(U0);

        for k=1:n
            Ut = 10^log10Ut;
            vn = stochastic_simulation_one_trait(N,si,Ut,Ut,steps);
            
            if (rate_flag == 0)     % find target U for fixed R
                log10Ut = min(log10(0.5),log10Ut + (1/k)*log10(trgt_rate*si/vn));
            else                    % find target U for fixed v
                log10Ut = min(log10(0.5),log10Ut + (1/k)*log10(trgt_rate/vn));
            end
        end
    
    end
%%
    function Ui = initial_U(Ni,trgt_rate_i,si,rate_flag)
    % function to get obtain initial estimate    
      
        if (rate_flag == 0)     % find initial U for fixed R 
            U1 = trgt_rate_i/(2*Ni*((1-exp(-si))/(1-exp(-Ni*si))));
            U2 = si*exp(-si/(2*trgt_rate_i)*(sqrt(1+8*trgt_rate_i*log(Ni*si)/si)-1));
            Ui = max(U1,U2);
        else                    % find initial U for fixed v
            U1 = trgt_rate_i/(2*Ni*si*((1-exp(-si))/(1-exp(-Ni*si))));
            U2 = si*exp(-si^2/(2*trgt_rate_i)*(sqrt(1+8*trgt_rate_i*log(Ni*si)/si^2)-1));
            Ui = max(U1,U2);
        end
    end
end
