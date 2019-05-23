function [Uarry,Uthry] = get_U_estimates(N,v,sarry,steps,n)
% main function for estimating U's to achieve target v
% inputs:
% N - pop size
% v - target rate of adaptation
% sarry - array with selection coefficients (column vector)
% steps - number of generations to use in simulation
% n - number of iterations to use in solving for U
% sim_num - string corresponding to simulation number
% 
% outputs:
% Uarry - corresponding U values (column vector)

digits(16)
rng(7);     % set seed for random number generator

Uarry = ones(size(sarry));
Uthry = ones(size(sarry));
data_pts = length(sarry);

for i=1:data_pts
    U0 = initial_U(N,v,sarry(i));
    Uthry(i) = U0;
    Uarry(i) = approximate_U(N,v,sarry(i),U0,steps,n);
end

%%    
    function Ut = approximate_U(N,vt,si,U0,steps,n)
    % Script attempts to find Ut such that vt = v(N,s,Ut), given 
    % N, s, and a target rate of adaptation vt, and starting at U0
    % algorithm based on robbins-monro approach
    % 
    % Inputs:
    % N - population size
    % s - selection coefficient
    % vt - target rate of adaptation
    % U0 - starting mutation rate
    % n - number of iterations to take
    % 
    % Outputs:
    % Ut - estimate of U that provides target rate of adaptation vt

        Ut = U0;
        log10Ut = log10(U0);

        for k=1:n
            Ut = 10^log10Ut;
            vn = stochastic_simulation_one_trait(N,si,Ut,Ut,steps);
            log10Ut = min(log10(0.5),log10Ut + (1/k)*log10(vt/vn));
        end
    
    end
%%
    function Ui = initial_U(Ni,vi,si)
    % function to get obtain initial estimate    
      
        U1 = vi/(2*Ni*si*((1-exp(-si))/(1-exp(-Ni*si))));
        U2 = si*exp(-si^2/(2*vi)*(sqrt(1+8*vi*log(Ni*si)/si^2)-1));

        Ui = max(U1,U2);
    end
end
