cd ~/Documents/mutBiasCI/

%% setup array for parameters N,s1,s2,u1,u2
N = 1e9;
s = 1e-2;
U = 1e-5;
data_pts = 50;
rng(7);                                                     % set seed for random number generator
v = s^2*(2*log(N*s)-log(s/U))/(log(s/U)^2);                 % extend this to a range of v

sarry = (1e-3)*(1e-1/1e-3).^((0:1:data_pts)./data_pts);     % use range specified in Gomez et al 2019
Uarry = ones(size(sarry));

% analytical Uarry contains numerical error
% Uarry = sarry.*exp((-0.5*sarry.^2/v).*(sqrt(8*v*(log(N*sarry)./(sarry.^2)))-1));
% sarry.^2.*(2*log(N*sarry)-log(sarry./Uarry))./(log(sarry./Uarry).^2)

% numerically solve for U given s and fixed value of v(N,s,U)
syms t;
for i=1:length(sarry)
    Uarry(i) = real(vpasolve(log(sarry(i)/t)^2+sarry(i)^2/v*log(sarry(i)/t)-2*log(N*sarry(i))*sarry(i)^2/v == 0,t));
end
clear t;

%% test simulations for initial figure

steps = 1e3;
start_time = 5e2;                     % collect data on distribution at start time
end_time = 1e3;                       % collect data on distribution at end time
outputfile = '~/Documents/mutBiasCI/data/mutBiasCI_data_for_2d_distribution_ml-01'; 
number_of_sims = 0.5*data_pts*(data_pts+1);
collect_distribution_data = [1;1;zeros(number_of_sims-2,1)];
% collect_distribution_data = rand(number_of_sims,1)>.95;
indx_of_collected_data = [];

NsU = zeros(number_of_sims,5);          % array that stores the parameters [N,s1,u1,s2,u2]
sim_data = zeros(number_of_sims,6);     % data collected [v,v1,v2,varx,vary,cov]

tic
for i=1:data_pts
    for j=i:data_pts
        indx = round(j+(i-1)*data_pts-0.5*i*(i-1));
        NsU(indx,:)=[N sarry(i) Uarry(i) sarry(j) Uarry(j)];
        sim_data(indx,:) = stochastic_simulation_two_traits(N,sarry(i),Uarry(i),sarry(j),Uarry(j),steps,collect_distribution_data(indx),start_time,end_time,[outputfile '-' num2str(indx)]);
        if(collect_distribution_data(indx))
            indx_of_collected_data = [indx_of_collected_data; indx];
        end
        indx
    end
end
toc

csvwrite('~/Documents/mutBiasCI/data/mutBiasCI_data_all_simulation_parameters_ml-01-0.dat',NsU);
csvwrite('~/Documents/mutBiasCI/data/mutBiasCI_data_all_simulation_grand_means_ml-01-1.dat',sim_data);
csvwrite('~/Documents/mutBiasCI/data/mutBiasCI_indx_of_collected_data_ml-01-2.dat',indx_of_collected_data);
