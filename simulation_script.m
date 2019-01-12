%% setup array for parameters N,s1,s2,u1,u2
N = 1e9;
s = 1e-2;
U = 1e-5;
v = s^2*(2*log(N*s)-log(s/U))/(log(s/U)^2);                 % extend this to a range of v
digits(16)

data_pts = 12;
rng(7);                                                     % set seed for random number generator

sarry = (1e-3)*(2e-2/1e-3).^((0:1:data_pts)./data_pts);     % range for possible s & U values
Uarry = ones(size(sarry));

for i=1:length(sarry)
    si = sarry(i);
    Ui = exp( (0.5*si^2/v) * ( 1 + 2*v*log(si)/si^2 - sqrt(1 + 8*v*log(N*si)/si^2) ) );
    Uarry(i) = Ui; 
    varry(i) = si.^2.*(2*log(N*si)-log(si./Ui))./(log(si./Ui).^2);      % checking that Ui is correct solution
    qarry(i) = 2*log(N*si)./log(si./Ui);
end

% Examine whether s,U combinations give the traveling wave regime
% [(1:51)' log10(sarry') log10(Uarry'./sarry') log10(N*Uarry') log10(N*Uarry'.*log(N*sarry')) qarry'] 

steps = 1e6;
start_time = 5e4;                     % collect data on distribution at start time
end_time = 6.5e4;                       % collect data on distribution at end time
outputfile = '~/Documents/mutBiasCI/data/mutBiasCI_data_for_2d_distribution_ml-01'; 
number_of_sims = 0.5*data_pts*(data_pts+1);

collect_distribution_data = zeros(number_of_sims,1);
collect_distribution_data([6 12 30 57])=1;
indx_of_collected_data = [];

NsU = zeros(number_of_sims,5);          % array that stores the parameters [N,s1,u1,s2,u2]
sim_data = zeros(number_of_sims,6);     % data collected [v,v1,v2,varx,vary,cov]

tic
for i=1:data_pts
    for j=i:data_pts
        indx = round(j+(i-1)*data_pts-0.5*i*(i-1));
        NsU(indx,:)=[N sarry(i) Uarry(i) sarry(j) Uarry(j)];
        [sim_data(indx,1),sim_data(indx,2),sim_data(indx,3),sim_data(indx,4),sim_data(indx,5),sim_data(indx,6)] ...
            = stochastic_simulation_two_traits(N,sarry(i),Uarry(i),sarry(j),Uarry(j),steps,collect_distribution_data(indx),start_time,end_time,[outputfile '-' num2str(indx)]);
        if(collect_distribution_data(indx))
            indx_of_collected_data = [indx_of_collected_data; indx];
        end
    end
end
toc

dlmwrite('~/Documents/mutBiasCI/data/mutBiasCI_data_all_simulation_parameters_ml-01-0.dat',NsU,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/mutBiasCI_data_all_simulation_grand_means_ml-01-1.dat',sim_data,'delimiter',',','precision',16);
dlmwrite('~/Documents/mutBiasCI/data/mutBiasCI_data_all_simulation_indx_of_collected_data_ml-01-2.dat',indx_of_collected_data,'delimiter',',','precision',16);
