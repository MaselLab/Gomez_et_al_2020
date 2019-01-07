cd ~
cd MATLAB/

%% setup array with parameters to simulate
N = 1e9;
s = 1e-2;
U = 1e-5;
data_pts = 50;

v = s^2*(2*log(N*s)-log(s/U))/(log(s/U)^2);     % extend this to a range of v

sarry = (1e-3)*(1e-1/1e-3).^((0:1:data_pts)./data_pts);     % use range specified in Gomez et al 2019
Uarry = ones(size(sarry));

% Analytical constructions of Uarry contain too much numerical error for use
% Uarry = sarry.*exp((-0.5*sarry.^2/v).*(sqrt(8*v*(log(N*sarry)./(sarry.^2)))-1));
% Uarry = (1e-6)*(1e-4/1e-6).^((0:1:data_pts)./data_pts);
% sarry.^2.*(2*log(N*sarry)-log(sarry./Uarry))./(log(sarry./Uarry).^2)

syms t;
for i=1:length(sarry)
    Uarry(i) = real(vpasolve(log(sarry(i)/t)^2+sarry(i)^2/v*log(sarry(i)/t)-2*log(N*sarry(i))*sarry(i)^2/v == 0,t));
end
clear t;

%% Simulations for main figure of new favor

steps = 2e4;
start_time = 1;                     % cpllect data on distribution at start time
end_time = 1;                       % collect data on distribution at end time
collect_data = false;               % flag indicating if data on distribution will be collected
outputfile = 'mutBiasCI_data_time_series_and_stats_ml-'; 
number_of_sims = 

sarry = (1e-3)*(1e-1/1e-3).^((0:1:data_pts)./data_pts);
Uarry = (1e-6)*(1e-4/1e-6).^((0:1:data_pts)./data_pts);

NsU = [Narry' s*x U*x; N*x sarry' U*x; N*x s*x Uarry'];
v = zeros(size(NsU,1),1);
vDF = zeros(size(NsU,1),1);
v1 = zeros(size(NsU,1),1);
v2 = zeros(size(NsU,1),1);
varx = zeros(size(NsU,1),1);
vary = zeros(size(NsU,1),1);
cov = zeros(size(NsU,1),1);

tic
for i=1:data_pts
    for j=i:data_pts
        [v(i),v1(i),v2(i),varx(i),vary(i),cov(i)] = stochastic_simulation2(NsU(i,1),NsU(i,2),NsU(i,3),collect_data,steps,start_time,end_time,outputfile);

    end
end
toc

csvwrite('~/Documents/kgrel2d/data/2dwave_data_time_avg_stats_ml-01-0.dat',NsU);
csvwrite('~/Documents/kgrel2d/data/2dwave_data_time_avg_stats_ml-01-1.dat',[v v1 v2 varx vary cov]);
