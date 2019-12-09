% use Robbins-Monro Algorithm to determine U achieving v, given s

Narry = [1e9,1e9,1e9,1e7,1e9,1e11];

% use this line if you want to set a target v (comment out if not used)
% For the paper, you need these target rates.
trgt_rate_arry = [5.3080436677064124e-06,5.3080436677064124e-05,... 
        5.3080436677064124e-04,5.3080436677064124e-05,... 
        5.3080436677064124e-05,5.3080436677064124e-05];

sarry = logspace(-3.3,-0.5,15)';
Uarry = ones(size(sarry));

steps = 5e4;
n = 30;

outputfile = 'data/SAapprox/mutBiasCI_estimate_U_ml-101-';     % name of output file

for i=1:6
    tic
    i
    [Uarry,Uthry,est_rate] = get_U_estimates(Narry(i),trgt_rate_arry(i),sarry,steps,n,i);
    dlmwrite([outputfile num2str(i) '-0.dat'],[Narry(i) trgt_rate_arry(i)],'delimiter',',','precision',16);
    dlmwrite([outputfile num2str(i) '-1.dat'],[sarry Uarry Uthry est_rate],'delimiter',',','precision',16);
    toc
end