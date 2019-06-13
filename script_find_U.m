% use Robbins-Monro Algorithm to determine U achieving v, given s

Narry = [1e9,1e9,1e9,1e7,1e9,1e11];

% use this line if you want to set a target v (comment out if not used)
trgt_rate_arry = [2.4127471216847332e-05,5.3080436677064124e-05,... 
        1.3993933305771449e-04,5.3080436677064124e-05,... 
        5.3080436677064124e-05,5.3080436677064124e-05];

if (rate_flag == 0)
    trgt_rate_arry = trgt_rate_arry*1e+2;
end

sarry = logspace(-3.3,-0.5,15)';
Uarry = ones(size(sarry));

steps = 1e5;
n = 35;

outputfile = 'data/mutBiasCI_estimate_U_ml-7-';     % name of output file

for i=1:6
    tic
    [Uarry,Uthry] = get_U_estimates(Narry(i),trgt_rate_arry(i),sarry,steps,n,rate_flag);
    dlmwrite([outputfile num2str(i) '-0.dat'],[Narry(i) trgt_rate_arry(i)],'delimiter',',','precision',16);
    dlmwrite([outputfile num2str(i) '-1.dat'],[sarry Uarry Uthry],'delimiter',',','precision',16);
    toc
end