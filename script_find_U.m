% use Robbins-Monro Algorithm to determine U achieving v, given s

Narry = [1e9,1e9,1e9,1e7,1e9,1e11];
varry = [2.4127471216847332e-05,5.3080436677064124e-05,... 
        1.3993933305771449e-04,5.3080436677064124e-05,... 
        5.3080436677064124e-05,5.3080436677064124e-05];
    
sarry = logspace(-3.3,-0.5,15)';
Uarry = ones(size(sarry));

steps = 1e5;
n = 50;

outputfile = 'data/mutBiasCI_estimate_U_ml-6-'; 

for i=1:6
    tic
    [Uarry,Uthry] = get_U_estimates(Narry(i),varry(i),sarry,steps,n);
    dlmwrite([outputfile num2str(i) '-0.dat'],[Narry(i) varry(i)],'delimiter',',','precision',16);
    dlmwrite([outputfile num2str(i) '-1.dat'],[sarry Uarry Uthry],'delimiter',',','precision',16);
    toc
end