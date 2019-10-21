% use Robbins-Monro Algorithm to determine U achieving v, given s

Narry = [1e9,1e9,1e9,1e7,1e9,1e11];
rate_flag = 1; % 0 - find R(U,s) contour lines, 1 - find v(U,s) countour lines

% use this line if you want to set a target v (comment out if not used)
% For the paper, you need these target rates.
trgt_rate_arry = [5.3080436677064124e-06,5.3080436677064124e-05,... 
        5.3080436677064124e-04,5.3080436677064124e-04,... 
        5.3080436677064124e-04,5.3080436677064124e-04];

if (rate_flag == 0)
    % The target rate is calibrated as v in the origin fixation regime
    % where v=2NUs^2 and R = v2NUs. Since s=1e+2. Multiplying by 1e2 
    % scales v to get R
    trgt_rate_arry = trgt_rate_arry*1e+2; 
end

sarry = logspace(-3.3,-0.5,15)';
Uarry = ones(size(sarry));

steps = 5e4;
n = 50;

outputfile = 'data/mutBiasCI_estimate_U_ml-22-';     % name of output file

for i=1:6
    tic
    [Uarry,Uthry,est_rate] = get_U_estimates(Narry(i),trgt_rate_arry(i),sarry,steps,n,rate_flag);
    dlmwrite([outputfile num2str(i) '-0.dat'],[Narry(i) trgt_rate_arry(i)],'delimiter',',','precision',16);
    dlmwrite([outputfile num2str(i) '-1.dat'],[sarry Uarry Uthry est_rate],'delimiter',',','precision',16);
    toc
end