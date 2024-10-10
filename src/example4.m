 
% example script to use compare.m

% a direct comparison for input pair TC, pH
% CO2SYS and QUODcarb should get the exact same parameter values
%   with small expected differences in calculated ucnertainties

load fakedata.mat; % variable name is 'data'

% populate opt structure
opt.K1K2 = 10; % option for K1K2 formulation
opt.KSO4 = 1;   % option for KSO4 formulation
opt.KF   = 2;   % option for KF formulation
opt.TB   = 2;   % option for TB formulation
opt.phscale     = 1; % 1 = tot, 2 = sws, 3 = free, 4 = nbs
opt.printcsv    = 1; % print output to CSV? (1=on, 0=off)
opt.fname       = 'QUODcarb.csv'; % not used if opt.printcsv = 0
opt.printmes    = 0; % print screen messages? (1=on, 0=off)
opt.co2press    = 1; % pressure correction for p2f and K0 (1=on, 0=off)
opt.Revelle     = 1; % calculate Revelle factor? (1=on, 0=off)

nD = length(data); % number of datapoints

for i = 1:nD
    % temperature and pressure INdependent
    obs(i).sal  = data(i,1); % salinity (PSU) 1st column of data.mat
    obs(i).esal = data(i,2); % salinity error 2nd, 1 sigma
    obs(i).TC   = data(i,3); % TC (umol/kg) 3rd column of data.mat
    obs(i).eTC  = data(i,4); % TC error is in the 4th column, 1 sigma
 
    obs(i).TP   = data(i,10); % total phosphate (umol/kg)
    obs(i).eTP  = data(i,11); % error TP, 1 sigma
    obs(i).TSi  = data(i,12); % total silicate (umol/kg)
    obs(i).eTSi = data(i,13); % error TSi, 1 sigma

    % temperature and pressure dependent 2nd system, for ph & co3 25 degC
    obs(i).tp(1).T  = 25.0; % temperature for ph in Celsius
    obs(i).tp(1).eT = 0.01; % ph temperature error, 1 sigma
    obs(i).tp(1).P  = 0.0; % pressure (dbar)
    obs(i).tp(1).eP = 1.0; % pressure error, 1 sigma
    obs(i).tp(1).ph     = data(i,7); % ph measured on total scale
    obs(i).tp(1).eph    = 0.004; % ph uncertainty, 1 sigma
end

[est,obs,sys,iflag,opt] = QUODcarb(obs,opt);

% now, compare output to CO2SYS
fid = 'compare_TC_ph.csv';
tp = 1;
A = compare(obs,est,opt,tp,2,fid);
                % note this '2' tells compare.m which pair you are
                % considering, see compare.m file for more options

% the excel file will have values of parameters of interest
% for an exactly-determined calculation, they should be the same values

