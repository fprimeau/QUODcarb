
% QUODcarb example 3
% this example will not work alone because it requires a data.mat file 

% measured 5 parameters n times saved in data.mat

load data.mat;
data = input_data; % save matrix to one variable

% populate opt structure
opt.K1K2 = 10; % option for K1K2 formulation
opt.KSO4 = 1;   % option for KSO4 formulation
opt.KF   = 2;   % option for KF formulation
opt.TB   = 2;   % option for TB formulation
opt.phscale     = 1; % 1 = tot, 2 = sws, 3 = free, 4 = nbs
opt.printcsv    = 0; % print output to CSV? (1=on, 0=off)
opt.fname       = 'QUODcarb.csv'; % not used if opt.printcsv = 0
opt.printmes    = 0; % print screen messages? (1=on, 0=off)
opt.co2press    = 1; % pressure correction for p2f and K0 (1=on, 0=off)
opt.Revelle     = 1; % calculate Revelle factor? (1=on, 0=off)

nD = length(data); % number of datapoints

for i = 1:nD
    % temperature and pressure INdependent
    obs(i).sal  = data(i,1); % salinity (PSU) 5th column
    obs(i).esal = data(i,1); % salinity error 6th column, 1 sigma
    obs(i).TC   = data(i,3); % TC (umol/kg) is the first column
    obs(i).eTC  = data(i,4); % TC error is in the 2nd column, 1 sigma
    obs(i).TA   = data(i,5); % TA (umol/kg) is in the 3rd column
    obs(i).eTA  = data(i,6); % TA error is in the 4th column, 1 sigma
    obs(i).TP   = data(i,7); % total phosphate (umol/kg)
    obs(i).eTP  = data(i,8); % error TP, 1 sigma
    obs(i).TSi  = data(i,9); % total silicate (umol/kg)
    obs(i).eTSi = data(i,10); % error TSi, 1 sigma

    % temperature and pressure dependent 1st system, in-situ
    obs(i).tp(1).T  = data(i,11); % in-situ temperature (Celsius)
    obs(i).tp(1).eT = data(i,12); % in-situ temperature error, 1 sigma
    obs(i).tp(1).P  = data(i,13); % in-situ pressure (dbar)
    obs(i).tp(1).eP = data(i,14); % in-situ pressure error (dbar), 1 sigma

    % temperature and pressure dependent 2nd system, for ph & co3 25 degC
    obs(i).tp(2).T  = 25.0 ; % temperature for ph in Celsius
    obs(i).tp(2).eT = data(i,15); % ph temperature error, 1 sigma
    obs(i).tp(2).P  = data(i,16); % pressure (dbar)
    obs(i).tp(2).eP = data(i,17); % pressure error, 1 sigma
    obs(i).tp(2).ph     = data(i,18); % ph measured on total scale
    obs(i).tp(2).eph    = data(i,19); % ph uncertainty, 1 sigma
    obs(i).tp(2).co3    = data(i,20); % total carbonate ion (umol/kg)
    obs(i).tp(2).eco3   = data(i,21); % CO3 error (umol/kg), 1 sigma

    % temperature and pressure dependent 3rd system, for pco2 20 degC
    obs(i).tp(3).T  = 20.0; % temperature for pCO2 in Celsius
    obs(i).tp(3).eT = 0.002; % pCO2 temperature error, 1 sigma
    obs(i).tp(3).P  = 0.0; % pressure in dbar for pCO2
    obs(i).tp(3).eP = 0.003; % dbar, pressure for pCO2 error, 1 sigma
    obs(i).tp(3).pco2   = data(i,22); % partial pressure CO2 (pCO2) (uatm)
    obs(i).tp(3).epco2  = data(i,23); % pCO2 error (uatm), 1 sigma

    % temperature and pressure dependend 4th system
    obs(i).tp(4).T  = 10.0; % want an output at 10 deg C
    obs(i).tp(4).eT = 0.01; 
    obs(i).tp(4).P  = 10.0; % want output at 10 dbar
    obs(i).tp(4).eP = 0.003;

end

[est,obs,sys,iflag] = QUODcarb(obs,opt);

save example3.mat est; % save output estimate as a mat file



