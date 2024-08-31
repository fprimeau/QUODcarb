
% QUODcarb example 3

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
    obs(i).TA   = data(i,5); % TA (umol/kg) is in the 5th column
    obs(i).eTA  = data(i,6); % TA error is in the 6th column, 1 sigma
    
    obs(i).TP   = data(i,10); % total phosphate (umol/kg)
    obs(i).eTP  = data(i,11); % error TP, 1 sigma
    obs(i).TSi  = data(i,12); % total silicate (umol/kg)
    obs(i).eTSi = data(i,13); % error TSi, 1 sigma

    % temperature and pressure dependent 1st system, in-situ
    obs(i).tp(1).T  = data(i,14); % in-situ temperature (Celsius)
    obs(i).tp(1).eT = 0.01; % in-situ temperature error, 1 sigma
    obs(i).tp(1).P  = data(i,15); % in-situ pressure (dbar)
    obs(i).tp(1).eP = 1.0; % in-situ pressure error (dbar), 1 sigma

    % temperature and pressure dependent 2nd system, for ph & co3 25 degC
    obs(i).tp(2).T  = 25.0; % temperature for ph in Celsius
    obs(i).tp(2).eT = 0.01; % ph temperature error, 1 sigma
    obs(i).tp(2).P  = 0.0; % pressure (dbar)
    obs(i).tp(2).eP = 1.0; % pressure error, 1 sigma
    obs(i).tp(2).ph     = data(i,7); % ph measured on total scale
    obs(i).tp(2).eph    = 0.004; % ph uncertainty, 1 sigma
    obs(i).tp(2).co3    = data(i,9); % total carbonate ion (umol/kg)
    obs(i).tp(2).eco3   = 2; % CO3 error (umol/kg), 1 sigma

    % temperature and pressure dependent 3rd system, for pco2 20 degC
    obs(i).tp(3).T  = 20.0; % temperature for pCO2 in Celsius
    obs(i).tp(3).eT = 0.002; % pCO2 temperature error, 1 sigma
    obs(i).tp(3).P  = 0.0; % pressure in dbar for pCO2
    obs(i).tp(3).eP = 0.003; % dbar, pressure for pCO2 error, 1 sigma
    obs(i).tp(3).pco2   = data(i,8); % partial pressure CO2 (pCO2) (uatm)
    obs(i).tp(3).epco2  = 2; % pCO2 error (uatm), 1 sigma

    % temperature and pressure dependend 4th system
    obs(i).tp(4).T  = 10.0; % want an output at 10 deg C
    obs(i).tp(4).eT = 0.01; 
    obs(i).tp(4).P  = 10.0; % want output at 10 dbar
    obs(i).tp(4).eP = 0.003;

end

[est,obs,sys,iflag,opt] = QUODcarb(obs,opt);

save example3.mat est; % save output estimate as a mat file

% output est(1) values should be:
% est(1).TC         = 2246.4;       % posterior TC
% est(1).TA         = 2418.8;       % posterior TA
% est(1).tp(1).T    = 20.3847;      % in-situ temperature
% est(1).tp(1).P    = 0.6964;       % in-situ pressure
% est(1).tp(2).T    = 24.9959;      % temperature at which ph and co3 were measured
% est(1).tp(2).ph   = 7.9240;       % posterior ph
% est(1).tp(2).co3  = 118.5302;     % posterior co3
% est(1).tp(3).T    = 20.0000;      % temperature at which pco2 was measured
% est(1).tp(3).pco2 = 711.5998;     % posterior pco2
% est(1).tp(4).T    = 10;           % output system @10degC
% est(1).tp(4).P    = 10;
% est(1).tp(2).OmegaAr = 1.9113;
% est(1).tp(2).OmegaCa = 2.9162;
% est(1).tp(2).Revelle = 13.6838;








