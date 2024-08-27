
% QUODcarb example 2

load fakedata.mat; % variable name is 'data'
n = length(data);

% populate opt structure
opt.K1K2 = 4;       % option for K1K2 formulation
opt.KSO4 = 1;       % option for KSO4 formulation
opt.KF   = 2;       % option for KF formulation
opt.TB   = 2;       % option for TB formulation
opt.phscale = 1;    % 1 = tot, 2 = sws, 3 = free, 4 = nbs
opt.printcsv = 1;   % print output to CSV? 1 = on, 0 = off
opt.fname = 'example2.csv'; % file name of desired CSV file
opt.printmes = 1;   % print screen messages? 1 = on, 0 = off
opt.co2press = 0;   % pressure correction for P2F (aka FugFac) and K0

for i = 1:n
    % temperature and pressure INdependent
    obs(i).sal  = data(i,1); % salinity (PSU)
    obs(i).esal = data(i,2); % salinity error, 1 sigma
    obs(i).TC   = data(i,3); % TC (umol/kg)
    obs(i).eTC  = data(i,4); % TC error, 1 sigma
    obs(i).TA   = data(i,5); % TA (umol/kg)
    obs(i).eTA  = data(i,6); % TA error, 1 sigma

    % temperature and pressure DEpendent 1st system
    obs(i).tp(1).T      = 25.0;     % temperature for ph in Celsius
    obs(i).tp(1).eT     = 0.01;     % temperature error for ph bath, 1 sigma
    obs(i).tp(1).P      = 0.0;      % pressure (dbar)
    obs(i).tp(1).eP     = 1.0;      % pressure error, 1 sigma
    obs(i).tp(1).ph     = data(i,7); % ph measured on total scale
    obs(i).tp(1).eph    = 0.004;     % ph uncertainty, 1 sigma

    % temperature and pressure DEpendend 2nd system
    obs(i).tp(2).T      = 20.0;     % temperature for pCO2 in Celsius
    obs(i).tp(2).eT     = 0.01;     % pCO2 temperature error, 1 sigma
    obs(i).tp(2).P      = 0.0;      % pressure in dbar for pCO2
    obs(i).tp(2).eP     = 1.0;      % dbar, pressure error, 1 sigma
    obs(i).tp(2).pco2   = data(i,8); % partial pressure CO2 (pCO2) (uatm)
    obs(i).tp(2).epco2  = data(i,9); % pCO2 error (uatm), 1 sigma
end

[est,obs,sys,iflag,opt] = QUODcarb(obs,opt); 

save example2.mat est; % save output est as mat file

% output est(1) values should be:
% est(1).TC         = 2238.6;       est(1).eTC          = 3.6729;
% est(1).TA         = 2427.3;       est(1).eTA          = 3.9231;
% est(1).tp(1).ph   = 7.9335;       est(1).tp(1).eph    = 0.0038;
% est(1).tp(1).pco2 = 619.4102;     est(1).tp(1).epco2  = 10.0685;
% est(1).tp(1).pK1  = 5.8638;       est(1).tp(1).epK1   = 0.0055;
% est(1).tp(1).pK2  = 9.1288;       est(1).tp(1).epK2   = 0.0085;
% est(1).tp(2).ph   = 7.9194;       est(1).tp(2).eph    = 0.0087;
% est(1).tp(2).pco2 = 608.1174;     est(1).tp(2).epco2  = 14.2019;
% est(1).tp(2).pK1  = 5.9024;       est(1).tp(2).epK1   = 0.0055;
% est(1).tp(2).pK2  = 9.0719;       est(1).tp(2).epK2   = 0.0100;







