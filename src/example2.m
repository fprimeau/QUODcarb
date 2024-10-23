
% QUODcarb example 2
% to show multiple datapoints and tp(1) and tp(2)

load fakedata.mat;
% variable name is 'fakedata'
% sal  in fakedata(1,:);
% TC   in fakedata(2,:);
% TA   in fakedata(3,:);
% ph   in fakedata(4,:); at 25 degC
% pco2 in fakedata(5,:); at 20 degC
% TP   in fakedata(7,:);
% TSi  in fakedata(8,i);

nD = length(fakedata);

% populate opt structure
opt.K1K2 = 4;       % option for K1K2 formulation
opt.KSO4 = 1;       % option for KSO4 formulation
opt.KF   = 2;       % option for KF formulation
opt.TB   = 2;       % option for TB formulation
opt.phscale     = 1;   % 1 = tot, 2 = sws, 3 = free, 4 = nbs
opt.printcsv    = 1;   % print output to CSV? 1 = on, 0 = off
opt.fname       = 'example2.csv'; % file name of desired CSV file
opt.printmes    = 1;   % print screen messages? 1 = on, 0 = off
opt.co2press    = 0;   % pressure correction for P2F (aka FugFac) and K0

for i = 1:nD
    % temperature and pressure INdependent
    obs(i).sal  = fakedata(1,i); % salinity
    obs(i).usal = 0.001; % salinity uncertainty, 1 sigma
    obs(i).TC   = fakedata(2,i);  % TC (umol/kg) 
    obs(i).uTC  = 2.00; % TC undertainty, 1 sigma
    obs(i).TA   = fakedata(3,i); % TA (umol/kg)
    obs(i).uTA  = 2.00; % TA uncertainty, 1 sigma
    
    obs(i).TP   = fakedata(7,i); % total phosphate (umol/kg)
    obs(i).uTP  = 0.01*fakedata(7,i); % 1% uncertainty TP, 1 sigma
    obs(i).TSi  = fakedata(8,i); % total silicate (umol/kg)
    obs(i).uTSi = 0.01*fakedata(8,i); % 1% uncertainty TSi, 1 sigma

    % temperature and pressure system
    obs(i).tp(1).T  = 25.0; % temperature for ph in Celsius
    obs(i).tp(1).uT = 0.1; % ph temperature uncertainty, 1 sigma
    obs(i).tp(1).P  = 0.0; % ph pressure (dbar)
    obs(i).tp(1).uP = 0.1; % ph pressure uncertainty, 1 sigma
    obs(i).tp(1).ph     = fakedata(4,i); % ph measured on total scale
    obs(i).tp(1).uph    = 0.01; % ph uncertainty, 1 sigma

    % temperature and pressure DEpendend 2nd system
    obs(i).tp(2).T      = 20.0; % temperature for pCO2 in Celsius
    obs(i).tp(2).uT     = 0.1; % pCO2 temperature uncertainty, 1 sigma
    obs(i).tp(2).P      = 0.0; % pressure in dbar for pCO2
    obs(i).tp(2).uP     = 0.1; % dbar, pressure uncertainty, 1 sigma
    obs(i).tp(2).pco2   = fakedata(5,i); % partial pressure CO2 (pCO2) (uatm)
    obs(i).tp(2).upco2  = 0.01*fakedata(5,i); % 1% pCO2 uncertainty (uatm), 1 sigma
end

[est,obs,sys,iflag,opt] = QUODcarb(obs,opt); 

% save example2.mat est; % save output est as mat file

% output est(1) values should be:
% est(1).TC         = 2196.0;       est(1).uTC          = 1.6555;
% est(1).TA         = 2306.8;       est(1).uTA          = 1.7115;
% est(1).tp(1).ph   = 7.6373;       est(1).tp(1).uph    = 0.0061;
% est(1).tp(1).pco2 = 1184.2;       est(1).tp(1).upco2  = 21.491;
% est(1).tp(1).pK1  = 5.8471;       est(1).tp(1).upK1   = 0.0055;
% est(1).tp(1).pK2  = 8.9666;       est(1).tp(1).upK2   = 0.0090;
% est(1).tp(2).ph   = 7.7147;       est(1).tp(2).uph    = 0.0056;
% est(1).tp(2).pco2 = 942.75;       est(1).tp(2).upco2  = 8.6015;
% est(1).tp(2).pK1  = 5.8827;       est(1).tp(2).upK1   = 0.0052;
% est(1).tp(2).pK2  = 9.0635;       est(1).tp(2).upK2   = 0.0087;







