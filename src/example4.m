 
% QUODcarb example 4
% example script to use compare.m

% a direct comparison for input pair TC, pH
% CO2SYS and QUODcarb should get the exact same parameter values
%   with small expected differences in calculated uncertainties

load fakedata.mat;
% variable name is 'fakedata'
% sal  in fakedata(1,:);
% TC   in fakedata(2,:);
% ph   in fakedata(4,:); at 25 degC
% TP   in fakedata(7,:);
% TSi  in fakedata(8,i);

% populate opt structure
opt.K1K2 = 10; % option for K1K2 formulation
opt.KSO4 = 1;   % option for KSO4 formulation
opt.KF   = 2;   % option for KF formulation
opt.TB   = 2;   % option for TB formulation
opt.phscale     = 1; % 1 = tot, 2 = sws, 3 = free, 4 = nbs
opt.printcsv    = 1; % print output to CSV? (1=on, 0=off)
opt.fname       = 'example4.csv'; % not used if opt.printcsv = 0
opt.printmes    = 0; % print screen messages? (1=on, 0=off)
opt.co2press    = 1; % pressure correction for p2f and K0 (1=on, 0=off)
opt.Revelle     = 1; % calculate Revelle factor? (1=on, 0=off)

nD = length(fakedata); % number of datapoints

for i = 1:nD
    % temperature and pressure INdependent
    obs(i).sal  = fakedata(1,i); % salinity
    obs(i).usal = 0.001; % salinity uncertainty, 1 sigma
    obs(i).TC   = fakedata(2,i);  % TC (umol/kg) 
    obs(i).uTC  = 2.00; % TC undertainty, 1 sigma
    
    obs(i).TP   = fakedata(7,i); % total phosphate (umol/kg)
    obs(i).uTP  = 0.01*fakedata(7,i);   % 1% uncertainty TP, 1 sigma
    obs(i).TSi  = fakedata(8,i); % total silicate (umol/kg)
    obs(i).uTSi = 0.01*fakedata(8,i);  % 1% uncertainty TSi, 1 sigma

    % temperature and pressure system
    obs(i).tp(1).T  = 25.0; % temperature for ph in Celsius
    obs(i).tp(1).uT = 0.05; % ph temperature uncertainty, 1 sigma
    obs(i).tp(1).P  = 0.0; % ph pressure (dbar)
    obs(i).tp(1).uP = 0.07; % ph pressure uncertainty, 1 sigma
    obs(i).tp(1).ph     = fakedata(4,i); % ph measured on total scale
    obs(i).tp(1).uph    = 0.01; % ph uncertainty, 1 sigma
end

[est,obs,sys,iflag,opt] = QUODcarb(obs,opt);

% now, compare output to CO2SYS
fid = 'compare_TC_ph.csv';
tp = 1;
out = compare(obs,est,opt,tp,2,3,fid);
                % note this '2' tells compare.m that parameter 1 is TC
                % the '3' tells compare.m that parameter 2 is ph
                % see compare.m file for more options

% the excel file will have values of parameters of interest
% for an exactly-determined calculation, they should be the same 
% parameter values with expected small differences between uncertainties

