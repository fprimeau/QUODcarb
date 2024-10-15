
% QUODcarb example 5
% to show example syntax of turning off pK1 

load fakedata.mat; 
% variable name is 'fakedata'
% sal  in fakedata(1,:);
% TC   in fakedata(2,:);
% TA   in fakedata(3,:);
% ph   in fakedata(4,:); at 25 degC
% pco2 in fakedata(5,:); at 20 degC
% co3  in fakedata(6,:); at 25 degC
% TP   in fakedata(7,:);
% TSi  in fakedata(8,i);
% T insitu in fakedata(9,:);
% P insitu in fakedata(10,:);

% populate opt structure
opt.K1K2 = 10;  % option for K1K2 formulation
opt.KSO4 = 1;   % option for KSO4 formulation
opt.KF   = 2;   % option for KF formulation
opt.TB   = 2;   % option for TB formulation
opt.phscale     = 1; % 1 = tot, 2 = sws, 3 = free, 4 = nbs
opt.printcsv    = 1; % print output to CSV? (1=on, 0=off)
opt.fname       = 'QUODcarb.csv'; % not used if opt.printcsv = 0
opt.printmes    = 0; % print screen messages? (1=on, 0=off)
opt.co2press    = 1; % pressure correction for p2f and K0 (1=on, 0=off)
opt.Revelle     = 0; % calculate Revelle factor? (1=on, 0=off)

opt.turnoff.TB  = 0; % Use TB formulation? (0=yes, 1=no)
opt.turnoff.pK1 = 0; % Use pK1 formulation? (0=yes, 1=no)
opt.turnoff.pK2 = 0; % Use pK2 formulation? (0=yes, 1=no)
opt.pKalpha     = 0; % Include organic alkalinity alpha? (1=on, 0=off)
opt.pKbeta      = 0; % Include organic alkalinity beta? (1=on, 0=off)

nD = length(fakedata); % number of datapoints

for i = 1:nD
    % temperature and pressure INdependent
    obs(i).sal  = fakedata(1,i); % salinity
    obs(i).usal = 0.001; % salinity uncertainty, 1 sigma
    obs(i).TC   = fakedata(2,i);  % TC (umol/kg) 
    obs(i).uTC  = 2.00; % TC undertainty, 1 sigma
    obs(i).TA   = fakedata(3,i); % TA (umol/kg)
    obs(i).uTA  = 2.00; % TA uncertainty 1 sigma
    
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


% with opt.turnoff.pK1 = 1, output est(1) values should be:
% est(1).TC = 2200.6;
% est(1).TA = 2304.4;
% est(1).tp(1).T = 25;
% est(1).tp(1).P = 0;
% est(1).tp(1).ph = 7.6617;
% est(1).tp(1).pco2 = 1714.4;
% est(1).tp(1).pK1 = 6.0354;
% est(1).tp(1).pK2 = 8.9640;

% with opt.turnoff.pK2 = 1, output est(1) values should be:
% est(1).TC = 2200.6;
% est(1).TA = 2304.4;
% est(1).tp(1).T = 25;
% est(1).tp(1).P = 0;
% est(1).tp(1).ph = 7.6617;
% est(1).tp(1).pco2 = 1127.6;
% est(1).tp(1).pK1 = 5.8465;
% est(1).tp(1).pK2 = 9.0480;

% note, you cannot turn off pK1 and pK2 for the same calculation

% for reference, with opt.turnoff.pK1 and opt.turnoff.pK2 both = 0,
% output est(1) values should be:
% est(1).TC = 2198.1;
% est(1).TA = 2306.9;
% est(1).tp(1).T = 24.9918;
% est(1).tp(1).P = -4.1633e-5; 
% est(1).tp(1).ph = 7.6377;
% est(1).tp(1).pco2 = 1187.9;
% est(1).tp(1).pK1 = 5.8480;
% est(1).tp(1).pK2 = 8.9773;
