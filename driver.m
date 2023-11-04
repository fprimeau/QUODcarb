

% driver to go with QUODcarb

% updated opt (options) structure

% g = GOMECC data

%format long

load datag.mat;
[in] = datag;
nD = 5; % 1186; % until UW start, or before: %length(in);

% choose options for opt structure
opt.K1K2 = 10; % option for K1K2 formulation
opt.KSO4 = 1;  % option for KSO4 formulation
opt.KF   = 2;  % option for KF formulation
opt.TB   = 2;  % option for TB formulation
opt.phscale  = 1;  % 1 = tot, 2 = sws, 3 = free, 4 = NBS
opt.printcsv = 0;  % print est to CSV? 1 = on , 0 = off
% opt.fname    = 'QUODcarb_output.csv'; % don't need it if printcsv is off
% opt.fname    = 'output_csv/Q5_Nov3.csv';
opt.co2press = 0; % 1 = on, 0 = off
opt.Revelle  = 0; % 1 = on, 0 = off 
opt.printmes = 0; % 1 = on, 0 = off

% read in GOMECC data and put into obs structure

for i = 1:nD
    
    % measurements that are independent of (T,P)
    obs(i).TC    = in(i,5); % (umol/kg)
    obs(i).eTC   = 2.01;    % TC error ±2.01 umol/kg
    obs(i).TA    = in(i,7);
    obs(i).eTA   = 1.78; % TA error ±2.01 umol/kg
    obs(i).sal   = in(i,3); % PSU
    obs(i).esal  = 0.002; % new = 0.001
    % nutrients P and Si also independent of (T,P)
    obs(i).TP    = in(i,24);
    obs(i).eTP   = in(i,24)*0.001; % 0.01% meas uncertainty
    obs(i).TSi   = in(i,18);
    obs(i).eTSi  = in(i,18)*0.001; % 0.01% meas uncertainty

    % first (T,P)-dependent measurement
    obs(i).tp(1).T  = in(i,2); % deg C, CTD temp
    obs(i).tp(1).eT = 0.02; % ±0.02 degC
    obs(i).tp(1).P  = in(i,1); % dbar
    obs(i).tp(1).eP = 0.63; % (max) ± 0.63 dbar

    % second(T,P)-dependent measurement
    obs(i).tp(2).T    = 25 ; % degC
    obs(i).tp(2).eT   = 0.05 ; % from cruise report
    obs(i).tp(2).P    = 0.0 ; %in(i+ad,1); % NOT in situ
    obs(i).tp(2).eP   = 0.63 ;
    obs(i).tp(2).ph   = in(i,9); % total scale
    obs(i).tp(2).eph  = 0.0004 ;
    obs(i).tp(2).co3  = in(i,15); % (µmol/kg)
    obs(i).tp(2).eco3 = 2.0 ;  % std ±2µmol/kg

    % third (T,P)-dependent measurement
    obs(i).tp(3).T     = 20 ; %degC
    obs(i).tp(3).eT    = 0.03 ; % from cruise report
    obs(i).tp(3).P     = 0.0 ; % dbar (surface pressure for pco2)
    obs(i).tp(3).eP    = 0.63 ;
    obs(i).tp(3).pco2  = in(i,12); % (µatm)
    obs(i).tp(3).epco2 = in(i,12)*0.0021; % 0.21% relative std error (avg)
end

obs_backup = obs;

%% Q5: All five input
% % CT AT pH pCO2 CO3 (Q5) (fid5)
obs = obs_backup;
[est,obs,sys,iflag] = QUODcarb(obs,opt);
est05 = est;


%% Q2: Input pairs

% % TA TC (Q2) (fid2)
for i = 1:nD
    obs(i).tp(2).ph = nan;   obs(i).tp(2).eph = nan;
    obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan;
    obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan;
end
[est,obs, ~, ~] = QUODcarb(obs,opt); % [est, obs, sys, iflag]
est02  = est;
% fid2   = 'compare_outs/compare_TC_TA.csv'; 
% tp     = 2; % second tp system for ph in there
% A      = compare(obs,est,opt,tp,1,fid2); % 1 for input pair TA TC

% % TA ph (Q2) (fid3)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;         obs(i).eTC = nan;
    obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan;
    obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan;
end
[est,obs, ~, ~] = QUODcarb(obs,opt);
est03   = est;
% tp      = 2;
% fid3    = 'compare_outs/compare_TA_ph.csv';
% [A]     = compare(obs,est,opt,tp,3,fid3);



