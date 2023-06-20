
% drivergV6 to go with QUODcarbV6

% updated opt (options) structure

% g = GOMECC data

format long

load datag.mat;
[in] = datag;
nD = length(in);

% choose options for opt structure
opt.K1K2 = 4; % option for K1K2 formulation
opt.KSO4 = 1; % option for KSO4 formulation
opt.KF   = 2; % option for KF formulation
opt.TB   = 2; % option for TB formulation
opt.phscale  = 1;  % 1 = tot, 2 = sws, 3 = free, 4 = NBS
opt.printcsv = 1; % print est to CSV? 1 = on , 0 = off
opt.fid      = 'CT_AT_K04.csv';
opt.printmes = 1; % print screen messages? 1 = on, 0 = off
opt.abr = 'all'; % option for which acid/base reactions to include
    % opt.abr = {'all'}; % same as below
    % opt.abr = {'borate','sulfate','fluoride','phosphate', ...
    %               'silicate','ammonia','sulfide','solubility'};
opt.co2press = 1; % 1 = on, 0 = off (default off)

ad = 1200; % for debugging, so can start at i = 1
% read in GOMECC data and put into obs structure

for i = 1:10 % need to initialize with i = 1
    
    % measurements that are independent of (T,P)
    obs(i).TC    = in(i+ad,5); % (umol/kg)
    obs(i).eTC   = 2.01;    % TC error ±2.01 umol/kg
    if i >= 1187
        obs(i).eTC = 2.24; % TC error ±2.24 for UW measurements
    end
    obs(i).TA    = in(i+ad,7);
    obs(i).eTA   = 1.78; % TA error ±2.01 umol/kg
    obs(i).sal   = in(i+ad,3); % PSU
    obs(i).esal  = 0.002;
    % nutrients P and Si also independent of (T,P)
    obs(i).TP    = in(i+ad,24);
    obs(i).eTP   = in(i+ad,24)*0.001; % 0.01% meas uncertainty
    obs(i).TSi   = in(i+ad,18);
    obs(i).eTSi  = in(i+ad,18)*0.001; % 0.01% meas uncertainty

    % first (T,P)-dependent measurement
    obs(i).m(1).T = in(i+ad,2); % deg C, CTD temp
    obs(i).m(1).eT = 0.02; % ±0.02 degC
    obs(i).m(1).P = in(i+ad,1); % dbar
    obs(i).m(1).eP = 0.63; % (max) ± 0.63 dbar

    % second(T,P)-dependent measurement
    obs(i).m(2).T    = 25 ; % degC
    obs(i).m(2).eT   = 0.05 ; % from cruise report
    obs(i).m(2).P    = in(i+ad,1); 
    obs(i).m(2).eP   = 0.63 ;
    obs(i).m(2).ph   = in(i+ad,9); % total scale
    obs(i).m(2).eph  = 0.0004 ;
    if i >= 1187
        obs(i).m(2).eph = 0.05; % for UW measurements
    end
    obs(i).m(2).co3  = in(i+ad,15); % (µmol/kg)
    obs(i).m(2).eco3 = 2.0 ;  % std ±2µmol/kg

    % third (T,P)-dependent measurement
    obs(i).m(3).T   = 20 ; %degC
    obs(i).m(3).eT  = 0.03 ; % from cruise report
    obs(i).m(3).P   = 0 ; % dbar (surface pressure for pco2)
    obs(i).m(3).eP  = 0.63 ;
    obs(i).m(3).pco2  = in(i+ad,12); % (µatm)
    obs(i).m(3).epco2 = in(i+ad,12)*0.0021; % 0.21% relative std error (avg)

end

obs_backup = obs;

% % CT AT (Q2) (fid2)
% for i = 1:nD
%     obs(i).m(2).ph = nan;    obs(i).m(2).eph = nan;
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
%     obs(i).m(2).co3 = nan;   obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,opt);
% est02 = est;
% 
% fid2 = 'compare_CT_AT_K04.csv';
% [A] = compare2(obs,est,opt,2,1,fid2);
% 
% % AT pH (Q2) (fid3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;         obs(i).eTC = nan;
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
%     obs(i).m(2).co3 = nan;   obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,opt); % need to change opt.fid for each
% est03 = est;
% 
% fid2 = 'compare_AT_pH_K04.csv';
% [A] = compare2(obs,est,opt,2,3,fid2);

% % TC pH (Q2) (fid4)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan;         obs(i).eTA = nan;
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
%     obs(i).m(2).co3 = nan;   obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,opt); % need to change opt.fid for each
% est04 = est;
% 
% fid2 = 'compare_TC_pH_K04.csv';
% [A] = compare2(obs,est,opt,2,2,fid2);

% % pH pCO2 (Q2) (fid10)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;            obs(i).eTC = nan;
%     obs(i).TA = nan;            obs(i).eTA = nan;
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,opt); % need to change opt.fid for each
% est10 = est;
% 
% fid2 = 'compare_pH_pCO2_K04.csv';
% [A] = compare2(obs,est,opt,2,4,fid2);

% CT AT pH pCO2 CO3 (Q5) (fid5)
obs = obs_backup;
[est,obs,iflag] = QUODcarbV6(obs,opt);
est05 = est;
















