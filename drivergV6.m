
% drivergV6 to go with QUODcarbV6

% updated opt (options) structure

% g = GOMECC data

format long

load datag.mat;
[in] = datag;
nD = length(in);
% nD = length(obs);

%nD = 12;

% choose options for opt structure
opt.K1K2 = 10; % option for K1K2 formulation
opt.KSO4 = 1; % option for KSO4 formulation
opt.KF   = 2; % option for KF formulation
opt.TB   = 2; % option for TB formulation
opt.phscale  = 1;  % 1 = tot, 2 = sws, 3 = free, 4 = NBS
opt.printcsv = 0; % print est to CSV? 1 = on , 0 (else) = off
opt.fid      = 'CT_AT_K10v6.csv';
opt.printmes = 1; % print screen messages? 1 = on, 0 = off
opt.abr = 'all'; % option for which acid/base reactions to include
    % opt.abr = {'all'}; % same as below
    % opt.abr =
    % {'borate','sulfate','fluoride','phosphate','silicate','ammonia'};
    % sys.abr = {'borate','sulfate','fluoride', ...
    %    'phosphate','silicate','ammonia','sulfide','solubility'};
opt.co2press = 1; % 1 = on, 0 = off (default)

fid2 = fopen('compare_TC_TA_K10v6.csv','w');

% read in GOMECC data and put into obs structure

for i = 1200:1210 %1:nD
    
    % measurements that are independent of (T,P)
    obs.TC    = in(i,5); % (umol/kg)
    obs.eTC   = 2.01;    % TC error ±2.01 umol/kg
    if i >= 1187
        obs.eTC = 2.24; % TC error ±2.24 for UW measurements
    end
    obs.TA    = in(i,7);
    obs.eTA   = 1.78; % TA error ±2.01 umol/kg
    obs.sal   = in(i,3); % PSU
    obs.esal  = 0.002;
    % nutrients P and Si also independent of (T,P)
    obs.TP    = in(i,24);
    obs.eTP   = in(i,24)*0.001; % 0.1% meas uncertainty
    obs.TSi   = in(i,18);
    obs.eTSi  = in(i,18)*0.001; % 0.1% meas uncertainty

    % first (T,P)-dependent measurement
    obs.m(1).T = in(i,2); % deg C, CTD temp
    obs.m(1).eT = 0.02; % ±0.02 degC
    obs.m(1).P = in(i,1); % dbar
    obs.m(1).eP = 0.63; % (max) ± 0.63 dbar

    % second(T,P)-dependent measurement
    obs.m(2).T    = 25 ; % degC
    obs.m(2).eT   = 0.05 ; % from cruise report
    obs.m(2).P    = in(i,1); 
    obs.m(2).eP   = 0.63 ;
    obs.m(2).ph   = in(i,9); % total scale
    obs.m(2).eph  = 0.0004 ;
    if i >= 1187
        obs.m(2).eph = 0.05; % for UW measurements
    end
    obs.m(2).co3  = in(i,15); % (µmol/kg)
    obs.m(2).eco3 = 2.0 ;  % std ±2µmol/kg

    % third (T,P)-dependent measurement
    obs.m(3).T   = 20 ; %degC
    obs.m(3).eT  = 0.03 ; % from cruise report
    obs.m(3).P   = 0 ; % dbar (surface pressure for pco2)
    obs.m(3).eP  = 0.63 ;
    obs.m(3).pco2  = in(i,12); % (µatm)
    obs.m(3).epco2 = in(i,12)*0.0021; % 0.21% relative std error (avg)

    obs_backup = obs;

    % CT AT (Q2) (fid2)
    obs.m(2).ph = nan;    obs.m(2).eph = nan;
    obs.m(3).pco2 = nan;  obs.m(3).epco2 = nan;
    obs.m(2).co3 = nan;   obs.m(2).eco3 = nan;
    [est,obs,iflag] = QUODcarbV6(obs,opt);
    est02(i) = est;

    if (length(est02) == 1 )
        compare2(fid2); % make columns on first it
    end
    compare2(obs,est02(i),opt,2,1,fid2); % 2 = m(2), 1 = pair = AT,CT

end

% fid2 = 'compare_TC_pH_K10v6.csv';
% [A] = compare2(obs,est,opt,2,2,fid2);

% % AT pH (Q2) (fid3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;         obs(i).eTC = nan;
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
%     obs(i).m(2).co3 = nan;   obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,opt); % need to change opt.fid for each
% est03 = est;


% % TC pH (Q2) (fid4)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan;         obs(i).eTA = nan;
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
%     obs(i).m(2).co3 = nan;   obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,opt); % need to change opt.fid for each
% est04 = est;

% fid2 = 'compare_TC_pH_K10v6.csv';
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


% % pCO2 CO3 (Q2) (fid11)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan; obs(i).eTC = nan;
%     obs(i).TA = nan; obs(i).eTA = nan;
%     obs(i).m(2).ph = nan; obs(i).m(2).eph = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est11 = est;


% % TC CO3 (fid12)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan; obs(i).eTA = nan;
%     obs(i).m(2).ph = nan; obs(i).m(2).eph = nan;
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est12 = est;


% % TC pCO2 (fid13)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan; obs(i).eTC = nan;
%     obs(i).m(2).ph = nan; obs(i).m(2).eph = nan;
%     obs(i).m(2).co3 = nan;   obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est13 = est;

% % TA pCO2 (fid14)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan; obs(i).eTC = nan;
%     obs(i).m(2).ph = nan; obs(i).m(2).eph = nan;
%     obs(i).m(2).co3 = nan;   obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est14 = est;

% % TA CO3 (fid15)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan; obs(i).eTC = nan;
%     obs(i).m(2).ph = nan; obs(i).m(2).eph = nan;
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est15 = est;


% % pH CO3 (fid 16)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan; obs(i).eTC = nan;
%     obs(i).TA = nan; obs(i).eTA = nan;
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est16 = est;


% % TC TA pH (fid20)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
%     obs(i).m(2).co3 = nan;   obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est20 = est;


% % TC TA pCO2 (fid21)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(2).ph = nan; obs(i).m(2).eph = nan;
%     obs(i).m(2).co3 = nan;   obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est21 = est;


% % TC TA CO3 (fid22)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(2).ph = nan; obs(i).m(2).eph = nan;
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est22 = est;


% % TC pH pCO2 (fid23)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan; obs(i).eTA = nan;
%     obs(i).m(2).co3 = nan;   obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est23 = est;


% % TC pH CO3 (fid24)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan; obs(i).eTA = nan;
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est24 = est;


% % TC pCO2 CO3 (fid25)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan; obs(i).eTA = nan;
%     obs(i).m(2).ph = nan; obs(i).m(2).eph = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est25 = est;


% % TA pH pCO2 (fid26)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan; obs(i).eTC = nan;
%     obs(i).m(2).co3 = nan;   obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est26 = est;


% % TA pH CO3 (fid27)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan; obs(i).eTC = nan;
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est27 = est;


% % TA pCO2 CO3 (fid28)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan; obs(i).eTC = nan;
%     obs(i).m(2).ph = nan; obs(i).m(2).eph = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est28 = est;


% % pH pCO2 CO3 (fid29)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan; obs(i).eTC = nan;
%     obs(i).TA = nan; obs(i).eTA = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est29 = est;


% % TC TA pH pCO2 (fid30)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(2).co3 = nan;   obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est30 = est;


% % TC TA pH CO3 (fid31)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(3).pco2 = nan;  obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,opt);
% est31 = est;


% % TC TA pCO2 CO3 (fid32)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(2).ph = nan; obs(i).m(2).eph = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est32 = est;


% % TC pH pCO2 CO3 (fid33)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan; obs(i).eTA = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est33 = est;


% % TA pH pCO2 CO3 (fid34)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan; obs(i).eTC = nan;
% end
% [est,obs,iflag] = QUODcarbV6(obs,sys);
% est34 = est;


% % CT AT pH pCO2 CO3 (Q5) (fid5)
% obs = obs_backup;
% [est,obs,iflag] = QUODcarbV6(obs,opt);
% est05 = est;


% save estK10v6/est02K10v6.mat est02;
% save estK10v6/est03K10v6.mat est03;
% save estK10v6/est04K10v6.mat est04;
% save estK10v6/est05K10v6.mat est05;
% save estK10v6/est10K10v6.mat est10;
% save estK10v6/est11K10v6.mat est11;
% save estK10v6/est12K10v6.mat est12;
% save estK10v6/est13K10v6.mat est13;
% save estK10v6/est14K10v6.mat est14;
% save estK10v6/est15K10v6.mat est15;
% save estK10v6/est16K10v6.mat est16;
% save estK10v6/est20K10v6.mat est20;
% save estK10v6/est21K10v6.mat est21;
% save estK10v6/est22K10v6.mat est22;
% save estK10v6/est23K10v6.mat est23;
% save estK10v6/est24K10v6.mat est24;
% save estK10v6/est25K10v6.mat est25;
% save estK10v6/est26K10v6.mat est26;
% save estK10v6/est27K10v6.mat est27;
% save estK10v6/est28K10v6.mat est28;
% save estK10v6/est29K10v6.mat est29;
% save estK10v6/est30K10v6.mat est30;
% save estK10v6/est31K10v6.mat est31;
% save estK10v6/est32K10v6.mat est32;
% save estK10v6/est33K10v6.mat est33;
% save estK10v6/est34K10v6.mat est34;



% % to make Orr-style plots, divide uncertainties by 1000 % try 100 next
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(1).epK1 = 0.0055/1000; % Lueker's uncertainty/1000
%     obs(i).m(2).epK1 = 0.0055/1000;
%     obs(i).m(3).epK1 = 0.0055/1000;
% end
% [est,obs,iflag] = QUODcarbV6(obs,opt);
% est05pK1 = est;
% 
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(1).epK2 = 0.0100/1000; % Lueker's uncertainty/1000
%     obs(i).m(2).epK2 = 0.0100/1000;
%     obs(i).m(3).epK2 = 0.0100/1000;
% end
% [est,obs,iflag] = QUODcarbV6(obs,opt);
% est05pK2 = est;
% 
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(1).epKar = 0.039/1000;
%     obs(i).m(2).epKar = 0.039/1000;
%     obs(i).m(3).epKar = 0.039/1000;
% end
% [est,obs,iflag] = QUODcarbV6(obs,opt);
% est05pKar = est;
% 
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(1).epK0 = 0.0055/1000;
%     obs(i).m(2).epK0 = 0.0055/1000;
%     obs(i).m(3).epK0 = 0.0055/1000;
% end
% [est,obs,iflag] = QUODcarbV6(obs,opt);
% est05pK0 = est;
% 
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(1).epKb = (0.004/log(10))/1000;
%     obs(i).m(2).epKb = (0.004/log(10))/1000;
%     obs(i).m(3).epKb = (0.004/log(10))/1000;
% end
% [est,obs,iflag] = QUODcarbV6(obs,opt);
% est05pKb = est;

% epTB I'll have to do in the code 

% epK0 I'll have to do in the code













