
% drivergV7 to go with QUODcarbV7

% updated opt (options) structure

% g = GOMECC data

format long

load datag.mat;
[in] = datag;
nD = 5; %1186; % until UW start, or before: %length(in);

% choose options for opt structure
opt.K1K2 = 10; % option for K1K2 formulation
opt.KSO4 = 1; % option for KSO4 formulation
opt.KF   = 2; % option for KF formulation
opt.TB   = 2; % option for TB formulation
opt.phscale  = 2;  % 1 = tot, 2 = sws, 3 = free, 4 = NBS
opt.printcsv = 1; % print est to CSV? 1 = on , 0 = off
opt.fid      = 'CT_ph.csv'; % don't need it if printcsv is off
opt.printmes = 1; % print screen messages? 1 = on, 0 = off
opt.abr = 'all'; % option for which acid/base reactions to include
    % opt.abr = {'all'}; % same as below
    % opt.abr = {'borate','sulfate','fluoride','phosphate', ...
    %               'silicate','ammonia','sulfide','solubility'};
opt.co2press = 1; % 1 = on, 0 = off (default off)

ad = 0; % for debugging, so can start at i = 1
% read in GOMECC data and put into obs structure

for i = 1:nD
    
    % measurements that are independent of (T,P)
    obs(i).TC    = in(i+ad,5); % (umol/kg)
    obs(i).eTC   = 2.01;    % TC error ±2.01 umol/kg
    if i >= 1187
        obs(i).eTC = 2.24; % TC error ±2.24 for UW measurements
    end
    obs(i).TA    = in(i+ad,7);
    obs(i).eTA   = 1.78; % TA error ±2.01 umol/kg
    obs(i).sal   = in(i+ad,3); % PSU
    obs(i).esal  = 0.002; % new = 0.001
    % nutrients P and Si also independent of (T,P)
    obs(i).TP    = in(i+ad,24);
    obs(i).eTP   = in(i+ad,24)*0.001; % 0.01% meas uncertainty
    obs(i).TSi   = in(i+ad,18);
    obs(i).eTSi  = in(i+ad,18)*0.001; % 0.01% meas uncertainty

    % first (T,P)-dependent measurement
    % obs(i).m(1).T = in(i+ad,2); % deg C, CTD temp
    % obs(i).m(1).eT = 0.02; % ±0.02 degC
    % obs(i).m(1).P = in(i+ad,1); % dbar
    % obs(i).m(1).eP = 0.63; % (max) ± 0.63 dbar

    % second(T,P)-dependent measurement
    obs(i).m(1).T    = 25 ; % degC
    obs(i).m(1).eT   = 0.05 ; % from cruise report
    obs(i).m(1).P    = 0.0 ; %in(i+ad,1); % NOT in situ
    obs(i).m(1).eP   = 0.63 ;
    obs(i).m(1).ph   = in(i+ad,9); % total scale
    obs(i).m(1).eph  = 0.0004 ;
    if i >= 1187 % DON"T WANT THIS ANYMORE
        obs(i).m(1).eph = 0.05; % for UW measurements
    end
    obs(i).m(1).co3  = in(i+ad,15); % (µmol/kg)
    obs(i).m(1).eco3 = 2.0 ;  % std ±2µmol/kg

    % third (T,P)-dependent measurement
    % obs(i).m(3).T   = 20 ; %degC
    % obs(i).m(3).eT  = 0.03 ; % from cruise report
    % obs(i).m(3).P   = 0.0 ; % dbar (surface pressure for pco2)
    % obs(i).m(3).eP  = 0.63 ;
    % obs(i).m(3).pco2  = in(i+ad,12); % (µatm)
    % obs(i).m(3).epco2 = in(i+ad,12)*0.0021; % 0.21% relative std error (avg)

end
fid2 = 'compare_CT_AT.csv';
% [A] = compare2(obs,est,opt,2,1,fid2);
obs_backup = obs;

% tp = 1;
% % TA TC (Q2) (fid2)
% for i = 1:nD
%     obs(i).m(1).ph = nan;   obs(i).m(1).eph = nan;
%     % obs(i).m(3).pco2 = nan; obs(i).m(3).epco2 = nan;
%     obs(i).m(1).co3 = nan;  obs(i).m(1).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est02   = est;
% [A] = compare2(obs,est,opt,tp,1,fid2);

% % TC ph (Q2) (fid4)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;        obs(i).eTA = nan;
    obs(i).m(1).pco2 = nan; obs(i).m(1).epco2 = nan; % m(3)
    obs(i).m(1).co3 = nan;  obs(i).m(1).eco3 = nan; % m(2)
end
[est,obs,iflag] = QUODcarbV7(obs,opt);
% est04   = est;
% tp = 1;
% fid3 = 'compare_CT_ph.csv';
% [A] = compare2(obs,est,opt,tp,2,fid3);


% CT AT pH pCO2 CO3 (Q5) (fid5)
% obs = obs_backup;
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est33 = est;


% 
% 
% % TA ph (Q2) (fid3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;        obs(i).eTC = nan;
%     obs(i).m(3).pco2 = nan; obs(i).m(3).epco2 = nan;
%     obs(i).m(2).co3 = nan;  obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est03   = est;
 
% % pH pCO2 (Q2) (fid10)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;        obs(i).eTC = nan;
%     obs(i).TA = nan;        obs(i).eTA = nan;
%     obs(i).m(2).co3 = nan;  obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est10   = est;
% 
% 
% % pCO2 CO3 (Q2) (fid11)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;        obs(i).eTC = nan;
%     obs(i).TA = nan;        obs(i).eTA = nan;
%     obs(i).m(2).ph = nan;   obs(i).m(2).eph = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est11   = est;
% 
% 
% % TC CO3 (Q2)(fid12)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan;        obs(i).eTA = nan;
%     obs(i).m(3).pco2 = nan; obs(i).m(3).epco2 = nan;
%     obs(i).m(2).ph = nan;   obs(i).m(2).eph = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est12   = est;
% 
% 
% % TC pCO2 (Q2)(fid13)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan;        obs(i).eTA = nan;
%     obs(i).m(2).ph = nan;   obs(i).m(2).eph = nan;
%     obs(i).m(2).co3 = nan;  obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est13   = est;
% 
% % TA pCO2 (Q2)(fid14)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;        obs(i).eTC = nan;
%     obs(i).m(2).ph = nan;   obs(i).m(2).eph = nan;
%     obs(i).m(2).co3 = nan;  obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est14   = est;
% 
% 
% % TA CO3 (Q2)(fid15)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;        obs(i).eTC = nan;
%     obs(i).m(2).ph = nan;   obs(i).m(2).eph = nan;
%     obs(i).m(3).pco2 = nan; obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est15   = est;
% 
% 
% % pH CO3 (Q2) (fid16)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan; obs(i).eTC = nan;
%     obs(i).TA = nan; obs(i).eTA = nan;
%     obs(i).m(3).pco2 = nan; obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est16   = est;
% 
% 
% % TC TA pH (Q3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(2).co3 = nan;  obs(i).m(2).eco3 = nan;
%     obs(i).m(3).pco2 = nan; obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est20   = est;
% 
% 
% % TC TA pco2 (Q3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(2).co3 = nan;  obs(i).m(2).eco3 = nan;
%     obs(i).m(2).ph = nan;   obs(i).m(2).eph = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est21   = est;
% 
% 
% % TC TA co3 (Q3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(3).pco2 = nan; obs(i).m(2).epco2 = nan;
%     obs(i).m(2).ph = nan;   obs(i).m(2).eph = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est22   = est;
% 
% 
% % TC pH pco2 (Q3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan;        obs(i).eTA = nan;
%     obs(i).m(2).co3 = nan;  obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est23   = est;
% 
% 
% % TC pH co3 (Q3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan;        obs(i).eTA = nan;
%     obs(i).m(3).pco2 = nan; obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est24   = est;
% 
% 
% % TC pco2 co3 (Q3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan;        obs(i).eTA = nan;
%     obs(i).m(2).ph = nan;   obs(i).m(2).eph = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est25   = est;
% 
% 
% % TA pH pco2 (Q3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;        obs(i).eTC = nan;
%     obs(i).m(2).co3 = nan;  obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est26   = est;
% 
% 
% % TA pH co3 (Q3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;        obs(i).eTC = nan;
%     obs(i).m(3).pco2 = nan; obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est27   = est;
% 
% 
% % TA pco2 co3 (Q3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;        obs(i).eTC = nan;
%     obs(i).m(2).ph = nan;   obs(i).m(2).eph = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est28   = est;
% 
% 
% % pH pco2 co3 (Q3)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan;    obs(i).eTA = nan;
%     obs(i).TC = nan;    obs(i).eTC = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est29   = est;
% 
% 
% % TC TA pH pco2 (Q4)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(2).co3 = nan; obs(i).m(2).eco3 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est30   = est;
% 
% 
% % TC TA pH co3 (Q4)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(3).pco2 = nan; obs(i).m(3).epco2 = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est31   = est;
% 
% 
% % TC TA pco2 co3 (Q4)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).m(2).ph = nan; obs(i).m(2).eph = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est32   = est;
% 
% 
% % TC pH pco2 co3 (Q4)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan; obs(i).eTA = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est33   = est;
% 
% 
% % TA pH pco2 co3 (Q4)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan; obs(i).eTC = nan;
% end
% [est,obs,iflag] = QUODcarbV7(obs,opt);
% est34  = est;



%save calc_f_files/est05V7K11.mat est05;
% 
% save estV7/est02.mat est02;
% save estV7/est03.mat est03;
% save estV7/est04.mat est04;
% save estV7/est05.mat est05;
% save estV7/est10.mat est10;
% save estV7/est11.mat est11;
% save estV7/est12.mat est12;
% save estV7/est13.mat est13;
% save estV7/est14.mat est14;
% save estV7/est15.mat est15;
% save estV7/est16.mat est16;
% save estV7/est20.mat est20;
% save estV7/est21.mat est21;
% save estV7/est22.mat est22;
% save estV7/est23.mat est23;
% save estV7/est24.mat est24;
% save estV7/est25.mat est25;
% save estV7/est26.mat est26;
% save estV7/est27.mat est27;
% save estV7/est28.mat est28;
% save estV7/est29.mat est29;
% save estV7v2/est30.mat est30;
% save estV7v2/est31.mat est31;
% save estV7v2/est32.mat est32;
% save estV7v2/est33.mat est33;
% save estV7v2/est34.mat est34;















