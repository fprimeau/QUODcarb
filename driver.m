
% driver to go with QUODcarb

load data.mat; % NEW as of Nov.11
[in] = data;
nD = 25; %length(in);

% choose options for opt structure
opt.K1K2 = 16; % option for K1K2 formulation
opt.KSO4 = 1;  % option for KSO4 formulation
opt.KF   = 2;  % option for KF formulation
opt.TB   = 2;  % option for TB formulation
opt.phscale  = 1;  % 1 = tot, 2 = sws, 3 = free, 4 = NBS
opt.printcsv = 0;  % print est to CSV? 1 = on , 0 = off
% opt.fname    = 'QUODcarb_output.csv'; % don't need it if printcsv is off
opt.fname    = 'output_csv/Q5K16.csv';
opt.co2press = 1; % 1 = on, 0 = off
opt.Revelle  = 0; % 1 = on, 0 = off 
opt.printmes = 0; % 1 = on, 0 = off


% read in GOMECC data and put into obs structure
for i = 1:nD
    
    % measurements that are independent of (T,P)
    obs(i).TC    = in(5,i); % (umol/kg)
    obs(i).eTC   = 2.01;    % TC error ±2.01 umol/kg
    obs(i).TA    = in(6,i);
    obs(i).eTA   = 1.78; % TA error ±2.01 umol/kg
    obs(i).sal   = in(1,i); % PSU
    obs(i).esal  = 0.001; % 1 new as of 1/23 old = 0.002
    % nutrients P and Si also independent of (T,P)
    obs(i).TP    = in(7,i);
    obs(i).eTP   = in(7,i)*0.003; % 0.30% meas precision NEW 4/17/24
    obs(i).TSi   = in(8,i);
    obs(i).eTSi  = in(8,i)*0.0031; % 0.31% meas uncertainty NEW 4/17/24

    % first (T,P)-dependent measurement
    obs(i).tp(1).T  = in(2,i); % deg C, CTD temp
    obs(i).tp(1).eT = 0.02; % ±0.02 degC
    obs(i).tp(1).P  = in(3,i); % dbar
    obs(i).tp(1).eP = 0.63; % (max) ± 0.63 dbar

    % second(T,P)-dependent measurement
    obs(i).tp(2).T    = 25 ; % degC
    obs(i).tp(2).eT   = 0.05 ; % from cruise report
    obs(i).tp(2).P    = 0.0 ; %in(i+ad,1); % NOT in situ
    obs(i).tp(2).eP   = 0.07 ;
    obs(i).tp(2).ph   = in(9,i); % total scale
    obs(i).tp(2).eph  = 0.0004 ;
    obs(i).tp(2).co3  = in(11,i); % (µmol/kg)
    obs(i).tp(2).eco3 = in(11,1)*0.02;  % 2% from Jon Sharp NEW 1/25/24

    % third (T,P)-dependent measurement
    obs(i).tp(3).T     = 20 ; %degC
    obs(i).tp(3).eT    = 0.03 ; % from cruise report
    obs(i).tp(3).P     = 0.0 ; % dbar (surface pressure for pco2)
    obs(i).tp(3).eP    = 0.07 ;
    obs(i).tp(3).pco2  = in(10,i); % (µatm)
    obs(i).tp(3).epco2 = in(10,i)*0.0021; % 0.21% relative std error (avg)
end

obs_backup = obs;

%% Q5: All five input
% CT AT pH pCO2 CO3 (Q5) (fid5)
[est,obs,sys,iflag] = QUODcarb(obs,opt); 


% K04
<<<<<<< HEAD
tic
opt.K1K2 = 4;
[est,~,sys,~,~] = QUODcarb(obs,opt); % [est, obs, sys, iflag]

estK04 = est;
% fprintf('K04 done \n')
toc

% K10
opt.K1K2 = 10;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK10 = est;

% K11
opt.K1K2 = 11;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK11 = est;

% K12
opt.K1K2 = 12;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK12 = est;

% K13
opt.K1K3 = 13;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK13 = est;

% K14
opt.K1K2 = 14;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK14 = est;

% K15
opt.K1K2 = 15;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK15 = est;

% K16
opt.K1K2 = 16;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK16 = est;

% K17
opt.K1K2 = 17;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK17 = est;
=======
<<<<<<< HEAD
<<<<<<< HEAD
%[est,~,~,~] = QUODcarb(obs,opt);
%estK04 = est;

% K10
%obs = obs_backup;
%opt.K1K2 = 10;
%[est,~,~,~] = QUODcarb(obs,opt);
%estK10 = est;

% K11
%obs = obs_backup;
%opt.K1K2 = 11;
%[est,~,~,~] = QUODcarb(obs,opt);
%estK11 = est;

% K12
%obs = obs_backup;
%opt.K1K2 = 12;
%[est,~,~,~] = QUODcarb(obs,opt);
%estK12 = est;

% K13
%obs = obs_backup;
%opt.K1K2 = 13;
%[est,~,~,~] = QUODcarb(obs,opt);
%estK13 = est;

% K14
%obs = obs_backup;
%opt.K1K2 = 14;
%[est,~,~,~] = QUODcarb(obs,opt);
%estK14 = est;

% K15
%obs = obs_backup;
%opt.K1K2 = 15;
%[est,~,~,~] = QUODcarb(obs,opt);
%estK15 = est;

% K16
%obs = obs_backup;
%opt.K1K2 = 16;
%[est,~,~,~] = QUODcarb(obs,opt);
%estK16 = est;

% K17
%obs = obs_backup;
%opt.K1K2 = 17;
%[est,~,~,~] = QUODcarb(obs,opt);
%estK17 = est;
>>>>>>> main

%save output_mat_files/pKi/estK10.mat estK10;
%save output_mat_files/pKi/estK04.mat estK04;
%save output_mat_files/pKi/estK11.mat estK11;
%save output_mat_files/pKi/estK12.mat estK12;
%save output_mat_files/pKi/estK13.mat estK13;
%save output_mat_files/pKi/estK14.mat estK14;
%save output_mat_files/pKi/estK15.mat estK15;
%save output_mat_files/pKi/estK16.mat estK16;
%save output_mat_files/pKi/estK17.mat estK17;
=======
% tic
% opt.K1K2 = 4;
% [est,~,~,~] = QUODcarb(obs,opt); % [est, obs, sys, iflag]
% estK04 = est;
% fprintf('K04 done \n')
% toc
% keyboard

% K10
% opt.K1K2 = 10;
% [est,~,~,~] = QUODcarb(obs,opt);
% estK10 = est;
>>>>>>> eef4a5cd5ee6bd54dfcc178ab32c16b76e173276

% K11
% opt.K1K2 = 11;
% [est,~,~,~] = QUODcarb(obs,opt);
% estK11 = est;

% K12
% opt.K1K2 = 12;
% [est,~,~,~] = QUODcarb(obs,opt);
% estK12 = est;

% K13
% opt.K1K3 = 13;
% [est,~,~,~] = QUODcarb(obs,opt);
% estK13 = est;

% K14
% opt.K1K2 = 14;
% [est,~,~,~] = QUODcarb(obs,opt);
% estK14 = est;

% K15
% opt.K1K2 = 15;
% [est,~,~,~] = QUODcarb(obs,opt);
% estK15 = est;

% K16
% opt.K1K2 = 16;
% [est,~,~,~] = QUODcarb(obs,opt);
% estK16 = est;

% K17
% opt.K1K2 = 17;
% [est,~,~,~] = QUODcarb(obs,opt);
% estK17 = est;

% for i = 1:nD
%     fK04(i) = estK04(i).f;
%     fK10(i) = estK10(i).f;
%     fK11(i) = estK11(i).f;
%     fK12(i) = estK12(i).f;
%     fK13(i) = estK13(i).f;
%     fK14(i) = estK14(i).f;
%     fK15(i) = estK15(i).f;
%     fK16(i) = estK16(i).f;
%     fK17(i) = estK17(i).f;
% end

% save output_mat_files/new_pK/fK04.mat fK04;
% save output_mat_files/new_pK/fK10.mat fK10;
% save output_mat_files/new_pK/fK11.mat fK11;
% save output_mat_files/new_pK/fK12.mat fK12;
% save output_mat_files/new_pK/fK13.mat fK13;
% save output_mat_files/new_pK/fK14.mat fK14;
% save output_mat_files/new_pK/fK15.mat fK15;
% save output_mat_files/new_pK/fK16.mat fK16;
% save output_mat_files/new_pK/fK17.mat fK17;

% get other stuff from driver_all.m, updated 1/25/24

% 
%% Q5: All five input
% CT AT pH pCO2 CO3 (Q5) (fid26)
% obs = obs_backup;
% [est,~,~,~,~] = QUODcarb(obs,opt);
% est26 = est;
% 
% %% Q2: Input Pairs
% 
% % TC TA (Q2) (fid01)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).tp(2).ph = nan;   obs(i).tp(2).eph = nan;
%     obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan;
%     obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan;
% end
% [est,~,~,~,~] = QUODcarb(obs,opt);
% est01  = est;
<<<<<<< HEAD
% % [est,obs,~,~,~] = QUODcarb(obs,opt); % need obs for compare
% % fid01   = 'compare_outs/compare_TC_TA.csv'; 
% % tp     = 2; % second tp system for ph in there
% % A      = compare(obs,est,opt,tp,1,fid01); % 1 for input pair TC TA
% 
% % TC ph (Q2) (fid02)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan;         obs(i).eTA = nan;
%     obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan; % tp(3)
%     obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan; % tp(2)
% end
% [est,~,~,~,~] = QUODcarb(obs,opt);
% est02   = est;
% % [est,obs,~,~,~] = QUODcarb(obs,opt);
% % tp      = 2;
% % fid02   = 'compare_outs/compare_TC_ph.csv';
% % [A]     = compare(obs,est,opt,tp,2,fid02);
=======
% [est,obs,~,~] = QUODcarb(obs,opt); % need obs for compare
% fid01   = 'compare_outs/compare_TC_TA.csv'; 
% tp     = 2; % second tp system for ph in there
% A      = compare(obs,est,opt,tp,1,fid01); % 1 for input pair TC TA

% TC ph (Q2) (fid02)
%obs = obs_backup;
%for i = 1:nD
%    obs(i).TA = nan;         obs(i).eTA = nan;
%    obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan; % tp(3)
%    obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan; % tp(2)
%end
%[est,~,~,~] = QUODcarb(obs,opt);
%est02   = est;
% [est,obs, ~, ~] = QUODcarb(obs,opt);
% tp      = 2;
% fid02   = 'compare_outs/compare_TC_ph.csv';
% [A]     = compare(obs,est,opt,tp,2,fid02);
>>>>>>> main
% 
% % TC pCO2 (Q2)(fid03)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan;         obs(i).eTA = nan;
%     obs(i).tp(2).ph = nan;   obs(i).tp(2).eph = nan;
%     obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan;
% end
% [est,~,~,~,~] = QUODcarb(obs,opt);
% est03   = est;
% % [est,obs,~,~,~] = QUODcarb(obs,opt);
% % tp      = 3;
% % fid03   = 'compare_outs/compare_TC_pco2.csv';
% % [A]     = compare(obs,est,opt,tp,7,fid03);
% 
% % TC CO3 (Q2)(fid04)
<<<<<<< HEAD
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TA = nan;         obs(i).eTA = nan;
%     obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan;
%     obs(i).tp(2).ph = nan;   obs(i).tp(2).eph = nan;
% end
% [est,~,~,~,~] = QUODcarb(obs,opt);
% est04   = est;
% % [est,obs,~,~,~] = QUODcarb(obs,opt);
% % tp      = 2;
% % fid04   = 'compare_outs/compare_TC_co3.csv';
% % [A]     = compare(obs,est,opt,tp,6,fid04);
=======
%obs = obs_backup;
%for i = 1:nD
%    obs(i).TA = nan;         obs(i).eTA = nan;
%    obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan;
%    obs(i).tp(2).ph = nan;   obs(i).tp(2).eph = nan;
%end
%[est,~,~,~] = QUODcarb(obs,opt);
%est04   = est;
% [est,obs, ~, ~] = QUODcarb(obs,opt);
% tp      = 2;
% fid04   = 'compare_outs/compare_TC_co3.csv';
% [A]     = compare(obs,est,opt,tp,6,fid04);
>>>>>>> main
% 
% % TA ph (Q2) (fid05)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;         obs(i).eTC = nan;
%     obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan;
%     obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan;
% end
% [est,~,~,~,~] = QUODcarb(obs,opt);
% est05   = est;
% % [est,obs,~,~,~] = QUODcarb(obs,opt);
% % tp      = 2;
% % fid05   = 'compare_outs/compare_TA_ph.csv';
% % [A]     = compare(obs,est,opt,tp,3,fid05);
% 
% % TA pCO2 (Q2)(fid06)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;        obs(i).eTC = nan;
%     obs(i).tp(2).ph = nan;  obs(i).tp(2).eph = nan;
%     obs(i).tp(2).co3 = nan; obs(i).tp(2).eco3 = nan;
% end
% [est,~,~,~,~] = QUODcarb(obs,opt);
% est06   = est;
% % [est,obs,~,~,~] = QUODcarb(obs,opt);
% % tp      = 3;
% % fid06   = 'compare_outs/compare_TA_pco2.csv';
% % [A]     = compare(obs,est,opt,tp,8,fid06);
% 
% % TA CO3 (Q2)(fid07)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;          obs(i).eTC = nan;
%     obs(i).tp(2).ph = nan;    obs(i).tp(2).eph = nan;
%     obs(i).tp(3).pco2 = nan;  obs(i).tp(3).epco2 = nan;
% end
% [est,~,~,~,~] = QUODcarb(obs,opt);
% est07   = est;
% % [est,obs,~,~,~] = QUODcarb(obs,opt);
% % tp      = 2;
% % fid07   = 'compare_outs/compare_TA_co3.csv';
% % [A]     = compare(obs,est,opt,tp,9,fid07);
% 
% % % pH CO3 (Q2) (fid08)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;            obs(i).eTC = nan;
%     obs(i).TA = nan;            obs(i).eTA = nan;
%     obs(i).tp(3).pco2 = nan;    obs(i).tp(3).epco2 = nan;
% end
% [est,~,~,~,~] = QUODcarb(obs,opt);
% est08   = est;
% % [est,obs,~,~,~] = QUODcarb(obs,opt);
% % tp      = 2;
% % fid08   = 'compare_outs/compare_ph_CO3.csv';
% % [A]     = compare(obs,est,opt,tp,10,fid08);
% 
% % pH pCO2 (Q2) (fid09)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;        obs(i).eTC = nan;
%     obs(i).TA = nan;        obs(i).eTA = nan;
%     obs(i).tp(2).co3 = nan; obs(i).tp(2).eco3 = nan;
% end
% [est,~,~,~,~] = QUODcarb(obs,opt);
% est09 = est;
% % % NO compare! Wrong tp's
% 
% % pCO2 CO3 (Q2) (fid10)
% obs = obs_backup;
% for i = 1:nD
%     obs(i).TC = nan;        obs(i).eTC = nan;
%     obs(i).TA = nan;        obs(i).eTA = nan;
%     obs(i).tp(2).ph = nan;  obs(i).tp(2).eph = nan;
% end
% [est,~,~,~,~] = QUODcarb(obs,opt);
% est10   = est;
% % % NO compare! Wrong tp's
<<<<<<< HEAD
=======

>>>>>>> main
