
% driver to go with QUODcarb

% user must create an 'output_mat_files' folder to hold mat files
% plus subfolders 'all_combos', 'all_pKs', and 'parsed_f'
% or change the 'save' lines to path of user's preference

load data.mat; % from 'parse_gomecc3_data.m' script
[in] = data;
nD = length(in);

% choose options for opt structure
opt.K1K2 = 16; % option for K1K2 formulation
opt.KSO4 = 2;  % option for KSO4 formulation
opt.KF   = 1;  % option for KF formulation
opt.TB   = 1;  % option for TB formulation
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
    obs(i).tp(3).T     = 20 ; % degC
    obs(i).tp(3).eT    = 0.03 ; % from cruise report
    obs(i).tp(3).P     = 0.0 ; % dbar (surface pressure for pco2)
    obs(i).tp(3).eP    = 0.07 ;
    obs(i).tp(3).pco2  = in(10,i); % (µatm)
    obs(i).tp(3).epco2 = in(10,i)*0.0021; % 0.21% relative std error (avg)
end

obs_backup = obs;

% run through all nD datapoints 35 times -> 
% 26 combinations plus 9 pK formulations

%% Q5: All five input
% CT AT pH pCO2 CO3 (Q5) (fid26)
obs = obs_backup;
[est,~,~,~,~] = QUODcarb(obs,opt);
est26 = est;
save output_mat_files/all_combos/est26.mat est26; clear est; clear obs;

%% Q2: Input pairs
% TC TA (Q2) (fid01)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(2).ph = nan;   obs(i).tp(2).eph = nan;
    obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan;
    obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est01  = est;
save output_mat_files/all_combos/est01.mat est01; clear est; clear obs;

% % TC ph (Q2) (fid02)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;         obs(i).eTA = nan;
    obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan; % tp(3)
    obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan; % tp(2)
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est02   = est;
save output_mat_files/all_combos/est02.mat est02; clear est; clear obs;
 
% TC pCO2 (Q2)(fid03)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;         obs(i).eTA = nan;
    obs(i).tp(2).ph = nan;   obs(i).tp(2).eph = nan;
    obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est03   = est;
save output_mat_files/all_combos/est03.mat est03; clear est; clear obs;

% % TC CO3 (Q2)(fid04)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;         obs(i).eTA = nan;
    obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan;
    obs(i).tp(2).ph = nan;   obs(i).tp(2).eph = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est04   = est;
save output_mat_files/all_combos/est04.mat est04; clear est; clear obs;

% TA ph (Q2) (fid05)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;         obs(i).eTC = nan;
    obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan;
    obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est05   = est;
save output_mat_files/all_combos/est05.mat est05; clear est; clear obs;

% TA pCO2 (Q2)(fid06)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;        obs(i).eTC = nan;
    obs(i).tp(2).ph = nan;  obs(i).tp(2).eph = nan;
    obs(i).tp(2).co3 = nan; obs(i).tp(2).eco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est06   = est;
save output_mat_files/all_combos/est06.mat est06; clear est; clear obs;

% TA CO3 (Q2)(fid07)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;          obs(i).eTC = nan;
    obs(i).tp(2).ph = nan;    obs(i).tp(2).eph = nan;
    obs(i).tp(3).pco2 = nan;  obs(i).tp(3).epco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est07   = est;
save output_mat_files/all_combos/est07.mat est07; clear est; clear obs;

% pH CO3 (Q2) (fid08)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;            obs(i).eTC = nan;
    obs(i).TA = nan;            obs(i).eTA = nan;
    obs(i).tp(3).pco2 = nan;    obs(i).tp(3).epco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est08   = est;
save output_mat_files/all_combos/est08.mat est08; clear est; clear obs;

% pH pCO2 (Q2) (fid09)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;        obs(i).eTC = nan;
    obs(i).TA = nan;        obs(i).eTA = nan;
    obs(i).tp(2).co3 = nan; obs(i).tp(2).eco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est09 = est;
save output_mat_files/all_combos/est09.mat est09; clear est; clear obs;

% pCO2 CO3 (Q2) (fid10)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;        obs(i).eTC = nan;
    obs(i).TA = nan;        obs(i).eTA = nan;
    obs(i).tp(2).ph = nan;  obs(i).tp(2).eph = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est10   = est;
save output_mat_files/all_combos/est10.mat est10; clear est; clear obs;


%% Q4: Four input

% TC TA pH pco2 (Q4) (fid21)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(2).co3 = nan; obs(i).tp(2).eco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est21   = est;
save output_mat_files/all_combos/est21.mat est21; clear est; clear obs;

% % TC TA pH co3 (Q4) (fid22)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(3).pco2 = nan; obs(i).tp(3).epco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est22   = est;
save output_mat_files/all_combos/est22.mat est22; clear est; clear obs;

% TC TA pco2 co3 (Q4) (fid23)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(2).ph = nan; obs(i).tp(2).eph = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est23   = est;
save output_mat_files/all_combos/est23.mat est23; clear est; clear obs;

% TC pH pco2 co3 (Q4) (fid24)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan; obs(i).eTA = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est24   = est;
save output_mat_files/all_combos/est24.mat est24; clear est; clear obs;

% TA ph pco2 co3 (Q4) (fid25)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan; obs(i).eTC = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est25   = est;
save output_mat_files/all_combos/est25.mat est25; clear est; clear obs;


%% Q3: Three input

% % TC TA pH (Q3) (fid11)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(2).co3 = nan;   obs(i).tp(2).eco3 = nan;
    obs(i).tp(3).pco2 = nan;  obs(i).tp(3).epco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt); 
est11   = est;
save output_mat_files/all_combos/est11.mat est11; clear est; clear obs;

% TC TA pco2 (Q3) (fid12)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(2).co3 = nan;  obs(i).tp(2).eco3 = nan;
    obs(i).tp(2).ph = nan;   obs(i).tp(2).eph = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est12   = est;
save output_mat_files/all_combos/est12.mat est12; clear est; clear obs;

% TC TA co3 (Q3) (fid13)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(3).pco2 = nan; obs(i).tp(2).epco2 = nan;
    obs(i).tp(2).ph = nan;   obs(i).tp(2).eph = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est13   = est;
save output_mat_files/all_combos/est13.mat est13; clear est; clear obs;

% TC pH pco2 (Q3) (fid14)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;          obs(i).eTA = nan;
    obs(i).tp(2).co3 = nan;   obs(i).tp(2).eco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est14   = est;
save output_mat_files/all_combos/est14.mat est14; clear est; clear obs;

% TC pH co3 (Q3) (fid15)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;          obs(i).eTA = nan;
    obs(i).tp(3).pco2 = nan;  obs(i).tp(3).epco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est15   = est;
save output_mat_files/all_combos/est15.mat est15; clear est; clear obs;

% TC pco2 co3 (Q3) (fid16)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;          obs(i).eTA = nan;
    obs(i).tp(2).ph = nan;    obs(i).tp(2).eph = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est16   = est;
save output_mat_files/all_combos/est16.mat est16; clear est; clear obs;

% TA pH pco2 (Q3) (fid17)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;          obs(i).eTC = nan;
    obs(i).tp(2).co3 = nan;   obs(i).tp(2).eco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est17   = est;
save output_mat_files/all_combos/est17.mat est17; clear est; clear obs;

% TA pH co3 (Q3) (fid18)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;          obs(i).eTC = nan;
    obs(i).tp(3).pco2 = nan;  obs(i).tp(3).epco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est18 = est;
save output_mat_files/all_combos/est18.mat est18; clear est; clear obs;

% % TA pco2 co3 (Q3) (fid19)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;          obs(i).eTC = nan;
    obs(i).tp(2).ph = nan;    obs(i).tp(2).eph = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est19   = est;
save output_mat_files/all_combos/est19.mat est19; clear est; clear obs;

% % pH pco2 co3 (Q3) (fid20)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;    obs(i).eTA = nan;
    obs(i).TC = nan;    obs(i).eTC = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est20   = est;
save output_mat_files/all_combos/est20.mat est20; clear est; clear obs;


for i = 1:nD
    f_all26(1,i) = est01(i).f;
    f_all26(2,i) = est02(i).f;
    f_all26(3,i) = est03(i).f;
    f_all26(4,i) = est04(i).f;
    f_all26(5,i) = est05(i).f;
    f_all26(6,i) = est06(i).f;
    f_all26(7,i) = est07(i).f;
    f_all26(8,i) = est08(i).f;
    f_all26(9,i) = est09(i).f;
    f_all26(10,i) = est10(i).f;
    f_all26(11,i) = est11(i).f;
    f_all26(12,i) = est12(i).f;
    f_all26(13,i) = est13(i).f;
    f_all26(14,i) = est14(i).f;
    f_all26(15,i) = est15(i).f;
    f_all26(16,i) = est16(i).f;
    f_all26(17,i) = est17(i).f;
    f_all26(18,i) = est18(i).f;
    f_all26(19,i) = est19(i).f;
    f_all26(20,i) = est20(i).f;
    f_all26(21,i) = est21(i).f;
    f_all26(22,i) = est22(i).f;
    f_all26(23,i) = est23(i).f;
    f_all26(24,i) = est24(i).f;
    f_all26(25,i) = est25(i).f;
    f_all26(26,i) = est26(i).f;
end

save output_mat_files/parsed_f/f_all26.mat f_all26;


%% run all 9 pK formulations with Q5

% K04
obs = obs_backup;
opt.K1K2 = 4;
[est,~,sys,~,~] = QUODcarb(obs,opt); % [est, obs, sys, iflag]
estK04 = est; clear est;

% K10
opt.K1K2 = 10;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK10 = est; clear est;

% K11
opt.K1K2 = 11;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK11 = est; clear est;

% K12
opt.K1K2 = 12;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK12 = est; clear est;

% K13
opt.K1K2 = 13;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK13 = est; clear est;

% K14
opt.K1K2 = 14;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK14 = est; clear est;

% K15
opt.K1K2 = 15;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK15 = est; clear est;

% K16
opt.K1K2 = 16;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK16 = est; clear est;

% K17
opt.K1K2 = 17;
[est,~,~,~,~] = QUODcarb(obs,opt);
estK17 = est; clear est;

save output_mat_files/all_pKs/estK04.mat estK04;
save output_mat_files/all_pKs/estK10.mat estK10;
save output_mat_files/all_pKs/estK11.mat estK11;
save output_mat_files/all_pKs/estK12.mat estK12;
save output_mat_files/all_pKs/estK13.mat estK13;
save output_mat_files/all_pKs/estK14.mat estK14;
save output_mat_files/all_pKs/estK15.mat estK15;
save output_mat_files/all_pKs/estK16.mat estK16;
save output_mat_files/all_pKs/estK17.mat estK17;

for i = 1:nD
    f_allpK(1,i) = estK04(i).f;
    f_allpK(2,i) = estK10(i).f;
    f_allpK(3,i) = estK11(i).f;
    f_allpK(4,i) = estK12(i).f;
    f_allpK(5,i) = estK13(i).f;
    f_allpK(6,i) = estK14(i).f;
    f_allpK(7,i) = estK15(i).f;
    f_allpK(8,i) = estK16(i).f;
    f_allpK(9,i) = estK17(i).f;
end

save output_mat_files/parsed_f/f_allpK.mat f_allpK;















