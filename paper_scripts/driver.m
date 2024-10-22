
% driver to go with 'paper_scripts/QUODcarb'

% user must create an 'output_mat_files' folder to hold mat files
% plus subfolder 'all_combos'
% or change the 'save' lines to path of user's preference

load data.mat; % from 'parse_gomecc3_data.m' script
[in] = data;
nD = length(in);

% choose options for opt structure
opt.K1K2 = 10; % option for K1K2 formulation
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
    obs(i).uTC   = 2.00;
    obs(i).TA    = in(6,i);
    obs(i).uTA   = 2.00;
    obs(i).sal   = in(1,i); % PSU
    obs(i).usal  = 0.001; % 1 new as of 1/23 old = 0.002
    % nutrients P and Si also independent of (T,P)
    obs(i).TP    = in(7,i);
    obs(i).uTP   = in(7,i)*0.02; % 2%
    obs(i).TSi   = in(8,i);
    obs(i).uTSi  = in(8,i)*0.02; % 2%

    % first (T,P)-dependent measurement for pH
    obs(i).tp(1).T    = 25 ; % degC
    obs(i).tp(1).uT   = 0.05 ; % from cruise report
    obs(i).tp(1).P    = 0.0 ;  % NOT in situ
    obs(i).tp(1).uP   = 0.07 ;
    obs(i).tp(1).ph   = in(9,i); % total scale
    obs(i).tp(1).uph  = 0.010; 

    % second (T,P)-dependent measurement for pCO2
    obs(i).tp(2).T     = 20 ; % degC
    obs(i).tp(2).uT    = 0.03 ; % from cruise report
    obs(i).tp(2).P     = 0.0 ; % dbar (surface pressure for pco2)
    obs(i).tp(2).uP    = 0.07 ;
    obs(i).tp(2).pco2  = in(10,i); % (µatm)
    obs(i).tp(2).upco2 = in(10,i)*0.01; % 1%

    % third (T,P)-dependent measurement for CO32-T
    obs(i).tp(3).T    = 25 ; % degC
    obs(i).tp(3).uT   = 0.05 ; % from cruise report
    obs(i).tp(3).P    = 0.0 ;  % NOT in situ
    obs(i).tp(3).uP   = 0.07 ;
    obs(i).tp(3).co3  = in(11,i); % (µmol/kg)
    obs(i).tp(3).uco3 = in(11,i)*0.02;  % 2% from Jon Sharp

    % % fourth (T,P)-dependent measurement IN SITU
    % obs(i).tp(4).T  = in(2,i); % deg C, CTD temp
    % obs(i).tp(4).uT = 0.02; % ±0.02 degC
    % obs(i).tp(4).P  = in(3,i); % dbar
    % obs(i).tp(4).uP = 0.63; % (max) ± 0.63 dbar
end

obs_backup = obs;

% run through all nD datapoints 26 times -> 26 combinations

%% Q5: All five input
% CT AT pH pCO2 CO3 (Q5) (fid26)
obs = obs_backup;
[est,~,~,~,~] = QUODcarb(obs,opt);
est26 = est;
save output_mat_files/est26.mat est26; clear est; clear obs;

%% Q2: Input pairs
% TC TA (Q2) (fid01)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(1).ph = nan;   obs(i).tp(1).uph = nan;
    obs(i).tp(2).pco2 = nan; obs(i).tp(2).upco2 = nan;
    obs(i).tp(3).co3 = nan;  obs(i).tp(3).uco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est01  = est;
save output_mat_files/est01.mat est01; clear est; clear obs;

% TC ph (Q2) (fid02)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;         obs(i).uTA = nan;
    obs(i).tp(2).pco2 = nan; obs(i).tp(2).upco2 = nan; % tp(3)
    obs(i).tp(3).co3 = nan;  obs(i).tp(3).uco3 = nan; % tp(2)
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est02   = est;
save output_mat_files/est02.mat est02; clear est; clear obs;

% TC pCO2 (Q2)(fid03)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;         obs(i).uTA = nan;
    obs(i).tp(1).ph = nan;   obs(i).tp(1).uph = nan;
    obs(i).tp(3).co3 = nan;  obs(i).tp(3).uco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est03   = est;
save output_mat_files/est03.mat est03; clear est; clear obs;

% % TC CO3 (Q2)(fid04)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;         obs(i).uTA = nan;
    obs(i).tp(1).ph = nan;   obs(i).tp(1).uph = nan; 
    obs(i).tp(2).pco2 = nan; obs(i).tp(2).upco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est04   = est;
save output_mat_files/est04.mat est04; clear est; clear obs;

% TA ph (Q2) (fid05)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;         obs(i).uTC = nan;
    obs(i).tp(2).pco2 = nan; obs(i).tp(2).upco2 = nan;
    obs(i).tp(3).co3 = nan;  obs(i).tp(3).uco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est05   = est;
save output_mat_files/est05.mat est05; clear est; clear obs;

% TA pCO2 (Q2)(fid06)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;        obs(i).uTC = nan;
    obs(i).tp(1).ph = nan;  obs(i).tp(1).uph = nan;
    obs(i).tp(3).co3 = nan; obs(i).tp(3).uco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est06   = est;
save output_mat_files/est06.mat est06; clear est; clear obs;

% TA CO3 (Q2)(fid07)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;          obs(i).uTC = nan;
    obs(i).tp(1).ph = nan;    obs(i).tp(1).uph = nan;
    obs(i).tp(2).pco2 = nan;  obs(i).tp(2).upco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est07   = est;
save output_mat_files/est07.mat est07; clear est; clear obs;

% pH CO3 (Q2) (fid08)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;            obs(i).uTC = nan;
    obs(i).TA = nan;            obs(i).uTA = nan;
    obs(i).tp(2).pco2 = nan;    obs(i).tp(2).upco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est08   = est;
save output_mat_files/est08.mat est08; clear est; clear obs;

% pH pCO2 (Q2) (fid09)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;        obs(i).uTC = nan;
    obs(i).TA = nan;        obs(i).uTA = nan;
    obs(i).tp(3).co3 = nan; obs(i).tp(3).uco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est09 = est;
save output_mat_files/est09.mat est09; clear est; clear obs;

% pCO2 CO3 (Q2) (fid10)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;        obs(i).uTC = nan;
    obs(i).TA = nan;        obs(i).uTA = nan;
    obs(i).tp(1).ph = nan;  obs(i).tp(1).uph = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est10   = est;
save output_mat_files/est10.mat est10; clear est; clear obs;


%% Q4: Four input

% TC TA pH pco2 (Q4) (fid21)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(3).co3 = nan; obs(i).tp(3).uco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est21   = est;
save output_mat_files/est21.mat est21; clear est; clear obs;

% % TC TA pH co3 (Q4) (fid22)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(2).pco2 = nan; obs(i).tp(2).upco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est22   = est;
save output_mat_files/est22.mat est22; clear est; clear obs;

% TC TA pco2 co3 (Q4) (fid23)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(1).ph = nan; obs(i).tp(1).uph = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est23   = est;
save output_mat_files/est23.mat est23; clear est; clear obs;

% TC pH pco2 co3 (Q4) (fid24)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan; obs(i).uTA = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est24   = est;
save output_mat_files/est24.mat est24; clear est; clear obs;

% TA ph pco2 co3 (Q4) (fid25)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan; obs(i).uTC = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est25   = est;
save output_mat_files/est25.mat est25; clear est; clear obs;


%% Q3: Three input

% % TC TA pH (Q3) (fid11)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(2).pco2 = nan;  obs(i).tp(2).upco2 = nan;
    obs(i).tp(3).co3 = nan;   obs(i).tp(3).uco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt); 
est11   = est;
save output_mat_files/est11.mat est11; clear est; clear obs;

% TC TA pco2 (Q3) (fid12)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(1).ph = nan;   obs(i).tp(1).uph = nan;
    obs(i).tp(3).co3 = nan;  obs(i).tp(3).uco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est12   = est;
save output_mat_files/est12.mat est12; clear est; clear obs;

% TC TA co3 (Q3) (fid13)
obs = obs_backup;
for i = 1:nD
    obs(i).tp(1).ph = nan;   obs(i).tp(1).uph = nan;
    obs(i).tp(2).pco2 = nan; obs(i).tp(2).upco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est13   = est;
save output_mat_files/est13.mat est13; clear est; clear obs;

% TC pH pco2 (Q3) (fid14)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;          obs(i).uTA = nan;
    obs(i).tp(3).co3 = nan;   obs(i).tp(3).uco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est14   = est;
save output_mat_files/est14.mat est14; clear est; clear obs;

% TC pH co3 (Q3) (fid15)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;          obs(i).uTA = nan;
    obs(i).tp(2).pco2 = nan;  obs(i).tp(2).upco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est15   = est;
save output_mat_files/est15.mat est15; clear est; clear obs;

% TC pco2 co3 (Q3) (fid16)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;          obs(i).uTA = nan;
    obs(i).tp(1).ph = nan;    obs(i).tp(1).uph = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est16   = est;
save output_mat_files/est16.mat est16; clear est; clear obs;

% TA pH pco2 (Q3) (fid17)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;          obs(i).uTC = nan;
    obs(i).tp(3).co3 = nan;   obs(i).tp(3).uco3 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est17   = est;
save output_mat_files/est17.mat est17; clear est; clear obs;

% TA pH co3 (Q3) (fid18)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;          obs(i).uTC = nan;
    obs(i).tp(2).pco2 = nan;  obs(i).tp(2).upco2 = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est18 = est;
save output_mat_files/est18.mat est18; clear est; clear obs;

% % TA pco2 co3 (Q3) (fid19)
obs = obs_backup;
for i = 1:nD
    obs(i).TC = nan;          obs(i).uTC = nan;
    obs(i).tp(1).ph = nan;    obs(i).tp(1).uph = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est19   = est;
save output_mat_files/est19.mat est19; clear est; clear obs;

% % pH pco2 co3 (Q3) (fid20)
obs = obs_backup;
for i = 1:nD
    obs(i).TA = nan;    obs(i).uTA = nan;
    obs(i).TC = nan;    obs(i).uTC = nan;
end
[est,~,~,~,~] = QUODcarb(obs,opt);
est20   = est;
save output_mat_files/est20.mat est20; clear est; clear obs;
















