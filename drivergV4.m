% drivergV4 changed to SWS scale like CO2SYS does
% implementing CO2SYS v3 changes

% driver for GOMECC data
format long
% take GOMECC data and put into QUODcarbV4

LOG10 = log(10);
p = @(x) -log10(x); % use p as in "pH" to denote the negative log base 10 function
q = @(x) 10.^(-x);  % use q, i.e., a backward p
w = @(x,e) p(1+e./x).^(-2); % convert x+/-e into precision for p(x)

load datag.mat; % data_gomecc
[in] = datag;
% a = in(776,24); in(776,24) = -a; % for some reason this is negative and throws the code
% b = in(809,24); in(809,24) = -b;
% c = in(811,24); in(811,24) = -c;
% d = in(812,24); in(812,24) = -d; % I think I save the new datag with this
nD = length(in);

%
% convert all data to -log10(data) and compute the precision of the resulting data
%
for i = 1:nD

    % choose equilibrium constant formulation
    obs.cK1K2 = 4; % 'choose' K1K2
    obs.cKSO4 = 1; % 'choose' KSO4
    obs.cKF   = 2; % 'choose' KF
    obs.cTB   = 2; % 'choose' TB
    obs.phscale = 1 ; % 1 = tot, 2 = sws, 3 = free, 4 = NBS

    % measurements that are independent of (T,P)
    obs.TC    = in(i,5); % (µmol/kg)
    obs.eTC   = 2.01 ; % TC error std ± 2.01µmol/kg
    if i >= 1187
         obs.eTC = 2.24 ; % DIC error std ±2.24µmol/kg for UW measurements
    end
    obs.TA    = in(i,7); % (µmol/kg)
    obs.eTA   = 1.78 ; % TA error std ±1.78µmol/kg   
    obs.sal   = in(i,3);
    obs.esal  = 0.002 ;
    % all nutrients measured @ 37degC, doesn't matter?
    obs.TP   = in(i,24); % total phosphate, µmol/kg
    obs.eTP  = in(i,24)*0.001; % 0.1% meas uncertainty
    obs.TSi  = in(i,18); % total Si, µmol/kg
    obs.eTSi = in(i,18)*0.001;

    % first (T,P)-dependent measurement
    obs.m(1).T   = in(i,2); % degC, CTD temp
    obs.m(1).eT  = 0.02 ; % precision of temp % std ±0.02 degC, (std ±0.0007 degC below 1000m, not included rn)
    obs.m(1).P   = in(i,1); % dbar
    obs.m(1).eP  = 0.63 ; % (max) std ±0.63 dbar

    % second (T,P)-dependent measurement
    obs.m(2).T    = 25 ; % degC
    obs.m(2).eT   = 0.05 ; % from cruise report
    obs.m(2).P    = in(i,1); 
    obs.m(2).eP   = 0.63 ;
    obs.m(2).ph   = in(i,9); % total scale
    obs.m(2).eph  = 0.001 ;
    obs.m(2).co3  = in(i,15); % (µmol/kg)
    obs.m(2).eco3 = 2.0 ;  % std ±2µmol/kg

    % third (T,P)-dependent measurement
    obs.m(3).T   = 20 ; %degC
    obs.m(3).eT  = 0.03 ; % from cruise report
    obs.m(3).P   = 0 ; % dbar (surface pressure for pco2)
    obs.m(3).eP  = 0.63 ;
    obs.m(3).pco2  = in(i,12); % (µatm)
    obs.m(3).epco2 = in(i,12)*0.0021; % 0.21% relative std error (avg)

    % fourth (T,P)-dependent measurement
    %obs.m(4).T   = 37 ; % degC
    %obs.m(4).eT  = 0.02 ; % not stated, assumed 0.02
    %obs.m(4).P   = in(i,1); % dbar (ctd pressure for pco2 I think?)
    %obs.m(4).eP  = 0.63 ;
    
    % no2 and no3 not included yet
    %obs.m(4).no2  = p(in(i,20)*1e-6); % nitrite NO2^-, molg/kg
    %obs.m(4).wno2 = w(in(i,20),(in(i,20)*0.001));
    %obs.m(4).no3  = p(in(i,22)*1e-6); % nitrate NO3^2-, mol/kg
    %obs.m(4).wno3 = w(in(i,22),(in(i,22)*0.001));   

    % zero silicate is very unlikely; reset minimum to 1 nanomolar
    if ((obs.TSi)==0)
        obs.TSi  = 1e-3 ; % µmol/kg-SW
        obs.eTSi = 1e-3;    
    end
    % zero po4 is very unlikely; reset minimum to 1 nanomolar
    if ((obs.TP)==0)
        obs.TP  = 1e-3 ; % µmol/kg-SW
        obs.eTP = 1e-3 ; % µmol/kg
    end

    if i == 1 

        % initialize the QUODcarb CO2-system solver
        %
        %sys.abr = {'air-sea','carbonate'};
        %sys.abr = {'air-sea','carbonate','water'};
        %sys.abr = {'air-sea','carbonate','water','borate'};
        %sys.abr = {'air-sea','carbonate','water','borate','sulfate'};
        %sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride'};
        %sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride','phosphate'};
        %sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride','phosphate','silicate'};
        %sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride', ...
        % 'phosphate','silicate','ammonia'};
        sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride', ...
         'phosphate','silicate','ammonia','sulfide','solubility'};
                        
        sys = mksysV3(obs,sys);

%         fid2 = fopen('gomecc_TA_TC_v4.csv','w');
%         fid6 = fopen('compare_TA_TC_v4.csv','w');
%         [A] = compare(fid6); % make headers
% 
%         fid3 = fopen('gomecc_TA_ph_v4.csv','w');
%         fid7 = fopen('compare_TA_ph_v4.csv','w');
%         [A] = compare(fid7);
% 
%         fid4 = fopen('gomecc_TC_ph_v4.csv','w');
%         fid8 = fopen('compare_TC_ph_v4.csv','w');
%         [A] = compare(fid8);
% 
%         fid5 = fopen('gomecc_Q5_v4.csv','w');
%         fid9 = fopen('compare_Q5_cp_v4.csv','w');
%         [A] = compare(fid9);
% 
%         fid10 = fopen('gomecc_pH_pco2_v4.csv','w');
%         fid110 = fopen('compare_pH_pco2_v4.csv','w');
%         [A] = compare(fid110);
% 
%         fid11 = fopen('gomecc_pco2_co3_v4.csv','w');
%         fid111 = fopen('compare_pco2_co3_v4.csv','w');
%         [A] = compare(fid111);
% 
%         fid12 = fopen('gomecc_TC_co3_v4.csv','w');
%         fid112 = fopen('compare_TC_co3_v4.csv','w');
%         [A] = compare(fid112);
% 
%         fid13 = fopen('gomecc_TC_pco2_v4.csv','w');
%         fid113 = fopen('compare_TC_pco2_v4.csv','w');
%         [A] = compare(fid113);
% 
%         fid14 = fopen('gomecc_TA_pco2_v4.csv','w');
%         fid114 = fopen('compare_TA_pco2_v4.csv','w');
%         [A] = compare(fid114);
% 
%         fid15 = fopen('gomecc_TA_co3_v4.csv','w');
%         fid115 = fopen('compare_TA_co3_v4.csv','w');
%         [A] = compare(fid115);
% 
%         fid16 = fopen('gomecc_pH_co3_v4.csv','w');
%         fid116 = fopen('compare_pH_co3_v4.csv','w');
%         [A] = compare(fid116);
% 
%         fid20 = fopen('gomecc_TC_TA_pH.csv','w');
%         fid21 = fopen('gomecc_TC_TA_pco2.csv','w');
%         fid22 = fopen('gomecc_TC_TA_co3.csv','w');
%         fid23 = fopen('gomecc_TC_pH_pco2.csv','w');
%         fid24 = fopen('gomecc_TC_pH_co3.csv','w');
%         fid25 = fopen('gomecc_TC_pco2_co3.csv','w');
%         fid26 = fopen('gomecc_TA_pH_pco2.csv','w');
%         fid27 = fopen('gomecc_TA_pH_co3.csv','w');
%         fid28 = fopen('gomecc_TA_pco2_co3.csv','w');
%         fid29 = fopen('gomecc_pH_pco2_co3.csv','w');
% 
%         fid30 = fopen('gomecc_TCTApHpco2.csv','w');
%         fid31 = fopen('gomecc_TCTApHco3.csv','w');
%         fid32 = fopen('gomecc_TCTApco2co3.csv','w');
%         fid33 = fopen('gomecc_TCpHpco2co3.csv','w');
%         fid34 = fopen('gomecc_TApHpco2co3.csv','w');

    end

    %
    % Call QUODcarb for all the data points
    %

    obs_backup = obs;

    % TA TC (Q2) (fid2)
    obs.m(2).ph = nan; obs.m(2).eph = nan;
    obs.m(3).pco2 = nan; obs.m(3).epco2 = nan;
    obs.m(2).co3 = nan; obs.m(2).eco3 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est02K04(i) = est;
%     if i == 1
%         % make headers
%         PrintCSVv3(sys,est,fid2);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid2);
%     [A] = compare(obs,1,est,fid6);
    clear obs
   
    % TA ph (Q2) (fid3)
    obs = obs_backup;
    obs.TC = nan; obs.eTC = nan;
    obs.m(3).pco2 = nan; obs.m(3).epco2 = nan;
    obs.m(2).co3 = nan; obs.m(2).eco3 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est03K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid3);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid3);
%     [A] = compare(obs,3,est,fid7);
    clear obs
    
    % TC ph (Q2) (fid4)
    obs = obs_backup;
    obs.TA = nan; obs.eTA = nan;
    obs.m(3).pco2 = nan; obs.m(3).epco2 = nan;
    obs.m(2).co3 = nan; obs.m(2).eco3 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est04K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid4);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid4);
%     [A] = compare(obs,2,est,fid8);
    clear obs

    % pH pCO2 (Q2) (fid10)
    obs = obs_backup;
    obs.TC = nan; obs.eTC = nan;
    obs.TA = nan; obs.eTA = nan;
    obs.m(2).co3 = nan; obs.m(2).eco3 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est10K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid10);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid10);
%     [A] = compare(obs,4,est,fid110);
    clear obs

    % pCO2 CO3 (Q2) (fid11)
    obs = obs_backup;
    obs.TC = nan; obs.eTC = nan;
    obs.TA = nan; obs.eTA = nan;
    obs.m(2).ph = nan; obs.m(2).eph = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est11K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid11);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid11);
%     [A] = compare(obs,5,est,fid111);
    clear obs

    % TC CO3 (Q2)(fid12)
    obs = obs_backup;
    obs.TA = nan; obs.eTA = nan;
    obs.m(3).pco2 = nan; obs.m(3).epco2 = nan;
    obs.m(2).ph = nan; obs.m(2).eph = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est12K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid12);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid12);
%     [A] = compare(obs,6,est,fid112);
    clear obs

    % TC pCO2 (Q2)(fid13)
    obs = obs_backup;
    obs.TA = nan; obs.eTA = nan;
    obs.m(2).ph = nan; obs.m(2).eph = nan;
    obs.m(2).co3 = nan; obs.m(2).eco3 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est13K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid13);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid13);
%     [A] = compare(obs,7,est,fid113); 
    clear obs

    % TA pCO2 (Q2)(fid14)
    obs = obs_backup;
    obs.TC = nan; obs.eTC = nan;
    obs.m(2).ph = nan; obs.m(2).eph = nan;
    obs.m(2).co3 = nan; obs.m(2).eco3 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est14K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid14);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid14);
%     [A] = compare(obs,8,est,fid114); 
    clear obs

    % TA CO3 (Q2)(fid15)
    obs = obs_backup;
    obs.TC = nan; obs.eTC = nan;
    obs.m(2).ph = nan; obs.m(2).eph = nan;
    obs.m(3).pco2 = nan; obs.m(3).epco2 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est15K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid15);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid15);
%     [A] = compare(obs,9,est,fid115); 
    clear obs

    % pH CO3 (Q2) (fid16)
    obs = obs_backup;
    obs.TC = nan; obs.eTC = nan;
    obs.TA = nan; obs.eTA = nan;
    obs.m(3).pco2 = nan; obs.m(3).epco2 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est16K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid16);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid16);
%     [A] = compare(obs,10,est,fid116); 
    clear obs
    
    % TC TA pH (Q3)
    obs = obs_backup;
    obs.m(2).co3 = nan; obs.m(2).eco3 = nan;
    obs.m(3).pco2 = nan; obs.m(3).epco2 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est20K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid20);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid20);
    clear obs

    % TC TA pco2 (Q3)
    obs = obs_backup;
    obs.m(2).co3 = nan; obs.m(2).eco3 = nan;
    obs.m(2).ph = nan; obs.m(2).eph = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est21K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid21);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid21);
    clear obs

    % TC TA co3 (Q3)
    obs = obs_backup;
    obs.m(3).pco2 = nan; obs.m(2).epco2 = nan;
    obs.m(2).ph = nan; obs.m(2).eph = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est22K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid22);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid22);
    clear obs

    % TC pH pco2 (Q3)
    obs = obs_backup;
    obs.TA = nan; obs.eTA = nan;
    obs.m(2).co3 = nan; obs.m(2).eco3 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est23K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid23);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid23);
    clear obs

    % TC pH co3 (Q3)
    obs = obs_backup;
    obs.TA = nan; obs.eTA = nan;
    obs.m(3).pco2 = nan; obs.m(3).epco2 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est24K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid24);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid24);
    clear obs

    % TC pco2 co3 (Q3)
    obs = obs_backup;
    obs.TA = nan; obs.eTA = nan;
    obs.m(2).ph = nan; obs.m(2).eph = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est25K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid25);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid25);
    clear obs

    % TA pH pco2 (Q3)
    obs = obs_backup;
    obs.TC = nan; obs.eTC = nan;
    obs.m(2).co3 = nan; obs.m(2).eco3 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est26K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid26);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid26);
    clear obs

    % TA pH co3 (Q3)
    obs = obs_backup;
    obs.TC = nan; obs.eTC = nan;
    obs.m(3).pco2 = nan; obs.m(3).epco2 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est27K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid27);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid27);
    clear obs

    % TA pco2 co3 (Q3)
    obs = obs_backup;
    obs.TC = nan; obs.eTC = nan;
    obs.m(2).ph = nan; obs.m(2).eph = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est28K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid28);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid28);
    clear obs

    % pH pco2 co3 (Q3)
    obs = obs_backup;
    obs.TA = nan; obs.eTA = nan;
    obs.TC = nan; obs.eTC = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est29K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid29);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid29);
    clear obs

    % TC TA pH pco2 (Q4)
    obs = obs_backup;
    obs.m(2).co3 = nan; obs.m(2).eco3 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est30K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid30);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid30);
    clear obs

    % TC TA pH co3 (Q4)
    obs = obs_backup;
    obs.m(3).pco2 = nan; obs.m(3).epco2 = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est31K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid31);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid31);
    clear obs

    % TC TA pco2 co3 (Q4)
    obs = obs_backup;
    obs.m(2).ph = nan; obs.m(2).eph = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est32K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid32);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid32);
    clear obs

    % TC pH pco2 co3 (Q4)
    obs = obs_backup;
    obs.TA = nan; obs.eTA = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est33K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid33);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid33);
    clear obs

    % TA pH pco2 co3 (Q4)
    obs = obs_backup;
    obs.TC = nan; obs.eTC = nan;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est34K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid34);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid34);
    clear obs
    
    % TA TC ph pCO2 co3 (Q5)
    obs = obs_backup;
    [est,obs,iflag] = QUODcarbV4(obs,sys);
    est05K04(i) = est;
%     if i == 1
%         PrintCSVv3(sys,est,fid5);
%     end
%     PrintCSVv3(sys,est,obs,iflag,fid5);
%     [A] = compare(obs,2,est,fid9);
    clear obs
       
end

save estK04/est02K04.mat est02K04;
save estK04/est03K04.mat est03K04;
save estK04/est04K04.mat est04K04;
save estK04/est05K04.mat est05K04;
save estK04/est10K04.mat est10K04;
save estK04/est11K04.mat est11K04;
save estK04/est12K04.mat est12K04;
save estK04/est13K04.mat est13K04;
save estK04/est14K04.mat est14K04;
save estK04/est15K04.mat est15K04;
save estK04/est16K04.mat est16K04;
save estK04/est20K04.mat est20K04;
save estK04/est21K04.mat est21K04;
save estK04/est22K04.mat est22K04;
save estK04/est23K04.mat est23K04;
save estK04/est24K04.mat est24K04;
save estK04/est25K04.mat est25K04;
save estK04/est26K04.mat est26K04;
save estK04/est27K04.mat est27K04;
save estK04/est28K04.mat est28K04;
save estK04/est29K04.mat est29K04;
save estK04/est30K04.mat est30K04;
save estK04/est31K04.mat est31K04;
save estK04/est32K04.mat est32K04;
save estK04/est33K04.mat est33K04;
save estK04/est34K04.mat est34K04;




% fclose(fid2);
% fclose(fid3);
% fclose(fid4);
% fclose(fid5);
% fclose(fid6);
% fclose(fid7);
% fclose(fid8);
% fclose(fid9);
% fclose(fid10);
% fclose(fid11);
% fclose(fid12);
% fclose(fid13);
% fclose(fid14);
% fclose(fid15);
% fclose(fid16);
% fclose(fid110);
% fclose(fid111);
% fclose(fid112);
% fclose(fid113);
% fclose(fid114);
% fclose(fid115);
% fclose(fid116);


