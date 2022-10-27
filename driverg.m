% driver for GOMECC data

% take GOMECC data and put into QUODcarbV2

LOG10 = log(10);
p = @(x) -log10(x); % use p as in "pH" to denote the negative log base 10 function
q = @(x) 10.^(-x);  % use q, i.e., a backward p
w = @(x,e) p(1+e./x).^(-2); % convert x+/-e into precision for p(x)

load datag.mat; % data_gomecc
[in] = datag;
nD = length(in);

%
% convert all data to -log10(data) and compute the precision of the resulting data
%
for i = 1:nD

    % measurements that are independent of (T,P)
    obs.TC    = in(i,5); % (µmol/kg)
    obs.eTC   = 2.01 ; % TC error std ± 2.01µmol/kg
    obs.TA    = in(i,7); % (µmol/kg)
    obs.eTA   = 1.78 ; % TA error std ±1.78µmol/kg   
%     for j = 1187:1310
%         wpdic(j) = w(in(i,5),(2.24*1e-6)); % DIC error std ±2.24µmol/kg for UW measurements
%     end
    obs.sal   = in(i,3);
    obs.esal  = 0.002 ;

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
    obs.m(4).T   = 37 ; % degC
    obs.m(4).eT  = 0.02 ; % not stated, assumed 0.02
    obs.m(4).P   = in(i,1); % dbar (ctd pressure for pco2 I think?)
    obs.m(4).eP  = 0.63 ;
    % all nutrients measured @ 37degC
    obs.m(4).TP   = in(i,24); % total phosphate, µmol/kg
    obs.m(4).eTP  = in(i,24)*0.001; % 0.1% meas uncertainty
    obs.m(4).TSi  = in(i,18); % total Si, µmol/kg
    obs.m(4).eTSi = in(i,18)*0.001;
    % no2 and no3 not included yet
    %obs.m(4).no2  = p(in(i,20)*1e-6); % nitrite NO2^-, molg/kg
    %obs.m(4).wno2 = w(in(i,20),(in(i,20)*0.001));
    %obs.m(4).no3  = p(in(i,22)*1e-6); % nitrate NO3^2-, mol/kg
    %obs.m(4).wno3 = w(in(i,22),(in(i,22)*0.001));   

    % zero silicate is very unlikely; reset minimum to 1 nanomolar
    if (q(obs.m(4).TSi)==0)
        obs.m(4).TSi  = 1e-3 ; % µmol/kg-SW
        obs.m(4).eTSi = 1e-3;    
    end
    % zero po4 is very unlikely; reset minimum to 1 nanomolar
    if (q(obs.m(4).TP)==0)
        obs.m(4).TP  = 1e-3 ; % µmol/kg-SW
        obs.m(4).eTP = 1e-3 ; % µmol/kg
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
        sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride','phosphate','silicate'};
        %sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride', ...
        % 'phosphate','silicate','ammonia'};
        %sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride', ...
        % 'phosphate','silicate','ammonia','sulfide'};
                        
        sys = mksysV2(obs,sys);

        fid2 = fopen('gomecc_TA_TC_v2.csv','w');
        fid3 = fopen('gomecc_TA_ph_v2.csv','w');
        fid4 = fopen('gomecc_TC_ph_v2.csv','w');
        fid5 = fopen('gomecc_Q5_v2.csv','w');
        fid6 = fopen('compare_TA_TC.csv','w');
        fid7 = fopen('compare_TA_ph.csv','w');
        fid8 = fopen('compare_TC_ph.csv','w');
        fid9 = fopen('compare_Q5_cp.csv','w');
% 
        fprintf(fid6,'%s, %s, %s, %s, %s ', 'cTA', 'eTA', 'qTA', 'eTA' );
        fprintf(fid6,'%s, %s, %s, %s, %s ', 'cTC', 'eTC', 'qTC', 'eTC' );
        fprintf(fid6,'%s, %s, %s, %s, %s ', 'cph', 'eph', 'qph', 'eph' );
        fprintf(fid6,'%s, %s, %s, %s, %s ', 'cpco2', 'epco2', 'qpco2', 'epco2' );
        fprintf(fid6,'%s, %s, %s, %s, %s ', 'cpK0', 'epK0', 'qpK0', 'epK0' );
        fprintf(fid6,'%s, %s, %s, %s, %s ', 'cpK1', 'epK1', 'qpK1', 'epK1' );
        fprintf(fid6,'%s, %s, %s, %s, %s ', 'cpK2', 'epK2', 'qpK2', 'epK2' );
        fprintf(fid6,'%s, %s, %s, %s, %s ', 'cpKw', 'epKw', 'qpKw', 'epKw' );
        fprintf(fid6,'%s, %s, %s, %s, %s ', 'cpKb', 'epKb', 'qpKb', 'epKb' );
        fprintf(fid6,'%s, %s, %s, %s, %s ', 'cpKs', 'qpKs', 'cpKf', 'qpKf' );
        fprintf(fid6,'%s, %s, %s, %s, %s ', 'cpK1p', 'qpK1p', 'cpK2p', 'qpK2p' );
        fprintf(fid6,'%s, %s, %s, %s, %s ', 'cpK3p', 'qpK3p', 'cpKsi', 'qpKsi' );
        fprintf(fid6, '\n');

        fprintf(fid7,'%s, %s, %s, %s, %s ', 'cTA', 'eTA', 'qTA', 'eTA' );
        fprintf(fid7,'%s, %s, %s, %s, %s ', 'cTC', 'eTC', 'qTC', 'eTC' );
        fprintf(fid7,'%s, %s, %s, %s, %s ', 'cph', 'eph', 'qph', 'eph' );
        fprintf(fid7,'%s, %s, %s, %s, %s ', 'cpco2', 'epco2', 'qpco2', 'epco2' );
        fprintf(fid7,'%s, %s, %s, %s, %s ', 'cpK0', 'epK0', 'qpK0', 'epK0' );
        fprintf(fid7,'%s, %s, %s, %s, %s ', 'cpK1', 'epK1', 'qpK1', 'epK1' );
        fprintf(fid7,'%s, %s, %s, %s, %s ', 'cpK2', 'epK2', 'qpK2', 'epK2' );
        fprintf(fid7,'%s, %s, %s, %s, %s ', 'cpKw', 'epKw', 'qpKw', 'epKw' );
        fprintf(fid7,'%s, %s, %s, %s, %s ', 'cpKb', 'epKb', 'qpKb', 'epKb' );
        fprintf(fid7,'%s, %s, %s, %s, %s ', 'cpKs', 'qpKs', 'cpKf', 'qpKf' );
        fprintf(fid7,'%s, %s, %s, %s, %s ', 'cpK1p', 'qpK1p', 'cpK2p', 'qpK2p' );
        fprintf(fid7,'%s, %s, %s, %s, %s ', 'cpK3p', 'qpK3p', 'cpKsi', 'qpKsi' );
        fprintf(fid7, '\n');

        fprintf(fid8,'%s, %s, %s, %s, %s ', 'cTA', 'eTA', 'qTA', 'eTA' );
        fprintf(fid8,'%s, %s, %s, %s, %s ', 'cTC', 'eTC', 'qTC', 'eTC' );
        fprintf(fid8,'%s, %s, %s, %s, %s ', 'cph', 'eph', 'qph', 'eph' );
        fprintf(fid8,'%s, %s, %s, %s, %s ', 'cpco2', 'epco2', 'qpco2', 'epco2' );
        fprintf(fid8,'%s, %s, %s, %s, %s ', 'cpK0', 'epK0', 'qpK0', 'epK0' );
        fprintf(fid8,'%s, %s, %s, %s, %s ', 'cpK1', 'epK1', 'qpK1', 'epK1' );
        fprintf(fid8,'%s, %s, %s, %s, %s ', 'cpK2', 'epK2', 'qpK2', 'epK2' );
        fprintf(fid8,'%s, %s, %s, %s, %s ', 'cpKw', 'epKw', 'qpKw', 'epKw' );
        fprintf(fid8,'%s, %s, %s, %s, %s ', 'cpKb', 'epKb', 'qpKb', 'epKb' );
        fprintf(fid8,'%s, %s, %s, %s, %s ', 'cpKs', 'qpKs', 'cpKf', 'qpKf' );
        fprintf(fid8,'%s, %s, %s, %s, %s ', 'cpK1p', 'qpK1p', 'cpK2p', 'qpK2p' );
        fprintf(fid8,'%s, %s, %s, %s, %s ', 'cpK3p', 'qpK3p', 'cpKsi', 'qpKsi' );
        fprintf(fid8, '\n');

        fprintf(fid9,'%s, %s, %s, %s, %s ', 'cTA', 'eTA', 'qTA', 'eTA' );
        fprintf(fid9,'%s, %s, %s, %s, %s ', 'cTC', 'eTC', 'qTC', 'eTC' );
        fprintf(fid9,'%s, %s, %s, %s, %s ', 'cph', 'eph', 'qph', 'eph' );
        fprintf(fid9,'%s, %s, %s, %s, %s ', 'cpco2', 'epco2', 'qpco2', 'epco2' );
        fprintf(fid9,'%s, %s, %s, %s, %s ', 'cpK0', 'epK0', 'qpK0', 'epK0' );
        fprintf(fid9,'%s, %s, %s, %s, %s ', 'cpK1', 'epK1', 'qpK1', 'epK1' );
        fprintf(fid9,'%s, %s, %s, %s, %s ', 'cpK2', 'epK2', 'qpK2', 'epK2' );
        fprintf(fid9,'%s, %s, %s, %s, %s ', 'cpKw', 'epKw', 'qpKw', 'epKw' );
        fprintf(fid9,'%s, %s, %s, %s, %s ', 'cpKb', 'epKb', 'qpKb', 'epKb' );
        fprintf(fid9,'%s, %s, %s, %s, %s ', 'cpKs', 'qpKs', 'cpKf', 'qpKf' );
        fprintf(fid9,'%s, %s, %s, %s, %s ', 'cpK1p', 'qpK1p', 'cpK2p', 'qpK2p' );
        fprintf(fid9,'%s, %s, %s, %s, %s ', 'cpK3p', 'qpK3p', 'cpKsi', 'qpKsi' );
        fprintf(fid9, '\n');
        %fid6 = 'compare.xls'; fopen(fid6);
%         writematrix(A,fid6,'Sheet',1);
%         writematrix(A,fid6,'Sheet',2);
%         writematrix(A,fid6,'Sheet',3);
%         writematrix(A,fid6,'Sheet',4);

    end

    %
    % Call QUODcarb for all the data points
    %

    obs_backup = obs;

    % TA TC 
    obs.m(1).pH = nan; obs.m(1).wpH = nan;
    obs.m(2).pco2 = nan; obs.m(2).wpco2 = nan;
    [est,obs,iflag] = QUODcarbV2(obs,sys); 
    if i == 1
        % make headers
        PrintCSVv2(sys,est,fid2);
    end
    PrintCSVv2(sys,est,obs,iflag,fid2);
%     keyboard
    [A] = compare(obs,1,est,fid6);
    %writematrix(A,fid6,'Sheet',1,'WriteMode','append'); %48
    %keyboard
    % TA ph
    obs = obs_backup;
    obs.TC = nan; obs.wTC = nan;
    obs.m(2).pco2 = nan; obs.m(2).wpco2 = nan;
    [est,obs,iflag] = QUODcarbV2(obs,sys);
    if i == 1
        % make headers
        PrintCSVv2(sys,est,fid3);
    end
    PrintCSVv2(sys,est,obs,iflag,fid3);
    [A] = compare(obs,3,est,fid7);
    %writematrix(A,fid6,'Sheet',2,'WriteMode','append');

    % TC ph
    obs = obs_backup;
    obs.TA = nan; obs.wTA = nan;
    obs.m(2).pco2 = nan; obs.m(2).wpco2 = nan;
    [est,obs,iflag] = QUODcarbV2(obs,sys);
    if i == 1
        % make headers
        PrintCSVv2(sys,est,fid4);
    end
    PrintCSVv2(sys,est,obs,iflag,fid4);
    [A] = compare(obs,2,est,fid8);
    %writematrix(A,fid6,'Sheet',3,'WriteMode','append');

    % TA TC ph pCO2
    obs = obs_backup;
    [est,obs,iflag] = QUODcarbV2(obs,sys);
    if i == 1
        % make headers
        PrintCSVv2(sys,est,fid5);
    end
    PrintCSVv2(sys,est,obs,iflag,fid5);
    [A] = compare(obs,2,est,fid9);
    %writematrix(A,fid6,'Sheet',4,'WriteMode','append');
       
    %keyboard
end

fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);
fclose(fid7);
fclose(fid8);
fclose(fid9);



