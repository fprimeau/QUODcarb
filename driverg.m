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
    obs.TC    = (in(i,5)*1e-6); % p(mol/kg)
    obs.wTC   = w((in(i,5)*1e-6),(2.01*1e-6)); % TC error std ± 2.01µmol/kg
    obs.TA    = (in(i,7)*1e-6); % p(mol/kg)
    obs.wTA   = w((in(i,7)*1e-6),(1.78*1e-6)); % TA error std ±1.78µmol/kg   
%     for j = 1187:1310
%         wpdic(j) = w(in(i,5),(2.24*1e-6)); % DIC error std ±2.24µmol/kg for UW measurements
%     end
    obs.sal   = in(i,3);
    obs.wsal  = (0.002)^(-2);

    % first (T,P)-dependent measurement
    obs.m(1).T   = in(i,2); % degC, CTD temp
    obs.m(1).wT  = (0.02)^(-2); % precision of temp % std ±0.02 degC, (std ±0.0007 degC below 1000m, not included rn)
    obs.m(1).P   = in(i,1); % dbar
    obs.m(1).wP  = (0.63)^(-2); % (max) std ±0.63 dbar

    % second (T,P)-dependent measurement
    obs.m(2).T    = 25; % degC
    obs.m(2).wT   = (0.05)^(-2); % from cruise report
    obs.m(2).P    = in(i,1); 
    obs.m(2).wP   = (0.63)^(-2);
    obs.m(2).ph   = in(i,9); % total scale
    obs.m(2).wph  = (0.001)^(-2);
    obs.m(2).co3  = (in(i,15)*1e-6); % p(mol/kg)
    obs.m(2).wco3 = w(in(i,15)*1e-6,2e-6);  % std ±2µmol/kg

    % third (T,P)-dependent measurement
    obs.m(3).T   = 20; %degC
    obs.m(3).wT  = (0.03)^(-2); % from cruise report
    obs.m(3).P   = 0; % dbar (surface pressure for pco2)
    obs.m(3).wP  = (0.63)^(-2);
    obs.m(3).pco2  = (in(i,12)*1e-6); % p(atm)
    obs.m(3).wpco2 = w(in(i,12),(in(i,12)*0.0021)); % 0.21% relative std error (avg)

    % fourth (T,P)-dependent measurement
    obs.m(4).T   = 37; % degC
    obs.m(4).wT  = (0.02)^(-2); % not stated, assumed 0.02
    obs.m(4).P   = in(i,1); % dbar (ctd pressure for pco2 I think?)
    obs.m(4).wP  = (0.63)^(-2);
    % all nutrients measured @ 37degC
    obs.m(4).TP   = (in(i,24)*1e-6); % total phosphate, mol/kg
    obs.m(4).wTP  = w(in(i,24),(in(i,24)*0.001)); % 0.1% meas uncertainty
    obs.m(4).TSi  = (in(i,18)*1e-6); % total Si, mol/kg
    obs.m(4).wTSi = w(in(i,18),(in(i,18)*0.001));
    % no2 and no3 not included yet
    %obs.m(4).no2  = p(in(i,20)*1e-6); % nitrite NO2^-, molg/kg
    %obs.m(4).wno2 = w(in(i,20),(in(i,20)*0.001));
    %obs.m(4).no3  = p(in(i,22)*1e-6); % nitrate NO3^2-, mol/kg
    %obs.m(4).wno3 = w(in(i,22),(in(i,22)*0.001));   

    % zero silicate is very unlikely; reset minimum to 0.1 micromolar
    if (q(obs.m(4).TSi)==0)
        obs.m(4).TSi  = (1e-7);
        obs.m(4).wTSi = w( 1e-7, 1e-7);    
    end
    % zero po4 is very unlikely; reset minimum to 1 nanomolar
    if (q(obs.m(4).TP)==0)
        obs.m(4).TP  = (1e-9); %
        obs.m(4).wTP = w( 1e-9 , 1e-9 );
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

        fid2 = fopen('gomecc_TA_TC_v2.csv','w');    %PrintCSV(sys,est,fid2);
        fid3 = fopen('gomecc_TA_ph_v2.csv','w');    %PrintCSV(sys,est,fid3);
        fid4 = fopen('gomecc_TC_ph_v2.csv','w');    %PrintCSV(sys,est,fid4);
        fid5 = fopen('gomecc_Q5_v2.csv','w');       %PrintCSV(sys,est,fid5);
        fid6 = 'compare.xls'; fopen(fid6);
        A = zeros(1,48);
        writematrix(A,fid6);

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
    keyboard
    [A] = compare(obs,1,est);
    writematrix(A,fid6,'Sheet',1,'WriteMode','append'); %48

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
    [A] = compare(obs,3,est);
    writematrix(A,fid6,'Sheet',2,'WriteMode','append');

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
    [A] = compare(obs,2,est);
    writematrix(A,fid6,'Sheet',3,'WriteMode','append');

    % TA TC ph pCO2
    obs = obs_backup;
    [est,obs,iflag] = QUODcarbV2(obs,sys);
    if i == 1
        % make headers
        PrintCSVv2(sys,est,fid5);
    end
    PrintCSVv2(sys,est,obs,iflag,fid5);
    [A] = compare(obs,2,est);
    writematrix(A,fid6,'Sheet',4,'WriteMode','append');
       
    
end

fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);
fclose(fid6);




