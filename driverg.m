% driver for GOMECC data

% take GOMECC data and put into QUODcarbV2

LOG10 = log(10);
p = @(x) -log10(x); % use p as in "pH" to denote the negative log base 10 function
q = @(x) 10.^(-x);  % use q, i.e., a backward p
w = @(x,e) p(1+e./x).^(-2); % convert x+/-e into precision for p(x)

load datag.mat; % data_gomecc
[in] = datag;

%
% convert all data to -log10(data) and compute the precision of the resulting data
%
for i = 1  %:1310

    % measurements that are independent of (T,P)
    obs.TC(i)    = p(in(i,5)*1e-6); % p(mol/kg)
    obs.wTC(i)   = w((in(i,5)*1e-6),(2.01*1e-6)); % TC error std ± 2.01µmol/kg
    obs.TA(i)    = p(in(i,7)*1e-6); % p(mol/kg)
    obs.wTA(i)   = w((in(i,7)*1e-6),(1.78*1e-6)); % TA error std ±1.78µmol/kg   
%     for j = 1187:1310
%         wpdic(j) = w(in(i,5),(2.24*1e-6)); % DIC error std ±2.24µmol/kg for UW measurements
%     end
    obs.sal(i)   = in(i,3);
    obs.wsal(i)  = (0.002)^(-2);

    % first (T,P)-dependent measurement
    obs.m(1).T(i)   = in(i,2); % degC, CTD temp
    obs.m(1).wT(i)  = (0.02)^(-2); % precision of temp % std ±0.02 degC, (std ±0.0007 degC below 1000m, not included rn)
    obs.m(1).P(i)   = in(i,1); % dbar
    obs.m(1).wP(i)  = (0.63)^(-2); % (max) std ±0.63 dbar

    % second (T,P)-dependent measurement
    obs.m(2).T(i)    = 25; % degC
    obs.m(2).wT(i)   = (0.05)^(-2); % from cruise report
    obs.m(2).P(i)    = 10.1325; % 1atm = 10dbar
    obs.m(2).wP(i)   = (0.63)^(-2);
    obs.m(2).ph(i)   = in(i,9); % total scale
    obs.m(2).wph(i)  = (0.001)^(-2);
    obs.m(2).co3(i)  = p(in(i,15)*1e-6); % p(mol/kg)
    obs.m(2).wco3(i) = w(in(i,15)*1e-6,2e-6);  % std ±2µmol/kg

    % third (T,P)-dependent measurement
    obs.m(3).T(i)   = 20; %degC
    obs.m(3).wT(i)  = (0.03)^(-2); % from cruise report
    obs.m(3).P(i)   = in(i,1); % dbar (ctd pressure for pco2 I think?)
    obs.m(3).wP(i)  = (0.63)^(-2);
    obs.m(3).pco2(i)  = p(in(i,12)*1e-6); % p(atm)
    obs.m(3).wpco2(i) = w(in(i,12),(in(i,12)*0.0021)); % 0.21% relative std error (avg)

    % fourth (T,P)-dependent measurement
    obs.m(4).T(i)   = 37; % degC
    obs.m(4).wT(i) = (0.02)^(-2); % not stated, assumed 0.02
    % all nutrients measured @ 37degC
    obs.m(4).TP(i)   = p(in(i,24)*1e-6); % total phosphate, mol/kg
    obs.m(4).wTP(i)  = w(in(i,24),(in(i,24)*0.001)); % 0.1% meas uncertainty
    obs.m(4).TSi(i)  = p(in(i,18)*1e-6); % total Si, mol/kg
    obs.m(4).wTSi(i) = w(in(i,18),(in(i,18)*0.001));
    % no2 and no3 not included yet
    %obs.m(4).no2  = p(in(i,20)*1e-6); % nitrite NO2^-, molg/kg
    %obs.m(4).wno2 = w(in(i,20),(in(i,20)*0.001));
    %obs.m(4).no3  = p(in(i,22)*1e-6); % nitrate NO3^2-, mol/kg
    %obs.m(4).wno3 = w(in(i,22),(in(i,22)*0.001));   

    % zero silicate is very unlikely; reset minimum to 0.1 micromolar
    if (q(obs.m(4).TSi(i))==0)
        obs.m(4).TSi(i)  = p(1e-7);
        obs.m(4).wTSi(i) = w( 1e-7, 1e-7);    
    end
    % zero po4 is very unlikely; reset minimum to 1 nanomolar
    if (q(obs.m(4).TP(i))==0)
        obs.m(4).TP(i) = p(1e-9); %
        obs.m(4).wTP(i) = w( 1e-9 , 1e-9 );
    end

end


                        
%
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
%PrintTable(sys);
%n = length(sys.variables);

%
% Call QUODcarb for all the data points
%
keyboard
fid2 = fopen('gomecc_TA_TC_v2.csv','w');    PrintCSV(sys,fid2);
fid3 = fopen('gomecc_TA_ph_v2.csv','w');    PrintCSV(sys,fid3);
fid4 = fopen('gomecc_TC_ph_v2.csv','w');    PrintCSV(sys,fid4);
fid5 = fopen('gomecc_TA_TC_ph_v2.csv','w'); PrintCSV(sys,fid5);

for i = 1  %1310

    yobs_backup = yobs;
    wobs_backup = wobs;

    % TA TC 
    [est,obs,iflag] = QUODcarbV2(yobs,sys); 
    %w = sigy.^(-2);
    PrintCSV(y,w,yobs_out,wobs_out,sys,iflag,fid2);


    % TA ph
    yobs = yobs_backup; wobs = wobs_backup;
    yobs(sys.iTC) = nan; wobs(sys.iTC) = nan;
    [y,sigy,yobs_out,wobs_out,iflag] = QUODcarb(yobs,wobs,temp(i),sal(i),pres(i),sys); 
    w = sigy.^(-2);
    PrintCSV(y,w,yobs_out,wobs_out,sys,iflag,fid3);

    % TC ph
    yobs = yobs_backup; wobs = wobs_backup;
    yobs(sys.iTA) = nan; wobs(sys.iTA) = nan;
    [y,sigy,yobs_out,wobs_out,iflag] = QUODcarb(yobs,wobs,temp(i),sal(i),pres(i),sys); 
    w = sigy.^(-2);
    PrintCSV(y,w,yobs_out,wobs_out,sys,iflag,fid4);
   
    % TA TC ph
    yobs = yobs_backup; wobs = wobs_backup;
    [y,sigy,yobs_out,wobs_out,iflag] = QUODcarb(yobs,wobs,temp(i),sal(i),pres(i),sys); 
    w = sigy.^(-2); 
    PrintCSV(y,w,yobs_out,wobs_out,sys,iflag,fid5);

    % try with other combos? also co3 and pco2
end
fclose(fid2);
fclose(fid3);
fclose(fid4);
fclose(fid5);




