% driver for GOMECC data

% take GOMECC data and put into QUODcarbV2

LOG10 = log(10);
p = @(x) -log10(x); % use p as in "pH" to denote the negative log base 10 function
q = @(x) 10.^(-x);  % use q, i.e., a backward p
w = @(x,e) p(1+e./x).^(-2); % convert x+/-e into precision for p(x)

load datag.mat; % data_gomecc
[in] = datag;

palk  = zeros(1310,1); wpalk  = zeros(1310,1); 
pdic  = zeros(1310,1); wpdic  = zeros(1310,1); 
pH    = zeros(1310,1); wpH    = zeros(1310,1);
ppco2 = zeros(1310,1); wppco2 = zeros(1310,1); % pCO2
pco3  = zeros(1310,1); wpco3  = zeros(1310,1); % CO3^2-
psil  = zeros(1310,1); wpsil  = zeros(1310,1); % silicate
pno2  = zeros(1310,1); wpno2  = zeros(1310,1); % nitrite
pno3  = zeros(1310,1); wpno3  = zeros(1310,1); % nitrate
ppo4  = zeros(1310,1); wppo4  = zeros(1310,1); % phosphate
ptb   = zeros(1310,1); wptb   = zeros(1310,1);
ptf   = zeros(1310,1); wptf   = zeros(1310,1);
pts   = zeros(1310,1); wpts   = zeros(1310,1);
ptnh3 = zeros(1310,1); wptnh4 = zeros(1310,1); % ammonia
pth2s = zeros(1310,1); wpth2s = zeros(1310,1); % hydrogen sulfide

sal  = zeros(1310,1);  esal   = zeros(1310,1); wsal  = zeros(1310,1);
temp = zeros(1310,1);  etemp  = zeros(1310,1); wtemp = zeros(1310,1);
pres = zeros(1310,1);  epres  = zeros(1310,1); wpres = zeros(1310,1);

%
% convert all data to -log10(data) and compute the precision of the resulting data
%
for i = 1:1310
    
    palk(i) = p(in(i,7)*1e-6); % mol/kg
    wpalk(i) = w(in(i,7),(1.78*1e-6)); % ALK error std ±1.78µmol/kg

    pdic(i) = p(in(i,5)*1e-6); % mol/kg
    wpdic(i) = w(in(i,5),(2.01*1e-6)); % DIC error std ±2.01µmol/kg
    
    for j = 1187:1310
        wpdic(j) = w(in(i,5),(2.24*1e-6)); % DIC error std ±2.24µmol/kg for UW measurements
    end

    pH(i) = in(i,9); % total pH scale
    wpH(i)  = (0.001).^(-2); % pH error
    tpH = 25; % temp pH @ 25 degC

    ppco2(i) = p( in(i,12) * 1e-6 ); % convert µatm to atm
    wppco2(i) = w( in(i,12), (in(i,12)*0.0021) ); % 0.21% relative std error on avg
    tpco2 = 20; % temp pCO2 @ 20 degC

    pco3(i) = p( in(i,15) * 1e-6 ); % convert µmol/kg to mol/kg
    wpco3(i) = w( in(i,15), 2); % std ±2µmol/kg
    tpco3 = 25; % temp pCO3 @ 25 degC

    sal(i) = in(i,3); % salinity
    esal(i) = (0.001); % want to calculate this ourselves w/duplicate meas
    wsal(i) = (esal(i)).^(-2);

    temp(i) = in(i,2); % deg C
    etemp(i) = (0.02); % std ±0.02 degC, (std ±0.0007 degC below 1000m, not included rn)
    wtemp(i) = (etemp(i)).^(-2);

    pres(i) = in(i,1); % dbar
    epres = 0.63; % (max) std ±0.63 dbar
    wpres(i) = (epres).^(-2);
    
    psil(i)  = p( in(i,18) * 1e-6 ); % total Si mol/kg
    wpsil(i) = w( in(i,18), (in(i,18)*0.001) ); % 0.1% meas uncertainty
    
    % zero silicate is very unlikely; reset minimum to 0.1 micromolar
    izero = find(q(psil)==0);
    psil(izero) = p(1e-7);
    wpsil(izero) = w( 1e-7, 1e-7);

    pno2(i) = p( in(i,20) * 1e-6); % nitrite NO2^-
    wpno2(i) = w( in(i,20), (in(i,20)*0.001) ); % 0.1% meas uncertainty

    pno3(i) = p( in(i,22) * 1e-6 ); % nitrate NO3^-
    wpno3(i) = w( in(i,22), (in(i,22)*0.001) ); % 0.1% meas uncertainty
    
    ppo4(i)  = p( in(i,24) * 1e-6 ); % total phosphate mol/kg
    wppo4(i) = w( in(i,24), (in(i,24)*0.001) ); % 0.1% meas uncertainty
    
    % zero po4 is very unlikely; reset minimum to 1 nanomolar
    izero = find(q(ppo4)==0);
    ppo4(izero) = p(1e-9); %
    wppo4(izero) = w( 1e-9 , 1e-9 );
    
    % tb(i) = 0.0004157.*sal(i)./35; % mol/kg-SW Uppstrom (1974)  
    % calculate sigma with crank 3 times method 
    % 0.232 ± 0.005 B(mg/kg)/Cl(per mil)
    tb  = ((0.000232+0e-6)/10.811)*(sal(i)/1.80655); % same as above but easier to compare to the error 
    tbu = ((0.000232+5e-6)/10.811)*(sal(i)/1.80655);
    tbl = ((0.000232-5e-6)/10.811)*(sal(i)/1.80655);
    etb = (tbu - tbl)/2;
    ptb(i) = p(tb);
    wptb(i) = w(tb,etb);

    tf = ((0.000067+0e-6)/18.998).*(sal(i)/1.80655); % in mol/kg-SW Riley (1965)
    % 6.7 ± 0.1 x10^5
    tfu   = ((0.000067+1e-6)/18.998).*(sal(i)/1.80655);
    tfl   = ((0.000067-1e-6)/18.998).*(sal(i)/1.80655);
    etf = (tfu - tfl)/2;
    ptf(i) = p(tf);
    wptf(i) = w(tf,etf);
    
    ts = ((0.14+0.0e-4)/96.062).*(sal(i)/1.80655); % in mol/kg-SW Morris (1966)
    % 0.14 ± 0.00023
    tsu   = ((0.14+2.3e-4)/96.062).*(sal(i)/1.80655);
    tsl   = ((0.14-2.3e-4)/96.062).*(sal(i)/1.80655);
    ets = (tsu - tsl)/2;
    pts(i) = p(ts);                              
    wpts(i) = w(ts,ets);%                                     
    
    % need to add NH4 and H2S
    
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
sys = mksys(sys);
PrintTable(sys);
n = length(sys.variables);
% 
% Call QUODcarb for all the data points
%

fid2 = fopen('hoppe_TA_TC_v2.csv','w');    PrintCSV(sys,fid2);
fid3 = fopen('hoppe_TA_ph_v2.csv','w');    PrintCSV(sys,fid3);
fid4 = fopen('hoppe_TC_ph_v2.csv','w');    PrintCSV(sys,fid4);
fid5 = fopen('hoppe_TA_TC_ph_v2.csv','w'); PrintCSV(sys,fid5);

for i = 1:1310
    
    yobs = nan(n,1); wobs = nan(n,1);
    yobs(sys.iT)  = temp(i); wobs(sys.iT)  = wtemp(i); 
    yobs(sys.iS)  = sal(i);  wobs(sys.iS)  = wsal(i);
    yobs(sys.iP)  = pres(i); wobs(sys.iP)  = wpres(i);

    yobs(sys.iTC) = pdic(i); wobs(sys.iTC) = wpdic(i);
    yobs(sys.iTA) = palk(i); wobs(sys.iTA) = wpalk(i);
    yobs(sys.ih)  = pH(i);   wobs(sys.ih)  = wpH(i);
    % add pCO2, pCO3, pNO2, pNO3
    yobs(sys.ipco2) = ppco2(i); wobs(sys.ipco2) = wppco2(i); % maybe named sys.ippco2
    yobs(sys.ipco3) = pco3(i);  wobs(sys.ipco3) = wpco3(i);
    % no2 and no3 not in mksys
    
    if (ismember('Kb',sys.variables))
        yobs(sys.iTB) = ptb(i);
        wobs(sys.iTB) = wptb(i);
    end
    if (ismember('Ks',sys.variables))
        yobs(sys.iTS) = pts(i);
        wobs(sys.iTS) = wpts(i);
    end
    if (ismember('Kf',sys.variables))
        yobs(sys.iTF) = ptf(i);
        wobs(sys.iTF) = wptf(i);
    end
    if (ismember('K1p',sys.variables))
        yobs(sys.iTP) = ppo4(i);
        wobs(sys.iTP) = wppo4(i);
    end
    if (ismember('Ksi',sys.variables))
        yobs(sys.iTSi) = psil(i);
        wobs(sys.iTSi) = wpsil(i);
    end
    % add NH3 and H2S
    if (ismember('Knh3',sys.variables))
        yobs(sys.iTNH3) = pnh3(i);
        wobs(sys.iTNH3) = wpnh3(i);
    end
    if (ismember('Kh2s',sys.variables))
        yobs(sys.iTH2S) = ph2s(i);
        wobs(sys.iTH2S) = wph2s(i);
    end

    yobs_backup = yobs;
    wobs_backup = wobs;

    % TA TC 
    yobs = yobs_backup; wobs = wobs_backup;
    yobs(sys.ih) = nan; wobs(sys.ih) = nan;
    [y,sigy,yobs_out,wobs_out,iflag] = QUODcarb(yobs,wobs,temp(i),sal(i),pres(i),sys); 
    w = sigy.^(-2);
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




