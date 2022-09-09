% take Hoppe data and put it into new QUODcarb

LOG10 = log(10);
p = @(x) -log10(x); % use p as in "pH" to denote the negative log base 10 function
q = @(x) 10.^(-x);  % use q, i.e., a backward p
w = @(x,e) p(1+e./x).^(-2); % convert x+/-e into precision for p(x)

load input.mat;
[in] = input;

palk  = zeros(102,1); wpalk  = zeros(102,1); 
pdic  = zeros(102,1); wpdic  = zeros(102,1); 
pH    = zeros(102,1); wpH    = zeros(102,1); 
psil  = zeros(102,1); wpsil  = zeros(102,1);
ppo4  = zeros(102,1); wppo4  = zeros(102,1);
ptb   = zeros(102,1); wptb   = zeros(102,1);
ptf   = zeros(102,1); wptf   = zeros(102,1);
pts   = zeros(102,1); wpts   = zeros(102,1);
ptnh3 = zeros(102,1); wptnh4 = zeros(102,1); % ammonia
pth2s = zeros(102,1); wpth2s = zeros(102,1); % hydrogen sulfide

sal  = zeros(102,1);  esal   = zeros(102,1); wsal  = zeros(102,1);
temp = zeros(102,1);  etemp  = zeros(102,1); wtemp = zeros(102,1);
pres = zeros(102,1);  epres  = zeros(102,1); wpres = zeros(102,1);



%
% convert all data to -log10(data) and compute the precision of the resulting data
%
for i = 1:102
    
    palk(i) = p(in(i,6)*1e-6); % mol/kg
    wpalk(i) = w(in(i,6),in(i,14));

    pdic(i) = p(in(i,7)*1e-6); % mol/kg
    wpdic(i) = w(in(i,7),in(i,15)); % DIC error

    pH(i) = in(i,8); % total pH scale
    wpH(i)  = in(i,16).^(-2); % pH error

    sal(i) = in(i,1); % salinity
    esal(i) = in(i,10);
    wsal(i) = (esal(i)).^(-2);

    temp(i) = in(i,2); % deg C
    etemp(i) = in(i,11); 
    wtemp(i) = (etemp(i)).^(-2);

    pres(i) = in(i,3); % dbar
    %epres = 0.03
    wpres(i) = (0.03).^(-2);
    
    psil(i)  = p( in(i,5) * 1e-6 ); % total Si mol/kg
    wpsil(i) = w( in(i,5), in(i,13) );
    
    % zero silicate is very unlikely; reset minimum to 0.1 micromolar
    izero = find(q(psil)==0);
    psil(izero) = p(1e-7);
    wpsil(izero) = w( 1e-7, 1e-7);
    
  
    ppo4(i)  = p( in(i,4) * 1e-6 ); % total phosphate mol/kg
    wppo4(i) = w( in(i,4), in(i,12) );
    
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
sys.abr = {'air-sea','carbonate'};
%sys.abr = {'air-sea','carbonate','water'};
%sys.abr = {'air-sea','carbonate','water','borate'};
%sys.abr = {'air-sea','carbonate','water','borate','sulfate'};
%sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride'};
%sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride','phosphate'};
%sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride','phosphate','silicate'};
%sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride', ...
% 'phosphate','silicate','ammonia'};
%sys.abr = {'air-sea','carbonate','water','borate','sulfate','fluoride', ...
% 'phosphate','silicate','ammonia','sulfide'};
sys = mksys(sys);
PrintTable(sys);
n = length(sys.variables);
% 
% Call QUODcarb for all the data points
%  1  2  3  4  5  6  7   8  9  10 11  12   13    14   15    16
%  T eT  S eS  P eP TC wTC TA wTA ph wph pco2 wpco2 hco3 whco3
D25 = zeros(102,16);
D20 = zeros(102,16);
D15 = zeros(102,16);
D10 = zeros(102,16);
D05 = zeros(102,16);

for i = 1:102
    
    yobs = nan(n,1); wobs = nan(n,1);
    yobs(sys.iT)  = temp(i); wobs(sys.iT)  = wtemp(i); 
    yobs(sys.iS)  = sal(i);  wobs(sys.iS)  = wsal(i);
    yobs(sys.iP)  = pres(i); wobs(sys.iP)  = wpres(i);

    yobs(sys.iTC) = pdic(i); wobs(sys.iTC) = wpdic(i);
    yobs(sys.iTA) = palk(i); wobs(sys.iTA) = wpalk(i);
    yobs(sys.ih)  = pH(i);   wobs(sys.ih)  = wpH(i);
    
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

    yobs(sys.iT) = 25;
    [y,sigy,yobs_out,wobs_out,iflag] = QUODcarb(yobs,wobs,temp(i),sal(i),pres(i),sys); 
    w = sigy.^(-2);
    D25(i,1) = y(sys.iT);
    D25(i,2) = w(sys.iT);
    D25(i,3) = y(sys.iS);
    D25(i,4) = w(sys.iS);
    D25(i,5) = y(sys.iP);
    D25(i,6) = w(sys.iP);
    D25(i,7) = q(y(sys.iTC));
    D25(i,8) = w(sys.iTC);
    D25(i,9) = q(y(sys.iTA));
    D25(i,10) = w(sys.iTA);
    D25(i,11) = y(sys.ih);
    D25(i,12) = w(sys.ih);
    D25(i,13) = q(y(sys.ipco2));
    D25(i,14) = w(sys.ipco2);
    D25(i,15) = q(y(sys.ihco3));
    D25(i,16) = w(sys.ihco3);


    yobs(sys.iT) = 20;
    [y,sigy,yobs_out,wobs_out,iflag] = QUODcarb(yobs,wobs,temp(i),sal(i),pres(i),sys); 
    w = sigy.^(-2);
    D20(i,1) = y(sys.iT);
    D20(i,2) = w(sys.iT);
    D20(i,3) = y(sys.iS);
    D20(i,4) = w(sys.iS);
    D20(i,5) = y(sys.iP);
    D20(i,6) = w(sys.iP);
    D20(i,7) = q(y(sys.iTC));
    D20(i,8) = w(sys.iTC);
    D20(i,9) = q(y(sys.iTA));
    D20(i,10) = w(sys.iTA);
    D20(i,11) = y(sys.ih);
    D20(i,12) = w(sys.ih);
    D20(i,13) = q(y(sys.ipco2));
    D20(i,14) = w(sys.ipco2);
    D20(i,15) = q(y(sys.ihco3));
    D20(i,16) = w(sys.ihco3);

    yobs(sys.iT) = 15;
    [y,sigy,yobs_out,wobs_out,iflag] = QUODcarb(yobs,wobs,temp(i),sal(i),pres(i),sys); 
    w = sigy.^(-2);
    D15(i,1) = y(sys.iT);
    D15(i,2) = w(sys.iT);
    D15(i,3) = y(sys.iS);
    D15(i,4) = w(sys.iS);
    D15(i,5) = y(sys.iP);
    D15(i,6) = w(sys.iP);
    D15(i,7) = q(y(sys.iTC));
    D15(i,8) = w(sys.iTC);
    D15(i,9) = q(y(sys.iTA));
    D15(i,10) = w(sys.iTA);
    D15(i,11) = y(sys.ih);
    D15(i,12) = w(sys.ih);
    D15(i,13) = q(y(sys.ipco2));
    D15(i,14) = w(sys.ipco2);
    D15(i,15) = q(y(sys.ihco3));
    D15(i,16) = w(sys.ihco3);

    yobs(sys.iT) = 10;
    [y,sigy,yobs_out,wobs_out,iflag] = QUODcarb(yobs,wobs,temp(i),sal(i),pres(i),sys); 
    w = sigy.^(-2);
    D10(i,1) = y(sys.iT);
    D10(i,2) = w(sys.iT);
    D10(i,3) = y(sys.iS);
    D10(i,4) = w(sys.iS);
    D10(i,5) = y(sys.iP);
    D10(i,6) = w(sys.iP);
    D10(i,7) = q(y(sys.iTC));
    D10(i,8) = w(sys.iTC);
    D10(i,9) = q(y(sys.iTA));
    D10(i,10) = w(sys.iTA);
    D10(i,11) = y(sys.ih);
    D10(i,12) = w(sys.ih);
    D10(i,13) = q(y(sys.ipco2));
    D10(i,14) = w(sys.ipco2);
    D10(i,15) = q(y(sys.ihco3));
    D10(i,16) = w(sys.ihco3);

    
    yobs(sys.iT) = 5;
    [y,sigy,yobs_out,wobs_out,iflag] = QUODcarb(yobs,wobs,temp(i),sal(i),pres(i),sys); 
    w = sigy.^(-2);
    D05(i,1) = y(sys.iT);
    D05(i,2) = w(sys.iT);
    D05(i,3) = y(sys.iS);
    D05(i,4) = w(sys.iS);
    D05(i,5) = y(sys.iP);
    D05(i,6) = w(sys.iP);
    D05(i,7) = q(y(sys.iTC));
    D05(i,8) = w(sys.iTC);
    D05(i,9) = q(y(sys.iTA));
    D05(i,10) = w(sys.iTA);
    D05(i,11) = y(sys.ih);
    D05(i,12) = w(sys.ih);
    D05(i,13) = q(y(sys.ipco2));
    D05(i,14) = w(sys.ipco2);
    D05(i,15) = q(y(sys.ihco3));
    D05(i,16) = w(sys.ihco3);

end

save fakedata.mat D25 D20 D15 D10 D05




