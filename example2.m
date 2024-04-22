
% QUODcarb example 2
% this example will not work alone because it requires a data.mat file 

% measured 4 parameters n times saved in data.mat

load data.mat;
data = input_data; % save matrix to one variable

% populate opt structure
opt.K1K2 = 4; % option for K1K2 formulation
opt.KSO4 = 1; % option for KSO4 formulation
opt.KF   = 2; % option for KF formulation
opt.TB   = 2; % option for TB formulation
opt.phscale = 1; % 1 = tot, 2 = sws, 3 = free, 4 = nbs
opt.printcsv = 1; % print output to CSV? 1 = on, 0 = off
opt.fname = 'example2.csv'; % file name of desired CSV file
opt.printmes = 1; % print screen messages? 1 = on, 0 = off
opt.co2press = 0; % pressure correction for P2F (aka FugFac) and K0

n = 100; % number of datapoints

for i = 1:n
    % temperature and pressure INdependent
    obs(i).TC = data(i,1); % TC (umol/kg) is the first column
    obs(i).eTC = data(i,2); % TC error is in the 2nd column, 1 sigma
    obs(i).TA = data(i,3); % TA (umol/kg) is in the 3rd column
    obs(i).eTA = data(i,4); % TA error is in the 4th column, 1 sigma
    obs(i).sal = data(i,5); % salinity (PSU) 5th column
    obs(i).esal = data(i,6); % salinity error 6th column, 1 sigma

    % temperature and pressure DEpendent 1st system
    obs(i).tp(1).T = 25.0 ; % temperature for ph in Celsius
    obs(i).tp(1).eT = data(i,8); % ph temperature error, 1 sigma
    obs(i).tp(1).P = data(i,9); % pressure (dbar)
    obs(i).tp(1).eP = data(i,10); % pressure error, 1 sigma
    obs(i).tp(1).ph = data(i,11); % ph measured on total scale
    obs(i).tp(1).eph = data(i,12); % ph uncertainty, 1 sigma

    % temperature and pressure DEpendend 2nd system
    obs(i).tp(2).T = 20.0; % temperature for pCO2 in Celsius
    obs(i).tp(2).eT = data(i,13); % pCO2 temperature error, 1 sigma
    obs(i).tp(2).P = 0.0; % pressure in dbar for pCO2
    obs(i).tp(2).eP = 0.003; % dbar, pressure for pCO2 error, 1 sigma
    obs(i).tp(2).pco2 = data(i,14); % partial pressure CO2 (pCO2) (uatm)
    obs(i).tp(2).epco2 = data(i,15); % pCO2 error (uatm), 1 sigma
end


[est,obs,sys,iflag] = QUODcarb(obs,opt);

save example2.mat est; % save output est as mat file

% uncomment if you want to compare output to CO2SYS's output
% compare(obs,est,opt,1,1,'compare.csv'); % to compare to TC;TA pair
% compare(obs,est,opt,1,2,'compare.csv'); % to compare to TC;ph pair
% compare(obs,est,opt,2,7,'compare.csv'); % to compare to TC;pCO2 pair



