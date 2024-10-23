
% QUODcarb example 1
clear all
% measured TC and TA only

% populate opt structure, see readme for default settings
opt.phscale = 1;        % which phscale are we interested in
opt.printcsv = 1;       % print output to CSV? 1 = on, 0 = off
opt.fname = 'example1.csv'; % file name of desired CSV file


% populate obs structure
obs.TC = 2150; % umol/kg-SW
obs.uTC = 5; % ± 5 umol/kg-SW, 1 sigma
obs.TA = 2300; % umol/kg-SW
obs.uTA = 5; % ± 5 umol/kg-SW, 1 sigma
obs.sal = 32.7; % PSU
obs.usal = 0.02; % ± 0.02 PSU, 1 sigma
obs.tp.T = 20; % deg Celsius
obs.tp.uT = 0.001; % deg Celsius, 1 sigma
obs.tp.P = 0; % dbar
obs.tp.uP = 0.005; % ± 0.005 dbar, 1 sigma

[est,obs,sys,iflag,opt] = QUODcarb(obs,opt);

% output est values should be:
% est.TC = 2150                 est.uTC = 4.9942
% est.TA = 2300                 est.uTA = 4.9946
% est.tp.ph = 7.8403            est.tp.uph = 0.0180
% est.tp.pco2 = 704.7821        est.tp.upco2 = 32.3876
% est.tp.co3 = 118.3386         est.tp.uco3 = 4.3372


% if want to compare output to CO2SYS's output...
% add 'compare.m' to current directory
% and uncomment below:
% compare(obs,est,opt,1,1,2,'compare.csv');


