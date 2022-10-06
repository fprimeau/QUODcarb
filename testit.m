% take Hoppe data and put it into new QUODcarb

LOG10 = log(10);
p = @(x) -log10(x); % use p as in "pH" to denote the negative log base 10 function
q = @(x) 10.^(-x);  % use q, i.e., a backward p
w = @(x,e) p(1+e./x).^(-2); % convert x+/-e into precision for p(x)

load fakedata.mat
i = 1;
% measurements that are independent of (T,P)
obs.TC   = D25(i,7) + nan; % nan means not measured
obs.wTC  = D25(i,8) + nan; % nan means not measured
obs.TA   = D25(i,9) + nan; % nan means not measured
obs.wTA  = D25(i,10)+ nan; % nan means not measured
obs.sal  = D25(i,3);
obs.wsal = D25(i,4);
% first (T,P)-dependent measurement
obs.m(1).T     = D25(i,1); % temp 1
obs.m(1).wT    = D25(i,2); % precision of temp
obs.m(1).P     = D25(i,5); % pressure 1
obs.m(1).wP    = D25(i,6);
obs.m(1).ph    = D25(i,11);
obs.m(1).wph   = D25(i,12);
% second (T,P)-dependent measurement
obs.m(2).T     = D15(i,1); % temp 2
obs.m(2).wT    = D15(i,2);
obs.m(2).P     = D15(i,5); % press 2
obs.m(2).wP    = D15(i,6);
obs.m(2).pco2  = D15(i,13);
obs.m(2).wpco2 = D15(i,14);
% third (T,P)-dependent measurement
obs.m(3).T     = D10(i,1); % temp 3 (could be a wanted output temp)
obs.m(3).wT    = D10(i,2);
obs.m(3).P     = D10(i,5); % pres 3
obs.m(3).wP    = D10(i,6);
%obs.m(3).hco3  = D10(i,15);
%obs.m(3).whco3 = D10(i,16);

%
% initialize the QUODcarb CO2-system solver
%
%sys.abr = {'carbonate'};
%sys.abr = {'carbonate','water'};
%sys.abr = {'carbonate','water','borate'};
%sys.abr = {'carbonate','water','borate','sulfate'};
%sys.abr = {'carbonate','water','borate','sulfate','fluoride'};
%sys.abr = {'carbonate','water','borate','sulfate','fluoride','phosphate'};
sys.abr = {'carbonate','water','borate','sulfate','fluoride','phosphate','silicate'};
%sys.abr = {'carbonate','water','borate','sulfate','fluoride', ...
%           'phosphate','silicate','ammonia','sulfide'};
sys = mksysV2(obs,sys);

[est,obs,iflag] = QUODcarbV2(obs,sys);
