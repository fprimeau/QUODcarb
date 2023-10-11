

% QUODcarbV8
% saved in FP_QUODcarb/v7

function [est,obs,iflag] = QUODcarbV8(obs,opt)
%
% OUTPUT:
%   est := posterior estimates of co2-system variables, equilibrium constants, and precisions
%   obs := same as input except that the pK's have been added
% iflag := 0 if solver converged to specified accuracy 
%          1 after reaching maximum number of iterations without converging
%          2 if it was one order of magnitude away from converging at
%               maximum iteration
%   CSV := optional CSV output, if turned on in options
%
% INPUT:
%   obs  := co2-system measured quantities with precisions 
%   opt  := solver options input by user 
%
% -------------------------------------------------------------------------
%
% SYNTAX example:
%   obs.sal         = salinity;     (PSU)           
%   obs.esal        = sal_error;    (±sigma)        
%   obs.TC          = total_c;      (umol/kg-SW)    
%   obs.eTC         = TC_error;     (±sigma)        
%   obs.TA          = total_alk;    (umol/kg-SW)    
%   obs.eTA         = alk_error;    (±sigma)        
%   obs.tp(1).T     = temp;         (deg C)         
%   obs.tp(1).eT    = temp_error;   (±sigma)        
%   obs.tp(1).P     = pressure;     (dbar)          
%   obs.tp(1).eP    = pres_error;   (±sigma)     
%   obs.tp(1).ph    = ph_meas;      
%   obs.tp(1).eph   = ph_error;     (±sigma)
%
%   opt.K1K2        = 10;           % (Lueker et al 2000)
%   opt.KSO4        = 1;            % (Dickson et al 1990a) 
%   opt.KF          = 2;            % (Perez and Fraga 1987
%   opt.TB          = 2;            % (Lee et al. 2010)
%   opt.phscale     = 1;            % (1=tot, 2=free, 3=sws, 4=nbs)
%   opt.printcsv    = 1;            % (1=on, 0=off)
%   opt.fid         = 'output.csv'; % (CSV filename)
%   opt.printmes    = 1;            % (1=on, 0=off)
%   opt.co2press    = 1;            % (1=on, 0=off)
%   opt.abr         = 'all';        % {'borate','sulfate','fluoride',...
%                                       'phosphate','silicate','ammonia', ...
%                                       'sulfide','solubility'}
%
%--------------------------------------------------------------------------
% 
% INPUT OPTIONS:
%   opt.K1K2  -> choice of K1 and K2 formulation
%           1 = Roy et al, 1993
%           2 = Goyet & Poisson, 1989
%           3 = Hansson, 1973          REFIT by Dickson and Millero, 1987
%           4 = Mehrbach et al., 1973  REFIT by Dickson and Millero, 1987
%           5 = Hansson, 1973 and Mehrbach, 1973 
%                                      REFIT by Dickson and Millero, 1987
%        x(6) = x(GEOSECS)            ~NOT AVAILABLE IN QUODcarb~
%        x(7) = x(Peng)               ~NOT AVAILABLE IN QUODcarb~
%        x(8) = x(Millero, 1979)      ~NOT AVAILABLE IN QUODcarb~
%           9 = Cai and Wang, 1998
%          10 = Lueker et al., 2000    (DEFAULT)
%          11 = Mojica Prieto and Millero, 2002
%          12 = Millero et al., 2000
%          13 = Millero et al., 2002
%          14 = Millero et al., 2006
%          15 = Waters, Millero, and Woosley, 2014
%          16 = Sulpis et al., 2020
%          17 = Schockman and Byrne, 2021
%
%   opt.KSO4  -> choice of KSO4 formulation
%           1 = Dickson (1990a)         (DEFAULT)
%           2 = Khoo et al., 1977
%           3 = Waters and Millero, 2013
%
%   opt.KF    -> choice of KF formulation
%           1 = Dickson and Riley, 1979
%           2 = Perez and Fraga, 1987  (DEFAULT)
%
%   opt.TB    -> choice of total borate formulation
%           1 = Uppstrom, 1979
%           2 = Lee et al., 2010       (DEFAULT)
%
%--------------------------------------------------------------------------
%
% OUTPUT:
%   est ->  'est' structure with best estimate contains:
%               1. p(value) and p(error) where p(x) = -log10(x)
%               2. value and average error about the value in 'q'
%                   where q(x) = x^(-10)
%               3. upper and lower bounds in 'q' space, not symmetric
%                   about the value in 'q' space
%   csv ->  csv file with most of est populated in a spreadsheet, 
%                 contains column headers with labels and units                 
%                    -does not include upper and lower errors
%
%--------------------------------------------------------------------------
%
% Changes? -> the only things you may want to change are: 
%               1. tolerance level of Newton solver -> line 120
%               2. Max Iteration number -> MAXIT in newtn.m
%
%--------------------------------------------------------------------------


    opt = check_opt(opt); % check opt structure

    sys = mksysV8(obs(1),opt.abr,opt.phscale); % create indexing

    nD = length(obs);  nv = size(sys.K,2);

    % populate obs, yobs, wobs at each datapoint
    [obs,yobs,wobs] = parse_input(obs,sys,opt,nD);
  
    for i = 1:nD

        z0 = init(yobs(i,:),sys,opt); % initialize
        tol = 1e-8;

        gun = @(z) grad_limpco2(z,yobs(i,:),wobs(i,:),sys,opt);
        [z,J,iflag(i)] = newtn(z0,gun,tol);
        if (iflag(i) ~=0) && (opt.printmes ~= 0)
            fprintf('Newton''s method iflag = %i at i = %i \n',iflag(i),i);
        end
       
        % posterior uncertainty
        warning('off');
        C = inv(J);
        warning('on');
        C = C(1:nv,1:nv);
        y = z(1:nv);
        sigy = sqrt(diag(C));

        if (sum(isnan(sigy)) > 0) && (opt.printmes ~= 0)
            fprintf('NaN found in output means faulty run. i = %i\n',i)
        end
        % populate est
        [est(i)] = parse_output(z,sigy,opt,sys);

        if opt.printcsv == 1
            % make QUODcarb.CSV if requested
            if i == 1
                fid = fopen(opt.fid,'w');
                PrintCSVv8(est(i),fid,opt); % make column headers
            end
             % fill with one row of data
            PrintCSVv8(opt,est(i),obs(i),iflag(i),fid);
        end
    end
end

% -------------------------------------------------------------------------

function [g,H] = grad_limpco2(z,y,w,sys,opt)
    I = eye(length(z));
    H = zeros(length(z),length(z));
    im = sqrt(-1);
    %test_g = zeros(length(z),1);
    for k = 1:length(z)
        [f,g] = limpco2( z + im * eps^3 * I(:,k), y, w, sys, opt);
        %   test_g(k) = imag(f)/eps^3;
        H(k,:) = imag( g(:) ) / eps^3 ;
    end
    g = real( g(:) );
    % [ g, test_g ]
    % keyboard
end

% ---------------------------------------------------------------------------------

function [f,g] = limpco2(z,y,w,sys,opt)
% [f,g] = limpco2(z,y,w,tc,s,P)
%
% Negative log probability for the co2 system  a.k.a. log improbability i.e., limp!
%
% INPUT:
%
%   z  := system state including equilibrium constants and lagrange multipliers
%   y  := measured components, non-measured components set to NaN
%   w  := measurement precisions, same size and order as y, non-measured components set to anything, they are ignored
% gpK  := first order derivatives of pK wrt T,S,P
% ggpK := second order derivatives of pK wrt T,S,P

% OUTPUT:
%
%   f := limp
%   g := grad f w.r.t. x
%   h := hessian of f w.r.t. x

    p   = sys.p;    
    q   = sys.q;
    M   = sys.M;    
    K   = sys.K;
    nrk     = size(K,1);
    nTP     = length(sys.tp); 
    nv      = size(M,2);
    nlam    = size(M,1) + size(K,1) + (nTP); 
            % one extra lagrange multiplier for each 
            % (T,P)-dependent ph_free
    x       =  z(1:nv);      % measureable variables
    lam     =  z(nv+1:end);  % Lagrange multipliers 
    
    % Make a vector of measured quantities    
    i   = find(~isnan(y));
    y   = y(i)';
    
    % Make a precision matrix
    W   = diag(w(i));
    
    % Build a matrix that Picks out the measured components of x
    I   = eye(nv); % for chain rule
    PP  = I(i,:); % picking/pick out the measured ones
    e   = PP*x - y; % calculated - measured (minus)
    
    % fill zpK and zgpK with associated calculated pK and gpK values
    % [zpK, zgpK, rph_tot, rph_sws, rph_free, rph_nbs] = parse_zpK(x,sys,opt);
    [zpK, zgpK, rph_free] = parse_zpK(x,sys,opt);

    % constraint equations
    c   = [  M * q( x ); ...
            (-K * x) + zpK;...
            rph_free] ;

    % keyboard
    f   = 0.5 *  e.' * W * e  + lam.' * c ;  % limp, method of lagrange multipliers    
    % -(-1/2 sum of squares) + constraint eqns, minimize f => grad(f) = 0
    
    if ( nargout > 1 ) % compute the gradient
        grph_free = zeros(nTP,nv);       
        for i = 1:nTP
            grph_free(i,[ sys.iTS           , sys.tp(i).iKs    , ... % (d/dx) _pTS,_pKs
                          sys.iTF           , sys.tp(i).iKf    , ... % (d/dx) _pTF,_pKf
                          sys.tp(i).ipfH     , sys.tp(i).iph    , ... % (d/dx) _pfH,_ph
                          sys.tp(i).iph_free  ])                 ... % (d/dx) _ph_free
                                            = sys.tp(i).grph_free(x);
        end
        % constraint eqns wrt -log10(concentrations)
        dcdx = [ M * diag( sys.dqdx( x ) ); ...
                (-K + zgpK)    ;            ... 
                grph_free    ] ;

        g    = [ e.' * W * PP +  lam.' * dcdx ,  c.' ];
        % keyboard
    end

    if ( nargout > 2 ) % compute the Hessian, not used as of v7
        ddq     =  diag( sys.d2qdx2( x ) ); % q"
        [nr,nc] = size(M);
        gg      = zeros(nc,1);
        for row = (1:nr)
            gg  = gg + lam(row)*diag(M(row,:))*ddq;
        end
        for row = (nr+1):(nr+nrk)
            gg  = gg + lam(row)*(-zggpK((row-nr),:,:)); % ggpK
        end
        dhfdx2  = zeros(nc,nc);
        ii      = [sys.iKs,sys.iTS];
        dhfdx2(ii,ii) = sys.ggf2t(x);
        gg      = gg + lam(nr+1)*dhfdx2;
        H       = [  PP.'*W*PP + gg , dcdx.'  ; ...
                        dcdx        , zeros(nlam)  ];
    end
    
end

% -----------------------------------------------------------------------------------

function [opt] = check_opt(opt)
    % check opt input
    isbad = @(thing) (isempty(thing) & sum(isnan(thing)));

    % opt.printmes
    if ~isfield(opt,'printmes') || isbad(opt.printmes)
        opt.printmes = 1; % default on
    end
    % opt.K1K2
    if ~isfield(opt,'K1K2') || isbad(opt.K1K2)
        if opt.printmes ~= 0
            fprintf('No K1K2 formulation chosen. Assuming opt.K1K2 = 4\n');
        end
        opt.K1K2 = 4; % default K1K2 setting
    elseif opt.K1K2 > 18 || opt.K1K2 < 1 || ...
                opt.K1K2 == 6 || opt.K1K2 == 7 || opt.K1K2 == 8
        if opt.printmes ~= 0
            fprintf(['Invalid K1K2 formulation chosen. ' ...
                'Assuming opt.K1K2 = 4\n']);
        end
        opt.K1K2 = 4; % default K1K2 setting
    end
    % opt.KSO4
    if ~isfield(opt,'KSO4') || isbad(opt.KSO4)
        if opt.printmes ~= 0
            fprintf('No KSO4 formulation chosen. Assuming opt.KSO4 = 1\n');
        end
        opt.KSO4 = 1; % default opt.KSO4 setting
    elseif opt.KSO4 > 3 || opt.KSO4 < 1
        if opt.printmes ~= 0
            fprintf(['Invalid KSO4 formulation chosen. ' ...
                'Assuming opt.KSO4 = 1\n']);
        end
        opt.KSO4 = 1; % default opt.KSO4 setting
    end
    % opt.KF
    if ~isfield(opt,'KF') || isbad(opt.KF)
        if opt.printmes ~= 0
            fprintf('No KF formulation chosen. Assuming opt.KF = 2 \n');
        end
        opt.KF = 2; % default KF
    elseif opt.KF > 2 || opt.KF < 1 
        if opt.printmes ~= 0
            fprintf(['Invalid KF formulation chosen. ' ...
                'Assuming opt.KF = 2 \n']);
        end
        opt.KF = 2;
    end
    % opt.TB
    if ~isfield(opt,'TB') || isbad(opt.TB)
        if opt.printmes ~= 0
            fprintf('No TB formulation chosen. Assuming opt.TB = 2\n');
        end
        opt.TB = 2;
    elseif opt.TB > 2 || opt.TB < 1
        if opt.printmes ~= 0
            fprintf(['Invalid TB formulation chosen. ' ...
                'Assuming opt.TB = 2\n']);
        end
        opt.TB = 2;
    end
    % opt.phscale
    if ~isfield(opt,'phscale') || isbad(opt.phscale)
            error(['No opt.phscale chosen, must choose 1 = tot, ' ...
                '2 = sws, 3 = free, 4 = NBS \n'])
    elseif opt.phscale > 4 || opt.phscale < 1
            eror(['Invalid opt.phscale chosen, must choose 1 = tot, ' ...
                '2 = sws, 3 = free, 4 = NBS \n'])
    end
    % opt.printcsv and opt.fid
    if ~isfield(opt,'printcsv') || isbad(opt.printcsv)
        opt.printcsv = 0; % default off
    elseif opt.printcsv > 1 || opt.printcsv < 0
        if opt.printmes ~= 0
            fprintf('Invalid CSV opt chosen. Assuming opt.csv = 1\n');
        end
    else
        if ~isfield(opt,'fid') || isbad(opt.fid)
            opt.fid = 'QUODcarb_output.csv';
        end
    end
    % opt.co2press
    if (~isfield(opt,'co2press')) || isbad(opt.co2press)
        opt.co2press = 0; % default co2press
    elseif opt.co2press > 1 || opt.co2press < 0
        if opt.printmes ~= 0
            fprintf('Invalid opt.co2press input, reset to default.\n')
        end
        opt.co2press = 0; % default co2press
    end
    % opt.abr (Acid/Base Reactions or systems to include)
    if (~isfield(opt,'abr')) || (strcmp(opt.abr,""))
        if opt.printmes ~= 0
            fprintf(['No acid base system chosen. ' ...
                'Assuming opt.abr = {''all''}\n'])
        end
        opt.abr = {'all'}; % default acid/base system = 'all'
    end
    if (strcmp(opt.abr,'all'))
        opt.abr = {'phosphate','silicate','ammonia','sulfide','solubility'};
    end
end


% ------------------------------------------------------------------------

function [obs,yobs,wobs] = parse_input(obs,sys,opt,nD)

    isgood  = @(thing) (~isempty(thing) & ~sum(isnan(thing)));
    p       = sys.p;
    q       = sys.q;
    % convert x+/-e into precision for p(x)
    w       = @(x,e) abs( p(1 + e./x) ).^(-2);

    nv      = size(sys.K,2);
    yobs    = nan(nD,nv);
    wobs    = nan(nD,nv);

    for i = 1:nD
        % make sure all the required fields in the obs struct exist
        if (~isfield(obs(i), 'sal'))
            if opt.printmes ~= 0
                error('Need to provide salinity measurement.');
            end
        else
            yobs(i,sys.isal) = obs(i).sal;
        end
        if (~isfield(obs(i), 'esal'))
            obs(i).esal         = 0.002; % std = 0.002 PSU
            wobs(i,sys.isal)    = (obs(i).esal)^(-2);
            if opt.printmes ~= 0
                fprintf('Warning: Assuming salinity uncertainty is 0.002 PSU');
            end
        else
            wobs(i,sys.isal)    = (obs(i).esal)^(-2); % std e -> w
        end
        if (~isfield(obs(i),'TC')) || (~isgood(obs(i).TC))
            obs(i).TC       = nan; %[]
            yobs(i,sys.iTC) = nan;
        else
            yobs(i,sys.iTC) = p((obs(i).TC)*1e-6); % convt to mol/kg
        end
        if (~isfield(obs(i),'eTC')) || (~isgood(obs(i).eTC))
            obs(i).eTC      = nan;
            wobs(i,sys.iTC) = nan;
        else
            wobs(i,sys.iTC) = w(obs(i).TC,obs(i).eTC); % std e -> w
        end
        if(~isfield(obs(i),'TA'))  || (~isgood(obs(i).TA))
            obs(i).TA       = nan; %[]
            yobs(i,sys.iTA) = nan;
        else
            yobs(i,sys.iTA) = p((obs(i).TA)*1e-6); % convt to mol/kg
        end
        if (~isfield(obs(i),'eTA'))  || (~isgood(obs(i).eTA))
            obs(i).eTA      = nan;
            wobs(i,sys.iTA) = nan;
        else
            wobs(i,sys.iTA) = w(obs(i).TA,obs(i).eTA); % std e -> w
        end
        if (~isfield(obs,'tp'))
            if opt.printmes ~= 0
                error('Need to provide temperature and pressure measurement.')
            end
        end

        % create obs structure with fieldnames
        if (~isfield(obs(i).tp(1),'epK0'))
            obs(i).tp(1).epK0   = [];
        end
        if (~isfield(obs(i).tp(1),'pK0'))
            obs(i).tp(1).pK0    = [];
        end
        if (~isfield(obs(i).tp(1),'epK1'))
            obs(i).tp(1).epK1   = [];
        end
        if (~isfield(obs(i).tp(1),'pK1'))
            obs(i).tp(1).pK1    = [];
        end
        if (~isfield(obs(i).tp(1),'epK2'))
            obs(i).tp(1).epK2   = [];
        end
        if (~isfield(obs(i).tp(1),'pK2'))
            obs(i).tp(1).pK2    = [];
        end
        if (~isfield(obs(i).tp(1),'epp2f')) || ...
                (~isgood(obs(i).tp(1).epp2f))
            obs(i).tp(1).epp2f  = [];
        end
        if (~isfield(obs(i).tp(1),'pp2f')) || ...
                (~isgood(obs(i).tp(1).pp2f))
            obs(i).tp(1).pp2f   = [];
        end
        if (~isfield(obs(i).tp(1),'efco2')) || ...
                (~isgood(obs(i).tp(1).efco2))
            obs(i).tp(1).efco2  = [];
        end
        if (~isfield(obs(i).tp(1),'fco2')) || ...
                (~isgood(obs(i).tp(1).fco2))
            obs(i).tp(1).fco2   = [];
        end
        if (~isfield(obs(i).tp(1),'epco2')) || ...
                (~isgood(obs(i).tp(1).epco2))
            obs(i).tp(1).epco2  = [];
        end
        if (~isfield(obs(i).tp(1),'pco2')) || ...
                (~isgood(obs(i).tp(1).pco2))
            obs(i).tp(1).pco2   = [];
        end
        if (~isfield(obs(i).tp(1),'eco2st')) || ...
                (~isgood(obs(i).tp(1).eco2st))
            obs(i).tp(1).eco2st = [];
        end
        if (~isfield(obs(i).tp(1),'co2st')) || ...
                (~isgood(obs(i).tp(1).co2st))
            obs(i).tp(1).co2st  = [];
        end
        if (~isfield(obs(i).tp(1),'ehco3')) || ...
                (~isgood(obs(i).tp(1).ehco3))
            obs(i).tp(1).ehco3  = [];
        end
        if (~isfield(obs(i).tp(1),'hco3')) || ...
                (~isgood(obs(i).tp(1).hco3))
            obs(i).tp(1).hco3   = [];
        end
        if (~isfield(obs(i).tp(1),'eco3')) || ...
                (~isgood(obs(i).tp(1).eco3))
            obs(i).tp(1).eco3   = [];
        end
        if (~isfield(obs(i).tp(1),'co3')) || (~isgood(obs(i).tp(1).co3))
            obs(i).tp(1).co3    = [];
        end
        if (~isfield(obs(i).tp(1),'eph')) || (~isgood(obs(i).tp(1).eph))
            obs(i).tp(1).eph    = [];
        end
        if (~isfield(obs(i).tp(1),'ph')) || (~isgood(obs(i).tp(1).ph))
            obs(i).tp(1).ph     = [];
        end
        if (~isfield(obs(i).tp(1),'eph_free')) || ...
                (~isgood(obs(i).tp(1).eph_free))
            obs(i).tp(1).eph_free = [];
        end
        if (~isfield(obs(i).tp(1),'ph_free')) || ...
                (~isgood(obs(i).tp(1).ph_free))
            obs(i).tp(1).ph_free = [];
        end
        if (~isfield(obs(i).tp(1),'epfH')) || ...
                (~isgood(obs(i).tp(1).epfH))
            obs(i).tp(1).epfH   = [];
        end
        if (~isfield(obs(i).tp(1),'pfH')) || (~isgood(obs(i).tp(1).pfH))
            obs(i).tp(1).pfH    = [];
        end


        if (~isfield(obs(i).tp(1), 'epKw'))
            obs(i).tp(1).epKw   = [];
        end
        if (~isfield(obs(i).tp(1), 'pKw'))
            obs(i).tp(1).pKw    = [];
        end
        if (~isfield(obs(i).tp(1),'eoh')) || (~isgood(obs(i).tp(1).eoh))
            obs(i).tp(1).eoh    = [];
        end
        if (~isfield(obs(i).tp(1),'oh')) || (~isgood(obs(i).tp(1).oh))
            obs(i).tp(1).oh     = [];
        end

        % borate system
        if (~isfield(obs(i).tp(1),'epKb'))
            obs(i).tp(1).epKb   = [];
        end
        if (~isfield(obs(i).tp(1),'pKb'))
            obs(i).tp(1).pKb    = [];
        end
        if (~isfield(obs(i),'TB')) || (~isgood(obs(i).TB))
            if opt.TB == 1
                % Uppstrom, L., Deep-Sea Research 21:161-162, 1974
                % ( copied from Orr's code )
                % TB = ( 0.000232/ 10.811) * (sal/1.80655)
                obs(i).TB       = (0.0004157 * obs(i).sal / 35) * (1e6) ; % convt to µmol/kg
                yobs(i,sys.iTB) = p(obs(i).TB*1e-6);
            elseif opt.TB == 2
                % Lee, Kim, Myrne, Millero, Feely, Yong-Ming Liu. 2010.
                % Geochemica Et Cosmochimica Acta 74 (6): 1801-1811.
                % ( copied from Sharp's code )
                % TB = (0.0002414/ 10.811) * (sal/1.80655)
                obs(i).TB       = (0.0004326 * obs(i).sal / 35)* (1e6) ; % umol/kg-SW
                yobs(i,sys.iTB) = p(obs(i).TB*1e-6);
            elseif opt.K1K2 == 6 || opt.K1K2 == 7
                % Culkin, F., in Chemical Oceanography,
                % ed. Riley and Skirrow, 1965: GEOSECS references this
                % (copied from Orr's code)
                obs(i).TB       = (0.0004106 * obs(i).sal / 35) * (1e6) ; % umol/kg
                yobs(i,sys.iTB) = p(obs(i).TB*1e-6);
            end
        else
            if ((obs(i).TB) == 0)
                obs(i).TB   = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.iTB) = p(obs(i).TB*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs(i), 'eTB'))  || (~isgood(obs(i).eTB))
            if opt.TB == 1 % Uppstrom 1974
                % std 5e-6 on avg 2.32e-4
                TBu = ( ( (2.32e-4 + 5e-6)/10.811) * obs(i).sal/1.80655 );
                TBl = ( ( (2.32e-4 - 5e-6)/10.811) * obs(i).sal/1.80655 );
                eTB = (TBu - TBl) /2 ;
                obs(i).eTB = (eTB) * 1e6 ; % µmol/kg
                wobs(i,sys.iTB) = w(obs(i).TB,obs(i).eTB);
            elseif opt.TB == 2 % Lee et al 2010
                % std 9e-6 on avg 2.414e-4
                TBu = ( ( (2.414e-4 + 9e-6)/10.811) * obs(i).sal/1.80655);
                TBl = ( ( (2.414e-4 - 9e-6)/10.811) * obs(i).sal/1.80655);
                eTB = (TBu - TBl) /2;
                obs(i).eTB = (eTB)*1e6; % umol/kg
                wobs(i,sys.iTB) = w(obs(i).TB,obs(i).eTB);
            end
        else
            wobs(i,sys.iTB)     = w(obs(i).TB,obs(i).eTB);
        end
        if (~isfield(obs(i).tp(1),'eboh4')) || ...
                (~isgood(obs(i).tp(1).eboh4))
            obs(i).tp(1).eboh4  = [];
        end
        if (~isfield(obs(i).tp(1),'boh4')) || ...
                (~isgood(obs(i).tp(1).boh4))
            obs(i).tp(1).boh4   = [];
        end
        if (~isfield(obs(i).tp(1),'eboh3')) || ...
                (~isgood(obs(i).tp(1).eboh3))
            obs(i).tp(1).eboh3  = [];
        end
        if (~isfield(obs(i).tp(1),'boh3')) || ...
                (~isgood(obs(i).tp(1).boh3))
            obs(i).tp(1).boh3   = [];
        end

        % sulfate system
        if (~isfield(obs(i).tp(1), 'epKs'))
            obs(i).tp(1).epKs   = [];
        end
        if (~isfield(obs(i).tp(1), 'pKs'))
            obs(i).tp(1).pKs    = [];
        end
        if (~isfield(obs(i), 'TS'))  || (~isgood(obs(i).TS))
            % Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
            % copied from Orr's code
            obs(i).TS = ( 0.14 / 96.062 ) * ( obs(i).sal / 1.80655 )*1e6; % µmol/kg
            yobs(i,sys.iTS)     = p(obs(i).TS*1e-6);
        else
            if ((obs(i).TS) == 0)
                obs(i).TS = 1e-3; % µmol/kg reset minimum to 1 nanomolar
            end
            yobs(i,sys.iTS)     = p(obs(i).TS*1e-6);
        end
        if (~isfield(obs(i), 'eTS'))  || (~isgood(obs(i).eTS))
            % 0.14000 ± 0.00023
            TSu = ( ( (0.14+0.00023)/96.062 ) * obs(i).sal/ 1.80655 );
            TSl = ( ( (0.14-0.00023)/96.062 ) * obs(i).sal/ 1.80655 );
            eTS                 = (TSu - TSl) / 2;
            obs(i).eTS          = (eTS*1e6) ;
            wobs(i,sys.iTS)     = w(obs(i).TS,obs(i).eTS);
        else
            wobs(i,sys.iTS)     = w(obs(i).TS,obs(i).eTS);
        end
        if (~isfield(obs(i).tp(1), 'so4')) || ...
                (~isgood(obs(i).tp(1).so4))
            obs(i).tp(1).so4    = [];
        end
        if (~isfield(obs(i).tp(1), 'eso4')) || ...
                (~isgood(obs(i).tp(1).eso4))
            obs(i).tp(1).eso4   = [];
        end
        if (~isfield(obs(i).tp(1), 'hso4')) || ...
                (~isgood(obs(i).tp(1).hso4))
            obs(i).tp(1).hso4   = [];
        end
        if (~isfield(obs(i).tp(1), 'ehso4')) || ...
                (~isgood(obs(i).tp(1).ehso4))
            obs(i).tp(1).ehso4  = [];
        end

        % fluoride system
        if (~isfield(obs(i).tp(1), 'epKf'))
            obs(i).tp(1).epKf   = [];
        end
        if (~isfield(obs(i).tp(1), 'pKf'))
            obs(i).tp(1).pKf    = [];
        end
        if (~isfield(obs(i), 'TF'))  || (~isgood(obs(i).TF))
            % Riley, J. P., Deep-Sea Research 12:219-220, 1965:
            % this is .000068.*Sali./35. = .00000195.*Sali
            obs(i).TF   = ( 0.000067 / 18.998 ) * ...
                ( obs(i).sal / 1.80655 )*1e6; % convt to µmol/kg-SW
            yobs(i,sys.iTF)     = p(obs(i).TF*1e-6);
        else
            if ((obs(i).TF) == 0)
                obs(i).TF       = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.iTF)     = p(obs(i).TF*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs(i), 'eTF'))  || (~isgood(obs(i).eTF))
            % 6.7 ± 0.1 e-5
            TFu = ( ( (6.7e-5 + 0.1e-5)/18.998) * obs(i).sal/1.80655 );
            TFl = ( ( (6.7e-5 - 0.1e-5)/18.998) * obs(i).sal/1.80655 );
            eTF                 = (TFu - TFl) / 2;
            obs(i).eTF          = (eTF)*1e6;
            wobs(i,sys.iTF)     = w(obs(i).TF,obs(i).eTF);
        else
            wobs(i,sys.iTF)     = w(obs(i).TF,obs(i).eTF);
        end
        if (~isfield(obs(i).tp(1), 'F')) || (~isgood(obs(i).tp(1).F))
            obs(i).tp(1).F      = [];
        end
        if (~isfield(obs(i).tp(1), 'eF')) || (~isgood(obs(i).tp(1).eF))
            obs(i).tp(1).eF     = [];
        end
        if (~isfield(obs(i).tp(1), 'HF')) || (~isgood(obs(i).tp(1).HF))
            obs(i).tp(1).HF     = [];
        end
        if (~isfield(obs(i).tp(1), 'eHF')) || (~isgood(obs(i).tp(1).eHF))
            obs(i).tp(1).eHF    = [];
        end

        if (ismember('phosphate',opt.abr))
            if (~isfield(obs(i).tp(1), 'epK1p'))
                obs(i).tp(1).epK1p  = [];
            end
            if (~isfield(obs(i).tp(1), 'pK1p'))
                obs(i).tp(1).pK1p   = [];
            end
            if (~isfield(obs(i).tp(1), 'epK2p'))
                obs(i).tp(1).epK2p  = [];
            end
            if (~isfield(obs(i).tp(1), 'pK2p'))
                obs(i).tp(1).pK2p   = [];
            end
            if (~isfield(obs(i).tp(1), 'epK3p'))
                obs(i).tp(1).epK3p  = [];
            end
            if (~isfield(obs(i).tp(1), 'pK3p'))
                obs(i).tp(1).pK3p   = [];
            end
            if (~isfield(obs(i), 'TP'))  || (~isgood(obs(i).TP))
                obs(i).TP = 1e-3; % µmol/kg
                yobs(i,sys.iTP)     = p(obs(i).TP*1e-6); % convt µmol/kg to mol/kg
            else
                if ((obs(i).TP) == 0) % zero po4 is very unlikely and breaks the code
                    obs(i).TP       = 1e-3; % umol/kg-SW, reset minimum to 1 nanomolar
                end
                yobs(i,sys.iTP)     = p(obs(i).TP*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i), 'eTP'))  || (~isgood(obs(i).eTP))
                obs(i).eTP          = 1e-3; % µmol/kg
                wobs(i,sys.iTP)     = w(obs(i).TP,obs(i).eTP);
            else
                if ((obs(i).eTP) == 0)
                    obs(i).eTP      = 1e-3; % umol/kg, reset minimum if zero
                end
                wobs(i,sys.iTP)     = w(obs(i).TP,obs(i).eTP);
            end
            if (~isfield(obs(i).tp(1), 'h3po4')) || ...
                    (~isgood(obs(i).tp(1).h3po4))
                obs(i).tp(1).h3po4  = [];
            end
            if (~isfield(obs(i).tp(1), 'eh3po4')) || ...
                    (~isgood(obs(i).tp(1).eh3po4))
                obs(i).tp(1).eh3po4 = [];
            end
            if (~isfield(obs(i).tp(1), 'h2po4')) || ...
                    (~isgood(obs(i).tp(1).h2po4))
                obs(i).tp(1).h2po4  = [];
            end
            if (~isfield(obs(i).tp(1), 'eh2po4')) || ...
                    (~isgood(obs(i).tp(1).eh2po4))
                obs(i).tp(1).eh2po4 = [];
            end
            if (~isfield(obs(i).tp(1), 'hpo4')) || ...
                    (~isgood(obs(i).tp(1).hpo4))
                obs(i).tp(1).hpo4   = [];
            end
            if (~isfield(obs(i).tp(1), 'ehpo4')) || ...
                    (~isgood(obs(i).tp(1).ehpo4))
                obs(i).tp(1).ehpo4  = [];
            end
            if (~isfield(obs(i).tp(1), 'po4')) || ...
                    (~isgood(obs(i).tp(1).po4))
                obs(i).tp(1).po4    = [];
            end
            if (~isfield(obs(i).tp(1), 'epo4')) || ...
                    (~isgood(obs(i).tp(1).epo4))
                obs(i).tp(1).epo4   = [];
            end
        end

        if (ismember('silicate',opt.abr))
            if (~isfield(obs(i).tp(1), 'epKsi'))
                obs(i).tp(1).epKsi  = [];
            end
            if (~isfield(obs(i).tp(1), 'pKsi'))
                obs(i).tp(1).pKsi   = [];
            end
            if (~isfield(obs(i), 'TSi'))  || (~isgood(obs(i).TSi))
                obs(i).TSi          = 1e-3; % µmol/kg
                yobs(i,sys.iTSi)    = p(obs(i).TSi*1e-6); % convt µmol/kg to mol/kg
            else
                if ((obs(i).TSi) == 0) % zero silicate very unlikely and breaks code
                    obs(i).TSi      = 1e-3; % umol/kg, reset minimum to 1 nanomolar
                end
                yobs(i,sys.iTSi)    = p(obs(i).TSi*1e-6);
            end
            if (~isfield(obs(i), 'eTSi'))  || (~isgood(obs(i).eTSi))
                obs(i).eTSi         = 1e-3; % µmol/kg
                wobs(i,sys.iTSi)    = w(obs(i).TSi,obs(i).eTSi);
            else
                if ((obs(i).eTSi) == 0)
                    obs(i).eTSi     = 1e-3; % umol/kg, reset minimum to 1 nanomolar
                end
                wobs(i,sys.iTSi)    = w(obs(i).TSi,obs(i).eTSi);
            end
            if (~isfield(obs(i).tp(1), 'sioh4')) || ...
                    (~isgood(obs(i).tp(1).sioh4))
                obs(i).tp(1).sioh4      = [];
            end
            if (~isfield(obs(i).tp(1), 'esioh4')) || ...
                    (~isgood(obs(i).tp(1).esioh4))
                obs(i).tp(1).esioh4     = [];
            end
            if (~isfield(obs(i).tp(1), 'siooh3')) || ...
                    (~isgood(obs(i).tp(1).siooh3))
                obs(i).tp(1).siooh3     = [];
            end
            if (~isfield(obs(i).tp(1), 'esiooh3')) || ...
                    (~isgood(obs(i).tp(1).esiooh3))
                obs(i).tp(1).esiooh3    = [];
            end
        end

        if (ismember('ammonia',opt.abr))
            if (~isfield(obs(i).tp(1), 'epKnh4'))
                obs(i).tp(1).epKnh4     = [];
            end
            if (~isfield(obs(i).tp(1), 'pKnh4'))
                obs(i).tp(1).pKnh4      = [];
            end
            if (~isfield(obs(i), 'TNH3'))  || (~isgood(obs(i).TNH3))
                obs(i).TNH3         = 1e-3; % µmol/kg
                yobs(i,sys.iTNH3)   = p(obs(i).TNH3*1e-6); % convt µmol/kg to mol/kg
            else
                if ((obs(i).TNH3) == 0)
                    obs(i).TNH3     = 1e-3; % umol/kg, reset minimum to 1 nanomolar
                end
                yobs(i,sys.iTNH3)   = p(obs(i).TNH3*1e-6);
            end
            if (~isfield(obs(i), 'eTNH3'))  || (~isgood(obs(i).eTNH3))
                obs(i).eTNH3        = 5e-4; % µmol/kg
                wobs(i,sys.iTNH3)   = w(obs(i).TNH3,obs(i).eTNH3);
            else
                wobs(i,sys.iTNH3)   = w(obs(i).TNH3,obs(i).eTNH3);
            end
            if (~isfield(obs(i).tp(1), 'nh3')) || ...
                    (~isgood(obs(i).tp(1).nh3))
                obs(i).tp(1).nh3    = [];
            end
            if (~isfield(obs(i).tp(1), 'enh3')) || ...
                    (~isgood(obs(i).tp(1).enh3))
                obs(i).tp(1).enh3   = [];
            end
            if (~isfield(obs(i).tp(1), 'nh4')) || ...
                    (~isgood(obs(i).tp(1).nh4))
                obs(i).tp(1).nh4    = [];
            end
            if (~isfield(obs(i).tp(1), 'enh4')) || ...
                    (~isgood(obs(i).tp(1).enh4))
                obs(i).tp(1).enh4   = [];
            end
        end

        if (ismember('sulfide',opt.abr))
            if (~isfield(obs(i).tp(1), 'epKh2s'))
                obs(i).tp(1).epKh2s = [];
            end
            if (~isfield(obs(i).tp(1), 'pKh2s'))
                obs(i).tp(1).pKh2s  = [];
            end
            if (~isfield(obs(i), 'TH2S'))  || (~isgood(obs(i).TH2S))
                obs(i).TH2S         = 1e-3; % µmol/kg
                yobs(i,sys.iTH2S)   = p(obs(i).TH2S*1e-6); % convt µmol/kg to mol/kg
            else
                if ((obs(i).TH2S) == 0)
                    obs(i).TH2S     = 1e-3; % umol/kg, reset minimum to 1 nanomolar
                end
                yobs(i,sys.iTH2S)   = p(obs(i).TH2S*1e-6);
            end
            if (~isfield(obs(i), 'eTH2S'))  || (~isgood(obs(i).eTH2S))
                obs(i).eTH2S        = 5e-4; % µmol/kg
                wobs(i,sys.iTH2S)   = w(obs(i).TH2S,obs(i).eTH2S);
            else
                wobs(i,sys.iTH2S)   = w(obs(i).TH2S,obs(i).eTH2S);
            end
            if (~isfield(obs(i).tp(1), 'hs')) || ...
                    (~isgood(obs(i).tp(1).hs))
                obs(i).tp(1).hs     = [];
            end
            if (~isfield(obs(i).tp(1), 'ehs')) || ...
                    (~isgood(obs(i).tp(1).ehs))
                obs(i).tp(1).ehs    = [];
            end
            if (~isfield(obs(i).tp(1), 'h2s')) || ...
                    (~isgood(obs(i).tp(1).h2s))
                obs(i).tp(1).h2s    = [];
            end
            if (~isfield(obs(i).tp(1), 'eh2s')) || ...
                    (~isgood(obs(i).tp(1).eh2s))
                obs(i).tp(1).eh2s   = [];
            end
        end

        if (ismember('solubility',opt.abr))
            if (~isfield(obs(i).tp(1), 'epKar'))
                obs(i).tp(1).epKar = [];
            end
            if (~isfield(obs(i).tp(1), 'pKar'))
                obs(i).tp(1).pKar   = [];
            end
            if (~isfield(obs(i), 'TCal'))  || (~isgood(obs(i).TCal))
                if opt.K1K2 == 6 || opt.K1K2 == 7
                    % Calculate Ca for GEOSECS, Riley and Skirrow 1965
                    obs(i).TCal     = (0.01026 .* obs(i).sal ./ 35) ;
                    yobs(i,sys.iTCal)   = p(obs(i).TCal); % mol/kg
                else
                    % Calculate Ca, Riley and Tongdui 1967
                    % this is 0.010285.*obs.sal./35;
                    obs(i).TCal = (0.02128./40.087.*(obs(i).sal./1.80655)) ; % stay on mol/kg
                    yobs(i,sys.iTCal)   = p(obs(i).TCal); % mol/kg
                end
            else
                if ((obs(i).TCal) == 0)
                    obs(i).TCal     = 1e-3; % mol/kg, reset minimum to 1 nanomolar
                end
                yobs(i,sys.iTCal)   = p(obs(i).TCal); % assume user input of mol/kg
            end
            if (~isfield(obs(i), 'eTCal'))  || (~isgood(obs(i).eTCal))
                obs(i).eTCal        = (6e-5); % mol/kg, from Riley and Tongdui 1967
                wobs(i,sys.iTCal)   = w(obs(i).TCal,obs(i).eTCal);
            else
                wobs(i,sys.iTCal)   = w(obs(i).TCal,obs(i).eTCal);
            end
            if (~isfield(obs(i).tp(1), 'ca')) || ...
                    (~isgood(obs(i).tp(1).ca))
                obs(i).tp(1).ca     = [];
            end
            if (~isfield(obs(i).tp(1), 'eca')) || ...
                    (~isgood(obs(i).tp(1).eca))
                obs(i).tp(1).eca    = [];
            end
            if (~isfield(obs(i).tp(1), 'OmegaAr')) || ...
                    (~isgood(obs(i).tp(1).OmegaAr))
                obs(i).tp(1).OmegaAr    = [];
            end
            if (~isfield(obs(i).tp(1), 'eOmegaAr')) || ...
                    (~isgood(obs(i).tp(1).eOmegaAr))
                obs(i).tp(1).eOmegaAr   = [];
            end
            if (~isfield(obs(i).tp(1), 'pKca'))
                obs(i).tp(1).pKca       = [];
            end
            if (~isfield(obs(i).tp(1), 'epKca'))
                obs(i).tp(1).epKca      = [];
            end
            if (~isfield(obs(i).tp(1), 'OmegaCa')) || ...
                    (~isgood(obs(i).tp(1).OmegaCa))
                obs(i).tp(1).OmegaCa    = [];
            end
            if (~isfield(obs(i).tp(1), 'eOmegaCa')) || ...
                    (~isgood(obs(i).tp(1).eOmegaCa))
                obs(i).tp(1).eOmegaCa   = [];
            end
        end

        nTP = length(obs(i).tp);

        for ii = 1:nTP % loop over all the pressure and temperature sensitive components
            yobs(i,sys.tp(ii).iT)   = obs(i).tp(ii).T;
            yobs(i,sys.tp(ii).iP)   = obs(i).tp(ii).P;

            wobs(i,sys.tp(ii).iT)   = (obs(i).tp(ii).eT)^(-2);
            wobs(i,sys.tp(ii).iP)   = (obs(i).tp(ii).eP)^(-2);

            [pK,gpK,epK]            = calc_pK(opt, obs(i).tp(ii).T, ...
                                        obs(i).sal, obs(i).tp(ii).P );
            
            pK0   = pK(1);      pK1  = pK(2);     pK2   = pK(3);  
            pKb   = pK(4);      pKw  = pK(5);     pKs   = pK(6);  
            pKf   = pK(7);      pK1p = pK(8);     pK2p  = pK(9);  
            pK3p  = pK(10);     pKsi = pK(11);    pKnh4 = pK(12);
            pKh2s = pK(13);     pp2f = pK(14);    pKar  = pK(15); 
            pKca  = pK(16);     pfH  = pK(17);

            epK0   = epK(1);    epK1  = epK(2);   epK2   = epK(3);  
            epKb   = epK(4);    epKw  = epK(5);   epKs   = epK(6);  
            epKf   = epK(7);    epK1p = epK(8);   epK2p  = epK(9);  
            epK3p  = epK(10);   epKsi = epK(11);  epKnh4 = epK(12);
            epKh2s = epK(13);   epp2f = epK(14);  epKar  = epK(15); 
            epKca  = epK(16);   epfH  = epK(17);
            
            % add "observations" for the equilibrium constants
            % and transfer from obs struct to yobs and wobs
            
            if (isgood(obs(i).tp(ii).epK0))
                wobs(i,sys.tp(ii).iK0)  = (obs(i).tp(ii).epK0)^(-2);
            else
                obs(i).tp(ii).epK0      = epK0;
                wobs(i,sys.tp(ii).iK0)  = (obs(i).tp(ii).epK0)^(-2);
            end
            if (isgood(obs(i).tp(ii).pK0))
                yobs(i,sys.tp(ii).iK0)  = obs(i).tp(ii).pK0;
            else
                yobs(i,sys.tp(ii).iK0)  = pK0;
                obs(i).tp(ii).pK0       = pK0;
            end

            if (isgood(obs(i).tp(ii).epK1))
                wobs(i,sys.tp(ii).iK1)  = (obs(i).tp(ii).epK1)^(-2);
            else
                obs(i).tp(ii).epK1      = epK1;
                wobs(i,sys.tp(ii).iK1)  = (obs(i).tp(ii).epK1)^(-2);
            end
            if (isgood(obs(i).tp(ii).pK1))
                yobs(i,sys.tp(ii).iK1)  = obs(i).tp(ii).pK1;
            else
                yobs(i,sys.tp(ii).iK1)  = pK1;
                obs(i).tp(ii).pK1       = pK1;
            end

            if (isgood(obs(i).tp(ii).epK2))
                wobs(i,sys.tp(ii).iK2)  = (obs(i).tp(ii).epK2)^(-2);
            else
                obs(i).tp(ii).epK2      = epK2;
                wobs(i,sys.tp(ii).iK2)  = (obs(i).tp(ii).epK2)^(-2);
            end
            if (isgood(obs(i).tp(ii).pK2))
                yobs(i,sys.tp(ii).iK2)  = obs(i).tp(ii).pK2;
            else
                yobs(i,sys.tp(ii).iK2)  = pK2;
                obs(i).tp(ii).pK2       = pK2;
            end

            if (isgood(obs(i).tp(ii).epp2f))
                wobs(i,sys.tp(ii).ip2f) = (obs(i).tp(ii).epp2f)^(-2);
            else
                obs(i).tp(ii).epp2f     = epp2f;
                wobs(i,sys.tp(ii).ip2f) = (obs(i).tp(ii).epp2f)^(-2);
            end
            if (isgood(obs(i).tp(ii).pp2f))
                yobs(i,sys.tp(ii).ip2f) = obs(i).tp(ii).pp2f;
            else
                yobs(i,sys.tp(ii).ip2f) = pp2f;
                obs(i).tp(ii).pp2f      = pp2f;
            end

            if (isgood(obs(i).tp(ii).co2st))
                yobs(i,sys.tp(ii).ico2st) = p(obs(i).tp(ii).co2st*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ico2st)   = nan;
                obs(i).tp(ii).co2st         = nan;
            end
            if (isgood(obs(i).tp(ii).eco2st))
                wobs(i,sys.tp(ii).ico2st)   = w(obs(i).tp(ii).co2st, ...
                                                obs(i).tp(ii).eco2st);
            else
                wobs(i,sys.tp(ii).ico2st)   = nan;
                obs(i).tp(ii).eco2st        = nan;
            end

            if (isgood(obs(i).tp(ii).hco3))
                yobs(i,sys.tp(ii).ihco3)    = p(obs(i).tp(ii).hco3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ihco3)    = nan;
                obs(i).tp(ii).hco3          = nan;
            end
            if (isgood(obs(i).tp(ii).ehco3))
                wobs(i,sys.tp(ii).ihco3)    = w(obs(i).tp(ii).hco3, ...
                                                obs(i).tp(ii).ehco3);
            else
                wobs(i,sys.tp(ii).ihco3)    = nan;
                obs(i).tp(ii).ehco3         = nan;
            end

            if (isgood(obs(i).tp(ii).co3))
                yobs(i,sys.tp(ii).ico3)     = p(obs(i).tp(ii).co3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ico3)     = nan;
                obs(i).tp(ii).co3           = nan;
            end
            if (isgood(obs(i).tp(ii).eco3))
                wobs(i,sys.tp(ii).ico3)     = w(obs(i).tp(ii).co3, ...
                                                obs(i).tp(ii).eco3);
            else
                wobs(i,sys.tp(ii).ico3)     = nan;
                obs(i).tp(ii).eco3          = nan;
            end

            if (isgood(obs(i).tp(ii).pco2))
                yobs(i,sys.tp(ii).ipco2)    = p(obs(i).tp(ii).pco2*1e-6); % convt µatm to atm
            else
                yobs(i,sys.tp(ii).ipco2)    = nan;
                obs(i).tp(ii).pco2          = nan;
            end
            if (isgood(obs(i).tp(ii).epco2))
                wobs(i,sys.tp(ii).ipco2)    = w(obs(i).tp(ii).pco2, ...
                                                obs(i).tp(ii).epco2);
            else
                wobs(i,sys.tp(ii).ipco2)    = nan;
                obs(i).tp(ii).epco2         = nan;
            end

            if (isgood(obs(i).tp(ii).fco2))
                yobs(i,sys.tp(ii).ifco2)    = p(obs(i).tp(ii).fco2*1e-6); % convt µatm to atm
            else
                yobs(i,sys.tp(ii).ifco2)    = nan;
                obs(i).tp(ii).fco2          = nan;
            end
            if (isgood(obs(i).tp(ii).efco2))
                wobs(i,sys.tp(ii).ifco2)    = w(obs(i).tp(ii).fco2, ...
                                                obs(i).tp(ii).efco2);
            else
                wobs(i,sys.tp(ii).ifco2)    = nan;
                obs(i).tp(ii).efco2         = nan;
            end
         
            if (isgood(obs(i).tp(ii).ph))
                yobs(i,sys.tp(ii).iph)      = obs(i).tp(ii).ph ;
            else
                yobs(i,sys.tp(ii).iph)      = nan;
                obs(i).tp(ii).ph            = nan;
            end
            if (isgood(obs(i).tp(ii).eph))
                wobs(i,sys.tp(ii).iph)      = (obs(i).tp(ii).eph).^(-2);
            else
                wobs(i,sys.tp(ii).iph)      = nan;
                obs(i).tp(ii).eph           = nan;
            end

            if (isgood(obs(i).tp(ii).ph_free)) % ph_free (ph on free scale)
                yobs(i,sys.tp(ii).iph_free) = obs(i).tp(ii).ph_free ;
            else
                yobs(i,sys.tp(ii).iph_free) = nan;
                obs(i).tp(ii).ph_free       = nan;
            end
            if (isgood(obs(i).tp(ii).eph_free)) 
                wobs(i,sys.tp(ii).iph_free) = obs(i).tp(ii).eph_free.^(-2);
            else
                wobs(i,sys.tp(ii).iph_free) = nan;
                obs(i).tp(ii).eph_free      = nan;
            end

            if (isgood(obs(i).tp(ii).ph))
                % calculate ph_tot, ph_sws, & ph_nbs via 'phscales' function
                ph_out      = phscales(obs(i).tp(ii).ph, opt.phscale, ...
                                p(obs(i).TS),pKs, p(obs(i).TF), pKf, pfH );
                % ph_out = phscales(phin,scalein,pTS,pKs,pTF,pKf,pfH);                
                obs(i).tp(ii).ph_tot    = ph_out(1);
                obs(i).tp(ii).ph_sws    = ph_out(2);
                obs(i).tp(ii).ph_nbs    = ph_out(3);
            else 
                obs(i).tp(ii).ph_tot    = nan;
                obs(i).tp(ii).ph_sws    = nan;
                obs(i).tp(ii).ph_nbs    = nan;
            end

            if (isgood(obs(i).tp(ii).eph))
                obs(i).tp(ii).eph_tot   = obs(i).tp(ii).eph;
                obs(i).tp(ii).eph_sws   = obs(i).tp(ii).eph;
                obs(i).tp(ii).eph_nbs   = obs(i).tp(ii).eph;
            else
                obs(i).tp(ii).eph_tot   = nan;
                obs(i).tp(ii).eph_sws   = nan;
                obs(i).tp(ii).eph_nbs   = nan;
            end

            if (isgood(obs(i).tp(ii).pfH)) % pfH activity coefficient
                yobs(i,sys.tp(ii).ipfH) = obs(i).tp(ii).pfH ;
            else
                yobs(i,sys.tp(ii).ipfH) = pfH;
                obs(i).tp(ii).pfH       = pfH;
            end
            if (isgood(obs(i).tp(ii).epfH)) % pfH activity coefficient
                wobs(i,sys.tp(ii).ipfH) = (obs(i).tp(ii).epfH).^(-2) ;
            else
                wobs(i,sys.tp(ii).ipfH) = (epfH).^(-2);
                obs(i).tp(ii).epfH      = epfH;
            end


            if (isgood(obs(i).tp(ii).epKw))
                wobs(i,sys.tp(ii).iKw)  = (obs(i).tp(ii).epKw).^(-2);
            else
                obs(i).tp(ii).epKw = epKw;
                wobs(i,sys.tp(ii).iKw)  = (obs(i).tp(ii).epKw).^(-2);
            end
            if (isgood(obs(i).tp(ii).pKw))
                yobs(i,sys.tp(ii).iKw)  = obs(i).tp(ii).pKw;
            else
                yobs(i,sys.tp(ii).iKw)  = pKw;
                obs(i).tp(ii).pKw       = pKw;
            end

            if (isgood(obs(i).tp(ii).oh))
                yobs(i,sys.tp(ii).ioh)  = p(obs(i).tp(ii).oh*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ioh)  = nan;
                obs(i).tp(ii).oh        = nan;
            end
            if (isgood(obs(i).tp(ii).eoh))
                wobs(i,sys.tp(ii).ioh)  = w(obs(i).tp(ii).oh, ...
                                            obs(i).tp(ii).eoh);
            else
                wobs(i,sys.tp(ii).ioh)  = nan;
                obs(i).tp(ii).eoh       = nan;
            end

            % borate system
            if (isgood(obs(i).tp(ii).epKb))
                wobs(i,sys.tp(ii).iKb)  = (obs(i).tp(ii).epKb).^(-2);
            else
                obs(i).tp(ii).epKb      = epKb;
                wobs(i,sys.tp(ii).iKb)  = (obs(i).tp(ii).epKb).^(-2);  % wKb = 1/(1 + (0.01/pKsys(4)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pKb))
                yobs(i,sys.tp(ii).iKb)  = obs(i).tp(ii).pKb;
            else
                yobs(i,sys.tp(ii).iKb)  = pKb; % from local_pK
                obs(i).tp(ii).pKb       = pKb;
            end

            if (isgood(obs(i).tp(ii).boh3))
                yobs(i,sys.tp(ii).iboh3) = p(obs(i).tp(ii).boh3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iboh3)    = nan;
                obs(i).tp(ii).boh3          = nan;
            end
            if (isgood(obs(i).tp(ii).eboh3))
                wobs(i,sys.tp(ii).iboh3)    = w(obs(i).tp(ii).boh3, ...
                                                obs(i).tp(ii).eboh3);
            else
                wobs(i,sys.tp(ii).iboh3)    = nan;
                obs(i).tp(ii).eboh3         = nan;
            end

            if (isgood(obs(i).tp(ii).boh4))
                yobs(i,sys.tp(ii).iboh4)    = p(obs(i).tp(ii).boh4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iboh4)    = nan;
                obs(i).tp(ii).boh4          = nan;
            end
            if (isgood(obs(i).tp(ii).eboh4))
                wobs(i,sys.tp(ii).iboh4)    = w(obs(i).tp(ii).boh4, ...
                                                    obs(i).tp(ii).eboh4);
            else
                wobs(i,sys.tp(ii).iboh4)    = nan;
                obs(i).tp(ii).eboh4         = nan;
            end

            % sulfate system
            if (isgood(obs(i).tp(ii).epKs))
                wobs(i,sys.tp(ii).iKs)      = (obs(i).tp(ii).epKs).^(-2);
            else
                obs(i).tp(ii).epKs          = epKs;
                wobs(i,sys.tp(ii).iKs)      = (obs(i).tp(ii).epKs).^(-2); % wKs = 1/(1 + (0.0021/pKsys(6)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pKs))
                yobs(i,sys.tp(ii).iKs)      = obs(i).tp(ii).pKs;
            else
                yobs(i,sys.tp(ii).iKs)      = pKs;
                obs(i).tp(ii).pKs           = pKs;
            end

            if (isgood(obs(i).tp(ii).hso4))
                yobs(i,sys.tp(ii).ihso4)    = p(obs(i).tp(ii).hso4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ihso4)    = nan;
                obs(i).tp(ii).hso4          = nan;
            end
            if (isgood(obs(i).tp(ii).ehso4))
                wobs(i,sys.tp(ii).ihso4)    = w(obs(i).tp(ii).hso4, ...
                                                obs(i).tp(ii).ehso4);
            else
                wobs(i,sys.tp(ii).ihso4)    = nan;
                obs(i).tp(ii).ehso4         = nan;
            end

            if (isgood(obs(i).tp(ii).so4))
                yobs(i,sys.tp(ii).iso4)     = p(obs(i).tp(ii).so4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iso4)     = nan;
                obs(i).tp(ii).so4           = nan;
            end
            if (isgood(obs(i).tp(ii).eso4))
                wobs(i,sys.tp(ii).iso4)     = w(obs(i).tp(ii).so4, ...
                                                    obs(i).tp(ii).eso4);
            else
                wobs(i,sys.tp(ii).iso4)     = nan;
                obs(i).tp(ii).eso4          = nan;
            end

            % fluoride system
            if (isgood(obs(i).tp(ii).epKf))
                wobs(i,sys.tp(ii).iKf)      = (obs(i).tp(ii).epKf).^(-2);
            else
                obs(i).tp(ii).epKf          = epKf;
                wobs(i,sys.tp(ii).iKf)      = (obs(i).tp(ii).epKf).^(-2);   % wKF = 1/(p(1 + 0.02/KF))^2 ; % 2% relative error
            end
            if (isgood(obs(i).tp(ii).pKf))
                yobs(i,sys.tp(ii).iKf)      = obs(i).tp(ii).pKf;
            else
                yobs(i,sys.tp(ii).iKf)      = pKf;
                obs(i).tp(ii).pKf           = pKf;
            end

            if (isgood(obs(i).tp(ii).F))
                yobs(i,sys.tp(ii).iF)       = p(obs(i).tp(ii).F*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iF)       = nan;
                obs(i).tp(ii).F             = nan;
            end
            if (isgood(obs(i).tp(ii).eF))
                wobs(i,sys.tp(ii).iF)       = w(obs(i).tp(ii).F, ...
                                                obs(i).tp(ii).eF);
            else
                wobs(i,sys.tp(ii).iF)       = nan;
                obs(i).tp(ii).eF            = nan;
            end

            if (isgood(obs(i).tp(ii).HF))
                yobs(i,sys.tp(ii).iHF)      = p(obs(i).tp(ii).HF*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iHF)      = nan;
                obs(i).tp(ii).HF            = nan;
            end
            if (isgood(obs(i).tp(ii).eHF))
                wobs(i,sys.tp(ii).iHF)      = w(obs(i).tp(ii).HF, ...
                                                obs(i).tp(ii).eHF);
            else
                wobs(i,sys.tp(ii).iHF)      = nan;
                obs(i).tp(ii).eHF           = nan;
            end

            if (ismember('phosphate',opt.abr))
                if (isgood(obs(i).tp(ii).epK1p))
                    wobs(i,sys.tp(ii).iK1p) = (obs(i).tp(ii).epK1p).^(-2);
                else
                    obs(i).tp(ii).epK1p     = epK1p;
                    wobs(i,sys.tp(ii).iK1p) = (obs(i).tp(ii).epK1p).^(-2);  % wK1p = 1/(1 + (0.09/pKsys(8)))^2 ;
                end
                if (isgood(obs(i).tp(ii).pK1p))
                    yobs(i,sys.tp(ii).iK1p) = obs(i).tp(ii).pK1p;
                else
                    yobs(i,sys.tp(ii).iK1p) = pK1p;
                    obs(i).tp(ii).pK1p      = pK1p;
                end

                if (isgood(obs(i).tp(ii).epK2p))
                    wobs(i,sys.tp(ii).iK2p) = (obs(i).tp(ii).epK2p).^(-2);
                else
                    obs(i).tp(ii).epK2p     = epK2p;
                    wobs(i,sys.tp(ii).iK2p) = (obs(i).tp(ii).epK2p).^(-2);  % wK2p = 1/(1 + (0.03/pKsys(9)))^2 ;
                end
                if (isgood(obs(i).tp(ii).pK2p))
                    yobs(i,sys.tp(ii).iK2p) = obs(i).tp(ii).pK2p;
                else
                    yobs(i,sys.tp(ii).iK2p) = pK2p;
                    obs(i).tp(ii).pK2p      = pK2p;
                end

                if (isgood(obs(i).tp(ii).epK3p))
                    wobs(i,sys.tp(ii).iK3p) = (obs(i).tp(ii).epK3p).^(-2);
                else
                    obs(i).tp(ii).epK3p     = epK3p;
                    wobs(i,sys.tp(ii).iK3p) = (obs(i).tp(ii).epK3p).^(-2);  % wK3p = 1/(1 + (0.02/pKsys(10)))^2 ;
                end
                if (isgood(obs(i).tp(ii).pK3p))
                    yobs(i,sys.tp(ii).iK3p) = obs(i).tp(ii).pK3p;
                else
                    yobs(i,sys.tp(ii).iK3p) = pK3p;
                    obs(i).tp(ii).pK3p      = pK3p;
                end

                if (isgood(obs(i).tp(ii).h3po4))
                    yobs(i,sys.tp(ii).ih3po4) = p(obs(i).tp(ii).h3po4*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.tp(ii).ih3po4)   = nan;
                    obs(i).tp(ii).h3po4         = nan;
                end
                if (isgood(obs(i).tp(ii).eh3po4))
                    wobs(i,sys.tp(ii).ih3po4)   = w(obs(i).tp(ii).h3po4, ...
                                                    obs(i).tp(ii).eh3po4);
                else
                    wobs(i,sys.tp(ii).ih3po4)   = nan;
                    obs(i).tp(ii).eh3po4        = nan;
                end

                if (isgood(obs(i).tp(ii).h2po4))
                    yobs(i,sys.tp(ii).ih2po4) = p(obs(i).tp(ii).h2po4*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.tp(ii).ih2po4)   = nan;
                    obs(i).tp(ii).h2po4         = nan;
                end
                if (isgood(obs(i).tp(ii).eh2po4))
                    wobs(i,sys.tp(ii).ih2po4)   = w(obs(i).tp(ii).h2po4, ...
                                                    obs(i).tp(ii).eh2po4);
                else
                    wobs(i,sys.tp(ii).ih2po4)   = nan;
                    obs(i).tp(ii).eh2po4        = nan;
                end

                if (isgood(obs(i).tp(ii).hpo4))
                    yobs(i,sys.tp(ii).ihpo4) = p(obs(i).tp(ii).hpo4*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.tp(ii).ihpo4)    = nan;
                    obs(i).tp(ii).hpo4          = nan;
                end
                if (isgood(obs(i).tp(ii).ehpo4))
                    wobs(i,sys.tp(ii).ihpo4)    = w(obs(i).tp(ii).hpo4, ...
                                                    obs(i).tp(ii).ehpo4);
                else
                    wobs(i,sys.tp(ii).ihpo4)    = nan;
                    obs(i).tp(ii).ehpo4         = nan;
                end

                if (isgood(obs(i).tp(ii).po4))
                    yobs(i,sys.tp(ii).ipo4) = p(obs(i).tp(ii).po4*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.tp(ii).ipo4)     = nan;
                    obs(i).tp(ii).po4           = nan;
                end
                if (isgood(obs(i).tp(ii).epo4))
                    wobs(i,sys.tp(ii).ipo4)     = w(obs(i).tp(ii).po4, ...
                                                    obs(i).tp(ii).epo4);
                else
                    wobs(i,sys.tp(ii).ipo4)     = nan;
                    obs(i).tp(ii).epo4          = nan;
                end

            end

            if (ismember('silicate',opt.abr))
                if (isgood(obs(i).tp(ii).epKsi))
                    wobs(i,sys.tp(ii).iKsi) = (obs(i).tp(ii).epKsi).^(-2);
                else
                    obs(i).tp(ii).epKsi     = epKsi;
                    wobs(i,sys.tp(ii).iKsi) = (obs(i).tp(ii).epKsi).^(-2);  % wKSi = 1/(1 + (0.02/pKsys(11)))^2 ;
                end
                if (isgood(obs(i).tp(ii).pKsi))
                    yobs(i,sys.tp(ii).iKsi) = obs(i).tp(ii).pKsi;
                else
                    yobs(i,sys.tp(ii).iKsi) = pKsi;
                    obs(i).tp(ii).pKsi      = pKsi;
                end

                if (isgood(obs(i).tp(ii).sioh4))
                    yobs(i,sys.tp(ii).isioh4) = p(obs(i).tp(ii).sioh4*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.tp(ii).isioh4)   = nan;
                    obs(i).tp(ii).sioh4         = nan;
                end
                if (isgood(obs(i).tp(ii).esioh4))
                    wobs(i,sys.tp(ii).isioh4)   = w(obs(i).tp(ii).sioh4, ...
                                                    obs(i).tp(ii).esioh4);
                else
                    wobs(i,sys.tp(ii).isioh4)   = nan;
                    obs(i).tp(ii).esioh4        = nan;
                end

                if (isgood(obs(i).tp(ii).siooh3))
                    yobs(i,sys.tp(ii).isiooh3) = p(obs(i).tp(ii).siooh3*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.tp(ii).isiooh3)  = nan;
                    obs(i).tp(ii).siooh3        = nan;
                end
                if (isgood(obs(i).tp(ii).esiooh3))
                    wobs(i,sys.tp(ii).isiooh3)  = w(obs(i).tp(ii).siooh3, ...
                                                    obs(i).tp(ii).esiooh3);
                else
                    wobs(i,sys.tp(ii).isiooh3)  = nan;
                    obs(i).tp(ii).esiooh3       = nan;
                end

            end

            if (ismember('ammonia',opt.abr))
                if (isgood(obs(i).tp(ii).epKnh4))
                    wobs(i,sys.tp(ii).iKnh4) = (obs(i).tp(ii).epKnh4).^(-2);
                else
                    obs(i).tp(ii).epKnh4     = epKnh4;
                    wobs(i,sys.tp(ii).iKnh4) = (obs(i).tp(ii).epKnh4).^(-2);  % wKnh4 = 1/(1 + (0.00017/pKsys(11)))^2 ;
                end
                if (isgood(obs(i).tp(ii).pKnh4))
                    yobs(i,sys.tp(ii).iKnh4) = obs(i).tp(ii).pKnh4;
                else
                    yobs(i,sys.tp(ii).iKnh4)    = pKnh4;
                    obs(i).tp(ii).pKnh4         = pKnh4;
                end

                if (isgood(obs(i).tp(ii).nh4))
                    yobs(i,sys.tp(ii).inh4) = p(obs(i).tp(ii).nh4*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.tp(ii).inh4)     = nan;
                    obs(i).tp(ii).nh4           = nan;
                end
                if (isgood(obs(i).tp(ii).enh4))
                    wobs(i,sys.tp(ii).inh4)     = w(obs(i).tp(ii).nh4, ...
                                                    obs(i).tp(ii).enh4);
                else
                    wobs(i,sys.tp(ii).inh4)     = nan;
                    obs(i).tp(ii).enh4          = nan;
                end

                if (isgood(obs(i).tp(ii).nh3))
                    yobs(i,sys.tp(ii).inh3) = p(obs(i).tp(ii).nh3*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.tp(ii).inh3)     = nan;
                    obs(i).tp(ii).nh3           = nan;
                end
                if (isgood(obs(i).tp(ii).enh3))
                    wobs(i,sys.tp(ii).inh3)     = w(obs(i).tp(ii).nh3, ...
                                                    obs(i).tp(ii).enh3);
                else
                    wobs(i,sys.tp(ii).inh3)     = nan;
                    obs(i).tp(ii).enh3          = nan;
                end
            end

            if (ismember('sulfide',opt.abr))
                if (isgood(obs(i).tp(ii).epKh2s))
                    wobs(i,sys.tp(ii).iKh2s) = (obs(i).tp(ii).epKh2s).^(-2);
                else
                    obs(i).tp(ii).epKh2s     = epKh2s;
                    wobs(i,sys.tp(ii).iKh2s) = (obs(i).tp(ii).epKh2s).^(-2);  % wKh2s = 1/(1 + (0.033/pKsys(11)))^2 ;
                end
                if (isgood(obs(i).tp(ii).pKh2s))
                    yobs(i,sys.tp(ii).iKh2s)    = obs(i).tp(ii).pKh2s;
                else
                    yobs(i,sys.tp(ii).iKh2s)    = pKh2s;
                    obs(i).tp(ii).pKh2s         = pKh2s;
                end

                if (isgood(obs(i).tp(ii).h2s))
                    yobs(i,sys.tp(ii).ih2s) = p(obs(i).tp(ii).h2s*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.tp(ii).ih2s) = nan;
                    obs(i).tp(ii).h2s       = nan;
                end
                if (isgood(obs(i).tp(ii).eh2s))
                    wobs(i,sys.tp(ii).ih2s) = w(obs(i).tp(ii).h2s, ...
                                                obs(i).tp(ii).eh2s);
                else
                    wobs(i,sys.tp(ii).ih2s) = nan;
                    obs(i).tp(ii).eh2s      = nan;
                end

                if (isgood(obs(i).tp(ii).hs))
                    yobs(i,sys.tp(ii).ihs)  = p(obs(i).tp(ii).hs*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.tp(ii).ihs)  = nan;
                    obs(i).tp(ii).hs        = nan;
                end
                if (isgood(obs(i).tp(ii).ehs))
                    wobs(i,sys.tp(ii).ihs)  = w(obs(i).tp(ii).hs, ...
                                                obs(i).tp(ii).ehs);
                else
                    wobs(i,sys.tp(ii).ihs)  = nan;
                    obs(i).tp(ii).ehs       = nan;
                end
            end

            if (ismember('solubility',opt.abr))
                if (isgood(obs(i).tp(ii).epKar))
                    wobs(i,sys.tp(ii).iKar) = (obs(i).tp(ii).epKar).^(-2);
                else
                    obs(i).tp(ii).epKar     = epKar;
                    wobs(i,sys.tp(ii).iKar) = (obs(i).tp(ii).epKar).^(-2);
                end
                if (isgood(obs(i).tp(ii).pKar))
                    yobs(i,sys.tp(ii).iKar) = obs(i).tp(ii).pKar;
                else
                    yobs(i,sys.tp(ii).iKar) = pKar;
                    obs(i).tp(ii).pKar      = pKar;
                end

                if (isgood(obs(i).tp(ii).ca))
                    yobs(i,sys.tp(ii).ica)  = p(obs(i).tp(ii).ca*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.tp(ii).ica)  = nan;
                    obs(i).tp(ii).ca        = nan;
                end
                if (isgood(obs(i).tp(ii).eca))
                    wobs(i,sys.tp(ii).ica)  = w(obs(i).tp(ii).ca, ...
                                                obs(i).tp(ii).eca);
                else
                    wobs(i,sys.tp(ii).ica)  = nan;
                    obs(i).tp(ii).eca       = nan;
                end

                if (isgood(obs(i).tp(ii).OmegaAr))
                    yobs(i,sys.tp(ii).iOmegaAr) = p(obs(i).tp(ii).OmegaAr); % Omega is dimensionless
                else
                    yobs(i,sys.tp(ii).iOmegaAr) = nan;
                    obs(i).tp(ii).OmegaAr       = nan;
                end
                if (isgood(obs(i).tp(ii).eOmegaAr))
                    wobs(i,sys.tp(ii).iOmegaAr) = w(obs(i).tp(ii).OmegaAr, ...
                                                    obs(i).tp(ii).eOmegaAr);
                else
                    wobs(i,sys.tp(ii).iOmegaAr) = nan;
                    obs(i).tp(ii).eOmegaAr      = nan;
                end
                
                if (isgood(obs(i).tp(ii).epKca))
                    wobs(i,sys.tp(ii).iKca) = (obs(i).tp(ii).epKca).^(-2);
                else
                    obs(i).tp(ii).epKca     = epKca;
                    wobs(i,sys.tp(ii).iKca) = (obs(i).tp(ii).epKca).^(-2);
                end
                if (isgood(obs(i).tp(ii).pKca))
                    yobs(i,sys.tp(ii).iKca) = obs(i).tp(ii).pKca;
                else
                    yobs(i,sys.tp(ii).iKca) = pKca;
                    obs(i).tp(ii).pKca      = pKca;
                end

                if (isgood(obs(i).tp(ii).OmegaCa))
                    yobs(i,sys.tp(ii).iOmegaCa) = p(obs(i).tp(ii).OmegaCa); % Omega is dimensionless
                else
                    yobs(i,sys.tp(ii).iOmegaCa) = nan;
                    obs(i).tp(ii).OmegaCa       = nan;
                end
                if (isgood(obs(i).tp(ii).eOmegaCa))
                    wobs(i,sys.tp(ii).iOmegaCa) = w(obs(i).tp(ii).OmegaCa, ...
                                                    obs(i).tp(ii).eOmegaCa);
                else
                    wobs(i,sys.tp(ii).iOmegaCa) = nan;
                    obs(i).tp(ii).eOmegaCa      = nan;
                end
            end
        end
    end
end

% --------------------------------------------------------------------------------

function [est] = parse_output(z,sigy,opt,sys)
    % populate est, output structure with best estimates
    %
    % INPUT:
    %   z  := system state including equilibrium constants and lagrange multipliers
    
    p       = sys.p;
    q       = sys.q;

    ebar    = @(j) (0.5 * ( q( z(j) - sigy(j) ) - q( z(j) + sigy(j) ) ) );
    ebar_l  = @(j) ( q( -sigy(j) ) ); % lower sigma
    ebar_u  = @(j) ( q( sigy(j) ) ); % upper sigma
        
    % populate 'est' structure with best estimate:
    %   1. p(value) and p(error) where p(x) = -log10(x)
    %   2. value and average error about the value in 'q' 
    %           where q(x) = x^(-10)
    %   3. upper and lower bounds in 'q' space, not symmetric
    %           about the value in 'q' space

    est.sal     = z(sys.isal);
    est.esal    = sigy(sys.isal);
    
    % TC
    est.pTC     = z(sys.iTC);               
    est.epTC    = sigy(sys.iTC);    
    est.TC      = q(z(sys.iTC))*1e6; % 1e6  converts mol/kg to µmol/kg
    est.eTC     = ebar(sys.iTC)*1e6;      
    est.eTC_l   = ebar_l(sys.iTC)*1e6;    
    est.eTC_u   = ebar_u(sys.iTC)*1e6;
    
    % TA
    est.pTA     = z(sys.iTA);               
    est.epTA    = sigy(sys.iTA);
    est.TA      = q(z(sys.iTA))*1e6; % convt 
    est.eTA     = ebar(sys.iTA)*1e6;      
    est.eTA_l   = ebar_l(sys.iTA)*1e6;    
    est.eTA_u   = ebar_u(sys.iTA)*1e6;
    
    % TB borate
    est.pTB     = z(sys.iTB);
    est.epTB    = sigy(sys.iTB);
    est.TB      = q(z(sys.iTB))*1e6; % convt mol/kg to µmol/kg
    est.eTB     = ebar(sys.iTB)*1e6;
    est.eTB_l   = ebar_l(sys.iTB)*1e6;
    est.eTB_u   = ebar_u(sys.iTB)*1e6;

    % TS sulfate
    est.pTS     = z(sys.iTS);
    est.epTS    = sigy(sys.iTS);
    est.TS      = q(z(sys.iTS))*1e6;  % convt mol/kg to µmol/kg
    est.eTS     = ebar(sys.iTS)*1e6;
    est.eTS_l   = ebar_l(sys.iTS)*1e6;
    est.eTS_u   = ebar_u(sys.iTS)*1e6;

    % TF fluoride
    est.pTF     = z(sys.iTF);
    est.epTF        = sigy(sys.iTF);
    est.TF      = q(z(sys.iTF))*1e6; % convt mol/kg to µmol/kg
    est.eTF     = ebar(sys.iTF)*1e6;
    est.eTF_l   = ebar_l(sys.iTF)*1e6;
    est.eTF_u   = ebar_u(sys.iTF)*1e6;
    
    if ismember('phosphate', opt.abr)
        % TP
        est.pTP     = z(sys.iTP);           
        est.epTP    = sigy(sys.iTP);
        est.TP      = q(z(sys.iTP))*1e6;   % convt mol/kg to µmol/kg   
        est.eTP     = ebar(sys.iTP)*1e6;
        est.eTP_l   = ebar_l(sys.iTP)*1e6; 
        est.eTP_u   = ebar_u(sys.iTP)*1e6;
    end
    if ismember('silicate', opt.abr)
        % TSi
        est.pTSi    = z(sys.iTSi);         
        est.epTSi   = sigy(sys.iTSi);
        est.TSi     = q(z(sys.iTSi))*1e6;  % convt mol/kg to µmol/kg 
        est.eTSi    = ebar(sys.iTSi)*1e6;
        est.eTSi_l  = ebar_l(sys.iTSi)*1e6; 
        est.eTSi_u  = ebar_u(sys.iTSi)*1e6;
    end
    if ismember('ammonia', opt.abr)
        % TNH3
        est.pTNH3       = z(sys.iTNH3);       
        est.epTNH3      = sigy(sys.iTNH3);
        est.TNH3        = q(z(sys.iTNH3))*1e6;   % convt mol/kg to µmol/kg
        est.eTNH3       = ebar(sys.iTNH3)*1e6;
        est.eTNH3_l     = ebar_l(sys.iTNH3)*1e6; 
        est.eTNH3_u     = ebar_u(sys.iTNH3)*1e6;
    end
    if ismember('sulfide', opt.abr)
        % TH2S
        est.pTH2S       = z(sys.iTH2S);       
        est.epTH2S      = sigy(sys.iTH2S);
        est.TH2S        = q(z(sys.iTH2S))*1e6; % convt mol/kg to µmol/kg
        est.eTH2S       = ebar(sys.iTH2S)*1e6;
        est.eTH2S_l     = ebar_l(sys.iTH2S)*1e6; 
        est.eTH2S_u     = ebar_u(sys.iTH2S)*1e6;
    end
    if ismember('solubility', opt.abr)
        % TCal
        est.pTCal       = z(sys.iTCal);       
        est.epTCal      = sigy(sys.iTCal);
        est.TCal        = q(z(sys.iTCal)); % no conversion, mol/kg
        est.eTCal       = ebar(sys.iTCal);
        est.eTCal_l     = ebar_l(sys.iTCal); 
        est.eTCal_u     = ebar_u(sys.iTCal);
    end
    
    nTP = length(sys.tp);
    for i = 1:nTP
        % temp (deg C)
        est.tp(i).T     = z(sys.tp(i).iT);
        est.tp(i).eT    = sigy(sys.tp(i).iT);
        est.tp(i).eT_l  = z(sys.tp(i).iT)-sigy(sys.tp(i).iT);
        est.tp(i).eT_u  = z(sys.tp(i).iT)+sigy(sys.tp(i).iT);
        % pressure (dbar)
        est.tp(i).P     = z(sys.tp(i).iP);        
        est.tp(i).eP    = sigy(sys.tp(i).iP);
        est.tp(i).eP_l  = z(sys.tp(i).iP)-sigy(sys.tp(i).iP);
        est.tp(i).eP_u  = z(sys.tp(i).iP)+sigy(sys.tp(i).iP);
       
        % pH
        est.tp(i).ph    = z(sys.tp(i).iph);
        est.tp(i).eph   = sigy(sys.tp(i).iph);
        est.tp(i).h_    = q(z(sys.tp(i).iph)) * 1e6;
        est.tp(i).eh_   = ebar(sys.tp(i).iph) * 1e6;
        est.tp(i).eh_l  = ebar_l(sys.tp(i).iph) * 1e6;
        est.tp(i).eh_u  = ebar_u(sys.tp(i).iph) * 1e6;
        % ph_free
        est.tp(i).ph_free   = z(sys.tp(i).iph_free);
        est.tp(i).eph_free  = sigy(sys.tp(i).iph_free);
        est.tp(i).h_free    = q(z(sys.tp(i).iph_free)) * 1e6; % H (free) = q(ph_free)
        est.tp(i).eh_free   = ebar(sys.tp(i).iph_free) * 1e6;
        est.tp(i).eh_free_l = ebar_l(sys.tp(i).iph_free) * 1e6;
        est.tp(i).eh_free_u = ebar_u(sys.tp(i).iph_free) * 1e6;
        % ph_tot, ph_sws, & ph_nbs calculate via 'phscales' function
        ph_out  = phscales(est.tp(i).ph, opt.phscale, est.pTS, ...
                            z(sys.tp(i).iKs), est.pTF, ...
                            z(sys.tp(i).iKf), z(sys.tp(i).ipfH) );
            % ph_out = phscales(phin,scalein,pTS,pKs,pTF,pKf,pfH)
            % ph_out = [ph_tot, ph_sws, ph_nbs];
        % ph_tot
        est.tp(i).ph_tot    = ph_out(1);
        est.tp(i).eph_tot   = sigy(sys.tp(i).iph); % same as iph
        est.tp(i).h_tot     = q(ph_out(1)) * 1e6; % H (tot) = q(ph_tot)
        est.tp(i).eh_tot    = ebar(sys.tp(i).iph) * 1e6; % same as iph
        est.tp(i).eh_tot_l  = ebar_l(sys.tp(i).iph) * 1e6;
        est.tp(i).eh_tot_u  = ebar_u(sys.tp(i).iph) * 1e6;
        % ph_sws
        est.tp(i).ph_sws    = ph_out(2);
        est.tp(i).eph_sws   = sigy(sys.tp(i).iph); % same as iph
        est.tp(i).h_sws     = q(ph_out(2)) * 1e6; % H (sws) = q(ph_sws)
        est.tp(i).eh_sws    = ebar(sys.tp(i).iph) * 1e6; % same as iph
        est.tp(i).eh_sws_l  = ebar_l(sys.tp(i).iph) * 1e6;
        est.tp(i).eh_sws_u  = ebar_u(sys.tp(i).iph) * 1e6;
        % ph_nbs
        est.tp(i).ph_nbs    = ph_out(3);
        est.tp(i).eph_nbs   = sigy(sys.tp(i).iph); % same as iph
        est.tp(i).h_nbs     = q(ph_out(3)) * 1e6; % H (nbs) = q(ph_nbs)
        est.tp(i).eh_nbs    = ebar(sys.tp(i).iph) * 1e6; % same as iph
        est.tp(i).eh_nbs_l  = ebar_l(sys.tp(i).iph) * 1e6;
        est.tp(i).eh_nbs_u  = ebar_u(sys.tp(i).iph) * 1e6;

        % fH = activity coefficient
        est.tp(i).pfH       = z(sys.tp(i).ipfH);
        est.tp(i).epfH      = sigy(sys.tp(i).ipfH);
        est.tp(i).fH        = q(z(sys.tp(i).ipfH))*1e6;
        est.tp(i).efH       = ebar(sys.tp(i).ipfH)*1e6;
        est.tp(i).efH_l     = ebar_l(sys.tp(i).ipfH)*1e6;
        est.tp(i).efH_u     = ebar_u(sys.tp(i).ipfH)*1e6;
        % fCO2
        est.tp(i).fco2      = q(z(sys.tp(i).ifco2))*1e6; % convt atm to µatm
        est.tp(i).efco2     = ebar(sys.tp(i).ifco2)*1e6;
        est.tp(i).efco2_l   = ebar_l(sys.tp(i).ifco2)*1e6;
        est.tp(i).efco2_u   = ebar_u(sys.tp(i).ifco2)*1e6;
        est.tp(i).pfco2     = z(sys.tp(i).ifco2);
        est.tp(i).epfco2    = sigy(sys.tp(i).ifco2);
        % pCO2
        est.tp(i).pco2      = q(z(sys.tp(i).ipco2))*1e6; % convt atm to µatm
        est.tp(i).epco2     = ebar(sys.tp(i).ipco2)*1e6;
        est.tp(i).epco2_l   = ebar_l(sys.tp(i).ipco2)*1e6;
        est.tp(i).epco2_u   = ebar_u(sys.tp(i).ipco2)*1e6;
        est.tp(i).ppco2     = z(sys.tp(i).ipco2);
        est.tp(i).eppco2    = sigy(sys.tp(i).ipco2);
        % HCO3
        est.tp(i).hco3      = q(z(sys.tp(i).ihco3))*1e6; % convt mol/kg to µmol/kg
        est.tp(i).ehco3     = ebar(sys.tp(i).ihco3)*1e6;
        est.tp(i).ehco3_l   = ebar_l(sys.tp(i).ihco3)*1e6;
        est.tp(i).ehco3_u   = ebar_u(sys.tp(i).ihco3)*1e6;
        est.tp(i).phco3     = z(sys.tp(i).ihco3);
        est.tp(i).ephco3    = sigy(sys.tp(i).ihco3);
        % CO3
        est.tp(i).co3       = q(z(sys.tp(i).ico3))*1e6;
        est.tp(i).eco3      = ebar(sys.tp(i).ico3)*1e6;
        est.tp(i).eco3_l    = ebar_l(sys.tp(i).ico3)*1e6;
        est.tp(i).eco3_u    = ebar_u(sys.tp(i).ico3)*1e6;
        est.tp(i).pco3      = z(sys.tp(i).ico3);
        est.tp(i).epco3     = sigy(sys.tp(i).ico3);
        % pCO2*
        est.tp(i).pco2st    = z(sys.tp(i).ico2st);
        est.tp(i).epco2st   = sigy(sys.tp(i).ico2st);
        est.tp(i).co2st     = q(z(sys.tp(i).ico2st));
        est.tp(i).eco2st    = ebar(sys.tp(i).ico2st);
        est.tp(i).eco2st_l  = ebar_l(sys.tp(i).ico2st);
        est.tp(i).eco2st_u  = ebar_u(sys.tp(i).ico2st);
        % pP2F (aka FugFac)
        est.tp(i).pp2f      = z(sys.tp(i).ip2f);
        est.tp(i).epp2f     = sigy(sys.tp(i).ip2f);
        est.tp(i).p2f       = q(z(sys.tp(i).ip2f));
        est.tp(i).ep2f      = ebar(sys.tp(i).ip2f);
        est.tp(i).ep2f_l    = ebar_l(sys.tp(i).ip2f);
        est.tp(i).ep2f_u    = ebar_u(sys.tp(i).ip2f);

        % pK0 
        est.tp(i).pK0       = z(sys.tp(i).iK0);
        est.tp(i).epK0      = sigy(sys.tp(i).iK0);
        est.tp(i).K0        = q(z(sys.tp(i).iK0));
        est.tp(i).eK0       = ebar(sys.tp(i).iK0);
        est.tp(i).eK0_l     = ebar_l(sys.tp(i).iK0);
        est.tp(i).eK0_u     = ebar_u(sys.tp(i).iK0);
        % pK1
        est.tp(i).pK1   = z(sys.tp(i).iK1);
        est.tp(i).epK1  = sigy(sys.tp(i).iK1);
        est.tp(i).K1    = q(z(sys.tp(i).iK1));
        est.tp(i).eK1   = ebar(sys.tp(i).iK1);
        est.tp(i).eK1_l = ebar_l(sys.tp(i).iK1); 
        est.tp(i).eK1_u = ebar_u(sys.tp(i).iK1);
        % pK2
        est.tp(i).pK2   = z(sys.tp(i).iK2);
        est.tp(i).epK2  = sigy(sys.tp(i).iK2);
        est.tp(i).K2    = q(z(sys.tp(i).iK2));
        est.tp(i).eK2   = ebar(sys.tp(i).iK2);
        est.tp(i).eK2_l = ebar_l(sys.tp(i).iK2);
        est.tp(i).eK2_u = ebar_u(sys.tp(i).iK2);
        % OH 
        est.tp(i).oh    = q(z(sys.tp(i).ioh))*1e6; % convt
        est.tp(i).eoh   = ebar(sys.tp(i).ioh)*1e6;
        est.tp(i).eoh_l = ebar_l(sys.tp(i).ioh)*1e6;
        est.tp(i).eoh_u = ebar_u(sys.tp(i).ioh)*1e6;
        est.tp(i).poh   = z(sys.tp(i).ioh);
        est.tp(i).epoh  = sigy(sys.tp(i).ioh);
        % pKw 
        est.tp(i).pKw   = z(sys.tp(i).iKw);
        est.tp(i).epKw  = sigy(sys.tp(i).iKw);
        est.tp(i).Kw    = q(z(sys.tp(i).iKw));
        est.tp(i).eKw   = ebar(sys.tp(i).iKw);
        est.tp(i).eKw_l = ebar_l(sys.tp(i).iKw);
        est.tp(i).eKw_u = ebar_u(sys.tp(i).iKw);

        % BOH4 borate
        est.tp(i).boh4      = q(z(sys.tp(i).iboh4))*1e6; % convt mol/kg to µmol/kg
        est.tp(i).eboh4     = ebar(sys.tp(i).iboh4)*1e6;
        est.tp(i).eboh4_l   = ebar_l(sys.tp(i).iboh4)*1e6;
        est.tp(i).eboh4_u   = ebar_u(sys.tp(i).iboh4)*1e6;
        est.tp(i).pboh4     = z(sys.tp(i).iboh4);
        est.tp(i).epboh4    = sigy(sys.tp(i).iboh4);
        % BOH3
        est.tp(i).boh3      = q(z(sys.tp(i).iboh3))*1e6;
        est.tp(i).eboh3     = ebar(sys.tp(i).iboh3)*1e6;
        est.tp(i).eboh3_l   = ebar_l(sys.tp(i).iboh3)*1e6;
        est.tp(i).eboh3_u   = ebar_u(sys.tp(i).iboh3)*1e6;
        est.tp(i).pboh3     = z(sys.tp(i).iboh3);
        est.tp(i).epboh3    = sigy(sys.tp(i).iboh3);

        % pKb
        est.tp(i).pKb       = z(sys.tp(i).iKb);
        est.tp(i).epKb      = sigy(sys.tp(i).iKb);
        est.tp(i).Kb        = q(z(sys.tp(i).iKb));
        est.tp(i).eKb       = ebar(sys.tp(i).iKb);
        est.tp(i).eKb_l     = ebar_l(sys.tp(i).iKb);
        est.tp(i).eKb_u     = ebar_u(sys.tp(i).iKb);

        % SO4 sulfide
        est.tp(i).so4       = q(z(sys.tp(i).iso4)); % mol/kg
        est.tp(i).eso4      = ebar(sys.tp(i).iso4);
        est.tp(i).eso4_l    = ebar_l(sys.tp(i).iso4);
        est.tp(i).eso4_u    = ebar_u(sys.tp(i).iso4);
        est.tp(i).pso4      = z(sys.tp(i).iso4);
        est.tp(i).epso4     = sigy(sys.tp(i).iso4);
        % HSO4
        est.tp(i).hso4      = q(z(sys.tp(i).ihso4))*1e6;
        est.tp(i).ehso4     = ebar(sys.tp(i).ihso4)*1e6;
        est.tp(i).ehso4_l   = ebar_l(sys.tp(i).ihso4)*1e6;
        est.tp(i).ehso4_u   = ebar_u(sys.tp(i).ihso4)*1e6;
        est.tp(i).phso4     = z(sys.tp(i).ihso4);
        est.tp(i).ephso4    = sigy(sys.tp(i).ihso4);
        % pKs
        est.tp(i).pKs       = z(sys.tp(i).iKs);
        est.tp(i).epKs      = sigy(sys.tp(i).iKs);
        est.tp(i).Ks        = q(z(sys.tp(i).iKs));
        est.tp(i).eKs       = ebar(sys.tp(i).iKs);
        est.tp(i).eKs_l     = ebar_l(sys.tp(i).iKs);
        est.tp(i).eKs_u     = ebar_u(sys.tp(i).iKs);

        % F fluoride 
        est.tp(i).F         = q(z(sys.tp(i).iF))*1e6; % convt
        est.tp(i).eF        = ebar(sys.tp(i).iF)*1e6;
        est.tp(i).eF_l      = ebar_l(sys.tp(i).iF)*1e6;
        est.tp(i).ef_u      = ebar_u(sys.tp(i).iF)*1e6;
        est.tp(i).pF        = z(sys.tp(i).iF);
        est.tp(i).epF       = sigy(sys.tp(i).iF);
        % HF 
        est.tp(i).HF        = q(z(sys.tp(i).iHF))*1e6;
        est.tp(i).eHF       = ebar(sys.tp(i).iHF)*1e6;
        est.tp(i).eHF_l     = ebar_l(sys.tp(i).iHF)*1e6;
        est.tp(i).eHF_u     = ebar_u(sys.tp(i).iHF)*1e6;
        est.tp(i).pHF       = z(sys.tp(i).iHF);
        est.tp(i).epHF      = sigy(sys.tp(i).iHF);
        % pKf
        est.tp(i).pKf       = z(sys.tp(i).iKf);
        est.tp(i).epKf      = sigy(sys.tp(i).iKf);
        est.tp(i).Kf        = q(z(sys.tp(i).iKf));
        est.tp(i).eKf       = ebar(sys.tp(i).iKf);
        est.tp(i).eKf_l     = ebar_l(sys.tp(i).iKf);
        est.tp(i).eKf_u     = ebar_u(sys.tp(i).iKf);

        if ismember('phosphate', opt.abr)
            % PO4
            est.tp(i).po4       = q(z(sys.tp(i).ipo4))*1e6; % convt
            est.tp(i).epo4      = ebar(sys.tp(i).ipo4)*1e6;
            est.tp(i).epo4_l    = ebar_l(sys.tp(i).ipo4)*1e6;
            est.tp(i).epo4_u    = ebar_u(sys.tp(i).ipo4)*1e6;
            est.tp(i).ppo4      = z(sys.tp(i).ipo4);
            est.tp(i).eppo4     = sigy(sys.tp(i).ipo4);
            % HPO4
            est.tp(i).hpo4      = q(z(sys.tp(i).ihpo4))*1e6;
            est.tp(i).ehpo4     = ebar(sys.tp(i).ihpo4)*1e6;
            est.tp(i).ehpo4_l   = ebar_l(sys.tp(i).ihpo4)*1e6;
            est.tp(i).ehpo4_u   = ebar_u(sys.tp(i).ihpo4)*1e6;
            est.tp(i).phpo4     = z(sys.tp(i).ihpo4);
            est.tp(i).ephpo4    = sigy(sys.tp(i).ihpo4);
            % H2PO4
            est.tp(i).h2po4     = q(z(sys.tp(i).ih2po4))*1e6;
            est.tp(i).eh2po4    = ebar(sys.tp(i).ih2po4)*1e6;
            est.tp(i).eh2po4_l  = ebar_l(sys.tp(i).ih2po4)*1e6;
            est.tp(i).eh2po4_u  = ebar_u(sys.tp(i).ih2po4)*1e6;
            est.tp(i).ph2po4    = z(sys.tp(i).ih2po4);
            est.tp(i).eph2po4   = sigy(sys.tp(i).ih2po4);            
            % H3PO4
            est.tp(i).h3po4     = q(z(sys.tp(i).ih3po4))*1e6;
            est.tp(i).eh3po4    = ebar(sys.tp(i).ih3po4)*1e6;
            est.tp(i).eh3po4_l  = ebar_l(sys.tp(i).ih3po4)*1e6;
            est.tp(i).eh3po4_u  = ebar_u(sys.tp(i).ih3po4)*1e6;
            est.tp(i).ph3po4    = z(sys.tp(i).ih3po4);
            est.tp(i).eph3po4   = sigy(sys.tp(i).ih3po4);
            % pK1p
            est.tp(i).pK1p      = z(sys.tp(i).iK1p);
            est.tp(i).epK1p     = sigy(sys.tp(i).iK1p);
            est.tp(i).K1p       = q(z(sys.tp(i).iK1p));
            est.tp(i).eK1p      = ebar(sys.tp(i).iK1p);
            est.tp(i).eK1p_l    = ebar_l(sys.tp(i).iK1p);
            est.tp(i).eK1p_u    = ebar_u(sys.tp(i).iK1p);
            % pK2p
            est.tp(i).pK2p      = z(sys.tp(i).iK2p);
            est.tp(i).epK2p     = sigy(sys.tp(i).iK2p);
            est.tp(i).K2p       = q(z(sys.tp(i).iK2p));
            est.tp(i).eK2p      = ebar(sys.tp(i).iK2p);
            est.tp(i).pK2p_l    = ebar_l(sys.tp(i).iK2p);
            est.tp(i).eK2p_u    = ebar_u(sys.tp(i).iK2p);
            % pK3p
            est.tp(i).pK3p      = z(sys.tp(i).iK3p);
            est.tp(i).epK3p     = sigy(sys.tp(i).iK3p);
            est.tp(i).K3p       = q(z(sys.tp(i).iK3p));
            est.tp(i).eK3p      = ebar(sys.tp(i).iK3p);
            est.tp(i).eK3p_l    = ebar_l(sys.tp(i).iK3p);
            est.tp(i).eK3p_u    = ebar_u(sys.tp(i).iK3p);
        end

        if ismember('silicate', opt.abr)
            % SiOH4
            est.tp(i).sioh4     = q(z(sys.tp(i).isioh4))*1e6; % convt
            est.tp(i).esioh4    = ebar(sys.tp(i).isioh4)*1e6;
            est.tp(i).esioh4_l  = ebar_l(sys.tp(i).isioh4)*1e6;
            est.tp(i).esioh4_u  = ebar_u(sys.tp(i).isioh4)*1e6;
            est.tp(i).psioh4    = z(sys.tp(i).isioh4);
            est.tp(i).epsioh4   = sigy(sys.tp(i).isioh4);
            % SiOH3
            est.tp(i).siooh3    = q(z(sys.tp(i).isiooh3))*1e6;
            est.tp(i).esiooh3   = ebar(sys.tp(i).isiooh3)*1e6;
            est.tp(i).esiooh3_l = ebar_l(sys.tp(i).isiooh3)*1e6;
            est.tp(i).esiooh3_u = ebar_u(sys.tp(i).isiooh3)*1e6;
            est.tp(i).psiooh3   = z(sys.tp(i).isiooh3);
            est.tp(i).epsiooh3  = sigy(sys.tp(i).isiooh3);
            % pKsi
            est.tp(i).pKsi      = z(sys.tp(i).iKsi);
            est.tp(i).epKsi     = sigy(sys.tp(i).iKsi);
            est.tp(i).Ksi       = q(z(sys.tp(i).iKsi));
            est.tp(i).eKsi      = ebar(sys.tp(i).iKsi);
            est.tp(i).eKsi_l    = ebar_l(sys.tp(i).iKsi);
            est.tp(i).eKsi_u    = ebar_u(sys.tp(i).iKsi);
        end

        if ismember('ammonia', opt.abr)
            % NH3
            est.tp(i).nh3       = q(z(sys.tp(i).inh3))*1e6; % convt
            est.tp(i).enh3      = ebar(sys.tp(i).inh3)*1e6;
            est.tp(i).enh3_l    = ebar_l(sys.tp(i).inh3)*1e6;
            est.tp(i).enh3_u    = ebar_u(sys.tp(i).inh3)*1e6;
            est.tp(i).pnh3      = z(sys.tp(i).inh3);
            est.tp(i).epnh3     = sigy(sys.tp(i).inh3);
            % NH4
            est.tp(i).nh4       = q(z(sys.tp(i).inh4))*1e6;
            est.tp(i).enh4      = ebar(sys.tp(i).inh4)*1e6;
            est.tp(i).enh4_l    = ebar_l(sys.tp(i).inh4)*1e6;
            est.tp(i).enh4_u    = ebar_u(sys.tp(i).inh4)*1e6;
            est.tp(i).pnh4      = z(sys.tp(i).inh4);
            est.tp(i).epnh4     = sigy(sys.tp(i).inh4);
            % pKNH4
            est.tp(i).pKnh4     = z(sys.tp(i).iKnh4);
            est.tp(i).epKnh4    = sigy(sys.tp(i).iKnh4);
            est.tp(i).Knh4      = q(z(sys.tp(i).iKnh4));
            est.tp(i).eKnh4     = ebar(sys.tp(i).iKnh4);
            est.tp(i).eKnh4_l   = ebar_l(sys.tp(i).iKnh4);
            est.tp(i).eKnh4_u   = ebar_u(sys.tp(i).iKnh4);
        end

        if ismember('sulfide', opt.abr)
            % HS
            est.tp(i).hs        = q(z(sys.tp(i).ihs))*1e6; % convt
            est.tp(i).ehs       = ebar(sys.tp(i).ihs)*1e6;
            est.tp(i).ehs_l     = ebar_l(sys.tp(i).ihs)*1e6;
            est.tp(i).ehs_u     = ebar_u(sys.tp(i).ihs)*1e6;
            est.tp(i).phs       = z(sys.tp(i).ihs);
            est.tp(i).ephs      = sigy(sys.tp(i).ihs);
            % H2S
            est.tp(i).h2s       = q(z(sys.tp(i).ih2s))*1e6;
            est.tp(i).eh2s      = ebar(sys.tp(i).ih2s)*1e6;
            est.tp(i).eh2s_l    = ebar_l(sys.tp(i).ih2s)*1e6;
            est.tp(i).ehs2_u    = ebar_u(sys.tp(i).ih2s)*1e6;
            est.tp(i).ph2s      = z(sys.tp(i).ih2s);
            est.tp(i).eph2s     = sigy(sys.tp(i).ih2s);
            % pKh2s
            est.tp(i).pKh2s     = z(sys.tp(i).iKh2s);
            est.tp(i).epKh2s    = sigy(sys.tp(i).iKh2s);
            est.tp(i).Kh2s      = q(z(sys.tp(i).iKh2s));
            est.tp(i).eKh2s     = ebar(sys.tp(i).iKh2s);
            est.tp(i).eKh2s_l   = ebar_l(sys.tp(i).iKh2s);
            est.tp(i).eKh2s_u   = ebar_u(sys.tp(i).iKh2s);
        end

        if ismember('solubility',opt.abr)
            % Ca
            est.tp(i).ca        = q(z(sys.tp(i).ica))*1e6;
            est.tp(i).eca       = ebar(sys.tp(i).ica)*1e6;
            est.tp(i).eca_l     = ebar_l(sys.tp(i).ica)*1e6;
            est.tp(i).eca_u     = ebar_u(sys.tp(i).ica)*1e6;
            est.tp(i).pca       = z(sys.tp(i).ica);
            est.tp(i).epca      = sigy(sys.tp(i).ica);
            % Omega_Ar
            est.tp(i).OmegaAr    = q(z(sys.tp(i).iOmegaAr)); % unitless
            est.tp(i).eOmegaAr   = ebar(sys.tp(i).iOmegaAr);
            est.tp(i).eOmegaAr_l = ebar_l(sys.tp(i).iOmegaAr);
            est.tp(i).eOmegaAr_u = ebar_u(sys.tp(i).iOmegaAr);
            est.tp(i).pOmegaAr   = z(sys.tp(i).iOmegaAr);
            est.tp(i).epOmegaAr  = sigy(sys.tp(i).iOmegaAr);
            % Omega_Ca
            est.tp(i).OmegaCa    = q(z(sys.tp(i).iOmegaCa));
            est.tp(i).eOmegaCa   = ebar(sys.tp(i).iOmegaCa);
            est.tp(i).eOmegaCa_l = ebar_l(sys.tp(i).iOmegaCa);
            est.tp(i).eOmegaCa_u = ebar_u(sys.tp(i).iOmegaCa);
            est.tp(i).pOmegaCa   = z(sys.tp(i).iOmegaCa);
            est.tp(i).epOmegaCa  = sigy(sys.tp(i).iOmegaCa);
            % pKar
            est.tp(i).pKar      = z(sys.tp(i).iKar);
            est.tp(i).epKar     = sigy(sys.tp(i).iKar);
            est.tp(i).Kar       = q(z(sys.tp(i).iKar));
            est.tp(i).eKar      = ebar(sys.tp(i).iKar);
            est.tp(i).eKar_l    = ebar_l(sys.tp(i).iKar);
            est.tp(i).eKar_u    = ebar_u(sys.tp(i).iKar);
            % pKca
            est.tp(i).pKca      = z(sys.tp(i).iKca);
            est.tp(i).epKca     = sigy(sys.tp(i).iKca);
            est.tp(i).Kca       = q(z(sys.tp(i).iKca));
            est.tp(i).eKca      = ebar(sys.tp(i).iKca);
            est.tp(i).eKca_l    = ebar_l(sys.tp(i).iKca);
            est.tp(i).eKca_u    = ebar_u(sys.tp(i).iKca);
        end
    end
end

% -------------------------------------------------------------------------

function [zpK,zgpK,rph_free] = parse_zpK(x,sys,opt)
% assigning proper calculated values to zpK and zgpK

% zpK  := equilibrium constants, aka pK's
% zgpK := first derivative of equilibrium constants, aka gpK's
% rph_free := value of output of residual ph_free equation

    nTP         = length(sys.tp);
    M           = sys.M;
    K           = sys.K;
    nv          = size(M,2);
    nrk         =size(K,1);
    zpK         = zeros(nrk,1);
    zgpK        = zeros(nrk,nv);
    rph_free    = [];

    for i = 1:nTP
        [pK, gpK]   = calc_pK(opt, x(sys.tp(i).iT), x(sys.isal), ...
                                x(sys.tp(i).iP) );
        iTSP        = [ sys.tp(i).iT, sys.isal, sys.tp(i).iP];

        rph_free    = [ rph_free; sys.tp(i).rph_free(x)];

        zpK(sys.tp(i).kK0)          = pK(1); % sys.tp(i).iK0
        zgpK(sys.tp(i).kK0, iTSP )  = gpK(1,:); % ∂T, ∂S, ∂P

        zpK(sys.tp(i).kK1)          = pK(2);
        zgpK(sys.tp(i).kK1, iTSP )  = gpK(2,:); % ∂T, ∂S, ∂P

        zpK(sys.tp(i).kK2)          = pK(3);
        zgpK(sys.tp(i).kK2, iTSP )  = gpK(3,:); % ∂T, ∂S, ∂P

        zpK(sys.tp(i).kp2f)         = pK(14);
        zgpK(sys.tp(i).kp2f, iTSP ) = gpK(14,:); % ∂T, ∂S, ∂P

        zpK(sys.tp(i).kKb)          = pK(4);
        zgpK(sys.tp(i).kKb, iTSP )  = gpK(4,:); % ∂T, ∂S, ∂P

        zpK(sys.tp(i).kKw)          = pK(5);
        zgpK(sys.tp(i).kKw, iTSP )  = gpK(5,:); % ∂T, ∂S, ∂P

        zpK(sys.tp(i).kKs)          = pK(6);
        zgpK(sys.tp(i).kKs, iTSP )  = gpK(6,:); % ∂T, ∂S, ∂P

        zpK(sys.tp(i).kKf)          = pK(7);
        zgpK(sys.tp(i).kKf, iTSP )  = gpK(7,:); % ∂T, ∂S, ∂P
        
        if (ismember('phosphate',opt.abr))
            zpK(sys.tp(i).kK1p)         = pK(8); 
            zgpK(sys.tp(i).kK1p, iTSP ) = gpK(8,:); % ∂T, ∂S, ∂P 
            
            zpK(sys.tp(i).kK2p)         = pK(9); 
            zgpK(sys.tp(i).kK2p,iTSP )  = gpK(9,:); % ∂T, ∂S, ∂P 
            
            zpK(sys.tp(i).kK3p)         = pK(10); 
            zgpK(sys.tp(i).kK3p, iTSP ) = gpK(10,:); % ∂T, ∂S, ∂P 
        end

        if (ismember('silicate',opt.abr))
            zpK(sys.tp(i).kKsi)         = pK(11); 
            zgpK(sys.tp(i).kKsi, iTSP ) = gpK(11,:); % ∂T, ∂S, ∂P 
        end

        if (ismember('ammonia',opt.abr))
            zpK(sys.tp(i).kKnh4)         = pK(12); 
            zgpK(sys.tp(i).kKnh4, iTSP ) = gpK(12,:); % ∂T, ∂S, ∂P  
        end

        if (ismember('sulfide',opt.abr))
            zpK(sys.tp(i).kKh2s)         = pK(13); 
            zgpK(sys.tp(i).kKh2s, iTSP ) = gpK(13,:); % ∂T, ∂S, ∂P 
        end

        if (ismember('solubility',opt.abr))
            zpK(sys.tp(i).kKar)         = pK(15);
            zgpK(sys.tp(i).kKar, iTSP ) = gpK(15,:); % ∂T, ∂S, ∂P 
            zpK(sys.tp(i).kKca)         = pK(16);
            zgpK(sys.tp(i).kKca, iTSP ) = gpK(16,:); % ∂T, ∂S, ∂P  
        end
    end
end

% ---------------------------------------------------------------------------------

function ph_out = phscales(phin,scalein,pTS,pKs,pTF,pKf,pfH)
    % convert input ph to three scales (excluding ph_free)
    p   = @(x) -log10(x);
    q   = @(x) 10.^(-x);

    % FREE2tot  = ( 1 + TS/Ks ) = (TS + KS) / Ks
    pFREE2tot   = p(q(pTS) + q(pKs)) - pKs;
    % SWS2tot   = ( 1 + TS/Ks ) / ( 1 + TS/Ks + TF/Kf )
    pSWS2tot    = p(q(pTS) + q(pKs)) - pKs - ... 
                    p( 1 + (q(pTS)./q(pKs)) + (q(pTF)./q(pKf)) );
    % pSWS2free = pFREE2tot - pSWS2tot = -p(1 + TS/Ks + TF/Kf)
    pSWS2free   = -p( 1 + (q(pTS)/q(pKs)) + (q(pTF)/q(pKf)) );

    if scalein == 1 % tot
        ph_tot = phin;
        ph_sws = phin - pSWS2tot;
        ph_nbs = phin - pSWS2tot + pfH;        
    elseif scalein == 2 % sws
        ph_tot = phin + pSWS2tot;
        ph_sws = phin;
        ph_nbs = phin + pfH;
    elseif scalein == 3 % free
        ph_tot = phin + pFREE2tot;
        ph_sws = phin - pSWS2free;
        ph_nbs = phin - pSWS2free + pfH;
    elseif scalein == 4 % nbs
        ph_tot = phin + pSWS2tot - pfH;
        ph_sws = phin - pfH;
        ph_nbs = phin;
    end
    ph_out = [ph_tot, ph_sws, ph_nbs];
    % ph_free calculated and output from parse_output
end

% ----------------------------------------------------------------------------------

function z0 = init(yobs,sys,opt)
    q   = sys.q;      p = sys.p;
    
    y0  = yobs; 
    dic = q(yobs(sys.iTC));
    alk = q(yobs(sys.iTA));
    if (isnan(dic))
        dic         = 2200e-6;
        y0(sys.iTC) = p(dic);
    end
    if (isnan(alk))
        alk         = 2200e-6;
        y0(sys.iTA) = p(alk);
    end

    nTP = length(sys.tp);
    for i = 1:nTP
        % solve for the [H+] using only the carbonate alkalinity
        gam     = dic/alk;
        K0      = q(y0(sys.tp(i).iK0));
        K1      = q(y0(sys.tp(i).iK1));
        K2      = q(y0(sys.tp(i).iK2));
        p2f     = q(y0(sys.tp(i).ip2f));
        h       = 0.5*( ( gam - 1 ) * K1 + ( ( 1 - gam )^2 * ...
                    K1^2 - 4 * K1 * K2 * ( 1 - 2 * gam ) ).^0.5 ) ;
        hco3    =  h * alk / (h + 2 * K2 );
        co2st   = h * hco3 / K1 ;
        co3     = dic*K1*K2/(K1*h + h*h + K1*K2) ;
        fco2    = co2st/K0;
        pco2    = fco2/p2f;

        y0(sys.tp(i).iph)    = p(h);
        y0(sys.tp(i).ihco3)  = p(hco3);
        y0(sys.tp(i).ico2st) = p(co2st);
        y0(sys.tp(i).ico3)   = p(co3);
        y0(sys.tp(i).ifco2)  = p(fco2);
        y0(sys.tp(i).ipco2)  = p(pco2);

        Kw      = q(y0(sys.tp(i).iKw));
        oh      = Kw / h;
        y0(sys.tp(i).ioh)   = p(oh);

        Kb      = q(y0(sys.tp(i).iKb));
        TB      = q(yobs(sys.iTB)); 
        boh4    = TB * Kb / (Kb + h) ;
        boh3    = TB - boh4;
        y0(sys.iTB)         = p(TB);
        y0(sys.tp(i).iboh3) = p(boh3);
        y0(sys.tp(i).iboh4) = p(boh4);

        Ks      = q(y0(sys.tp(i).iKs));
        TS      = q(yobs(sys.iTS));
        h_free  = h / ( 1 + TS / Ks );
        hso4    = TS / ( 1 + Ks / h_free);
        so4     = Ks * hso4 / h_free;        
        y0(sys.iTS)             = p(TS);
        y0(sys.tp(i).iph_free)  = p(h_free); % ph_free
        y0(sys.tp(i).ihso4)     = p(hso4);
        y0(sys.tp(i).iso4)      = p(so4);

        Kf      = q(y0(sys.tp(i).iKf));
        TF      = q(yobs(sys.iTF));
        HF      = TF / ( 1 + Kf / h_free );
        F       = Kf * HF / h_free;
        y0(sys.iTF)         = p(TF);
        y0(sys.tp(i).iF)    = p(F);
        y0(sys.tp(i).iHF)   = p(HF);

        if (ismember('phosphate',opt.abr))
            K1p     = q(y0(sys.tp(i).iK1p));
            K2p     = q(y0(sys.tp(i).iK2p));
            K3p     = q(y0(sys.tp(i).iK3p));
            TP      = q(yobs(sys.iTP));
            d = ( h^3 + K1p * h^2 + K1p * K2p * h + K1p * K2p * K3p);
            h3po4   = TP * h^3 / d;
            h2po4   = TP * K1p * h^2 / d;
            hpo4    = TP * K1p * K2p * h / d;
            po4     = TP * K1p * K2p * K3p / d;
            y0(sys.iTP)             = p(TP);
            y0(sys.tp(i).ih3po4)    = p(h3po4);
            y0(sys.tp(i).ih2po4)    = p(h2po4);
            y0(sys.tp(i).ihpo4)     = p(hpo4);
            y0(sys.tp(i).ipo4)      = p(po4);
        end

        if (ismember('silicate',opt.abr))
            Ksi     = q(y0(sys.tp(i).iKsi));
            TSi     = q(yobs(sys.iTSi));
            siooh3  = TSi / ( 1 + h / Ksi );
            sioh4   = TSi - siooh3;
            y0(sys.iTSi)            = p(TSi);
            y0(sys.tp(i).isiooh3)   = p(siooh3);
            y0(sys.tp(i).isioh4)    = p(sioh4);
        end
        if (ismember('ammonia',opt.abr))
            Knh4    = q(y0(sys.tp(i).iKnh4));
            TNH3    = q(yobs(sys.iTNH3));
            nh3     = TNH3 / ( 1 + h / Knh4 );
            nh4     = TNH3 - nh3 ;
            y0(sys.iTNH3)       = p(TNH3);
            y0(sys.tp(i).inh3)  = p(nh3);
            y0(sys.tp(i).inh4)  = p(nh4);
        end
        if (ismember('sulfide',opt.abr))
            Kh2s    = q(y0(sys.tp(i).iKh2s));
            TH2S    = q(yobs(sys.iTH2S));
            hs      = TH2S / ( 1 + h / Kh2s );
            h2s     = TH2S - hs ;
            y0(sys.iTH2S)       = p(TH2S);
            y0(sys.tp(i).ihs)   = p(hs);
            y0(sys.tp(i).ih2s)  = p(h2s);
        end
        if (ismember('solubility',opt.abr))
            Kar     = q(y0(sys.tp(i).iKar));
            TCal    = q(yobs(sys.iTCal));
            OmegaAr = co3 * TCal / Kar;
            Kca     = q(y0(sys.tp(i).iKca));
            OmegaCa = co3 * TCal / Kca ;
            y0(sys.iTCal)           = p(TCal);
            y0(sys.tp(i).ica)       = p(TCal);
            y0(sys.tp(i).iOmegaAr)  = p(OmegaAr);
            y0(sys.tp(i).iOmegaCa)  = p(OmegaCa);
        end
    end
    % add the Lagrange multipliers
    nlam    = size(sys.M,1) + size(sys.K,1) + nTP;
    lam     = zeros(nlam,1);
    z0      = [y0(:);lam(:)];
end






