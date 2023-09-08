
% QUODcarbV7
% saved in FP_QUODcarb/v7

function [est,obs,iflag] = QUODcarbV7(obs,opt)
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
%   obs.m(1).T      = temp;         (deg C)         
%   obs.m(1).eT     = temp_error;   (±sigma)        
%   obs.m(1).P      = pressure;     (dbar)          
%   obs.m(1).eP     = pres_error;   (±sigma)     
%   obs.m(1).ph     = pH_meas;      
%   obs.m(1).eph    = ph_error;     (±sigma)
%
%   opt.K1K2        = 10;           (Lueker et al 2000)
%   opt.KSO4        = 1;            (Dickson et al 1990a) 
%   opt.KF          = 2;            (Perez and Fraga 1987
%   opt.TB          = 2;            (Lee et al. 2010)
%   opt.phscale     = 1;            (1=tot, 2=free, 3=sws, 4=nbs)
%   opt.printcsv    = 1;            (1=on, 0=off)
%   opt.fid         = 'output.csv'; (CSV filename)
%   opt.printmes    = 1;            (1=on, 0=off)
%   opt.co2press    = 1;            (1=on, 0=off)
%   opt.abr         = 'all';        {'borate','sulfate','fluoride',...
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
%           6 = GEOSECS                ~NOT AVAILABLE IN QUODCARB~
%           7 = Peng                   ~NOT AVAILABLE IN QUODCARB~
%           8 = Millero, 1979          ~NOT AVAILABLE IN QUODCARB~
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
%   est     ->   'est' structure with best estimate contains:
%               1. p(value) and p(error) where p(x) = -log10(x)
%               2. value and average error about the value in 'q'
%                   where q(x) = x^(-10)
%               3. upper and lower bounds in 'q' space, not symmetric
%                   about the value in 'q' space
%   csv     ->   csv file with most of est populated in a spreadsheet, 
%                 contains column headers with labels and units                 
%                    -does not include upper and lower errors
%
%--------------------------------------------------------------------------
%
% Changes? -> the only things you may want to change are: 
%               1. tolerance level of Newton solver -> line 117
%               2. Max Iteration number -> MAXIT in newtn.m
%
%--------------------------------------------------------------------------


    opt = check_opt(opt); % check opt structure

    sys = mksysV7(obs(1),opt.abr);

    nD = length(obs);
    nv = size(sys.K,2);

    % populate obs, yobs, wobs at each datapoint
    [obs,yobs,wobs] = parse_input(obs,sys,opt,nD);
  
    for i = 1:nD

        z0 = init(yobs(i,:),sys,opt);
        tol = 1e-7;

        gun = @(z) grad_limpco2(z,yobs(i,:),wobs(i,:),sys,opt);
        [z,J,iflag(i)] = newtn(z0,gun,tol);
        if (iflag(i) ~=0) && (opt.printmes ~= 0)
            fprintf('Newton''s method iflag = %i\n',iflag(i));
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
            % keyboard
        end
        % populate est
        [est(i)] = parse_output(z,sigy,opt,sys);

        if opt.printcsv == 1
            if i == 1
                fid = fopen(opt.fid,'w');
                PrintCSVv7(est(i),fid,opt); % make column headers
            end
            PrintCSVv7(opt,est(i),obs(i),iflag(i),fid); % fill with one row of data
        end

    end
  


end

% -----------------------------------------------------------------------------

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
  
    p = sys.p;
    q = sys.q;
    M = sys.M;
    K = sys.K;
    nrk = size(K,1);
    nTP = length(sys.m); 
    nv = size(M,2);
    nlam = size(M,1) + size(K,1) + (nTP); 
        % one extra lagrange multiplier for each 
        % (T,P)-dependent free to total ph conversions
    % nlam = size(M,1)+size(K,1)+(nTP*3);
        % got rid of f2t and replaced it with ph_sws, _free, _nbs (3 things)
    x   =  z(1:nv);      % measureable variables
    lam =  z(nv+1:end);  % Lagrange multipliers 
    
    % Make a vector of measured quantities    
    i = find(~isnan(w)); % was 'y', changed to w because MF put some defintions into obs w/o precisions
    y = y(i)';
    
    % Make a precision matrix
    W = diag(w(i));
    
    % Build a matrix that Picks out the measured components of x
    I = eye(nv); % for chain rule
    PP = I(i,:); % picking/pick out the measured ones
    e = PP*x - y; % calculated - measured (minus)
    
    % fill zpK and zgpK with associated calculated pK and gpK values
    % [zpK, zgpK, ph_all] = parse_zpK(x,sys,opt);
    [zpK, zgpK, ph_free] = parse_zpK(x,sys,opt);
    
    % constraint equations
    c = [  M * q( x ); ...
        (-K * x) + zpK;...
        ph_free]; %ph_all] ;
    
    f = 0.5 *  e.' * W * e  + lam.' * c ;  % limp, method of lagrange multipliers    
    % -(-1/2 sum of squares) + constraint eqns, minimize f => grad(f) = 0
    
    if ( nargout > 1 ) % compute the gradient
        % gph_sws = zeros(nTP,nv);
        gph_free = zeros(nTP,nv);
        % gph_nbs = zeros(nTP,nv);
        for i = 1:nTP
            % gph_sws(i,[ sys.iTS, sys.m(i).iKs, sys.iTF, sys.m(i).iKf, ...
            %     sys.m(i).iph, sys.m(i).ipfH]) = sys.m(i).gph_sws(x);
            gph_free(i,[ sys.iTS, sys.m(i).iKs, sys.iTF, sys.m(i).iKf, ...
                sys.m(i).iph, sys.m(i).ipfH]) = sys.m(i).gph_free(x);
            % gph_nbs(i,[ sys.iTS, sys.m(i).iKs, sys.iTF, sys.m(i).iKf, ...
            %     sys.m(i).iph, sys.m(i).ipfH]) = sys.m(i).gph_nbs(x);
        end
        dcdx = [ M * diag( sys.dqdx( x ) ); ...
            (-K + zgpK) ;... % constraint eqns wrt -log10(concentrations)
            gph_free];
            % gph_sws; gph_free; gph_nbs];

        % gf2t = zeros(nTP,nv);
        % for i = 1:nTP
        %     gf2t(i, [ sys.m(i).iph, sys.m(i).iKs, sys.iTS, ...
        %         sys.m(i).ipfH ] ) = sys.m(i).gf2t(x);
        % end
        % dcdx = [ M * diag( sys.dqdx( x ) ); ...
        %         (-K + zgpK) ;...
        %         gf2t ]; % constraint eqns wrt -log10(concentrations)
        g = [ e.' * W * PP +  lam.' * dcdx ,  c.' ];
    end
    %     
    if ( nargout > 2 ) % compute the Hessian
        ddq =  diag( sys.d2qdx2( x ) ); % q"
        [nr,nc] = size(M);
        gg = zeros(nc,1);
        for row = (1:nr)
            gg = gg + lam(row)*diag(M(row,:))*ddq;
        end
        for row = (nr+1):(nr+nrk)
            gg = gg + lam(row)*(-zggpK((row-nr),:,:)); % ggpK
        end
        dhfdx2 = zeros(nc,nc);
        ii = [sys.iKs,sys.iTS];
        dhfdx2(ii,ii) = sys.ggf2t(x);
        gg = gg + lam(nr+1)*dhfdx2;
        H = [  PP.'*W*PP + gg , dcdx.'  ; ...
               dcdx         , zeros(nlam)  ];
        keyboard
    end
    % keyboard
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
            fprintf('No acid base system chosen. Assuming opt.abr = {''all''}\n')
        end
        opt.abr = {'all'}; % default acid/base system = 'all'
    end
    if (strcmp(opt.abr,'all'))
        opt.abr = {'phosphate','silicate','ammonia','sulfide','solubility'};
%         opt.abr = {'borate','sulfate','fluoride','phosphate',...
%             'silicate','ammonia','sulfide','solubility'};
    end
end


% ------------------------------------------------------------------------

% NEED TO ADD pH SCALE CONVERTER ON INPUT !!!!!!!!!!!%%%%%%%%%%%%
function [obs,yobs,wobs] = parse_input(obs,sys,opt,nD)

    isgood = @(thing) (~isempty(thing) & ~sum(isnan(thing)));
    p = sys.p;
    q = sys.q;
    w = @(x,e) p(1+e./x).^(-2); % convert x+/-e into precision for p(x)

    nv = size(sys.K,2);
    yobs = nan(nD,nv);
    wobs = nan(nD,nv);

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
            obs(i).esal = 0.002; % std = 0.002 PSU
            wobs(i,sys.isal) = (obs(i).esal)^(-2);
            if opt.printmes ~= 0
                fprintf('Warning: Assuming salinity uncertainty is 0.002 PSU');
            end
        else
            wobs(i,sys.isal) = (obs(i).esal)^(-2); % std e -> w
        end
        if (~isfield(obs(i),'TC')) || (~isgood(obs(i).TC))
            obs(i).TC = [];
            yobs(i,sys.iTC) = nan;
        else
            yobs(i,sys.iTC) = p((obs(i).TC)*1e-6); % convt to mol/kg
        end
        if (~isfield(obs(i),'eTC')) || (~isgood(obs(i).eTC))
            obs(i).eTC = [];
            wobs(i,sys.iTC) = nan;
        else
            wobs(i,sys.iTC) = w(obs(i).TC,obs(i).eTC); % std e -> w
        end
        if(~isfield(obs(i),'TA'))  || (~isgood(obs(i).TA))
            obs(i).TA = [];
            yobs(i,sys.iTA) = nan;
        else
            yobs(i,sys.iTA) = p((obs(i).TA)*1e-6); % convt to mol/kg
        end
        if (~isfield(obs(i),'eTA'))  || (~isgood(obs(i).eTA))
            obs(i).eTA = [];
            wobs(i,sys.iTA) = nan;
        else
            wobs(i,sys.iTA) = w(obs(i).TA,obs(i).eTA); % std e -> w
        end
        if (~isfield(obs,'m'))
            if opt.printmes ~= 0
                error('Need to provide temperature and pressure measurement.')
            end
        end

        % create obs structure with fieldnames
        if (~isfield(obs(i).m(1),'epK0'))
            obs(i).m(1).epK0 = [];
        end
        if (~isfield(obs(i).m(1),'pK0'))
            obs(i).m(1).pK0 = [];
        end
        if (~isfield(obs(i).m(1),'epK1'))
            obs(i).m(1).epK1 = [];
        end
        if (~isfield(obs(i).m(1),'pK1'))
            obs(i).m(1).pK1 = [];
        end
        if (~isfield(obs(i).m(1),'epK2'))
            obs(i).m(1).epK2 = [];
        end
        if (~isfield(obs(i).m(1),'pK2'))
            obs(i).m(1).pK2 = [];
        end
        if (~isfield(obs(i).m(1),'epp2f')) || (~isgood(obs(i).m(1).epp2f))
            obs(i).m(1).epp2f = [];
        end
        if (~isfield(obs(i).m(1),'pp2f')) || (~isgood(obs(i).m(1).pp2f))
            obs(i).m(1).pp2f = [];
        end
        if (~isfield(obs(i).m(1),'efco2')) || (~isgood(obs(i).m(1).efco2))
            obs(i).m(1).efco2 = [];
        end
        if (~isfield(obs(i).m(1),'fco2')) || (~isgood(obs(i).m(1).fco2))
            obs(i).m(1).fco2 = [];
        end
        if (~isfield(obs(i).m(1),'epco2')) || (~isgood(obs(i).m(1).epco2))
            obs(i).m(1).epco2 = [];
        end
        if (~isfield(obs(i).m(1),'pco2')) || (~isgood(obs(i).m(1).pco2))
            obs(i).m(1).pco2 = [];
        end
        if (~isfield(obs(i).m(1),'eco2st')) || ...
                (~isgood(obs(i).m(1).eco2st))
            obs(i).m(1).eco2st = [];
        end
        if (~isfield(obs(i).m(1),'co2st')) || ...
                (~isgood(obs(i).m(1).co2st))
            obs(i).m(1).co2st = [];
        end
        if (~isfield(obs(i).m(1),'ehco3')) || (~isgood(obs(i).m(1).ehco3))
            obs(i).m(1).ehco3 = [];
        end
        if (~isfield(obs(i).m(1),'hco3')) || (~isgood(obs(i).m(1).hco3))
            obs(i).m(1).hco3 = [];
        end
        if (~isfield(obs(i).m(1),'eco3')) || (~isgood(obs(i).m(1).eco3))
            obs(i).m(1).eco3 = [];
        end
        if (~isfield(obs(i).m(1),'co3')) || (~isgood(obs(i).m(1).co3))
            obs(i).m(1).co3 = [];
        end
        if (~isfield(obs(i).m(1),'eph')) || (~isgood(obs(i).m(1).eph))
            obs(i).m(1).eph = [];
        end
        if (~isfield(obs(i).m(1),'ph')) || (~isgood(obs(i).m(1).ph))
            obs(i).m(1).ph = [];
        end
        % if (~isfield(obs(i).m(1),'ph_sws')) || (~isgood(obs(i).m(1).ph_sws))
        %     obs(i).m(1).ph_sws = [];
        % end
        if (~isfield(obs(i).m(1),'ph_free')) || (~isgood(obs(i).m(1).ph_free))
            obs(i).m(1).ph_free = [];
        end
        % if (~isfield(obs(i).m(1),'ph_nbs')) || (~isgood(obs(i).m(1).ph_nbs))
        %     obs(i).m(1).ph_nbs = [];
        % end
        if (~isfield(obs(i).m(1),'epfH')) || (~isgood(obs(i).m(1).epfH))
            obs(i).m(1).epfH = [];
        end
        if (~isfield(obs(i).m(1),'pfH')) || (~isgood(obs(i).m(1).pfH))
            obs(i).m(1).pfH = [];
        end


        if (~isfield(obs(i).m(1), 'epKw'))
            obs(i).m(1).epKw = [];
        end
        if (~isfield(obs(i).m(1), 'pKw'))
            obs(i).m(1).pKw = [];
        end
        if (~isfield(obs(i).m(1),'eoh')) || (~isgood(obs(i).m(1).eoh))
            obs(i).m(1).eoh = [];
        end
        if (~isfield(obs(i).m(1),'oh')) || (~isgood(obs(i).m(1).oh))
            obs(i).m(1).oh = [];
        end

        % borate system
        if (~isfield(obs(i).m(1),'epKb'))
            obs(i).m(1).epKb = [];
        end
        if (~isfield(obs(i).m(1),'pKb'))
            obs(i).m(1).pKb = [];
        end
        if (~isfield(obs(i),'TB')) || (~isgood(obs(i).TB))
            if opt.TB == 1
                % Uppstrom, L., Deep-Sea Research 21:161-162, 1974
                % ( copied from Orr's code )
                % TB = ( 0.000232/ 10.811) * (sal/1.80655)
                obs(i).TB = (0.0004157 * obs(i).sal / 35) * (1e6) ; % convt to µmol/kg
                yobs(i,sys.iTB) = p(obs(i).TB*1e-6);
            elseif opt.TB == 2
                % Lee, Kim, Myrne, Millero, Feely, Yong-Ming Liu. 2010.
                % Geochemica Et Cosmochimica Acta 74 (6): 1801-1811.
                % ( copied from Sharp's code )
                % TB = (0.0002414/ 10.811) * (sal/1.80655)
                obs(i).TB = (0.0004326 * obs(i).sal / 35)* (1e6) ; % umol/kg-SW
                yobs(i,sys.iTB) = p(obs(i).TB*1e-6);
            elseif opt.K1K2 == 6 || opt.K1K2 == 7
                % Culkin, F., in Chemical Oceanography,
                % ed. Riley and Skirrow, 1965: GEOSECS references this
                % (copied from Orr's code)
                obs(i).TB = (0.0004106 * obs(i).sal / 35) * (1e6) ; % umol/kg
                yobs(i,sys.iTB) = p(obs(i).TB*1e-6);
            end
        else
            if ((obs(i).TB) == 0)
                obs(i).TB = 1e-3; % umol/kg, reset minimum to 1 nanomolar
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
            wobs(i,sys.iTB) = w(obs(i).TB,obs(i).eTB);
        end
        if (~isfield(obs(i).m(1),'eboh4')) || ...
                (~isgood(obs(i).m(1).eboh4))
            obs(i).m(1).eboh4 = [];
        end
        if (~isfield(obs(i).m(1),'boh4')) || ...
                (~isgood(obs(i).m(1).boh4))
            obs(i).m(1).boh4 = [];
        end
        if (~isfield(obs(i).m(1),'eboh3')) || ...
                (~isgood(obs(i).m(1).eboh3))
            obs(i).m(1).eboh3 = [];
        end
        if (~isfield(obs(i).m(1),'boh3')) || ...
                (~isgood(obs(i).m(1).boh3))
            obs(i).m(1).boh3 = [];
        end

        % sulfate system
        if (~isfield(obs(i).m(1), 'epKs'))
            obs(i).m(1).epKs = [];
        end
        if (~isfield(obs(i).m(1), 'pKs'))
            obs(i).m(1).pKs = [];
        end
        if (~isfield(obs(i), 'TS'))  || (~isgood(obs(i).TS))
            % Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
            % copied from Orr's code
            obs(i).TS = ( 0.14 / 96.062 ) * ( obs(i).sal / 1.80655 ); % mol/kg
            yobs(i,sys.iTS) = p(obs(i).TS);
        else
            if ((obs(i).TS) == 0)
                obs(i).TS = 1e-9; % mol/kg reset minimum to 1 nanomolar
            end
            yobs(i,sys.iTS) = p(obs(i).TS);
        end
        if (~isfield(obs(i), 'eTS'))  || (~isgood(obs(i).eTS))
            % 0.14000 ± 0.00023
            TSu = ( ( (0.14+0.00023)/96.062 ) * obs(i).sal/ 1.80655 );
            TSl = ( ( (0.14-0.00023)/96.062 ) * obs(i).sal/ 1.80655 );
            eTS = (TSu - TSl) / 2;
            obs(i).eTS = (eTS) ;
            wobs(i,sys.iTS) = w(obs(i).TS,obs(i).eTS);
        else
            wobs(i,sys.iTS) = w(obs(i).TS,obs(i).eTS);
        end
        % if (~isfield(obs(i).m(1), 'phf')) || ... % changed to phfree in v6.5
        %         (~isgood(obs(i).m(1).phf)) % phfree
        %     obs(i).m(1).phf = []; % phfree
        % end
        % if (~isfield(obs(i).m(1), 'ephf')) || ... % ephfree
        %         (~isgood(obs(i).m(1).ephf)) % ephfree
        %     obs(i).m(1).ephf = []; % ephfree
        % end
        if (~isfield(obs(i).m(1), 'so4')) || ...
                (~isgood(obs(i).m(1).so4))
            obs(i).m(1).so4 = [];
        end
        if (~isfield(obs(i).m(1), 'eso4')) || ...
                (~isgood(obs(i).m(1).eso4))
            obs(i).m(1).eso4 = [];
        end
        if (~isfield(obs(i).m(1), 'hso4')) || ...
                (~isgood(obs(i).m(1).hso4))
            obs(i).m(1).hso4 = [];
        end
        if (~isfield(obs(i).m(1), 'ehso4')) || ...
                (~isgood(obs(i).m(1).ehso4))
            obs(i).m(1).ehso4 = [];
        end

        % fluoride system
        if (~isfield(obs(i).m(1), 'epKf'))
            obs(i).m(1).epKf = [];
        end
        if (~isfield(obs(i).m(1), 'pKf'))
            obs(i).m(1).pKf = [];
        end
        if (~isfield(obs(i), 'TF'))  || (~isgood(obs(i).TF))
            % Riley, J. P., Deep-Sea Research 12:219-220, 1965:
            % this is .000068.*Sali./35. = .00000195.*Sali
            obs(i).TF = ( 0.000067 / 18.998 ) * ( obs(i).sal / 1.80655 )*1e6; % convt to µmol/kg-SW
            yobs(i,sys.iTF) = p(obs(i).TF*1e-6);
        else
            if ((obs(i).TF) == 0)
                obs(i).TF = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.iTF) = p(obs(i).TF*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs(i), 'eTF'))  || (~isgood(obs(i).eTF))
            % 6.7 ± 0.1 e-5
            TFu = ( ( (6.7e-5 + 0.1e-5)/18.998) * obs(i).sal/1.80655 );
            TFl = ( ( (6.7e-5 - 0.1e-5)/18.998) * obs(i).sal/1.80655 );
            eTF = (TFu - TFl) / 2;
            obs(i).eTF = (eTF)*1e6;
            wobs(i,sys.iTF) = w(obs(i).TF,obs(i).eTF);
        else
            wobs(i,sys.iTF) = w(obs(i).TF,obs(i).eTF);
        end
        if (~isfield(obs(i).m(1), 'F')) || (~isgood(obs(i).m(1).F))
            obs(i).m(1).F = [];
        end
        if (~isfield(obs(i).m(1), 'eF')) || (~isgood(obs(i).m(1).eF))
            obs(i).m(1).eF = [];
        end
        if (~isfield(obs(i).m(1), 'HF')) || (~isgood(obs(i).m(1).HF))
            obs(i).m(1).HF = [];
        end
        if (~isfield(obs(i).m(1), 'eHF')) || (~isgood(obs(i).m(1).eHF))
            obs(i).m(1).eHF = [];
        end

        if (ismember('phosphate',opt.abr))
            if (~isfield(obs(i).m(1), 'epK1p'))
                obs(i).m(1).epK1p = [];
            end
            if (~isfield(obs(i).m(1), 'pK1p'))
                obs(i).m(1).pK1p = [];
            end
            if (~isfield(obs(i).m(1), 'epK2p'))
                obs(i).m(1).epK2p = [];
            end
            if (~isfield(obs(i).m(1), 'pK2p'))
                obs(i).m(1).pK2p = [];
            end
            if (~isfield(obs(i).m(1), 'epK3p'))
                obs(i).m(1).epK3p = [];
            end
            if (~isfield(obs(i).m(1), 'pK3p'))
                obs(i).m(1).pK3p = [];
            end
            if (~isfield(obs(i), 'TP'))  || (~isgood(obs(i).TP))
                obs(i).TP = 1e-3; % µmol/kg
                yobs(i,sys.iTP) = p(obs(i).TP*1e-6); % convt µmol/kg to mol/kg
            else
                if ((obs(i).TP) == 0) % zero po4 is very unlikely and breaks the code
                    obs(i).TP  = 1e-3; % umol/kg-SW, reset minimum to 1 nanomolar
                end
                yobs(i,sys.iTP) = p(obs(i).TP*1e-6); % convt µmol/kg to mol/kg
            end
            if (~isfield(obs(i), 'eTP'))  || (~isgood(obs(i).eTP))
                obs(i).eTP = 1e-3; % µmol/kg
                wobs(i,sys.iTP) = w(obs(i).TP,obs(i).eTP);
            else
                if ((obs(i).eTP) == 0)
                    obs(i).eTP = 1e-3; % umol/kg, reset minimum if zero
                end
                wobs(i,sys.iTP) = w(obs(i).TP,obs(i).eTP);
            end
            if (~isfield(obs(i).m(1), 'h3po4')) || ...
                    (~isgood(obs(i).m(1).h3po4))
                obs(i).m(1).h3po4 = [];
            end
            if (~isfield(obs(i).m(1), 'eh3po4')) || ...
                    (~isgood(obs(i).m(1).eh3po4))
                obs(i).m(1).eh3po4 = [];
            end
            if (~isfield(obs(i).m(1), 'h2po4')) || ...
                    (~isgood(obs(i).m(1).h2po4))
                obs(i).m(1).h2po4 = [];
            end
            if (~isfield(obs(i).m(1), 'eh2po4')) || ...
                    (~isgood(obs(i).m(1).eh2po4))
                obs(i).m(1).eh2po4 = [];
            end
            if (~isfield(obs(i).m(1), 'hpo4')) || ...
                    (~isgood(obs(i).m(1).hpo4))
                obs(i).m(1).hpo4 = [];
            end
            if (~isfield(obs(i).m(1), 'ehpo4')) || ...
                    (~isgood(obs(i).m(1).ehpo4))
                obs(i).m(1).ehpo4 = [];
            end
            if (~isfield(obs(i).m(1), 'po4')) || ...
                    (~isgood(obs(i).m(1).po4))
                obs(i).m(1).po4 = [];
            end
            if (~isfield(obs(i).m(1), 'epo4')) || ...
                    (~isgood(obs(i).m(1).epo4))
                obs(i).m(1).epo4 = [];
            end
        end

        if (ismember('silicate',opt.abr))
            if (~isfield(obs(i).m(1), 'epKsi'))
                obs(i).m(1).epKsi = [];
            end
            if (~isfield(obs(i).m(1), 'pKsi'))
                obs(i).m(1).pKsi = [];
            end
            if (~isfield(obs(i), 'TSi'))  || (~isgood(obs(i).TSi))
                obs(i).TSi = 1e-3; % µmol/kg
                yobs(i,sys.iTSi) = p(obs(i).TSi*1e-6); % convt µmol/kg to mol/kg
            else
                if ((obs(i).TSi) == 0) % zero silicate very unlikely and breaks code
                    obs(i).TSi = 1e-3; % umol/kg, reset minimum to 1 nanomolar
                end
                yobs(i,sys.iTSi) = p(obs(i).TSi*1e-6);
            end
            if (~isfield(obs(i), 'eTSi'))  || (~isgood(obs(i).eTSi))
                obs(i).eTSi = 1e-3; % µmol/kg
                wobs(i,sys.iTSi) = w(obs(i).TSi,obs(i).eTSi);
            else
                if ((obs(i).eTSi) == 0)
                    obs(i).eTSi = 1e-3; % umol/kg, reset minimum to 1 nanomolar
                end
                wobs(i,sys.iTSi) = w(obs(i).TSi,obs(i).eTSi);
            end
            if (~isfield(obs(i).m(1), 'sioh4')) || ...
                    (~isgood(obs(i).m(1).sioh4))
                obs(i).m(1).sioh4 = [];
            end
            if (~isfield(obs(i).m(1), 'esioh4')) || ...
                    (~isgood(obs(i).m(1).esioh4))
                obs(i).m(1).esioh4 = [];
            end
            if (~isfield(obs(i).m(1), 'siooh3')) || ...
                    (~isgood(obs(i).m(1).siooh3))
                obs(i).m(1).siooh3 = [];
            end
            if (~isfield(obs(i).m(1), 'esiooh3')) || ...
                    (~isgood(obs(i).m(1).esiooh3))
                obs(i).m(1).esiooh3 = [];
            end
        end

        if (ismember('ammonia',opt.abr))
            if (~isfield(obs(i).m(1), 'epKnh4'))
                obs(i).m(1).epKnh4 = [];
            end
            if (~isfield(obs(i).m(1), 'pKnh4'))
                obs(i).m(1).pKnh4 = [];
            end
            if (~isfield(obs(i), 'TNH3'))  || (~isgood(obs(i).TNH3))
                obs(i).TNH3 = 1e-3; % µmol/kg
                yobs(i,sys.iTNH3) = p(obs(i).TNH3*1e-6); % convt µmol/kg to mol/kg
            else
                if ((obs(i).TNH3) == 0)
                    obs(i).TNH3 = 1e-3; % umol/kg, reset minimum to 1 nanomolar
                end
                yobs(i,sys.iTNH3) = p(obs(i).TNH3*1e-6);
            end
            if (~isfield(obs(i), 'eTNH3'))  || (~isgood(obs(i).eTNH3))
                obs(i).eTNH3 = 1e-3; % µmol/kg
                wobs(i,sys.iTNH3) = w(obs(i).TNH3,obs(i).eTNH3);
            else
                wobs(i,sys.iTNH3) = w(obs(i).TNH3,obs(i).eTNH3);
            end
            if (~isfield(obs(i).m(1), 'nh3')) || ...
                    (~isgood(obs(i).m(1).nh3))
                obs(i).m(1).nh3 = [];
            end
            if (~isfield(obs(i).m(1), 'enh3')) || ...
                    (~isgood(obs(i).m(1).enh3))
                obs(i).m(1).enh3 = [];
            end
            if (~isfield(obs(i).m(1), 'nh4')) || ...
                    (~isgood(obs(i).m(1).nh4))
                obs(i).m(1).nh4 = [];
            end
            if (~isfield(obs(i).m(1), 'enh4')) || ...
                    (~isgood(obs(i).m(1).enh4))
                obs(i).m(1).enh4 = [];
            end
        end

        if (ismember('sulfide',opt.abr))
            if (~isfield(obs(i).m(1), 'epKh2s'))
                obs(i).m(1).epKh2s = [];
            end
            if (~isfield(obs(i).m(1), 'pKh2s'))
                obs(i).m(1).pKh2s = [];
            end
            if (~isfield(obs(i), 'TH2S'))  || (~isgood(obs(i).TH2S))
                obs(i).TH2S = 1e-3; % µmol/kg
                yobs(i,sys.iTH2S) = p(obs(i).TH2S*1e-6); % convt µmol/kg to mol/kg
            else
                if ((obs(i).TH2S) == 0)
                    obs(i).TH2S = 1e-3; % umol/kg, reset minimum to 1 nanomolar
                end
                yobs(i,sys.iTH2S) = p(obs(i).TH2S*1e-6);
            end
            if (~isfield(obs(i), 'eTH2S'))  || (~isgood(obs(i).eTH2S))
                obs(i).eTH2S = 1e-3; % µmol/kg
                wobs(i,sys.iTH2S) = w(obs(i).TH2S,obs(i).eTH2S);
            else
                wobs(i,sys.iTH2S) = w(obs(i).TH2S,obs(i).eTH2S);
            end
            if (~isfield(obs(i).m(1), 'hs')) || ...
                    (~isgood(obs(i).m(1).hs))
                obs(i).m(1).hs = [];
            end
            if (~isfield(obs(i).m(1), 'ehs')) || ...
                    (~isgood(obs(i).m(1).ehs))
                obs(i).m(1).ehs = [];
            end
            if (~isfield(obs(i).m(1), 'h2s')) || ...
                    (~isgood(obs(i).m(1).h2s))
                obs(i).m(1).h2s = [];
            end
            if (~isfield(obs(i).m(1), 'eh2s')) || ...
                    (~isgood(obs(i).m(1).eh2s))
                obs(i).m(1).eh2s = [];
            end
        end

        if (ismember('solubility',opt.abr))
            if (~isfield(obs(i).m(1), 'epKar'))
                obs(i).m(1).epKar = [];
            end
            if (~isfield(obs(i).m(1), 'pKar'))
                obs(i).m(1).pKar = [];
            end
            if (~isfield(obs(i), 'TCal'))  || (~isgood(obs(i).TCal))
                if opt.K1K2 == 6 || opt.K1K2 == 7
                    % Calculate Ca for GEOSECS, Riley and Skirrow 1965
                    obs(i).TCal = (0.01026 .* obs(i).sal ./ 35) * 1e6 ;
                    yobs(i,sys.iTCal) = p(obs(i).TCal*1e-6); % convt µmol/kg to mol/kg
                else
                    % Calculate Ca, Riley and Tongdui 1967
                    % this is 0.010285.*obs.sal./35;
                    obs(i).TCal = (0.02128./40.087.*(obs(i).sal./1.80655)) * 1e6 ; % convt to umol/kg
                    yobs(i,sys.iTCal) = p(obs(i).TCal*1e-6); % convert back to mol/kg
                end
            else
                if ((obs(i).TCal) == 0)
                    obs(i).TCal = 1e-3; % umol/kg, reset minimum to 1 nanomolar
                end
                yobs(i,sys.iTCal) = p(obs(i).TCal*1e-6); % assume user input of umol/kg
            end
            if (~isfield(obs(i), 'eTCal'))  || (~isgood(obs(i).eTCal))
                obs(i).eTCal = (6e-5)*1e6; % umol/kg, from Riley and Tongdui 1967
                wobs(i,sys.iTCal) = w(obs(i).TCal,obs(i).eTCal);
            else
                wobs(i,sys.iTCal) = w(obs(i).TCal,obs(i).eTCal);
            end
            if (~isfield(obs(i).m(1), 'ca')) || ...
                    (~isgood(obs(i).m(1).ca))
                obs(i).m(1).ca = [];
            end
            if (~isfield(obs(i).m(1), 'eca')) || ...
                    (~isgood(obs(i).m(1).eca))
                obs(i).m(1).eca = [];
            end
            if (~isfield(obs(i).m(1), 'OmegaAr')) || ...
                    (~isgood(obs(i).m(1).OmegaAr))
                obs(i).m(1).OmegaAr = [];
            end
            if (~isfield(obs(i).m(1), 'eOmegaAr')) || ...
                    (~isgood(obs(i).m(1).eOmegaAr))
                obs(i).m(1).eOmegaAr = [];
            end
            if (~isfield(obs(i).m(1), 'pKca'))
                obs(i).m(1).pKca = [];
            end
            if (~isfield(obs(i).m(1), 'epKca'))
                obs(i).m(1).epKca = [];
            end
            if (~isfield(obs(i).m(1), 'OmegaCa')) || ...
                    (~isgood(obs(i).m(1).OmegaCa))
                obs(i).m(1).OmegaCa = [];
            end
            if (~isfield(obs(i).m(1), 'eOmegaCa')) || ...
                    (~isgood(obs(i).m(1).eOmegaCa))
                obs(i).m(1).eOmegaCa = [];
            end
        end

        nTP = length(obs(i).m);

        for ii = 1:nTP % loop over all the pressure and temperature sensitive components
            yobs(i,sys.m(ii).iT) = obs(i).m(ii).T;
            yobs(i,sys.m(ii).iP) = obs(i).m(ii).P;

            wobs(i,sys.m(ii).iT) = obs(i).m(ii).eT;
            wobs(i,sys.m(ii).iP) = obs(i).m(ii).eP;

            [pK,gpK,epK] = calc_pK(opt, obs(i).m(ii).T, obs(i).sal, ...
                obs(i).m(ii).P );
            
            pK0   = pK(1);  pK1  = pK(2);  pK2  = pK(3);  pKb   = pK(4);
            pKw   = pK(5);  pKs  = pK(6);  pKf  = pK(7);  pK1p  = pK(8);
            pK2p  = pK(9);  pK3p = pK(10); pKsi = pK(11); pKnh4 = pK(12);
            pKh2s = pK(13); pp2f = pK(14); pKar = pK(15); pKca  = pK(16);
            pfH   = pK(17);

            epK0   = epK(1);  epK1  = epK(2);  epK2  = epK(3);  epKb   = epK(4);
            epKw   = epK(5);  epKs  = epK(6);  epKf  = epK(7);  epK1p  = epK(8);
            epK2p  = epK(9);  epK3p = epK(10); epKsi = epK(11); epKnh4 = epK(12);
            epKh2s = epK(13); epp2f = epK(14); epKar = epK(15); epKca  = epK(16);
            epfH = epK(17);
            %
            % add "observations" for the equilibrium constants
            % and transfer from obs struct to yobs and wobs vectors (or
            % matrixes if nD >1)
            %
            if (isgood(obs(i).m(ii).epK0))
                wobs(i,sys.m(ii).iK0) = (obs(i).m(ii).epK0)^(-2);
            else
                obs(i).m(ii).epK0 = epK0;
                wobs(i,sys.m(ii).iK0) = (obs(i).m(ii).epK0)^(-2);
            end
            if (isgood(obs(i).m(ii).pK0))
                yobs(i,sys.m(ii).iK0) = obs(i).m(ii).pK0;
            else
                yobs(i,sys.m(ii).iK0) = pK0;
                obs(i).m(ii).pK0 = pK0;
            end

            if (isgood(obs(i).m(ii).epK1))
                wobs(i,sys.m(ii).iK1) = (obs(i).m(ii).epK1)^(-2);
            else
                obs(i).m(ii).epK1 = epK1;
                wobs(i,sys.m(ii).iK1) = (obs(i).m(ii).epK1)^(-2);
            end
            if (isgood(obs(i).m(ii).pK1))
                yobs(i,sys.m(ii).iK1) = obs(i).m(ii).pK1;
            else
                yobs(i,sys.m(ii).iK1) = pK1;
                obs(i).m(ii).pK1 = pK1;
            end

            if (isgood(obs(i).m(ii).epK2))
                wobs(i,sys.m(ii).iK2) = (obs(i).m(ii).epK2)^(-2);
            else
                obs(i).m(ii).epK2 = epK2;
                wobs(i,sys.m(ii).iK2) = (obs(i).m(ii).epK2)^(-2);  % wK2 = 1/(1 + (0.02/pKsys(3)))^2 ;
            end

            if (isgood(obs(i).m(ii).pK2))
                yobs(i,sys.m(ii).iK2) = obs(i).m(ii).pK2;
            else
                yobs(i,sys.m(ii).iK2) = pK2;
                obs(i).m(ii).pK2 = pK2;
            end

            if (isgood(obs(i).m(ii).epp2f))
                wobs(i,sys.m(ii).ip2f) = (obs(i).m(ii).epp2f)^(-2);
            else
                obs(i).m(ii).epp2f = epp2f;
                wobs(i,sys.m(ii).ip2f) = (obs(i).m(ii).epp2f)^(-2);
            end

            if (isgood(obs(i).m(ii).pp2f))
                yobs(i,sys.m(ii).ip2f) = obs(i).m(ii).pp2f;
            else
                yobs(i,sys.m(ii).ip2f) = pp2f;
                obs(i).m(ii).pp2f = pp2f;
            end

            if (isgood(obs(i).m(ii).co2st))
                yobs(i,sys.m(ii).ico2st) = p(obs(i).m(ii).co2st*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.m(ii).ico2st) = nan;
                obs(i).m(ii).co2st = nan;
            end
            if (isgood(obs(i).m(ii).eco2st))
                wobs(i,sys.m(ii).ico2st) = w(obs(i).m(ii).co2st,obs(i).m(ii).eco2st);
            else
                wobs(i,sys.m(ii).ico2st) = nan;
                obs(i).m(ii).eco2st = nan;
            end
            if (isgood(obs(i).m(ii).hco3))
                yobs(i,sys.m(ii).ihco3) = p(obs(i).m(ii).hco3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.m(ii).ihco3) = nan;
                obs(i).m(ii).hco3 = nan;
            end
            if (isgood(obs(i).m(ii).ehco3))
                wobs(i,sys.m(ii).ihco3) = w(obs(i).m(ii).hco3,obs(i).m(ii).ehco3);
            else
                wobs(i,sys.m(ii).ihco3) = nan;
                obs(i).m(ii).ehco3 = nan;
            end
            if (isgood(obs(i).m(ii).co3))
                yobs(i,sys.m(ii).ico3) = p(obs(i).m(ii).co3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.m(ii).ico3) = nan;
                obs(i).m(ii).co3 = nan;
            end
            if (isgood(obs(i).m(ii).eco3))
                wobs(i,sys.m(ii).ico3) = w(obs(i).m(ii).co3,obs(i).m(ii).eco3);
            else
                wobs(i,sys.m(ii).ico3) = nan;
                obs(i).m(ii).eco3 = nan;
            end
            if (isgood(obs(i).m(ii).pco2))
                yobs(i,sys.m(ii).ipco2) = p(obs(i).m(ii).pco2*1e-6); % convt µatm to atm
            else
                yobs(i,sys.m(ii).ipco2) = nan;
                obs(i).m(ii).pco2 = nan;
            end
            if (isgood(obs(i).m(ii).epco2))
                wobs(i,sys.m(ii).ipco2) = w(obs(i).m(ii).pco2,obs(i).m(ii).epco2);
            else
                wobs(i,sys.m(ii).ipco2) = nan;
                obs(i).m(ii).epco2 = nan;
            end
            if (isgood(obs(i).m(ii).fco2))
                yobs(i,sys.m(ii).ifco2) = p(obs(i).m(ii).fco2*1e-6); % convt µatm to atm
            else
                yobs(i,sys.m(ii).ifco2) = nan;
                obs(i).m(ii).fco2 = nan;
            end
            if (isgood(obs(i).m(ii).efco2))
                wobs(i,sys.m(ii).ifco2) = w(obs(i).m(ii).fco2,obs(i).m(ii).efco2);
            % elseif (isgood(wobs(i,sys.m(ii).ifco2)))
            %     wobs(i,sys.m(ii).ifco2) = wobs(i,sys.m(ii).ifco2);
            else
                wobs(i,sys.m(ii).ifco2) = nan;
                obs(i).m(ii).efco2 = nan;
            end
            if (isgood(obs(i).m(ii).ph)) % ph_tot
                yobs(i,sys.m(ii).iph) = obs(i).m(ii).ph ;
            else
                yobs(i,sys.m(ii).iph) = nan;
                obs(i).m(ii).ph = nan;
            end
            if (isgood(obs(i).m(ii).eph)) % eph same for all ph scales
                wobs(i,sys.m(ii).iph) = (obs(i).m(ii).eph).^(-2);
            else
                wobs(i,sys.m(ii).iph) = nan;
                obs(i).m(ii).eph = nan;
            end
            % if (isgood(obs(i).m(ii).ph_sws)) % ph_sws
            %     yobs(i,sys.m(ii).iph_sws) = obs(i).m(ii).ph_sws ;
            % else
            %     yobs(i,sys.m(ii).iph_sws) = nan;
            %     obs(i).m(ii).ph_sws = nan;
            % end
            if (isgood(obs(i).m(ii).ph_free)) % ph_free (ph on free scale)
                yobs(i,sys.m(ii).iph_free) = obs(i).m(ii).ph_free ;
            else
                yobs(i,sys.m(ii).iph_free) = nan;
                obs(i).m(ii).ph_free = nan;
            end
            % if (isgood(obs(i).m(ii).ph_nbs)) % ph_nbs
            %     yobs(i,sys.m(ii).iph_nbs) = obs(i).m(ii).ph_nbs ;
            % else
            %     yobs(i,sys.m(ii).iph_nbs) = nan;
            %     obs(i).m(ii).ph_nbs = nan;
            % end
            if (isgood(obs(i).m(ii).pfH)) % pfH activity coefficient
                yobs(i,sys.m(ii).ipfH) = obs(i).m(ii).pfH ;
            else
                yobs(i,sys.m(ii).ipfH) = pfH;
                obs(i).m(ii).pfH = pfH;
            end
            if (isgood(obs(i).m(ii).epfH)) % pfH activity coefficient
                wobs(i,sys.m(ii).ipfH) = (obs(i).m(ii).epfH).^(-2) ;
            else
                wobs(i,sys.m(ii).ipfH) = (epfH).^(-2);
                obs(i).m(ii).pfH = epfH;
            end


            if (isgood(obs(i).m(ii).epKw))
                wobs(i,sys.m(ii).iKw) = (obs(i).m(ii).epKw).^(-2);
            else
                obs(i).m(ii).epKw = epKw;
                wobs(i,sys.m(ii).iKw) = (obs(i).m(ii).epKw).^(-2);  % wKw = 1/(1 + (0.01/pKsys(5)))^2 ;
            end
            if (isgood(obs(i).m(ii).pKw))
                yobs(i,sys.m(ii).iKw) = obs(i).m(ii).pKw;
            else
                yobs(i,sys.m(ii).iKw) = pKw;
                obs(i).m(ii).pKw = pKw;
            end
            if (isgood(obs(i).m(ii).oh))
                yobs(i,sys.m(ii).ioh) = p(obs(i).m(ii).oh*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.m(ii).ioh) = nan;
                obs(i).m(ii).oh = nan;
            end
            if (isgood(obs(i).m(ii).eoh))
                wobs(i,sys.m(ii).ioh) = w(obs(i).m(ii).oh,obs(i).m(ii).eoh);
            else
                wobs(i,sys.m(ii).ioh) = nan;
                obs(i).m(ii).eoh = nan;
            end

            % borate system
            if (isgood(obs(i).m(ii).epKb))
                wobs(i,sys.m(ii).iKb) = (obs(i).m(ii).epKb).^(-2);
            else
                obs(i).m(ii).epKb = epKb;
                wobs(i,sys.m(ii).iKb) = (obs(i).m(ii).epKb).^(-2);  % wKb = 1/(1 + (0.01/pKsys(4)))^2 ;
            end
            if (isgood(obs(i).m(ii).pKb))
                yobs(i,sys.m(ii).iKb) = obs(i).m(ii).pKb;
            else
                yobs(i,sys.m(ii).iKb) = pKb; % from local_pK
                obs(i).m(ii).pKb = pKb;
            end
            if (isgood(obs(i).m(ii).boh3))
                yobs(i,sys.m(ii).iboh3) = p(obs(i).m(ii).boh3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.m(ii).iboh3) = nan;
                obs(i).m(ii).boh3 = nan;
            end
            if (isgood(obs(i).m(ii).eboh3))
                wobs(i,sys.m(ii).iboh3) = w(obs(i).m(ii).boh3,obs(i).m(ii).eboh3);
            else
                wobs(i,sys.m(ii).iboh3) = nan;
                obs(i).m(ii).eboh3 = nan;
            end
            if (isgood(obs(i).m(ii).boh4))
                yobs(i,sys.m(ii).iboh4) = p(obs(i).m(ii).boh4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.m(ii).iboh4) = nan;
                obs(i).m(ii).boh4 = nan;
            end
            if (isgood(obs(i).m(ii).eboh4))
                wobs(i,sys.m(ii).iboh4) = w(obs(i).m(ii).boh4,obs(i).m(ii).eboh4);
            else
                wobs(i,sys.m(ii).iboh4) = nan;
                obs(i).m(ii).eboh4 = nan;
            end

            % sulfate system
            if (isgood(obs(i).m(ii).epKs))
                wobs(i,sys.m(ii).iKs) = (obs(i).m(ii).epKs).^(-2);
            else
                obs(i).m(ii).epKs = epKs;
                wobs(i,sys.m(ii).iKs) = (obs(i).m(ii).epKs).^(-2); % wKs = 1/(1 + (0.0021/pKsys(6)))^2 ;
            end
            if (isgood(obs(i).m(ii).pKs))
                yobs(i,sys.m(ii).iKs) = obs(i).m(ii).pKs;
            else
                yobs(i,sys.m(ii).iKs) = pKs;
                obs(i).m(ii).pKs = pKs;
            end
            if (isgood(obs(i).m(ii).hso4))
                yobs(i,sys.m(ii).ihso4) = p(obs(i).m(ii).hso4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.m(ii).ihso4) = nan;
                obs(i).m(ii).hso4 = nan;
            end
            if (isgood(obs(i).m(ii).ehso4))
                wobs(i,sys.m(ii).ihso4) = w(obs(i).m(ii).hso4,obs(i).m(ii).ehso4);
            else
                wobs(i,sys.m(ii).ihso4) = nan;
                obs(i).m(ii).ehso4 = nan;
            end
            if (isgood(obs(i).m(ii).so4))
                yobs(i,sys.m(ii).iso4) = p(obs(i).m(ii).so4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.m(ii).iso4) = nan;
                obs(i).m(ii).so4 = nan;
            end
            if (isgood(obs(i).m(ii).eso4))
                wobs(i,sys.m(ii).iso4) = w(obs(i).m(ii).so4,obs(i).m(ii).eso4);
            else
                wobs(i,sys.m(ii).iso4) = nan;
                obs(i).m(ii).eso4 = nan;
            end
            % if (isgood(obs(i).m(ii).phf)) % phfree
            %     yobs(i,sys.m(ii).iphf) = obs(i).m(ii).phf; % hydrogen free
            % else
            %     yobs(i,sys.m(ii).iphf) = nan; % from CO2SYS
            %     obs(i).m(ii).phf = nan; % phfree
            % end
            % if (isgood(obs(i).m(ii).ephf)) % ephfree
            %     wobs(i,sys.m(ii).iphf) = (obs(i).m(ii).ephf)^(-2);
            % else
            %     wobs(i,sys.m(ii).iphf) = nan;
            %     obs(i).m(ii).ephf = nan; % ephfree
            % end

            % fluoride system
            if (isgood(obs(i).m(ii).epKf))
                wobs(i,sys.m(ii).iKf) = (obs(i).m(ii).epKf).^(-2);
            else
                obs(i).m(ii).epKf = epKf;
                wobs(i,sys.m(ii).iKf) = (obs(i).m(ii).epKf).^(-2);   % wKF = 1/(p(1 + 0.02/KF))^2 ; % 2% relative error
            end
            if (isgood(obs(i).m(ii).pKf))
                yobs(i,sys.m(ii).iKf) = obs(i).m(ii).pKf;
            else
                yobs(i,sys.m(ii).iKf) = pKf;
                obs(i).m(ii).pKf = pKf;
            end
            if (isgood(obs(i).m(ii).F))
                yobs(i,sys.m(ii).iF) = p(obs(i).m(ii).F*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.m(ii).iF) = nan;
                obs(i).m(ii).F = nan;
            end
            if (isgood(obs(i).m(ii).eF))
                wobs(i,sys.m(ii).iF) = w(obs(i).m(ii).F,obs(i).m(ii).eF);
            else
                wobs(i,sys.m(ii).iF) = nan;
                obs(i).m(ii).eF = nan;
            end
            if (isgood(obs(i).m(ii).HF))
                yobs(i,sys.m(ii).iHF) = p(obs(i).m(ii).HF*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.m(ii).iHF) = nan;
                obs(i).m(ii).HF = nan;
            end
            if (isgood(obs(i).m(ii).eHF))
                wobs(i,sys.m(ii).iHF) = w(obs(i).m(ii).HF,obs(i).m(ii).eHF);
            else
                wobs(i,sys.m(ii).iHF) = nan;
                obs(i).m(ii).eHF = nan;
            end

            if (ismember('phosphate',opt.abr))
                if (isgood(obs(i).m(ii).epK1p))
                    wobs(i,sys.m(ii).iK1p) = (obs(i).m(ii).epK1p).^(-2);
                else
                    obs(i).m(ii).epK1p = epK1p;
                    wobs(i,sys.m(ii).iK1p) = (obs(i).m(ii).epK1p).^(-2);  % wK1p = 1/(1 + (0.09/pKsys(8)))^2 ;
                end
                if (isgood(obs(i).m(ii).pK1p))
                    yobs(i,sys.m(ii).iK1p) = obs(i).m(ii).pK1p;
                else
                    yobs(i,sys.m(ii).iK1p) = pK1p;
                    obs(i).m(ii).pK1p = pK1p;
                end
                if (isgood(obs(i).m(ii).epK2p))
                    wobs(i,sys.m(ii).iK2p) = (obs(i).m(ii).epK2p).^(-2);
                else
                    obs(i).m(ii).epK2p = epK2p;
                    wobs(i,sys.m(ii).iK2p) = (obs(i).m(ii).epK2p).^(-2);  % wK2p = 1/(1 + (0.03/pKsys(9)))^2 ;
                end
                if (isgood(obs(i).m(ii).pK2p))
                    yobs(i,sys.m(ii).iK2p) = obs(i).m(ii).pK2p;
                else
                    yobs(i,sys.m(ii).iK2p) = pK2p;
                    obs(i).m(ii).pK2p = pK2p;
                end
                if (isgood(obs(i).m(ii).epK3p))
                    wobs(i,sys.m(ii).iK3p) = (obs(i).m(ii).epK3p).^(-2);
                else
                    obs(i).m(ii).epK3p = epK3p;
                    wobs(i,sys.m(ii).iK3p) = (obs(i).m(ii).epK3p).^(-2);  % wK3p = 1/(1 + (0.02/pKsys(10)))^2 ;
                end
                if (isgood(obs(i).m(ii).pK3p))
                    yobs(i,sys.m(ii).iK3p) = obs(i).m(ii).pK3p;
                else
                    yobs(i,sys.m(ii).iK3p) = pK3p;
                    obs(i).m(ii).pK3p = pK3p;
                end
                if (isgood(obs(i).m(ii).h3po4))
                    yobs(i,sys.m(ii).ih3po4) = p(obs(i).m(ii).h3po4*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.m(ii).ih3po4) = nan;
                    obs(i).m(ii).h3po4 = nan;
                end
                if (isgood(obs(i).m(ii).eh3po4))
                    wobs(i,sys.m(ii).ih3po4) = w(obs(i).m(ii).h3po4,obs(i).m(ii).eh3po4);
                else
                    wobs(i,sys.m(ii).ih3po4) = nan;
                    obs(i).m(ii).eh3po4 = nan;
                end
                if (isgood(obs(i).m(ii).h2po4))
                    yobs(i,sys.m(ii).ih2po4) = p(obs(i).m(ii).h2po4*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.m(ii).ih2po4) = nan;
                    obs(i).m(ii).h2po4 = nan;
                end
                if (isgood(obs(i).m(ii).eh2po4))
                    wobs(i,sys.m(ii).ih2po4) = w(obs(i).m(ii).h2po4,obs(i).m(ii).eh2po4);
                else
                    wobs(i,sys.m(ii).ih2po4) = nan;
                    obs(i).m(ii).eh2po4 = nan;
                end
                if (isgood(obs(i).m(ii).hpo4))
                    yobs(i,sys.m(ii).ihpo4) = p(obs(i).m(ii).hpo4*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.m(ii).ihpo4) = nan;
                    obs(i).m(ii).hpo4 = nan;
                end
                if (isgood(obs(i).m(ii).ehpo4))
                    wobs(i,sys.m(ii).ihpo4) = w(obs(i).m(ii).hpo4,obs(i).m(ii).ehpo4);
                else
                    wobs(i,sys.m(ii).ihpo4) = nan;
                    obs(i).m(ii).ehpo4 = nan;
                end
                if (isgood(obs(i).m(ii).po4))
                    yobs(i,sys.m(ii).ipo4) = p(obs(i).m(ii).po4*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.m(ii).ipo4) = nan;
                    obs(i).m(ii).po4 = nan;
                end
                if (isgood(obs(i).m(ii).epo4))
                    wobs(i,sys.m(ii).ipo4) = w(obs(i).m(ii).po4,obs(i).m(ii).epo4);
                else
                    wobs(i,sys.m(ii).ipo4) = nan;
                    obs(i).m(ii).epo4 = nan;
                end

            end

            if (ismember('silicate',opt.abr))
                if (isgood(obs(i).m(ii).epKsi))
                    wobs(i,sys.m(ii).iKsi) = (obs(i).m(ii).epKsi).^(-2);
                else
                    obs(i).m(ii).epKsi = epKsi;
                    wobs(i,sys.m(ii).iKsi) = (obs(i).m(ii).epKsi).^(-2);  % wKSi = 1/(1 + (0.02/pKsys(11)))^2 ;
                end
                if (isgood(obs(i).m(ii).pKsi))
                    yobs(i,sys.m(ii).iKsi) = obs(i).m(ii).pKsi;
                else
                    yobs(i,sys.m(ii).iKsi) = pKsi;
                    obs(i).m(ii).pKsi = pKsi;
                end
                if (isgood(obs(i).m(ii).sioh4))
                    yobs(i,sys.m(ii).isioh4) = p(obs(i).m(ii).sioh4*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.m(ii).isioh4) = nan;
                    obs(i).m(ii).sioh4 = nan;
                end
                if (isgood(obs(i).m(ii).esioh4))
                    wobs(i,sys.m(ii).isioh4) = w(obs(i).m(ii).sioh4,obs(i).m(ii).esioh4);
                else
                    wobs(i,sys.m(ii).isioh4) = nan;
                    obs(i).m(ii).esioh4 = nan;
                end
                if (isgood(obs(i).m(ii).siooh3))
                    yobs(i,sys.m(ii).isiooh3) = p(obs(i).m(ii).siooh3*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.m(ii).isiooh3) = nan;
                    obs(i).m(ii).siooh3 = nan;
                end
                if (isgood(obs(i).m(ii).esiooh3))
                    wobs(i,sys.m(ii).isiooh3) = w(obs(i).m(ii).siooh3,obs(i).m(ii).esiooh3);
                else
                    wobs(i,sys.m(ii).isiooh3) = nan;
                    obs(i).m(ii).esiooh3 = nan;
                end

            end

            if (ismember('ammonia',opt.abr))
                if (isgood(obs(i).m(ii).epKnh4))
                    wobs(i,sys.m(ii).iKnh4) = (obs(i).m(ii).epKnh4).^(-2);
                else
                    obs(i).m(ii).epKnh4 = epKnh4;
                    wobs(i,sys.m(ii).iKnh4) = (obs(i).m(ii).epKnh4).^(-2);  % wKnh4 = 1/(1 + (0.00017/pKsys(11)))^2 ;
                end
                if (isgood(obs(i).m(ii).pKnh4))
                    yobs(i,sys.m(ii).iKnh4) = obs(i).m(ii).pKnh4;
                else
                    yobs(i,sys.m(ii).iKnh4) = pKnh4;
                    obs(i).m(ii).pKnh4 = pKnh4;
                end
                if (isgood(obs(i).m(ii).nh4))
                    yobs(i,sys.m(ii).inh4) = p(obs(i).m(ii).nh4*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.m(ii).inh4) = nan;
                    obs(i).m(ii).nh4 = nan;
                end
                if (isgood(obs(i).m(ii).enh4))
                    wobs(i,sys.m(ii).inh4) = w(obs(i).m(ii).nh4,obs(i).m(ii).enh4);
                else
                    wobs(i,sys.m(ii).inh4) = nan;
                    obs(i).m(ii).enh4 = nan;
                end
                if (isgood(obs(i).m(ii).nh3))
                    yobs(i,sys.m(ii).inh3) = p(obs(i).m(ii).nh3*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.m(ii).inh3) = nan;
                    obs(i).m(ii).nh3 = nan;
                end
                if (isgood(obs(i).m(ii).enh3))
                    wobs(i,sys.m(ii).inh3) = w(obs(i).m(ii).nh3,obs(i).m(ii).enh3);
                else
                    wobs(i,sys.m(ii).inh3) = nan;
                    obs(i).m(ii).enh3 = nan;
                end
            end

            if (ismember('sulfide',opt.abr))
                if (isgood(obs(i).m(ii).epKh2s))
                    wobs(i,sys.m(ii).iKh2s) = (obs(i).m(ii).epKh2s).^(-2);
                else
                    obs(i).m(ii).epKh2s = epKh2s;
                    wobs(i,sys.m(ii).iKh2s) = (obs(i).m(ii).epKh2s).^(-2);  % wKh2s = 1/(1 + (0.033/pKsys(11)))^2 ;
                end
                if (isgood(obs(i).m(ii).pKh2s))
                    yobs(i,sys.m(ii).iKh2s) = obs(i).m(ii).pKh2s;
                else
                    yobs(i,sys.m(ii).iKh2s) = pKh2s;
                    obs(i).m(ii).pKh2s = pKh2s;
                end
                if (isgood(obs(i).m(ii).h2s))
                    yobs(i,sys.m(ii).ih2s) = p(obs(i).m(ii).h2s*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.m(ii).ih2s) = nan;
                    obs(i).m(ii).h2s = nan;
                end
                if (isgood(obs(i).m(ii).eh2s))
                    wobs(i,sys.m(ii).ih2s) = w(obs(i).m(ii).h2s,obs(i).m(ii).eh2s);
                else
                    wobs(i,sys.m(ii).ih2s) = nan;
                    obs(i).m(ii).eh2s = nan;
                end
                if (isgood(obs(i).m(ii).hs))
                    yobs(i,sys.m(ii).ihs) = p(obs(i).m(ii).hs*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.m(ii).ihs) = nan;
                    obs(i).m(ii).hs = nan;
                end
                if (isgood(obs(i).m(ii).ehs))
                    wobs(i,sys.m(ii).ihs) = w(obs(i).m(ii).hs,obs(i).m(ii).ehs);
                else
                    wobs(i,sys.m(ii).ihs) = nan;
                    obs(i).m(ii).ehs = nan;
                end
            end

            if (ismember('solubility',opt.abr))
                if (isgood(obs(i).m(ii).epKar))
                    wobs(i,sys.m(ii).iKar) = (obs(i).m(ii).epKar).^(-2);
                else
                    obs(i).m(ii).epKar = epKar;
                    wobs(i,sys.m(ii).iKar) = (obs(i).m(ii).epKar).^(-2);
                end
                if (isgood(obs(i).m(ii).pKar))
                    yobs(i,sys.m(ii).iKar) = obs(i).m(ii).pKar;
                else
                    yobs(i,sys.m(ii).iKar) = pKar;
                    obs(i).m(ii).pKar = pKar;
                end
                if (isgood(obs(i).m(ii).ca))
                    yobs(i,sys.m(ii).ica) = p(obs(i).m(ii).ca*1e-6); % convt µmol/kg to mol/kg
                else
                    yobs(i,sys.m(ii).ica) = nan;
                    obs(i).m(ii).ca = nan;
                end
                if (isgood(obs(i).m(ii).eca))
                    wobs(i,sys.m(ii).ica) = w(obs(i).m(ii).ca,obs(i).m(ii).eca);
                else
                    wobs(i,sys.m(ii).ica) = nan;
                    obs(i).m(ii).eca = nan;
                end
                if (isgood(obs(i).m(ii).OmegaAr))
                    yobs(i,sys.m(ii).iOmegaAr) = p(obs(i).m(ii).OmegaAr); % Omega is dimensionless
                else
                    yobs(i,sys.m(ii).iOmegaAr) = nan;
                    obs(i).m(ii).OmegaAr = nan;
                end
                if (isgood(obs(i).m(ii).eOmegaAr))
                    wobs(i,sys.m(ii).iOmegaAr) = w(obs(i).m(ii).OmegaAr,obs(i).m(ii).eOmegaAr);
                else
                    wobs(i,sys.m(ii).iOmegaAr) = nan;
                    obs(i).m(ii).eOmegaAr = nan;
                end
                if (isgood(obs(i).m(ii).epKca))
                    wobs(i,sys.m(ii).iKca) = (obs(i).m(ii).epKca).^(-2);
                else
                    obs(i).m(ii).epKca = epKca;
                    wobs(i,sys.m(ii).iKca) = (obs(i).m(ii).epKca).^(-2);
                end
                if (isgood(obs(i).m(ii).pKca))
                    yobs(i,sys.m(ii).iKca) = obs(i).m(ii).pKca;
                else
                    yobs(i,sys.m(ii).iKca) = pKca;
                    obs(i).m(ii).pKca = pKca;
                end
                if (isgood(obs(i).m(ii).OmegaCa))
                    yobs(i,sys.m(ii).iOmegaCa) = p(obs(i).m(ii).OmegaCa); % Omega is dimensionless
                else
                    yobs(i,sys.m(ii).iOmegaCa) = nan;
                    obs(i).m(ii).OmegaCa = nan;
                end
                if (isgood(obs(i).m(ii).eOmegaCa))
                    wobs(i,sys.m(ii).iOmegaCa) = w(obs(i).m(ii).OmegaCa,obs(i).m(ii).eOmegaCa);
                else
                    wobs(i,sys.m(ii).iOmegaCa) = nan;
                    obs(i).m(ii).eOmegaCa = nan;
                end
            end
        end
    end
end

% --------------------------------------------------------------------------------

function [est] = parse_output(z,sigy,opt,sys)
    % populate est
    %
    % INPUT:
    %
    %   z  := system state including equilibrium constants and lagrange multipliers
    
    p = sys.p;
    q = sys.q;

    ebar = @(j) (0.5 * ( q( z(j) - sigy(j) ) - q( z(j) + sigy(j) ) ) );
    ebar_l = @(j) ( q( -sigy(j) ) ); % lower sigma
    ebar_u = @(j) ( q( sigy(j) ) ); % upper sigma
        
    % populate 'est' structure with best estimate:
    %   1. p(value) and p(error) where p(x) = -log10(x)
    %   2. value and average error about the value in 'q' 
    %           where q(x) = x^(-10)
    %   3. upper and lower bounds in 'q' space, not symmetric
    %           about the value in 'q' space

    est.sal = z(sys.isal);
    est.esal = sigy(sys.isal);
    
    % TC
    est.pTC = z(sys.iTC);               
    est.epTC = sigy(sys.iTC);    
    est.TC = q(z(sys.iTC))*1e6; % 1e6  converts mol/kg to µmol/kg
    est.eTC = ebar(sys.iTC)*1e6;      
    est.eTC_l = ebar_l(sys.iTC)*1e6;    
    est.eTC_u = ebar_u(sys.iTC)*1e6;
    
    % TA
    est.pTA = z(sys.iTA);               
    est.epTA = sigy(sys.iTA);
    est.TA = q(z(sys.iTA))*1e6; % convt 
    est.eTA = ebar(sys.iTA)*1e6;      
    est.eTA_l = ebar_l(sys.iTA)*1e6;    
    est.eTA_u = ebar_u(sys.iTA)*1e6;
    
    % TB borate
    est.pTB = z(sys.iTB);
    est.epTB = sigy(sys.iTB);
    est.TB = q(z(sys.iTB))*1e6; % convt mol/kg to µmol/kg
    est.eTB = ebar(sys.iTB)*1e6;
    est.eTB_l = ebar_l(sys.iTB)*1e6;
    est.eTB_u = ebar_u(sys.iTB)*1e6;

    % TS sulfate
    est.pTS = z(sys.iTS);
    est.epTS = sigy(sys.iTS);
    est.TS = q(z(sys.iTS));  % no conversion on TS
    est.eTS = ebar(sys.iTS);
    est.eTS_l = ebar_l(sys.iTS);
    est.eTS_u = ebar_u(sys.iTS);

    % TF fluoride
    est.pTF = z(sys.iTF);
    est.epTF = sigy(sys.iTF);
    est.TF = q(z(sys.iTF))*1e6; % convt mol/kg to µmol/kg
    est.eTF = ebar(sys.iTF)*1e6;
    est.eTF_l = ebar_l(sys.iTF)*1e6;
    est.eTF_u = ebar_u(sys.iTF)*1e6;
    
    if ismember('phosphate', opt.abr)
        % TP
        est.pTP = z(sys.iTP);           
        est.epTP = sigy(sys.iTP);
        est.TP = q(z(sys.iTP))*1e6;   % convt mol/kg to µmol/kg   
        est.eTP = ebar(sys.iTP)*1e6;
        est.eTP_l = ebar_l(sys.iTP)*1e6; 
        est.eTP_u = ebar_u(sys.iTP)*1e6;
    end
    if ismember('silicate', opt.abr)
        % TSi
        est.pTSi = z(sys.iTSi);         
        est.epTSi = sigy(sys.iTSi);
        est.TSi = q(z(sys.iTSi))*1e6;  % convt mol/kg to µmol/kg 
        est.eTSi = ebar(sys.iTSi)*1e6;
        est.eTSi_l = ebar_l(sys.iTSi)*1e6; 
        est.eTSi_u = ebar_u(sys.iTSi)*1e6;
    end
    if ismember('ammonia', opt.abr)
        % TNH3
        est.pTNH3 = z(sys.iTNH3);       
        est.epTNH3 = sigy(sys.iTNH3);
        est.TNH3 = q(z(sys.iTNH3))*1e6;   % convt mol/kg to µmol/kg
        est.eTNH3 = ebar(sys.iTNH3)*1e6;
        est.eTNH3_l = ebar_l(sys.iTNH3)*1e6; 
        est.eTNH3_u = ebar_u(sys.iTNH3)*1e6;
    end
    if ismember('sulfide', opt.abr)
        % TH2S
        est.pTH2S = z(sys.iTH2S);       
        est.epTH2S = sigy(sys.iTH2S);
        est.TH2S = q(z(sys.iTH2S))*1e6; % convt mol/kg to µmol/kg
        est.eTH2S = ebar(sys.iTH2S)*1e6;
        est.eTH2S_l = ebar_l(sys.iTH2S)*1e6; 
        est.eTH2S_u = ebar_u(sys.iTH2S)*1e6;
    end
    if ismember('solubility', opt.abr)
        % TCal
        est.pTCal = z(sys.iTCal);       
        est.epTCal = sigy(sys.iTCal);
        est.TCal = q(z(sys.iTCal)); % no conversion, mol/kg
        est.eTCal = ebar(sys.iTCal);
        est.eTCal_l = ebar_l(sys.iTCal); 
        est.eTCal_u = ebar_u(sys.iTCal);
    end
    
    nTP = length(sys.m);
    for i = 1:nTP
        % temp (deg C)
        est.m(i).T     = z(sys.m(i).iT);
        est.m(i).eT    = sigy(sys.m(i).iT);
        est.m(i).eT_l  = z(sys.m(i).iT)-sigy(sys.m(i).iT);
        est.m(i).eT_u  = z(sys.m(i).iT)+sigy(sys.m(i).iT);
        
        % pressure (dbar)
        est.m(i).P     = z(sys.m(i).iP);        
        est.m(i).eP    = sigy(sys.m(i).iP);
        est.m(i).eP_l  = z(sys.m(i).iP)-sigy(sys.m(i).iP);
        est.m(i).eP_u  = z(sys.m(i).iP)+sigy(sys.m(i).iP);
        
        % pH
        est.m(i).ph    = z(sys.m(i).iph);
        est.m(i).eph   = sigy(sys.m(i).iph);
        % output pH on all scales
        ph_all = phscales(est.m(i).ph, opt.phscale, ... % pH_in, pHscale_in
            est.TS, q(z(sys.m(i).iKs)), est.TF, ... % TS, Ks, TF
            q(z(sys.m(i).iKf)), z(sys.m(i).ipfH) ); % Kf, pfH
            % q(z(sys.m(i).iKf)), z(sys.m(i).iphf) ); % Kf, phf

        est.m(i).ph      = ph_all(1); % ph_tot is default
        est.m(i).ph_sws  = ph_all(2);
        est.m(i).ph_free = ph_all(3);
        est.m(i).ph_nbs  = ph_all(4);

    % ebar = @(j) (0.5 * ( q( z(j) - sigy(j) ) - q( z(j) + sigy(j) ) ) );

        % H (free) = q(ph_free)
        est.m(i).h_free    = q(z(sys.m(i).iph_free))*1e6;
        est.m(i).eh_free = (0.5 * (q (z(sys.m(i).iph_free) - ...
            sigy(sys.m(i).iph) ) - ( q( z(sys.m(i).iph_free) + ...
            sigy(sys.m(i).iph) ) ) ) )  * 1e6; % eph_free = eph
        est.m(i).eh_free_l = (q( -sigy(sys.m(i).iph) ) ) * 1e6;
        est.m(i).eh_free_u = (q( sigy(sys.m(i).iph) ) ) * 1e6;
        est.m(i).ph_free   = z(sys.m(i).iph_free);
        est.m(i).eph_free  = sigy(sys.m(i).iph);

        % fH = activity coefficient
        est.m(i).fH    = q(z(sys.m(i).ipfH))*1e6; % hfree & phfree
        est.m(i).efH   = ebar(sys.m(i).ipfH)*1e6;
        est.m(i).efH_l = ebar_l(sys.m(i).ipfH)*1e6;
        est.m(i).efH_u = ebar_u(sys.m(i).ipfH)*1e6;
        est.m(i).pfH   = z(sys.m(i).ipfH);
        est.m(i).epfH  = sigy(sys.m(i).ipfH); % hfree & phfree

        % fCO2
        est.m(i).fco2    = q(z(sys.m(i).ifco2))*1e6; % convt atm to µatm
        est.m(i).efco2   = ebar(sys.m(i).ifco2)*1e6;
        est.m(i).efco2_l = ebar_l(sys.m(i).ifco2)*1e6;
        est.m(i).efco2_u = ebar_u(sys.m(i).ifco2)*1e6;
        est.m(i).pfco2   = z(sys.m(i).ifco2);
        est.m(i).epfco2  = sigy(sys.m(i).ifco2);
        
        % pCO2
        est.m(i).pco2    = q(z(sys.m(i).ipco2))*1e6; % convt atm to µatm
        est.m(i).epco2   = ebar(sys.m(i).ipco2)*1e6;
        est.m(i).epco2_l = ebar_l(sys.m(i).ipco2)*1e6;
        est.m(i).epco2_u = ebar_u(sys.m(i).ipco2)*1e6;
        est.m(i).ppco2   = z(sys.m(i).ipco2);
        est.m(i).eppco2  = sigy(sys.m(i).ipco2);
        
        % HCO3
        est.m(i).hco3    = q(z(sys.m(i).ihco3))*1e6; % convt mol/kg to µmol/kg
        est.m(i).ehco3   = ebar(sys.m(i).ihco3)*1e6;
        est.m(i).ehco3_l = ebar_l(sys.m(i).ihco3)*1e6;
        est.m(i).ehco3_u = ebar_u(sys.m(i).ihco3)*1e6;
        est.m(i).phco3   = z(sys.m(i).ihco3);
        est.m(i).ephco3  = sigy(sys.m(i).ihco3);

        % CO3
        est.m(i).co3    = q(z(sys.m(i).ico3))*1e6;
        est.m(i).eco3   = ebar(sys.m(i).ico3)*1e6;
        est.m(i).eco3_l = ebar_l(sys.m(i).ico3)*1e6;
        est.m(i).eco3_u = ebar_u(sys.m(i).ico3)*1e6;
        est.m(i).pco3   = z(sys.m(i).ico3);
        est.m(i).epco3  = sigy(sys.m(i).ico3);

        % pCO2*
        est.m(i).pco2st   = z(sys.m(i).ico2st);
        est.m(i).epco2st  = sigy(sys.m(i).ico2st);
        est.m(i).co2st    = q(z(sys.m(i).ico2st));
        est.m(i).eco2st   = ebar(sys.m(i).ico2st);
        est.m(i).eco2st_l = ebar_l(sys.m(i).ico2st);
        est.m(i).eco2st_u = ebar_u(sys.m(i).ico2st);

        % pP2F (aka FugFac)
        est.m(i).pp2f   = z(sys.m(i).ip2f);
        est.m(i).epp2f  = sigy(sys.m(i).ip2f);
        est.m(i).p2f    = q(z(sys.m(i).ip2f));
        est.m(i).ep2f   = ebar(sys.m(i).ip2f);
        est.m(i).ep2f_l = ebar_l(sys.m(i).ip2f);
        est.m(i).ep2f_u = ebar_u(sys.m(i).ip2f);
        
        % pK0
        est.m(i).pK0   = z(sys.m(i).iK0);
        est.m(i).epK0  = sigy(sys.m(i).iK0);
        est.m(i).K0    = q(z(sys.m(i).iK0));
        est.m(i).eK0   = ebar(sys.m(i).iK0);
        est.m(i).eK0_l = ebar_l(sys.m(i).iK0);
        est.m(i).eK0_u = ebar_u(sys.m(i).iK0);
        
        % pK1
        est.m(i).pK1   = z(sys.m(i).iK1);
        est.m(i).epK1  = sigy(sys.m(i).iK1);
        est.m(i).K1    = q(z(sys.m(i).iK1));
        est.m(i).eK1   = ebar(sys.m(i).iK1);
        est.m(i).eK1_l = ebar_l(sys.m(i).iK1); 
        est.m(i).eK1_u = ebar_u(sys.m(i).iK1);
        
        % pK2
        est.m(i).pK2   = z(sys.m(i).iK2);
        est.m(i).epK2  = sigy(sys.m(i).iK2);
        est.m(i).K2    = q(z(sys.m(i).iK2));
        est.m(i).eK2   = ebar(sys.m(i).iK2);
        est.m(i).eK2_l = ebar_l(sys.m(i).iK2);
        est.m(i).eK2_u = ebar_u(sys.m(i).iK2);

        % OH
        est.m(i).oh    = q(z(sys.m(i).ioh))*1e6; % convt
        est.m(i).eoh   = ebar(sys.m(i).ioh)*1e6;
        est.m(i).eoh_l = ebar_l(sys.m(i).ioh)*1e6;
        est.m(i).eoh_u = ebar_u(sys.m(i).ioh)*1e6;
        est.m(i).poh   = z(sys.m(i).ioh);
        est.m(i).epoh  = sigy(sys.m(i).ioh);

        % pKw
        est.m(i).pKw   = z(sys.m(i).iKw);
        est.m(i).epKw  = sigy(sys.m(i).iKw);
        est.m(i).Kw    = q(z(sys.m(i).iKw));
        est.m(i).eKw   = ebar(sys.m(i).iKw);
        est.m(i).eKw_l = ebar_l(sys.m(i).iKw);
        est.m(i).eKw_u = ebar_u(sys.m(i).iKw);

        % BOH4 borate
        est.m(i).boh4    = q(z(sys.m(i).iboh4))*1e6; % convt mol/kg to µmol/kg
        est.m(i).eboh4   = ebar(sys.m(i).iboh4)*1e6;
        est.m(i).eboh4_l = ebar_l(sys.m(i).iboh4)*1e6;
        est.m(i).eboh4_u = ebar_u(sys.m(i).iboh4)*1e6;
        est.m(i).pboh4   = z(sys.m(i).iboh4);
        est.m(i).epboh4  = sigy(sys.m(i).iboh4);

        % BOH3
        est.m(i).boh3    = q(z(sys.m(i).iboh3))*1e6;
        est.m(i).eboh3   = ebar(sys.m(i).iboh3)*1e6;
        est.m(i).eboh3_l = ebar_l(sys.m(i).iboh3)*1e6;
        est.m(i).eboh3_u = ebar_u(sys.m(i).iboh3)*1e6;
        est.m(i).pboh3   = z(sys.m(i).iboh3);
        est.m(i).epboh3  = sigy(sys.m(i).iboh3);

        % pKb
        est.m(i).pKb   = z(sys.m(i).iKb);
        est.m(i).epKb  = sigy(sys.m(i).iKb);
        est.m(i).Kb    = q(z(sys.m(i).iKb));
        est.m(i).eKb   = ebar(sys.m(i).iKb);
        est.m(i).eKb_l = ebar_l(sys.m(i).iKb);
        est.m(i).eKb_u = ebar_u(sys.m(i).iKb);

        % SO4 sulfide
        est.m(i).so4    = q(z(sys.m(i).iso4))*1e6; % convt
        est.m(i).eso4   = ebar(sys.m(i).iso4)*1e6;
        est.m(i).eso4_l = ebar_l(sys.m(i).iso4)*1e6;
        est.m(i).eso4_u = ebar_u(sys.m(i).iso4)*1e6;
        est.m(i).pso4   = z(sys.m(i).iso4);
        est.m(i).epso4  = sigy(sys.m(i).iso4);

        % HSO4
        est.m(i).hso4    = q(z(sys.m(i).ihso4))*1e6;
        est.m(i).ehso4   = ebar(sys.m(i).ihso4)*1e6;
        est.m(i).ehso4_l = ebar_l(sys.m(i).ihso4)*1e6;
        est.m(i).ehso4_u = ebar_u(sys.m(i).ihso4)*1e6;
        est.m(i).phso4   = z(sys.m(i).ihso4);
        est.m(i).ephso4  = sigy(sys.m(i).ihso4);

        % pKs
        est.m(i).pKs   = z(sys.m(i).iKs);
        est.m(i).epKs  = sigy(sys.m(i).iKs);
        est.m(i).Ks    = q(z(sys.m(i).iKs));
        est.m(i).eKs   = ebar(sys.m(i).iKs);
        est.m(i).eKs_l = ebar_l(sys.m(i).iKs);
        est.m(i).eKs_u = ebar_u(sys.m(i).iKs);

        % F fluoride
        est.m(i).F    = q(z(sys.m(i).iF))*1e6; % convt
        est.m(i).eF   = ebar(sys.m(i).iF)*1e6;
        est.m(i).eF_l = ebar_l(sys.m(i).iF)*1e6;
        est.m(i).ef_u = ebar_u(sys.m(i).iF)*1e6;
        est.m(i).pF   = z(sys.m(i).iF);
        est.m(i).epF  = sigy(sys.m(i).iF);

        % HF
        est.m(i).HF    = q(z(sys.m(i).iHF))*1e6;
        est.m(i).eHF   = ebar(sys.m(i).iHF)*1e6;
        est.m(i).eHF_l = ebar_l(sys.m(i).iHF)*1e6;
        est.m(i).eHF_u = ebar_u(sys.m(i).iHF)*1e6;
        est.m(i).pHF   = z(sys.m(i).iHF);
        est.m(i).epHF  = sigy(sys.m(i).iHF);

        % pKf
        est.m(i).pKf   = z(sys.m(i).iKf);
        est.m(i).epKf  = sigy(sys.m(i).iKf);
        est.m(i).Kf    = q(z(sys.m(i).iKf));
        est.m(i).eKf   = ebar(sys.m(i).iKf);
        est.m(i).eKf_l = ebar_l(sys.m(i).iKf);
        est.m(i).eKf_u = ebar_u(sys.m(i).iKf);

        if ismember('phosphate', opt.abr)
            % PO4
            est.m(i).po4    = q(z(sys.m(i).ipo4))*1e6; % convt
            est.m(i).epo4   = ebar(sys.m(i).ipo4)*1e6;
            est.m(i).epo4_l = ebar_l(sys.m(i).ipo4)*1e6;
            est.m(i).epo4_u = ebar_u(sys.m(i).ipo4)*1e6;
            est.m(i).ppo4   = z(sys.m(i).ipo4);
            est.m(i).eppo4  = sigy(sys.m(i).ipo4);

            % HPO4
            est.m(i).hpo4    = q(z(sys.m(i).ihpo4))*1e6;
            est.m(i).ehpo4   = ebar(sys.m(i).ihpo4)*1e6;
            est.m(i).ehpo4_l = ebar_l(sys.m(i).ihpo4)*1e6;
            est.m(i).ehpo4_u = ebar_u(sys.m(i).ihpo4)*1e6;
            est.m(i).phpo4   = z(sys.m(i).ihpo4);
            est.m(i).ephpo4  = sigy(sys.m(i).ihpo4);

            % H2PO4
            est.m(i).h2po4    = q(z(sys.m(i).ih2po4))*1e6;
            est.m(i).eh2po4   = ebar(sys.m(i).ih2po4)*1e6;
            est.m(i).eh2po4_l = ebar_l(sys.m(i).ih2po4)*1e6;
            est.m(i).eh2po4_u = ebar_u(sys.m(i).ih2po4)*1e6;
            est.m(i).ph2po4   = z(sys.m(i).ih2po4);
            est.m(i).eph2po4  = sigy(sys.m(i).ih2po4);            

            % H3PO4
            est.m(i).h3po4    = q(z(sys.m(i).ih3po4))*1e6;
            est.m(i).eh3po4   = ebar(sys.m(i).ih3po4)*1e6;
            est.m(i).eh3po4_l = ebar_l(sys.m(i).ih3po4)*1e6;
            est.m(i).eh3po4_u = ebar_u(sys.m(i).ih3po4)*1e6;
            est.m(i).ph3po4   = z(sys.m(i).ih3po4);
            est.m(i).eph3po4  = sigy(sys.m(i).ih3po4);

            % pK1p
            est.m(i).pK1p   = z(sys.m(i).iK1p);
            est.m(i).epK1p  = sigy(sys.m(i).iK1p);
            est.m(i).K1p    = q(z(sys.m(i).iK1p));
            est.m(i).eK1p   = ebar(sys.m(i).iK1p);
            est.m(i).eK1p_l = ebar_l(sys.m(i).iK1p);
            est.m(i).eK1p_u = ebar_u(sys.m(i).iK1p);

            % pK2p
            est.m(i).pK2p   = z(sys.m(i).iK2p);
            est.m(i).epK2p  = sigy(sys.m(i).iK2p);
            est.m(i).K2p    = q(z(sys.m(i).iK2p));
            est.m(i).eK2p   = ebar(sys.m(i).iK2p);
            est.m(i).pK2p_l = ebar_l(sys.m(i).iK2p);
            est.m(i).eK2p_u = ebar_u(sys.m(i).iK2p);

            % pK3p
            est.m(i).pK3p   = z(sys.m(i).iK3p);
            est.m(i).epK3p  = sigy(sys.m(i).iK3p);
            est.m(i).K3p    = q(z(sys.m(i).iK3p));
            est.m(i).eK3p   = ebar(sys.m(i).iK3p);
            est.m(i).eK3p_l = ebar_l(sys.m(i).iK3p);
            est.m(i).eK3p_u = ebar_u(sys.m(i).iK3p);
        end

        if ismember('silicate', opt.abr)
            % SiOH4
            est.m(i).sioh4     = q(z(sys.m(i).isioh4))*1e6; % convt
            est.m(i).esioh4    = ebar(sys.m(i).isioh4)*1e6;
            est.m(i).esioh4_l  = ebar_l(sys.m(i).isioh4)*1e6;
            est.m(i).esioh4_u  = ebar_u(sys.m(i).isioh4)*1e6;
            est.m(i).psioh4    = z(sys.m(i).isioh4);
            est.m(i).epsioh4   = sigy(sys.m(i).isioh4);

            % SiOH3
            est.m(i).siooh3     = q(z(sys.m(i).isiooh3))*1e6;
            est.m(i).esiooh3    = ebar(sys.m(i).isiooh3)*1e6;
            est.m(i).esiooh3_l  = ebar_l(sys.m(i).isiooh3)*1e6;
            est.m(i).esiooh3_u  = ebar_u(sys.m(i).isiooh3)*1e6;
            est.m(i).psiooh3    = z(sys.m(i).isiooh3);
            est.m(i).epsiooh3   = sigy(sys.m(i).isiooh3);

            % pKsi
            est.m(i).pKsi    = z(sys.m(i).iKsi);
            est.m(i).epKsi   = sigy(sys.m(i).iKsi);
            est.m(i).Ksi     = q(z(sys.m(i).iKsi));
            est.m(i).eKsi    = ebar(sys.m(i).iKsi);
            est.m(i).eKsi_l  = ebar_l(sys.m(i).iKsi);
            est.m(i).eKsi_u  = ebar_u(sys.m(i).iKsi);
        end

        if ismember('ammonia', opt.abr)
            % NH3
            est.m(i).nh3     = q(z(sys.m(i).inh3))*1e6; % convt
            est.m(i).enh3    = ebar(sys.m(i).inh3)*1e6;
            est.m(i).enh3_l  = ebar_l(sys.m(i).inh3)*1e6;
            est.m(i).enh3_u  = ebar_u(sys.m(i).inh3)*1e6;
            est.m(i).pnh3    = z(sys.m(i).inh3);
            est.m(i).epnh3   = sigy(sys.m(i).inh3);

            % NH4
            est.m(i).nh4     = q(z(sys.m(i).inh4))*1e6;
            est.m(i).enh4    = ebar(sys.m(i).inh4)*1e6;
            est.m(i).enh4_l  = ebar_l(sys.m(i).inh4)*1e6;
            est.m(i).enh4_u  = ebar_u(sys.m(i).inh4)*1e6;
            est.m(i).pnh4    = z(sys.m(i).inh4);
            est.m(i).epnh4   = sigy(sys.m(i).inh4);

            % pKNH4
            est.m(i).pKnh4    = z(sys.m(i).iKnh4);
            est.m(i).epKnh4   = sigy(sys.m(i).iKnh4);
            est.m(i).Knh4     = q(z(sys.m(i).iKnh4));
            est.m(i).eKnh4    = ebar(sys.m(i).iKnh4);
            est.m(i).eKnh4_l  = ebar_l(sys.m(i).iKnh4);
            est.m(i).eKnh4_u  = ebar_u(sys.m(i).iKnh4);
        end

        if ismember('sulfide', opt.abr)
            % HS
            est.m(i).hs     = q(z(sys.m(i).ihs))*1e6; % convt
            est.m(i).ehs    = ebar(sys.m(i).ihs)*1e6;
            est.m(i).ehs_l  = ebar_l(sys.m(i).ihs)*1e6;
            est.m(i).ehs_u  = ebar_u(sys.m(i).ihs)*1e6;
            est.m(i).phs    = z(sys.m(i).ihs);
            est.m(i).ephs   = sigy(sys.m(i).ihs);

            % H2S
            est.m(i).h2s     = q(z(sys.m(i).ih2s))*1e6;
            est.m(i).eh2s    = ebar(sys.m(i).ih2s)*1e6;
            est.m(i).eh2s_l  = ebar_l(sys.m(i).ih2s)*1e6;
            est.m(i).ehs2_u  = ebar_u(sys.m(i).ih2s)*1e6;
            est.m(i).ph2s    = z(sys.m(i).ih2s);
            est.m(i).eph2s   = sigy(sys.m(i).ih2s);

            % pKh2s
            est.m(i).pKh2s    = z(sys.m(i).iKh2s);
            est.m(i).epKh2s   = sigy(sys.m(i).iKh2s);
            est.m(i).Kh2s     = q(z(sys.m(i).iKh2s));
            est.m(i).eKh2s    = ebar(sys.m(i).iKh2s);
            est.m(i).eKh2s_l  = ebar_l(sys.m(i).iKh2s);
            est.m(i).eKh2s_u  = ebar_u(sys.m(i).iKh2s);
        end

        if ismember('solubility',opt.abr)
            % Ca
            est.m(i).ca     = q(z(sys.m(i).ica))*1e6;
            est.m(i).eca    = ebar(sys.m(i).ica)*1e6;
            est.m(i).eca_l  = ebar_l(sys.m(i).ica)*1e6;
            est.m(i).eca_u  = ebar_u(sys.m(i).ica)*1e6;
            est.m(i).pca    = z(sys.m(i).ica);
            est.m(i).epca   = sigy(sys.m(i).ica);

            % Omega_Ar
            est.m(i).OmegaAr     = q(z(sys.m(i).iOmegaAr)); % unitless
            est.m(i).eOmegaAr    = ebar(sys.m(i).iOmegaAr);
            est.m(i).eOmegaAr_l  = ebar_l(sys.m(i).iOmegaAr);
            est.m(i).eOmegaAr_u  = ebar_u(sys.m(i).iOmegaAr);
            est.m(i).pOmegaAr    = z(sys.m(i).iOmegaAr);
            est.m(i).epOmegaAr   = sigy(sys.m(i).iOmegaAr);

            % Omega_Ca
            est.m(i).OmegaCa     = q(z(sys.m(i).iOmegaCa));
            est.m(i).eOmegaCa    = ebar(sys.m(i).iOmegaCa);
            est.m(i).eOmegaCa_l  = ebar_l(sys.m(i).iOmegaCa);
            est.m(i).eOmegaCa_u  = ebar_u(sys.m(i).iOmegaCa);
            est.m(i).pOmegaCa    = z(sys.m(i).iOmegaCa);
            est.m(i).epOmegaCa   = sigy(sys.m(i).iOmegaCa);

            % pKar
            est.m(i).pKar    = z(sys.m(i).iKar);
            est.m(i).epKar   = sigy(sys.m(i).iKar);
            est.m(i).Kar     = q(z(sys.m(i).iKar));
            est.m(i).eKar    = ebar(sys.m(i).iKar);
            est.m(i).eKar_l  = ebar_l(sys.m(i).iKar);
            est.m(i).eKar_u  = ebar_u(sys.m(i).iKar);

            % pKca
            est.m(i).pKca    = z(sys.m(i).iKca);
            est.m(i).epKca   = sigy(sys.m(i).iKca);
            est.m(i).Kca     = q(z(sys.m(i).iKca));
            est.m(i).eKca    = ebar(sys.m(i).iKca);
            est.m(i).eKca_l  = ebar_l(sys.m(i).iKca);
            est.m(i).eKca_u  = ebar_u(sys.m(i).iKca);
        end
    end
end

% -------------------------------------------------------------------------

function [zpK, zgpK, ph_free] = parse_zpK(x,sys,opt) % was f2t instead of ph_all
% assigning proper calculated values to zpK and zgpK

% zpK  := equilibrium constants, aka pK's
% zgpK := first derivative of equilibrium constants, aka gpK's

    nTP = length(sys.m);
    M = sys.M;
    K = sys.K;
    nv = size(M,2);
    nrk = size(K,1);
    zpK = zeros(nrk,1);
    zgpK = zeros(nrk,nv);
    % ph_all = [];
    ph_free = [];

    for i = 1:nTP
        [pK, gpK] = calc_pK(opt, x(sys.m(i).iT), x(sys.isal), x(sys.m(i).iP) );
        iTSP = [ sys.m(i).iT, sys.isal, sys.m(i).iP];

        zpK(sys.m(i).kK0)          = pK(1);
        zgpK(sys.m(i).kK0, iTSP )  = gpK(1,:); % ∂T, ∂S, ∂P

        zpK(sys.m(i).kK1)          = pK(2);
        zgpK(sys.m(i).kK1, iTSP )  = gpK(2,:); % ∂T, ∂S, ∂P

        zpK(sys.m(i).kK2)          = pK(3);
        zgpK(sys.m(i).kK2, iTSP )  = gpK(3,:); % ∂T, ∂S, ∂P

        zpK(sys.m(i).kp2f)         = pK(14);
        zgpK(sys.m(i).kp2f, iTSP ) = gpK(14,:); % ∂T, ∂S, ∂P

        zpK(sys.m(i).kKb)          = pK(4);
        zgpK(sys.m(i).kKb, iTSP )  = gpK(4,:); % ∂T, ∂S, ∂P

        zpK(sys.m(i).kKw)          = pK(5);
        zgpK(sys.m(i).kKw, iTSP )  = gpK(5,:); % ∂T, ∂S, ∂P

        zpK(sys.m(i).kKs)          = pK(6);
        zgpK(sys.m(i).kKs, iTSP )  = gpK(6,:); % ∂T, ∂S, ∂P
        % f2t = [f2t;sys.m(i).f2t(x)];
        ph_free = [ph_free; sys.m(i).ph_free(x)];
        % ph_all = [ph_all; sys.m(i).ph_nbs(x); sys.m(i).ph_sws(x); ...
        %     sys.m(i).ph_free(x)] ;

        zpK(sys.m(i).kKf)          = pK(7);
        zgpK(sys.m(i).kKf, iTSP )  = gpK(7,:); % ∂T, ∂S, ∂P
        
        if (ismember('phosphate',opt.abr))
            zpK(sys.m(i).kK1p)         = pK(8); 
            zgpK(sys.m(i).kK1p, iTSP ) = gpK(8,:); % ∂T, ∂S, ∂P 
            
            zpK(sys.m(i).kK2p)         = pK(9); 
            zgpK(sys.m(i).kK2p,iTSP )  = gpK(9,:); % ∂T, ∂S, ∂P 
            
            zpK(sys.m(i).kK3p)         = pK(10); 
            zgpK(sys.m(i).kK3p, iTSP ) = gpK(10,:); % ∂T, ∂S, ∂P 
        end

        if (ismember('silicate',opt.abr))
            zpK(sys.m(i).kKsi)         = pK(11); 
            zgpK(sys.m(i).kKsi, iTSP ) = gpK(11,:); % ∂T, ∂S, ∂P 
        end

        if (ismember('ammonia',opt.abr))
            zpK(sys.m(i).kKnh4)         = pK(12); 
            zgpK(sys.m(i).kKnh4, iTSP ) = gpK(12,:); % ∂T, ∂S, ∂P  
        end

        if (ismember('sulfide',opt.abr))
            zpK(sys.m(i).kKh2s)         = pK(13); 
            zgpK(sys.m(i).kKh2s, iTSP ) = gpK(13,:); % ∂T, ∂S, ∂P 
        end

        if (ismember('solubility',opt.abr))
            zpK(sys.m(i).kKar)         = pK(15);
            zgpK(sys.m(i).kKar, iTSP ) = gpK(15,:); % ∂T, ∂S, ∂P 
            zpK(sys.m(i).kKca)         = pK(16);
            zgpK(sys.m(i).kKca, iTSP ) = gpK(16,:); % ∂T, ∂S, ∂P  
        end
    end
end
% ---------------------------------------------------------------------------------

function ph_all = phscales(phin,scalein,TS,Ks,TF,Kf,pfH) % pfH
    % convert input pH to all scales
    q = @(x) 10.^(-x);
    % input TS, Ks, TF, Kf
    free2tot = (1 + TS./Ks);
    sws2tot  = (1 + TS./Ks)./(1 + TS./Ks + TF./Kf);
    fH = q(pfH);
    % fH = q(phf);
    % 1 = total scale, 2 = sea water scale, 3 = free scale, 4 = NBS
    if scalein == 1
        % total scale
        factor = 0;
    elseif scalein == 2
        % seawater scale
        factor = -log(sws2tot)./log(0.1);
    elseif scalein == 3
        % free scale
        factor = -log(free2tot)./log(0.1);
    elseif scalein == 4
        % NBS
        factor = -log(sws2tot)./log(0.1) + log(fH)./log(0.1);
    elseif scalein < 1 || scalein > 4
        fprintf('Warning: Incorrect pH scale factor used.\n');
    end
    ph_tot  = phin - factor;
    ph_nbs  = ph_tot - log(sws2tot)./log(0.1) + log(fH)/log(0.1);
    ph_free = ph_tot - log(free2tot)./log(0.1);
    ph_sws  = ph_tot - log(sws2tot)./log(0.1);

    ph_all = [ph_tot, ph_sws, ph_free, ph_nbs];

end

% ----------------------------------------------------------------------------------

function z0 = init(yobs,sys,opt)
    q = sys.q;
    p = sys.p;
    
    y0  = yobs;   
    dic = q(yobs(sys.iTC));
    alk = q(yobs(sys.iTA));
    if (isnan(dic))
        dic = 2200e-6;
        y0(sys.iTC) = p(dic);
    end
    if (isnan(alk))
        alk = 2200e-6;
        y0(sys.iTA) = p(alk);
    end

    nTP = length(sys.m);
    for i = 1:nTP
        % solve for the [H+] ion concentration using only the carbonate alkalinity
        gam = dic/alk;
        K0 = q(y0(sys.m(i).iK0));
        K1 = q(y0(sys.m(i).iK1));
        K2 = q(y0(sys.m(i).iK2));
        p2f = q(y0(sys.m(i).ip2f));
        h = 0.5*( ( gam - 1 ) * K1 + ( ( 1 - gam )^2 * K1^2 - 4 * K1 * ...
            K2 * ( 1 - 2 * gam ) ).^0.5 ) ;
        hco3 =  h * alk / (h + 2 * K2 );
        co2st = h * hco3 / K1 ;
        %co3 = 0.5 * ( alk - hco3 ) ;
        co3 = dic*K1*K2/(K1*h + h*h + K1*K2) ;
        fco2 = co2st/K0;
        pco2 = fco2/p2f;

        y0(sys.m(i).iph)    = p(h);
        y0(sys.m(i).ihco3)  = p(hco3);
        y0(sys.m(i).ico2st) = p(co2st);
        y0(sys.m(i).ico3)   = p(co3);
        y0(sys.m(i).ifco2)  = p(fco2);
        y0(sys.m(i).ipco2)  = p(pco2);

        Kw = q(y0(sys.m(i).iKw));
        oh = Kw / h;
        y0(sys.m(i).ioh) = p(oh);

        Kb = q(y0(sys.m(i).iKb));
        TB = q(yobs(sys.iTB));
        %boh4 = TB / ( 1 + h / Kb );
        boh4 = TB * Kb / (Kb + h) ;
        boh3 = TB - boh4;
        y0(sys.iTB)   = p(TB);
        y0(sys.m(i).iboh3) = p(boh3);
        y0(sys.m(i).iboh4) = p(boh4);

        Ks = q(y0(sys.m(i).iKs));
        TS = q(yobs(sys.iTS));
        fH = h / ( 1 + TS / Ks );
        hso4 = TS / ( 1 + Ks / fH);
        so4  = Ks * hso4 / fH;
        y0(sys.iTS)   = p(TS);
        % y0(sys.m(ii).iphfree)   = p(hfree);
        y0(sys.m(i).ipfH)   = p(fH);
        y0(sys.m(i).ihso4) = p(hso4);
        y0(sys.m(i).iso4)  = p(so4);

        Kf = q(y0(sys.m(i).iKf));
        TF = q(yobs(sys.iTF));
        HF = TF / ( 1 + Kf / h );
        F  = Kf * HF / h;
        y0(sys.iTF) = p(TF);
        y0(sys.m(i).iF)  = p(F);
        y0(sys.m(i).iHF) = p(HF);

        free2tot = (1 + TS./Ks);
        % sws2tot  = (1 + TS./Ks)./(1 + TS./Ks + TF./Kf);
        % fH = q(y0(sys.m(i).ipfH));

        ph_tot = p(h);
        % ph_nbs  = ph_tot - log(sws2tot)./log(0.1) + log(fH)/log(0.1);
        ph_free = ph_tot - log(free2tot)./log(0.1);
        % ph_sws  = ph_tot - log(sws2tot)./log(0.1);
        % y0(sys.m(i).iph_sws) = ph_sws;
        % y0(sys.m(i).iph_free) = ph_free;
        % y0(sys.m(i).iph_nbs) = ph_nbs;

        if (ismember('phosphate',opt.abr))
            K1p = q(y0(sys.m(i).iK1p));
            K2p = q(y0(sys.m(i).iK2p));
            K3p = q(y0(sys.m(i).iK3p));
            TP = q(yobs(sys.iTP));
            d = ( h^3 + K1p * h^2 + K1p * K2p * h + K1p * K2p * K3p);
            h3po4 = TP * h^3 / d;
            h2po4 = TP * K1p * h^2 / d;
            hpo4  = TP * K1p * K2p * h / d;
            po4   = TP * K1p * K2p * K3p / d;
            y0(sys.iTP)    = p(TP);
            y0(sys.m(i).ih3po4) = p(h3po4);
            y0(sys.m(i).ih2po4) = p(h2po4);
            y0(sys.m(i).ihpo4)  = p(hpo4);
            y0(sys.m(i).ipo4)   = p(po4);
        end

        if (ismember('silicate',opt.abr))
            Ksi = q(y0(sys.m(i).iKsi));
            TSi = q(yobs(sys.iTSi));
            siooh3 = TSi / ( 1 + h / Ksi );
            sioh4  = TSi - siooh3;
            y0(sys.iTSi)    = p(TSi);
            y0(sys.m(i).isiooh3) = p(siooh3);
            y0(sys.m(i).isioh4)  = p(sioh4);
        end
        if (ismember('ammonia',opt.abr))
            Knh4 = q(y0(sys.m(i).iKnh4));
            TNH3 = q(yobs(sys.iTNH3));
            nh3 = TNH3 / ( 1 + h / Knh4 );
            nh4 = TNH3 - nh3 ;
            y0(sys.iTNH3)      = p(TNH3);
            y0(sys.m(i).inh3)  = p(nh3);
            y0(sys.m(i).inh4)  = p(nh4);
        end
        if (ismember('sulfide',opt.abr))
            Kh2s = q(y0(sys.m(i).iKh2s));
            TH2S = q(yobs(sys.iTH2S));
            hs = TH2S / ( 1 + h / Kh2s );
            h2s = TH2S - hs ;
            y0(sys.iTH2S)     = p(TH2S);
            y0(sys.m(i).ihs)  = p(hs);
            y0(sys.m(i).ih2s) = p(h2s);
        end
        if (ismember('solubility',opt.abr))
            Kar = q(y0(sys.m(i).iKar));
            TCal = q(yobs(sys.iTCal));
            OmegaAr = co3 * TCal / Kar;
            Kca = q(y0(sys.m(i).iKca));
            OmegaCa = co3 * TCal / Kca ;
            y0(sys.iTCal)         = p(TCal);
            y0(sys.m(i).ica)      = p(TCal);
            y0(sys.m(i).iOmegaAr) = p(OmegaAr);
            y0(sys.m(i).iOmegaCa) = p(OmegaCa);
        end
    end
    % add the Lagrange multipliers
    nlam = size(sys.M,1) + size(sys.K,1) + nTP; % (old)
    % ^ + nTP was for each f2t (x3), now we also have ph_nbs, ph_free,
    % and ph_sws with f2t, so nTP * 4 things
    % nlam = size(sys.M,1) + size(sys.K,1) + (nTP*3);
    lam = zeros(nlam,1);
    z0 = [y0(:);lam(:)];
  
keyboard
    % q = sys.q;
    % p = sys.p;
    % 
    % nD = size(yobs,1);
    % y0  = yobs;
    % 
    % nTP = length(sys.m);
    % nz0 = length(yobs) + size(sys.M,1) + size(sys.K,1) + (nTP*4);
    % 
    % z0 = zeros(nTP,nz0); %initialize it
    % 
    % for i = 1:nD
    % 
    %     dic = q(yobs(i,sys.iTC));
    %     alk = q(yobs(i,sys.iTA));
    %     if (isnan(dic))
    %         dic = 2200e-6;
    %         y0(sys.iTC) = p(dic);
    %     end
    %     if (isnan(alk))
    %         alk = 2200e-6;
    %         y0(sys.iTA) = p(alk);
    %     end
    % 
    %     for ii = 1:nTP
    %         % solve for the [H+] ion concentration using only the carbonate alkalinity
    %         gam = dic/alk;
    %         K0 = q(y0(i,sys.m(ii).iK0));
    %         K1 = q(y0(i,sys.m(ii).iK1));
    %         K2 = q(y0(i,sys.m(ii).iK2));
    %         p2f = q(y0(i,sys.m(ii).ip2f));
    %         h = 0.5*( ( gam - 1 ) * K1 + ( ( 1 - gam )^2 * K1^2 - 4 * K1 * K2 * ( 1 - 2 * gam ) ).^0.5 ) ;
    %         hco3 =  h * alk / (h + 2 * K2 );
    %         co2st = h * hco3 / K1 ;
    %         %co3 = 0.5 * ( alk - hco3 ) ;
    %         co3 = dic*K1*K2/(K1*h + h*h + K1*K2) ;
    %         fco2 = co2st/K0;
    %         pco2 = fco2/p2f;
    % 
    %         y0(i,sys.m(ii).iph)    = p(h);
    %         y0(i,sys.m(ii).ihco3)  = p(hco3);
    %         y0(i,sys.m(ii).ico2st) = p(co2st);
    %         y0(i,sys.m(ii).ico3)   = p(co3);
    %         y0(i,sys.m(ii).ifco2)  = p(fco2);
    %         y0(i,sys.m(ii).ipco2)  = p(pco2);
    % 
    %         Kw = q(y0(i,sys.m(ii).iKw));
    %         oh = Kw / h;
    %         y0(i,sys.m(ii).ioh) = p(oh);
    % 
    %         Kb = q(y0(i,sys.m(ii).iKb));
    %         TB = q(yobs(i,sys.iTB));
    %         %boh4 = TB / ( 1 + h / Kb );
    %         boh4 = TB * Kb / (Kb + h) ;
    %         boh3 = TB - boh4;
    %         y0(i,sys.iTB)   = p(TB);
    %         y0(i,sys.m(ii).iboh3) = p(boh3);
    %         y0(i,sys.m(ii).iboh4) = p(boh4);
    % 
    %         Ks = q(y0(i,sys.m(ii).iKs));
    %         TS = q(yobs(i,sys.iTS));
    %         hfree = h / ( 1 + TS / Ks );
    %         hso4 = TS / ( 1 + Ks / hfree);
    %         so4  = Ks * hso4 / hfree;
    %         y0(i,sys.iTS)   = p(TS);
    %         y0(i,sys.m(ii).iphfree)   = p(hfree);
    %         y0(i,sys.m(ii).ihso4) = p(hso4);
    %         y0(i,sys.m(ii).iso4)  = p(so4);
    % 
    %         Kf = q(y0(i,sys.m(ii).iKf));
    %         TF = q(yobs(i,sys.iTF));
    %         HF = TF / ( 1 + Kf / h );
    %         F  = Kf * HF / h;
    %         y0(i,sys.iTF) = p(TF);
    %         y0(i,sys.m(ii).iF)  = p(F);
    %         y0(i,sys.m(ii).iHF) = p(HF);
    % 
    %         free2tot = (1 + TS./Ks);
    %         sws2tot  = (1 + TS./Ks)./(1 + TS./Ks + TF./Kf);
    %         fH = q(y0(i,sys.m(ii).ipfH));
    % 
    %         ph_tot = p(h);
    %         ph_nbs  = ph_tot - log(sws2tot)./log(0.1) + log(fH)/log(0.1);
    %         ph_free = ph_tot - log(free2tot)./log(0.1);
    %         ph_sws  = ph_tot - log(sws2tot)./log(0.1);
    %         y0(i,sys.m(ii).iph_sws) = ph_sws;
    %         y0(i,sys.m(ii).iph_free) = ph_free;
    %         y0(i,sys.m(ii).iph_nbs) = ph_nbs;
    % 
    %         if (ismember('phosphate',opt.abr))
    %             K1p = q(y0(i,sys.m(ii).iK1p));
    %             K2p = q(y0(i,sys.m(ii).iK2p));
    %             K3p = q(y0(i,sys.m(ii).iK3p));
    %             TP = q(yobs(i,sys.iTP));
    %             d = ( h^3 + K1p * h^2 + K1p * K2p * h + K1p * K2p * K3p);
    %             h3po4 = TP * h^3 / d;
    %             h2po4 = TP * K1p * h^2 / d;
    %             hpo4  = TP * K1p * K2p * h / d;
    %             po4   = TP * K1p * K2p * K3p / d;
    %             y0(i,sys.iTP)    = p(TP);
    %             y0(i,sys.m(ii).ih3po4) = p(h3po4);
    %             y0(i,sys.m(ii).ih2po4) = p(h2po4);
    %             y0(i,sys.m(ii).ihpo4)  = p(hpo4);
    %             y0(i,sys.m(ii).ipo4)   = p(po4);
    %         end
    % 
    %         if (ismember('silicate',opt.abr))
    %             Ksi = q(y0(i,sys.m(ii).iKsi));
    %             TSi = q(yobs(i,sys.iTSi));
    %             siooh3 = TSi / ( 1 + h / Ksi );
    %             sioh4  = TSi - siooh3;
    %             y0(i,sys.iTSi)    = p(TSi);
    %             y0(i,sys.m(ii).isiooh3) = p(siooh3);
    %             y0(i,sys.m(ii).isioh4)  = p(sioh4);
    %         end
    %         if (ismember('ammonia',opt.abr))
    %             %error('Need to implement initialization for ammonia');
    %             Knh4 = q(y0(i,sys.m(ii).iKnh4));
    %             TNH3 = q(yobs(i,sys.iTNH3));
    %             nh3 = TNH3 / ( 1 + h / Knh4 );
    %             nh4 = TNH3 - nh3 ;
    %             y0(i,sys.iTNH3)      = p(TNH3);
    %             y0(i,sys.m(ii).inh3)  = p(nh3);
    %             y0(i,sys.m(ii).inh4)  = p(nh4);
    %         end
    %         if (ismember('sulfide',opt.abr))
    %             %error('Need to implement initialization for sulfide');
    %             Kh2s = q(y0(i,sys.m(ii).iKh2s));
    %             TH2S = q(yobs(i,sys.iTH2S));
    %             hs = TH2S / ( 1 + h / Kh2s );
    %             h2s = TH2S - hs ;
    %             y0(i,sys.iTH2S)     = p(TH2S);
    %             y0(i,sys.m(ii).ihs)  = p(hs);
    %             y0(i,sys.m(ii).ih2s) = p(h2s);
    %         end
    %         if (ismember('solubility',opt.abr))
    %             Kar = q(y0(i,sys.m(ii).iKar));
    %             TCal = q(yobs(i,sys.iTCal));
    %             OmegaAr = co3 * TCal / Kar;
    %             Kca = q(y0(i,sys.m(ii).iKca));
    %             OmegaCa = co3 * TCal / Kca ;
    %             y0(i,sys.iTCal)         = p(TCal);
    %             y0(i,sys.m(ii).ica)      = p(TCal);
    %             y0(i,sys.m(ii).iOmegaAr) = p(OmegaAr);
    %             y0(i,sys.m(ii).iOmegaCa) = p(OmegaCa);
    %         end
    %     end
    %     % add the Lagrange multipliers
    %     % nlam = size(sys.M,1) + size(sys.K,1) + nTP; % (old)
    %     % ^ + nTP was for each f2t (x3), now we also have ph_nbs, ph_free,
    %     % and ph_sws with f2t, so nTP * 4 things
    %     nlam = size(sys.M,1) + size(sys.K,1) + (nTP*4);
    %     lam = zeros(nlam,1);
    %     z0(i,:) = [y0(i,:)';lam(:)];
    % end
end






