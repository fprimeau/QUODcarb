function [est,obs,sys,iflag] = QUODcarb(obs,opt)
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


    opt = check_opt(opt);                       % check opt structure

    sys = mksys(obs(1),opt.phscale);  % create indexing
    nD = length(obs);  nv = size(sys.K,2);

    % populate obs, yobs, wobs at each datapoint
    [obs,yobs,wobs] = parse_input(obs,sys,opt,nD);
    for i = 1:nD

        z0 = init(yobs(i,:),sys,opt);           % initialize
        tol = 1e-9;
        fun = @(z) limp(z,yobs(i,:),wobs(i,:),obs(i),sys,opt);

        [zhat,J,iflag(i)] = newtn(z0,fun,tol);

        if (iflag(i) ~=0) && (opt.printmes ~= 0)
            fprintf('Newton''s method iflag = %i at i = %i \n',iflag(i),i);
        end
       
        % calculate posterior uncertainty
        warning('off');
        C = inv(J);
        warning('on');
        C = C(1:nv,1:nv);
        %y = zhat(1:nv);
        sigx = sqrt(full(diag(C)));
        if (sum(isnan(sigx)) > 0) && (opt.printmes ~= 0)
            fprintf('NaN found in output means faulty run. i = %i\n',i)
        end
        [est(i)] = parse_output(zhat,sigx,opt,sys);    % populate est


        %
        % finite difference Revelle factor computation
        %
        dpTC = 100*sqrt(eps);
        yobs_pert_plus = yobs;
        yobs_pert_plus(i,sys.ipTC) = yobs(i,sys.ipTC) + dpTC;

        yobs_pert_minus = yobs;
        yobs_pert_minus(i,sys.ipTC) = yobs(i,sys.ipTC) - dpTC;

        fun_pert_plus = @(z) limp(z,yobs_pert_plus(i,:),wobs(i,:),obs(i),sys,opt);
        fun_pert_minus = @(z) limp(z,yobs_pert_minus(i,:),wobs(i,:),obs(i),sys,opt);
        Z0 = zhat; 
        [z_pert_plus,dJ,iflag(i)] = newtn(Z0,fun_pert_plus,tol/100);
        [z_pert_minus,dJ,iflag(i)] = newtn(Z0,fun_pert_minus,tol/100);
        for j = 1:length(sys.tp(:))
            dpfco2 = z_pert_plus(sys.tp(j).ipfco2)-z_pert_minus(sys.tp(j).ipfco2);
            est(i).tp(j).Revelle = dpfco2./(2*dpTC);
        end

        %fixed = intersect(sys.fixed,union(sys.tp(1).mc,sys.tp(1).kc));
        %free = setdiff(sys.tp(1).mc,fixed);

        %ML = sys.M(sys.tp(1).mr,:);
        %KL = sys.K(sys.tp(1).kr,:);
        %MR = sys.M(sys.tp(1).mr,:);
        %KR = sys.K(sys.tp(1).kr,:);
        
        %ML = ML(:,free);
        %KL = KL(:,free);
        %MR = MR(:,fixed);
        %KR = KR(:,fixed);

        %Jac = [ML*diag(sys.dqdx(zhat(free)));...
        %       -KL];
        %keyboard
        %rhs = [MR*sys.q(zhat(fixed));...
        %       -KR*zhat(fixed)];
        %keyboard
    end
end

function [g,H,f] = limp(z,y,w,obs,sys,opt)
% [g,H,f] = limp(z,y,w,obs,sys,opt)
%
% Negative log probability for the co2 system  a.k.a. log improbability i.e., limp!
%
% INPUT:
%
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
    x       =  z(1:nv);      % measured variables

    % fill ypK, gypK, ggypK, and with associated calculated pK, gpK, and ggpK values
    % update the pK values based on the new estimate of (T,P)

    [y, gy, ggy ] = update_y(y,x,obs,sys,opt);
    % Make a vector of measured quantities    
    id   = find(~isnan(y));
    y   = y(id).';
    gy = gy(id,:);
    ggy = ggy(id,:,:);
    
    % Make a precision matrix
    W   = diag(w(id));

    % Build a matrix that Picks out the measured components of x
    I   = eye(nv); % for chain rule
    PP  = I(id,:); % picking/pick out the measured ones
    e   = PP*x - y; % calculated - measured (minus)
    ge = PP - gy;
    
    nlam    = size(M,1) + size(K,1);
    % one extra lagrange multiplier for each 
    % (T,P)-dependent ph_free
    lam     =  z(nv+1:end);  % Lagrange multipliers 

    % constraint equations
    c   = [  M * q( x ); ...
             K * x ] ;
    f   = 0.5 * e.' * W * e  + lam.' * c ;  % limp, method of lagrange multipliers    
    % -(-1/2 sum of squares) + constraint eqns, minimize f => grad(f) = 0
    dcdx = [ M * diag( sys.dqdx( x ) ); ...
             K  ];

    
    %iTSP = [[sys.tp(:).iT], sys.isal, [sys.tp(:).iP]];
    g    = [ e.'*W*ge+lam.'*dcdx,  c.' ];

    ddq     =  diag( sys.d2qdx2( x ) ); % q"
        
    [nr,nc] = size(M);
    gg      = zeros(nc);
    for row = (1:nr)
        gg  = gg + lam(row)*diag(M(row,:))*ddq;
    end
    tmp = zeros(length(x));
    eW = e.'*W;
    for jj = 1:size(ggy,1);
        tmp = tmp + eW(jj)*(squeeze(ggy(jj,:,:)));
    end
    H   = [  ge.'*W*ge-tmp+gg,  dcdx.'    ; ...
             dcdx            ,  zeros(nlam)  ];
    g = g(:); % make sure g is returned as a column vector
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

    if (~isfield(obs,'tp'))
        if opt.printmes ~= 0
            error('Need to provide temperature and pressure measurement.')
        end
    end
    
    nTP = length(obs(1).tp);
    
    for i = 1:nD % loop over all the stations
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
            yobs(i,sys.ipTC) = nan;
        else
            yobs(i,sys.ipTC) = p((obs(i).TC)*1e-6); % convt to mol/kg
        end
        if (~isfield(obs(i),'eTC')) || (~isgood(obs(i).eTC))
            obs(i).eTC      = nan;
            wobs(i,sys.ipTC) = nan;
        else
            wobs(i,sys.ipTC) = w(obs(i).TC,obs(i).eTC); % std e -> w
        end
        if(~isfield(obs(i),'TA'))  || (~isgood(obs(i).TA))
            obs(i).TA       = nan; %[]
            yobs(i,sys.ipTA) = nan;
        else
            yobs(i,sys.ipTA) = p((obs(i).TA)*1e-6); % convt to mol/kg
        end
        if (~isfield(obs(i),'eTA'))  || (~isgood(obs(i).eTA))
            obs(i).eTA      = nan;
            wobs(i,sys.ipTA) = nan;
        else
            wobs(i,sys.ipTA) = w(obs(i).TA,obs(i).eTA); % std e -> w
        end
        % total borate
        if (~isfield(obs(i),'TB')) || (~isgood(obs(i).TB))
            if opt.TB == 1
                % Uppstrom, L., Deep-Sea Research 21:161-162, 1974
                % ( copied from Orr's code )
                % TB = ( 0.000232/ 10.811) * (sal/1.80655)
                obs(i).TB       = nan;  
                yobs(i,sys.ipTB) = p(0.0004157 * obs(i).sal / 35);
            elseif opt.TB == 2
                % Lee, Kim, Myrne, Millero, Feely, Yong-Ming Liu. 2010.
                % Geochemica Et Cosmochimica Acta 74 (6): 1801-1811.
                % ( copied from Sharp's code )
                % TB = (0.0002414/ 10.811) * (sal/1.80655)
                obs(i).TB       = nan; 
                yobs(i,sys.ipTB) = p( 0.0004326 * obs(i).sal / 35 )  ; 
            elseif opt.K1K2 == 6 || opt.K1K2 == 7
                % Culkin, F., in Chemical Oceanography,
                % ed. Riley and Skirrow, 1965: GEOSECS references this
                % (copied from Orr's code)
                obs(i).TB       = nan 
                yobs(i,sys.ipTB) = p( 0.0004106 * obs(i).sal / 35 );
            end
        else
            if ((obs(i).TB) == 0)
                obs(i).TB   = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTB) = p(obs(i).TB*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs(i), 'eTB'))  || (~isgood(obs(i).eTB))
            if opt.TB == 1 % Uppstrom 1974
                           % std 5e-6 on avg 2.32e-4
                TB = (2.34e-4/10.811) * obs(i).sal/1.80655;
                TBu = ( ( (2.32e-4 + 5e-6)/10.811) * obs(i).sal/1.80655 );
                TBl = ( ( (2.32e-4 - 5e-6)/10.811) * obs(i).sal/1.80655 );
                eTB = (TBu - TBl) /2 ;
                obs(i).eTB = (eTB) ; % µmol/kg
                wobs(i,sys.ipTB) = w(TB,obs(i).eTB);
            elseif opt.TB == 2 % Lee et al 2010
                               % std 9e-6 on avg 2.414e-4
                TB = (2.414e-4/10.811) * obs(i).sal/1.80655;
                TBu = ( ( (2.414e-4 + 9e-6)/10.811) * obs(i).sal/1.80655);
                TBl = ( ( (2.414e-4 - 9e-6)/10.811) * obs(i).sal/1.80655);
                eTB = (TBu - TBl) /2;
                obs(i).eTB = eTB; % umol/kg
                wobs(i,sys.ipTB) = w(TB,obs(i).eTB);
            end
        else
            TB = (2.414e-4/10.811) * obs(i).sal/1.80655; % arbitrary choice of Lee's 2010 formula
            wobs(i,sys.ipTB)     = w(obs(i).TB,obs(i).eTB);
        end
        
        % total sulfate
        if (~isfield(obs(i), 'TS'))  || (~isgood(obs(i).TS))
            % Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
            % copied from Orr's code
            obs(i).TS = nan;
            yobs(i,sys.ipTS)     = p(( 0.14 / 96.062 ) * ( obs(i).sal / 1.80655 ) );
        else
            if ((obs(i).TS) == 0)
                obs(i).TS = 1e-3; % µmol/kg reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTS)     = p(obs(i).TS*1e-6);
        end
        if (~isfield(obs(i), 'eTS'))  || (~isgood(obs(i).eTS))
            % 0.14000 ± 0.00023
            TS = (0.14/96.062)*obs(i).sal/1.80655;
            TSu = ( ( (0.14+0.00023)/96.062 ) * obs(i).sal/ 1.80655 );
            TSl = ( ( (0.14-0.00023)/96.062 ) * obs(i).sal/ 1.80655 );
            eTS                 = (TSu - TSl) / 2;
            obs(i).eTS          = (eTS) ;
            wobs(i,sys.ipTS)     = w(TS,obs(i).eTS);
        else
            TS = (0.14/96.062)*obs(i).sal/1.80655;
            wobs(i,sys.ipTS)     = w(TS,obs(i).eTS);
        end
        % total fluoride
        if (~isfield(obs(i), 'TF'))  || (~isgood(obs(i).TF))
            % Riley, J. P., Deep-Sea Research 12:219-220, 1965:
            % this is .000068.*Sali./35. = .00000195.*Sali
            obs(i).TF   = nan; ; % convt to µmol/kg-SW
            yobs(i,sys.ipTF) = p(( 0.000067 / 18.998 ) * ...
                                 ( obs(i).sal / 1.80655 ) );
        else
            if ((obs(i).TF) == 0)
                obs(i).TF       = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTF)     = p(obs(i).TF*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs(i), 'eTF'))  || (~isgood(obs(i).eTF))
            % 6.7 ± 0.1 e-5
            TF = (6.7e-5/18.998)*obs(i).sal/1.80655;
            TFu = ( ( (6.7e-5 + 0.1e-5)/18.998) * obs(i).sal/1.80655 );
            TFl = ( ( (6.7e-5 - 0.1e-5)/18.998) * obs(i).sal/1.80655 );
            eTF                 = (TFu - TFl) / 2;
            obs(i).eTF          = (eTF);
            wobs(i,sys.ipTF)     = w(TF,obs(i).eTF);
        else
            TF = (6.7e-5/18.998)*obs(i).sal/1.80655;
            wobs(i,sys.ipTF)     = w(TF,obs(i).eTF);
        end        
        
        % total phosphate
        if (~isfield(obs(i), 'TP'))  || (~isgood(obs(i).TP))
            obs(i).TP = 1e-3; % µmol/kg
            yobs(i,sys.ipTP)     = p(obs(i).TP*1e-6); % convt µmol/kg to mol/kg
        else
            if ((obs(i).TP) == 0) % zero po4 is very unlikely and breaks the code
                obs(i).TP       = 1e-3; % umol/kg-SW, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTP)     = p(obs(i).TP*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs(i), 'eTP'))  || (~isgood(obs(i).eTP))
            obs(i).eTP          = 1e-3; % µmol/kg
            wobs(i,sys.ipTP)     = w(1e-3,obs(i).eTP);
        else
            if ((obs(i).eTP) == 0)
                obs(i).eTP      = 1e-3; % umol/kg, reset minimum if zero
            end
            wobs(i,sys.ipTP)     = w(obs(i).TP,obs(i).eTP);
        end
        
        % total silicate
        if (~isfield(obs(i), 'TSi'))  || (~isgood(obs(i).TSi))
            obs(i).TSi          = 1e-3; % µmol/kg
            yobs(i,sys.ipTSi)    = p(obs(i).TSi*1e-6); % convt µmol/kg to mol/kg
        else
            if ((obs(i).TSi) == 0) % zero silicate very unlikely and breaks code
                obs(i).TSi      = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTSi)    = p(obs(i).TSi*1e-6);
        end
        if (~isfield(obs(i), 'eTSi'))  || (~isgood(obs(i).eTSi))
            obs(i).eTSi         = 1e-3; % µmol/kg
            wobs(i,sys.ipTSi)    = w(obs(i).TSi,obs(i).eTSi);
        else
            if ((obs(i).eTSi) == 0)
                obs(i).eTSi     = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            wobs(i,sys.ipTSi)    = w(obs(i).TSi,obs(i).eTSi);
        end
        % total amonia
        if (~isfield(obs(i), 'TNH3'))  || (~isgood(obs(i).TNH3))
            obs(i).TNH3         = 1e-3; % µmol/kg
            yobs(i,sys.ipTNH3)   = p(obs(i).TNH3*1e-6); % convt µmol/kg to mol/kg
        else
            if ((obs(i).TNH3) == 0)
                obs(i).TNH3     = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTNH3)   = p(obs(i).TNH3*1e-6);
        end
        if (~isfield(obs(i), 'eTNH3'))  || (~isgood(obs(i).eTNH3))
            obs(i).eTNH3        = 5e-4; % µmol/kg
            wobs(i,sys.ipTNH3)   = w(obs(i).TNH3,obs(i).eTNH3);
        else
            wobs(i,sys.ipTNH3)   = w(obs(i).TNH3,obs(i).eTNH3);
        end
        % total sulfide
        if (~isfield(obs(i), 'TH2S'))  || (~isgood(obs(i).TH2S))
            obs(i).TH2S         = 1e-3; % µmol/kg
            yobs(i,sys.ipTH2S)   = p(obs(i).TH2S*1e-6); % convt µmol/kg to mol/kg
        else
            if ((obs(i).TH2S) == 0)
                obs(i).TH2S     = 1e-3; % umol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTH2S)   = p(obs(i).TH2S*1e-6);
        end
        if (~isfield(obs(i), 'eTH2S'))  || (~isgood(obs(i).eTH2S))
            obs(i).eTH2S        = 5e-4; % µmol/kg
            wobs(i,sys.ipTH2S)   = w(obs(i).TH2S,obs(i).eTH2S);
        else
            wobs(i,sys.ipTH2S)   = w(obs(i).TH2S,obs(i).eTH2S);
        end
            
        % total calcium
        if (~isfield(obs(i), 'TCa'))  || (~isgood(obs(i).TCa))
            if opt.K1K2 == 6 || opt.K1K2 == 7
                % Calculate Ca for GEOSECS, Riley and Skirrow 1965
                obs(i).TCa     = nan ;
                yobs(i,sys.ipTCa)   = p( 0.01026 .* obs(i).sal ./ 35 );
            else
                % Calculate Ca, Riley and Tongdui 1967
                % this is 0.010285.*obs.sal./35;
                obs(i).TCa = nan ; % stay on mol/kg
                yobs(i,sys.ipTCa)   = p( 0.02128./40.087.*(obs(i).sal./1.80655) ); % mol/kg
            end
        else
            if ((obs(i).TCa) == 0)
                obs(i).TCa     = 1e-3; % mol/kg, reset minimum to 1 nanomolar
            end
            yobs(i,sys.ipTCa)   = p(obs(i).TCa); % assume user input of mol/kg
        end
        if (~isfield(obs(i), 'eTCa'))  || (~isgood(obs(i).eTCa))
            obs(i).eTCa        = (6e-5); % mol/kg, from Riley and Tongdui 1967
            TCa = 0.02128./40.087.*(obs(i).sal./1.80655) ;
            wobs(i,sys.ipTCa)   = w(TCa,obs(i).eTCa);
        else
            TCa = 0.02128./40.087.*(obs(i).sal./1.80655); % arbitrary choice of Riley and Tongdui 67 formula
            wobs(i,sys.ipTCa)   = w(TCa,obs(i).eTCa);
        end
    
        for j = 1:nTP % loop over (T,P) pairs
            
            % create obs structure with fieldnames
            if (~isfield(obs(i).tp(j),'epK0'))
                obs(i).tp(j).epK0   = [];
            end
            if (~isfield(obs(i).tp(j),'pK0'))
                obs(i).tp(j).pK0    = nan;
            end
            if (~isfield(obs(i).tp(j),'epK1'))
                obs(i).tp(j).epK1   = [];
            end
            if (~isfield(obs(i).tp(j),'pK1'))
                obs(i).tp(j).pK1    = nan;
            end
            if (~isfield(obs(i).tp(j),'epK2'))
                obs(i).tp(j).epK2   = [];
            end
            if (~isfield(obs(i).tp(j),'pK2'))
                obs(i).tp(j).pK2    = nan;
            end
            if (~isfield(obs(i).tp(j),'efco2')) || ...
                    (~isgood(obs(i).tp(j).efco2))
                obs(i).tp(j).efco2  = [];
            end
            if (~isfield(obs(i).tp(j),'fco2')) || ...
                    (~isgood(obs(i).tp(j).fco2))
                obs(i).tp(j).fco2   = nan;
            end
            if (~isfield(obs(i).tp(j),'eco2st')) || ...
                    (~isgood(obs(i).tp(j).eco2st))
                obs(i).tp(j).eco2st = [];
            end
            if (~isfield(obs(i).tp(j),'co2st')) || ...
                    (~isgood(obs(i).tp(j).co2st))
                obs(i).tp(j).co2st  = nan;
            end
            if (~isfield(obs(i).tp(j),'ehco3')) || ...
                    (~isgood(obs(i).tp(j).ehco3))
                obs(i).tp(j).ehco3  = [];
            end
            if (~isfield(obs(i).tp(j),'hco3')) || ...
                    (~isgood(obs(i).tp(j).hco3))
                obs(i).tp(j).hco3   = nan;
            end
            if (~isfield(obs(i).tp(j),'eco3')) || ...
                    (~isgood(obs(i).tp(j).eco3))
                obs(i).tp(j).eco3   = [];
            end
            if (~isfield(obs(i).tp(j),'co3')) || (~isgood(obs(i).tp(j).co3))
                obs(i).tp(j).co3    = nan;
            end
            if (~isfield(obs(i).tp(j),'eph')) || (~isgood(obs(i).tp(j).eph))
                obs(i).tp(j).eph    = [];
            end
            if (~isfield(obs(i).tp(j),'ph')) || (~isgood(obs(i).tp(j).ph))
                obs(i).tp(j).ph     = nan;
            end

            % borate system
            if (~isfield(obs(i).tp(j),'epKb'))
                obs(i).tp(j).epKb   = [];
            end
            if (~isfield(obs(i).tp(j),'pKb'))
                obs(i).tp(j).pKb    = nan;
            end
            if (~isfield(obs(i).tp(j),'eboh4')) || (~isgood(obs(i).tp(j).eboh4))
                obs(i).tp(j).eboh4  = [];
            end
            if (~isfield(obs(i).tp(j),'boh4')) || (~isgood(obs(i).tp(j).boh4))
                obs(i).tp(j).boh4   = [];
            end
            if (~isfield(obs(i).tp(j),'eboh3')) || (~isgood(obs(i).tp(j).eboh3))
                obs(i).tp(j).eboh3  = [];
            end
            if (~isfield(obs(i).tp(j),'boh3')) || (~isgood(obs(i).tp(j).boh3))
                obs(i).tp(j).boh3   = [];
            end

            % water 
            if (~isfield(obs(i).tp(j), 'epKw'))
                obs(i).tp(j).epKw   = [];
            end
            if (~isfield(obs(i).tp(j), 'pKw'))
                obs(i).tp(j).pKw    = nan;
            end
            if (~isfield(obs(i).tp(j),'eoh')) || (~isgood(obs(i).tp(j).eoh))
                obs(i).tp(j).eoh    = [];
            end
            if (~isfield(obs(i).tp(j),'oh')) || (~isgood(obs(i).tp(j).oh))
                obs(i).tp(j).oh     = [];
            end
            
            % sulfate system
            if (~isfield(obs(i).tp(j), 'epKs'))
                obs(i).tp(j).epKs   = [];
            end
            if (~isfield(obs(i).tp(j), 'pKs'))
                obs(i).tp(j).pKs    = nan;
            end
            if (~isfield(obs(i).tp(j), 'so4')) || (~isgood(obs(i).tp(j).so4))
                obs(i).tp(j).so4    = [];
            end
            if (~isfield(obs(i).tp(j), 'eso4')) || (~isgood(obs(i).tp(j).eso4))
                obs(i).tp(j).eso4   = [];
            end
            if (~isfield(obs(i).tp(j), 'hso4')) || (~isgood(obs(i).tp(j).hso4))
                obs(i).tp(j).hso4   = [];
            end
            if (~isfield(obs(i).tp(j), 'ehso4')) || (~isgood(obs(i).tp(j).ehso4))
                obs(i).tp(j).ehso4  = [];
            end
            % fluoride system
            if (~isfield(obs(i).tp(j), 'epKf'))
                obs(i).tp(j).epKf   = [];
            end
            if (~isfield(obs(i).tp(j), 'pKf'))
                obs(i).tp(j).pKf    = nan;
            end
            if (~isfield(obs(i).tp(j), 'F')) || (~isgood(obs(i).tp(j).F))
                obs(i).tp(j).F      = [];
            end
            if (~isfield(obs(i).tp(j), 'eF')) || (~isgood(obs(i).tp(j).eF))
                obs(i).tp(j).eF     = [];
            end
            if (~isfield(obs(i).tp(j), 'HF')) || (~isgood(obs(i).tp(j).HF))
                obs(i).tp(j).HF     = [];
            end
            if (~isfield(obs(i).tp(j), 'eHF')) || (~isgood(obs(i).tp(j).eHF))
                obs(i).tp(j).eHF    = [];
            end
            
            if (~isfield(obs(i).tp(j), 'epK1p'))
                obs(i).tp(j).epK1p  = [];
            end
            if (~isfield(obs(i).tp(j), 'pK1p'))
                obs(i).tp(j).pK1p   = nan;
            end
            if (~isfield(obs(i).tp(j), 'epK2p'))
                obs(i).tp(j).epK2p  = [];
            end
            if (~isfield(obs(i).tp(j), 'pK2p'))
                obs(i).tp(j).pK2p   = nan;
            end
            if (~isfield(obs(i).tp(j), 'epK3p'))
                obs(i).tp(j).epK3p  = [];
            end
            if (~isfield(obs(i).tp(j), 'pK3p'))
                obs(i).tp(j).pK3p   = nan;
            end
            if (~isfield(obs(i).tp(j), 'h3po4')) || (~isgood(obs(i).tp(j).h3po4))
                obs(i).tp(j).h3po4  = [];
            end
            if (~isfield(obs(i).tp(j), 'eh3po4')) || (~isgood(obs(i).tp(j).eh3po4))
                obs(i).tp(j).eh3po4 = [];
            end
            if (~isfield(obs(i).tp(j), 'h2po4')) || (~isgood(obs(i).tp(j).h2po4))
                obs(i).tp(j).h2po4  = [];
            end
            if (~isfield(obs(i).tp(j), 'eh2po4')) || (~isgood(obs(i).tp(j).eh2po4))
                obs(i).tp(j).eh2po4 = [];
            end
            if (~isfield(obs(i).tp(j), 'hpo4')) || (~isgood(obs(i).tp(j).hpo4))
                obs(i).tp(j).hpo4   = [];
            end
            if (~isfield(obs(i).tp(j), 'ehpo4')) || (~isgood(obs(i).tp(j).ehpo4))
                obs(i).tp(j).ehpo4  = [];
            end
            if (~isfield(obs(i).tp(j), 'po4')) || (~isgood(obs(i).tp(j).po4))
                obs(i).tp(j).po4    = [];
            end
            if (~isfield(obs(i).tp(j), 'epo4')) || (~isgood(obs(i).tp(j).epo4))
                obs(i).tp(j).epo4   = [];
            end
            
            if (~isfield(obs(i).tp(j), 'epKsi'))
                obs(i).tp(j).epKsi  = [];
            end
            if (~isfield(obs(i).tp(j), 'pKsi'))
                obs(i).tp(j).pKsi   = nan;
            end
            if (~isfield(obs(i).tp(j), 'sioh4')) || (~isgood(obs(i).tp(j).sioh4))
                obs(i).tp(j).sioh4      = [];
            end
            if (~isfield(obs(i).tp(j), 'esioh4')) || (~isgood(obs(i).tp(j).esioh4))
                obs(i).tp(j).esioh4     = [];
            end
            if (~isfield(obs(i).tp(j), 'siooh3')) || (~isgood(obs(i).tp(j).siooh3))
                obs(i).tp(j).siooh3     = [];
            end
            if (~isfield(obs(i).tp(j), 'esiooh3')) || (~isgood(obs(i).tp(j).esiooh3))
                obs(i).tp(j).esiooh3    = [];
            end
            if (~isfield(obs(i).tp(j), 'epKnh4'))
                obs(i).tp(j).epKnh4     = [];
            end
            if (~isfield(obs(i).tp(j), 'pKnh4'))
                obs(i).tp(j).pKnh4      = nan;
            end
            if (~isfield(obs(i).tp(j), 'nh3')) || (~isgood(obs(i).tp(j).nh3))
                obs(i).tp(j).nh3    = [];
            end
            if (~isfield(obs(i).tp(j), 'enh3')) || (~isgood(obs(i).tp(j).enh3))
                obs(i).tp(j).enh3   = [];
            end
            if (~isfield(obs(i).tp(j), 'nh4')) || (~isgood(obs(i).tp(j).nh4))
                obs(i).tp(j).nh4    = [];
            end
            if (~isfield(obs(i).tp(j), 'enh4')) || (~isgood(obs(i).tp(j).enh4))
                obs(i).tp(j).enh4   = [];
            end
            if (~isfield(obs(i).tp(j), 'epKh2s'))
                obs(i).tp(j).epKh2s = [];
            end
            if (~isfield(obs(i).tp(j), 'pKh2s'))
                obs(i).tp(j).pKh2s  = nan;
            end
            if (~isfield(obs(i).tp(j), 'HS')) || (~isgood(obs(i).tp(j).HS))
                obs(i).tp(j).HS     = [];
            end
            if (~isfield(obs(i).tp(j), 'eHS')) || (~isgood(obs(i).tp(j).eHS))
                obs(i).tp(j).eHS    = [];
            end
            if (~isfield(obs(i).tp(j), 'H2S')) || (~isgood(obs(i).tp(j).H2S))
                obs(i).tp(j).H2S    = [];
            end
            if (~isfield(obs(i).tp(j), 'eH2S')) || (~isgood(obs(i).tp(j).eH2S))
                obs(i).tp(j).eH2S   = [];
            end
            if (~isfield(obs(i).tp(j),'epco2')) || ...
                    (~isgood(obs(i).tp(j).epco2))
                obs(i).tp(j).epco2  = [];
            end
            if (~isfield(obs(i).tp(j),'epp2f')) || ...
                    (~isgood(obs(i).tp(j).epp2f))
                obs(i).tp(j).epp2f  = [];
            end
            if (~isfield(obs(i).tp(j),'pp2f')) || ...
                    (~isgood(obs(i).tp(j).pp2f))
                obs(i).tp(j).pp2f   = nan;
            end
            if (~isfield(obs(i).tp(j),'pco2')) || ...
                    (~isgood(obs(i).tp(j).pco2))
                obs(i).tp(j).pco2   = nan;
            end

            if (~isfield(obs(i).tp(j), 'epKar'))
                obs(i).tp(j).epKar = [];
            end
            if (~isfield(obs(i).tp(j), 'pKar'))
                obs(i).tp(j).pKar   = nan;
            end
            if (~isfield(obs(i).tp(j), 'ca')) || (~isgood(obs(i).tp(j).ca))
                obs(i).tp(j).ca     = [];
            end
            if (~isfield(obs(i).tp(j), 'eca')) || (~isgood(obs(i).tp(j).eca))
                obs(i).tp(j).eca    = [];
            end
            if (~isfield(obs(i).tp(j), 'OmegaAr')) || (~isgood(obs(i).tp(j).OmegaAr))
                obs(i).tp(j).OmegaAr    = [];
            end
            if (~isfield(obs(i).tp(j), 'eOmegaAr')) || (~isgood(obs(i).tp(j).eOmegaAr))
                obs(i).tp(j).eOmegaAr   = [];
            end 
            if (~isfield(obs(i).tp(j), 'pKar'))
                obs(i).tp(j).pKar       = nan;
            end
            if (~isfield(obs(i).tp(j), 'epKca'))
                obs(i).tp(j).epKar      = [];
            end
            if (~isfield(obs(i).tp(j), 'pKca'))
                obs(i).tp(j).pKca       = nan;
            end
            if (~isfield(obs(i).tp(j), 'epKca'))
                obs(i).tp(j).epKca      = [];
            end
            if (~isfield(obs(i).tp(j), 'OmegaCa')) || (~isgood(obs(i).tp(j).OmegaCa))
                obs(i).tp(j).OmegaCa    = [];
            end
            if (~isfield(obs(i).tp(j), 'eOmegaCa')) || (~isgood(obs(i).tp(j).eOmegaCa))
                obs(i).tp(j).eOmegaCa   = [];
            end
            
            if (~isfield(obs(i).tp(j),'eph_free')) || ...
                    (~isgood(obs(i).tp(j).eph_free))
                obs(i).tp(j).eph_free = [];
            end
            if (~isfield(obs(i).tp(j),'ph_free')) || ...
                    (~isgood(obs(i).tp(j).ph_free))
                obs(i).tp(j).ph_free = nan;
            end
            if (~isfield(obs(i).tp(j),'epfh')) || ...
                    (~isgood(obs(i).tp(j).epfh))
                obs(i).tp(j).epfh   = [];
            end
            if (~isfield(obs(i).tp(j),'pfh')) || (~isgood(obs(i).tp(j).pfh))
                obs(i).tp(j).pfh    = nan;
            end
        end

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
            pKca  = pK(16);     pfh  = pK(17);
            
            epK0   = epK(1);    epK1  = epK(2);   epK2   = epK(3);  
            epKb   = epK(4);    epKw  = epK(5);   epKs   = epK(6);  
            epKf   = epK(7);    epK1p = epK(8);   epK2p  = epK(9);  
            epK3p  = epK(10);   epKsi = epK(11);  epKnh4 = epK(12);
            epKh2s = epK(13);   epp2f = epK(14);  epKar  = epK(15); 
            epKca  = epK(16);   epfh  = epK(17);
            
            % add "observations" for the equilibrium constants
            % and transfer from obs struct to yobs and wobs
            
            %
            % co2 solubility and fugacity
            %
            if (isgood(obs(i).tp(ii).epK0))
                wobs(i,sys.tp(ii).ipK0) = (obs(i).tp(ii).epK0)^(-2);
            else
                obs(i).tp(ii).epK0 = epK0;
                wobs(i,sys.tp(ii).ipK0) = (obs(i).tp(ii).epK0)^(-2);
            end
            if (isgood(obs(i).tp(ii).pK0))
                yobs(i,sys.tp(ii).ipK0) = obs(i).tp(ii).pK0;
            else
                yobs(i,sys.tp(ii).ipK0) = pK0;
                obs(i).tp(ii).pK0 = nan;
            end
            if (isgood(obs(i).tp(ii).fco2))
                yobs(i,sys.tp(ii).ipfco2) = p(obs(i).tp(ii).fco2*1e-6); % convt µatm to atm
            else
                yobs(i,sys.tp(ii).ipfco2) = nan;
                obs(i).tp(ii).fco2 = nan;
            end
            if (isgood(obs(i).tp(ii).efco2))
                wobs(i,sys.tp(ii).ifco2) = w(obs(i).tp(ii).fco2, obs(i).tp(ii).efco2);
            else
                wobs(i,sys.tp(ii).ipfco2) = nan;
                obs(i).tp(ii).efco2 = nan;
            end
            %
            % carbonate system
            %
            if (isgood(obs(i).tp(ii).epK1))
                wobs(i,sys.tp(ii).ipK1) = (obs(i).tp(ii).epK1)^(-2);
            else
                obs(i).tp(ii).epK1 = epK1;
                wobs(i,sys.tp(ii).ipK1) = (obs(i).tp(ii).epK1)^(-2);
            end
            if (isgood(obs(i).tp(ii).pK1))
                yobs(i,sys.tp(ii).ipK1) = obs(i).tp(ii).pK1;
            else
                yobs(i,sys.tp(ii).ipK1) = pK1;
                obs(i).tp(ii).pK1  = nan;
            end
            if (isgood(obs(i).tp(ii).epK2))
                wobs(i,sys.tp(ii).ipK2) = (obs(i).tp(ii).epK2)^(-2);
            else
                obs(i).tp(ii).epK2 = epK2;
                wobs(i,sys.tp(ii).ipK2)  = (obs(i).tp(ii).epK2)^(-2);
            end
            if (isgood(obs(i).tp(ii).pK2))
                yobs(i,sys.tp(ii).ipK2) = obs(i).tp(ii).pK2;
            else
                yobs(i,sys.tp(ii).ipK2) = pK2;
                obs(i).tp(ii).pK2 = nan;
            end
            if (isgood(obs(i).tp(ii).co2st))
                yobs(i,sys.tp(ii).ipco2st) = p(obs(i).tp(ii).co2st*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipco2st) = nan;
                obs(i).tp(ii).co2st = nan;
            end
            if (isgood(obs(i).tp(ii).eco2st))
                wobs(i,sys.tp(ii).ipco2st) = w(obs(i).tp(ii).co2st, ...
                                               obs(i).tp(ii).eco2st);
            else
                wobs(i,sys.tp(ii).ipco2st) = nan;
                obs(i).tp(ii).eco2st = nan;
            end
            if (isgood(obs(i).tp(ii).hco3))
                yobs(i,sys.tp(ii).iphco3) = p(obs(i).tp(ii).hco3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iphco3) = nan;
                obs(i).tp(ii).hco3 = nan;
            end
            if (isgood(obs(i).tp(ii).ehco3))
                wobs(i,sys.tp(ii).iphco3) = w(obs(i).tp(ii).hco3, ...
                                              obs(i).tp(ii).ehco3);
            else
                wobs(i,sys.tp(ii).iphco3) = nan;
                obs(i).tp(ii).ehco3 = nan;
            end
            if (isgood(obs(i).tp(ii).co3))
                yobs(i,sys.tp(ii).ipco3) = p(obs(i).tp(ii).co3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipco3) = nan;
                obs(i).tp(ii).co3 = nan;
            end
            if (isgood(obs(i).tp(ii).eco3))
                wobs(i,sys.tp(ii).ipco3) = w(obs(i).tp(ii).co3, ...
                                             obs(i).tp(ii).eco3);
            else
                wobs(i,sys.tp(ii).ipco3) = nan;
                obs(i).tp(ii).eco3 = nan;
            end
            if (isgood(obs(i).tp(ii).ph))
                yobs(i,sys.tp(ii).iph) = obs(i).tp(ii).ph ;
            else
                yobs(i,sys.tp(ii).iph) = nan;
                obs(i).tp(ii).ph = nan;
            end
            if (isgood(obs(i).tp(ii).eph))
                wobs(i,sys.tp(ii).iph) = (obs(i).tp(ii).eph).^(-2);
            else
                wobs(i,sys.tp(ii).iph) = nan;
                obs(i).tp(ii).eph = nan;
            end
            % borate system 
            if (isgood(obs(i).tp(ii).epKb))
                wobs(i,sys.tp(ii).ipKb) = (obs(i).tp(ii).epKb).^(-2);
            else
                obs(i).tp(ii).epKb = epKb;
                wobs(i,sys.tp(ii).ipKb)  = (obs(i).tp(ii).epKb).^(-2);  % wKb = 1/(1 + (0.01/pKsys(4)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pKb))
                yobs(i,sys.tp(ii).ipKb) = obs(i).tp(ii).pKb;
            else
                yobs(i,sys.tp(ii).ipKb) = pKb; 
                obs(i).tp(ii).pKb = nan;
            end
            if (isgood(obs(i).tp(ii).boh3))
                yobs(i,sys.tp(ii).ipboh3) = p(obs(i).tp(ii).boh3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipboh3) = nan;
                obs(i).tp(ii).boh3 = nan;
            end
            if (isgood(obs(i).tp(ii).eboh3))
                wobs(i,sys.tp(ii).ipboh3) = w(obs(i).tp(ii).boh3, obs(i).tp(ii).eboh3);
            else
                wobs(i,sys.tp(ii).ipboh3) = nan;
                obs(i).tp(ii).eboh3 = nan;
            end
            if (isgood(obs(i).tp(ii).boh4))
                yobs(i,sys.tp(ii).ipboh4) = p(obs(i).tp(ii).boh4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipboh4) = nan;
                obs(i).tp(ii).boh4 = nan;
            end
            if (isgood(obs(i).tp(ii).eboh4))
                wobs(i,sys.tp(ii).ipboh4) = w(obs(i).tp(ii).boh4, obs(i).tp(ii).eboh4);
            else
                wobs(i,sys.tp(ii).ipboh4) = nan;
                obs(i).tp(ii).eboh4 = nan;
            end
            
            % water dissociation
            if (isgood(obs(i).tp(ii).epKw))
                wobs(i,sys.tp(ii).iKw) = (obs(i).tp(ii).epKw).^(-2);
            else
                    obs(i).tp(ii).epKw = epKw;
                    wobs(i,sys.tp(ii).ipKw) = (obs(i).tp(ii).epKw).^(-2);
            end
            if (isgood(obs(i).tp(ii).pKw))
                yobs(i,sys.tp(ii).ipKw) = obs(i).tp(ii).pKw;
            else
                yobs(i,sys.tp(ii).ipKw) = pKw;
                obs(i).tp(ii).pKw = nan;
            end
            if (isgood(obs(i).tp(ii).oh))
                yobs(i,sys.tp(ii).ipoh) = p(obs(i).tp(ii).oh*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipoh) = nan;
                obs(i).tp(ii).oh = nan;
            end
            if (isgood(obs(i).tp(ii).eoh))
                wobs(i,sys.tp(ii).ipoh) = w(obs(i).tp(ii).oh, obs(i).tp(ii).eoh);
            else
                wobs(i,sys.tp(ii).ipoh) = nan;
                obs(i).tp(ii).eoh = nan;
            end
            
            % sulfate system
            if (isgood(obs(i).tp(ii).epKs))
                wobs(i,sys.tp(ii).ipKs) = (obs(i).tp(ii).epKs).^(-2);
            else
                obs(i).tp(ii).epKs = epKs;
                wobs(i,sys.tp(ii).ipKs) = (obs(i).tp(ii).epKs).^(-2); % wKs = 1/(1 + (0.0021/pKsys(6)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pKs))
                yobs(i,sys.tp(ii).ipKs) = obs(i).tp(ii).pKs;
            else
                yobs(i,sys.tp(ii).ipKs) = pKs;
                obs(i).tp(ii).pKs = nan;
            end
            if (isgood(obs(i).tp(ii).hso4))
                yobs(i,sys.tp(ii).iphso4) = p(obs(i).tp(ii).hso4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iphso4) = nan;
                obs(i).tp(ii).hso4 = nan;
            end
            if (isgood(obs(i).tp(ii).ehso4))
                wobs(i,sys.tp(ii).iphso4) = w(obs(i).tp(ii).hso4, obs(i).tp(ii).ehso4);
            else
                wobs(i,sys.tp(ii).iphso4) = nan;
                obs(i).tp(ii).ehso4 = nan;
            end
            if (isgood(obs(i).tp(ii).so4))
                yobs(i,sys.tp(ii).ipso4) = p(obs(i).tp(ii).so4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipso4) = nan;
                obs(i).tp(ii).so4 = nan;
            end
            if (isgood(obs(i).tp(ii).eso4))
                wobs(i,sys.tp(ii).ipso4) = w(obs(i).tp(ii).so4, obs(i).tp(ii).eso4);
            else
                wobs(i,sys.tp(ii).ipso4) = nan;
                obs(i).tp(ii).eso4 = nan;
            end
            % fluoride system
            if (isgood(obs(i).tp(ii).epKf))
                wobs(i,sys.tp(ii).ipKf) = (obs(i).tp(ii).epKf).^(-2);
            else
                obs(i).tp(ii).epKf = epKf;
                wobs(i,sys.tp(ii).ipKf) = (obs(i).tp(ii).epKf).^(-2);   % wKF = 1/(p(1 + 0.02/KF))^2 ; % 2% relative error
            end
            if (isgood(obs(i).tp(ii).pKf))
                yobs(i,sys.tp(ii).ipKf) = obs(i).tp(ii).pKf;
            else
                yobs(i,sys.tp(ii).ipKf) = pKf;
                obs(i).tp(ii).pKf = nan;
            end
            if (isgood(obs(i).tp(ii).F))
                yobs(i,sys.tp(ii).ipF) = p(obs(i).tp(ii).F*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipF) = nan;
                obs(i).tp(ii).F = nan;
            end
            if (isgood(obs(i).tp(ii).eF))
                wobs(i,sys.tp(ii).ipF) = w(obs(i).tp(ii).F, obs(i).tp(ii).eF);
            else
                wobs(i,sys.tp(ii).ipF) = nan;
                obs(i).tp(ii).eF = nan;
            end
            if (isgood(obs(i).tp(ii).HF))
                yobs(i,sys.tp(ii).ipHF) = p(obs(i).tp(ii).HF*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipHF) = nan;
                obs(i).tp(ii).HF = nan;
            end
            if (isgood(obs(i).tp(ii).eHF))
                wobs(i,sys.tp(ii).ipHF) = w(obs(i).tp(ii).HF, obs(i).tp(ii).eHF);
            else
                wobs(i,sys.tp(ii).ipHF) = nan;
                obs(i).tp(ii).eHF = nan;
            end
            
            % phosphate system
            if (isgood(obs(i).tp(ii).epK1p))
                wobs(i,sys.tp(ii).ipK1p) = (obs(i).tp(ii).epK1p).^(-2);
            else
                obs(i).tp(ii).epK1p = epK1p;
                wobs(i,sys.tp(ii).ipK1p) = (obs(i).tp(ii).epK1p).^(-2);  % wK1p = 1/(1 + (0.09/pKsys(8)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pK1p))
                yobs(i,sys.tp(ii).ipK1p) = obs(i).tp(ii).pK1p;
            else
                yobs(i,sys.tp(ii).ipK1p) = pK1p;
                obs(i).tp(ii).pK1p = nan;
            end
            if (isgood(obs(i).tp(ii).epK2p))
                wobs(i,sys.tp(ii).ipK2p) = (obs(i).tp(ii).epK2p).^(-2);
            else
                obs(i).tp(ii).epK2p = epK2p;
                wobs(i,sys.tp(ii).ipK2p) = (obs(i).tp(ii).epK2p).^(-2);  % wK2p = 1/(1 + (0.03/pKsys(9)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pK2p))
                yobs(i,sys.tp(ii).ipK2p) = obs(i).tp(ii).pK2p;
            else
                yobs(i,sys.tp(ii).ipK2p) = pK2p;
                obs(i).tp(ii).pK2p = nan;
            end
            if (isgood(obs(i).tp(ii).epK3p))
                wobs(i,sys.tp(ii).ipK3p) = (obs(i).tp(ii).epK3p).^(-2);
            else
                obs(i).tp(ii).epK3p = epK3p;
                wobs(i,sys.tp(ii).ipK3p) = (obs(i).tp(ii).epK3p).^(-2);  % wK3p = 1/(1 + (0.02/pKsys(10)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pK3p))
                yobs(i,sys.tp(ii).ipK3p) = obs(i).tp(ii).pK3p;
            else
                yobs(i,sys.tp(ii).ipK3p) = pK3p;
                obs(i).tp(ii).pK3p = nan;
            end
            if (isgood(obs(i).tp(ii).h3po4))
                yobs(i,sys.tp(ii).iph3po4) = p(obs(i).tp(ii).h3po4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iph3po4) = nan;
                obs(i).tp(ii).h3po4 = nan;
            end
            if (isgood(obs(i).tp(ii).eh3po4))
                wobs(i,sys.tp(ii).iph3po4) = w(obs(i).tp(ii).h3po4, obs(i).tp(ii).eh3po4);
            else
                wobs(i,sys.tp(ii).iph3po4) = nan;
                obs(i).tp(ii).eh3po4 = nan;
            end
            
            if (isgood(obs(i).tp(ii).h2po4))
                yobs(i,sys.tp(ii).iph2po4) = p(obs(i).tp(ii).h2po4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iph2po4) = nan;
                obs(i).tp(ii).h2po4 = nan;
            end
            if (isgood(obs(i).tp(ii).eh2po4))
                wobs(i,sys.tp(ii).iph2po4) = w(obs(i).tp(ii).h2po4, obs(i).tp(ii).eh2po4);
            else
                wobs(i,sys.tp(ii).iph2po4) = nan;
                obs(i).tp(ii).eh2po4 = nan;
            end
            if (isgood(obs(i).tp(ii).hpo4))
                yobs(i,sys.tp(ii).iphpo4) = p(obs(i).tp(ii).hpo4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).iphpo4) = nan;
                obs(i).tp(ii).hpo4 = nan;
            end
            if (isgood(obs(i).tp(ii).ehpo4))
                wobs(i,sys.tp(ii).iphpo4) = w(obs(i).tp(ii).hpo4, obs(i).tp(ii).ehpo4);
            else
                wobs(i,sys.tp(ii).iphpo4) = nan;
                obs(i).tp(ii).ehpo4 = nan;
            end
            if (isgood(obs(i).tp(ii).po4))
                yobs(i,sys.tp(ii).ippo4) = p(obs(i).tp(ii).po4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ippo4) = nan;
                obs(i).tp(ii).po4 = nan;
            end
            if (isgood(obs(i).tp(ii).epo4))
                wobs(i,sys.tp(ii).ippo4) = w(obs(i).tp(ii).po4, obs(i).tp(ii).epo4);
            else
                wobs(i,sys.tp(ii).ippo4) = nan;
                obs(i).tp(ii).epo4 = nan;
            end
            % silicate system
            if (isgood(obs(i).tp(ii).epKsi))
                wobs(i,sys.tp(ii).ipKsi) = (obs(i).tp(ii).epKsi).^(-2);
            else
                obs(i).tp(ii).epKsi = epKsi;
                wobs(i,sys.tp(ii).ipKsi) = (obs(i).tp(ii).epKsi).^(-2);  % wKSi = 1/(1 + (0.02/pKsys(11)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pKsi))
                yobs(i,sys.tp(ii).ipKsi) = obs(i).tp(ii).pKsi;
            else
                yobs(i,sys.tp(ii).ipKsi) = pKsi;
                obs(i).tp(ii).pKsi = nan;
            end
            if (isgood(obs(i).tp(ii).sioh4))
                yobs(i,sys.tp(ii).ipsioh4) = p(obs(i).tp(ii).sioh4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipsioh4) = nan;
                obs(i).tp(ii).sioh4 = nan;
            end
            if (isgood(obs(i).tp(ii).esioh4))
                wobs(i,sys.tp(ii).ipsioh4) = w(obs(i).tp(ii).sioh4, obs(i).tp(ii).esioh4);
            else
                wobs(i,sys.tp(ii).ipsioh4) = nan;
                obs(i).tp(ii).esioh4 = nan;
            end
            if (isgood(obs(i).tp(ii).siooh3))
                yobs(i,sys.tp(ii).ipsiooh3) = p(obs(i).tp(ii).siooh3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipsiooh3) = nan;
                obs(i).tp(ii).siooh3 = nan;
            end
            if (isgood(obs(i).tp(ii).esiooh3))
                wobs(i,sys.tp(ii).ipsiooh3) = w(obs(i).tp(ii).siooh3, obs(i).tp(ii).esiooh3);
            else
                wobs(i,sys.tp(ii).ipsiooh3) = nan;
                obs(i).tp(ii).esiooh3 = nan;
            end
            % amonia system
            if (isgood(obs(i).tp(ii).epKnh4))
                wobs(i,sys.tp(ii).ipKnh4) = (obs(i).tp(ii).epKnh4).^(-2);
            else
                obs(i).tp(ii).epKnh4 = epKnh4;
                wobs(i,sys.tp(ii).ipKnh4) = (obs(i).tp(ii).epKnh4).^(-2);  % wKnh4 = 1/(1 + (0.00017/pKsys(11)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pKnh4))
                yobs(i,sys.tp(ii).ipKnh4) = obs(i).tp(ii).pKnh4;
            else
                yobs(i,sys.tp(ii).ipKnh4) = pKnh4;
                obs(i).tp(ii).pKnh4 = nan;
            end
            if (isgood(obs(i).tp(ii).nh4))
                yobs(i,sys.tp(ii).ipnh4) = p(obs(i).tp(ii).nh4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipnh4) = nan;
                obs(i).tp(ii).nh4 = nan;
            end
            if (isgood(obs(i).tp(ii).enh4))
                wobs(i,sys.tp(ii).ipnh4) = w(obs(i).tp(ii).nh4, obs(i).tp(ii).enh4);
            else
                wobs(i,sys.tp(ii).ipnh4) = nan;
                obs(i).tp(ii).enh4 = nan;
            end
            if (isgood(obs(i).tp(ii).nh3))
                yobs(i,sys.tp(ii).ipnh3) = p(obs(i).tp(ii).nh3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipnh3) = nan;
                obs(i).tp(ii).nh3 = nan;
            end
            if (isgood(obs(i).tp(ii).enh3))
                wobs(i,sys.tp(ii).ipnh3) = w(obs(i).tp(ii).nh3, obs(i).tp(ii).enh3);
            else
                wobs(i,sys.tp(ii).ipnh3) = nan;
                obs(i).tp(ii).enh3 = nan;
            end
            % sulfide system
            if (isgood(obs(i).tp(ii).epKh2s))
                wobs(i,sys.tp(ii).ipKh2s) = (obs(i).tp(ii).epKh2s).^(-2);
            else
                obs(i).tp(ii).epKh2s = epKh2s;
                wobs(i,sys.tp(ii).ipKh2s) = (obs(i).tp(ii).epKh2s).^(-2);  % wKh2s = 1/(1 + (0.033/pKsys(11)))^2 ;
            end
            if (isgood(obs(i).tp(ii).pKh2s))
                yobs(i,sys.tp(ii).ipKh2s) = obs(i).tp(ii).pKh2s;
            else
                yobs(i,sys.tp(ii).ipKh2s) = pKh2s;
                obs(i).tp(ii).pKh2s = nan;
            end
            if (isgood(obs(i).tp(ii).H2S))
                yobs(i,sys.tp(ii).ipH2S) = p(obs(i).tp(ii).H2S*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipH2S) = nan;
                obs(i).tp(ii).H2S = nan;
            end
            if (isgood(obs(i).tp(ii).eH2S))
                wobs(i,sys.tp(ii).ipH2S) = w(obs(i).tp(ii).H2S, ...
                                             obs(i).tp(ii).eH2S);
            else
                wobs(i,sys.tp(ii).ipH2S) = nan;
                obs(i).tp(ii).eH2S = nan;
            end
            if (isgood(obs(i).tp(ii).HS))
                yobs(i,sys.tp(ii).ipHS) = p(obs(i).tp(ii).HS*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipHS) = nan;
                obs(i).tp(ii).HS = nan;
            end
            if (isgood(obs(i).tp(ii).eHS))
                wobs(i,sys.tp(ii).ipHS) = w(obs(i).tp(ii).HS, obs(i).tp(ii).eHS);
            else
                wobs(i,sys.tp(ii).ipHS) = nan;
                obs(i).tp(ii).eHS = nan;
            end
            if (isgood(obs(i).tp(ii).epp2f))
                wobs(i,sys.tp(ii).ipp2f) = (obs(i).tp(ii).epp2f)^(-2);
            else
                obs(i).tp(ii).epp2f = epp2f;
                wobs(i,sys.tp(ii).ipp2f) = (obs(i).tp(ii).epp2f)^(-2);
            end
            if (isgood(obs(i).tp(ii).pp2f))
                yobs(i,sys.tp(ii).ipp2f) = obs(i).tp(ii).pp2f;
            else
                yobs(i,sys.tp(ii).ipp2f) = pp2f;
                obs(i).tp(ii).pp2f = nan;
            end
            if (isgood(obs(i).tp(ii).pco2))
                yobs(i,sys.tp(ii).ippco2) = p(obs(i).tp(ii).pco2*1e-6); % convt µatm to atm
            else
                yobs(i,sys.tp(ii).ippco2) = nan;
                obs(i).tp(ii).pco2 = nan;
            end
            if (isgood(obs(i).tp(ii).epco2))
                wobs(i,sys.tp(ii).ipco2) = w(obs(i).tp(ii).pco2, obs(i).tp(ii).epco2);
            else
                wobs(i,sys.tp(ii).ippco2) = nan;
                obs(i).tp(ii).epco2 = nan;
            end
            
            % calcium carbonate mineral solubility
            if (isgood(obs(i).tp(ii).epKar))
                wobs(i,sys.tp(ii).ipKar) = (obs(i).tp(ii).epKar).^(-2);
            else
                obs(i).tp(ii).epKar = epKar;
                wobs(i,sys.tp(ii).ipKar) = (obs(i).tp(ii).epKar).^(-2);
            end
            if (isgood(obs(i).tp(ii).pKar))
                yobs(i,sys.tp(ii).ipKar) = obs(i).tp(ii).pKar;
            else
                yobs(i,sys.tp(ii).ipKar) = pKar;
                obs(i).tp(ii).pKar = nan;
            end
            if (isgood(obs(i).tp(ii).ca))
                yobs(i,sys.tp(ii).ipca) = p(obs(i).tp(ii).ca*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(i,sys.tp(ii).ipca) = nan;
                obs(i).tp(ii).ca = nan;
            end
            if (isgood(obs(i).tp(ii).eca))
                wobs(i,sys.tp(ii).ipca)  = w(obs(i).tp(ii).ca, obs(i).tp(ii).eca);
            else
                wobs(i,sys.tp(ii).ipca) = nan;
                obs(i).tp(ii).eca = nan;
            end
            if (isgood(obs(i).tp(ii).OmegaAr))
                yobs(i,sys.tp(ii).ipOmegaAr) = p(obs(i).tp(ii).OmegaAr); % Omega is dimensionless
            else
                yobs(i,sys.tp(ii).ipOmegaAr) = nan;
                obs(i).tp(ii).OmegaAr = nan;
            end
            if (isgood(obs(i).tp(ii).eOmegaAr))
                wobs(i,sys.tp(ii).ipOmegaAr) = w(obs(i).tp(ii).OmegaAr, obs(i).tp(ii).eOmegaAr);
            else
                wobs(i,sys.tp(ii).ipOmegaAr) = nan;
                obs(i).tp(ii).eOmegaAr = nan;
            end
            if (isgood(obs(i).tp(ii).epKca))
                wobs(i,sys.tp(ii).ipKca) = (obs(i).tp(ii).epKca).^(-2);
            else
                obs(i).tp(ii).epKca = epKca;
                wobs(i,sys.tp(ii).ipKca) = (obs(i).tp(ii).epKca).^(-2);
            end
            if (isgood(obs(i).tp(ii).pKca))
                yobs(i,sys.tp(ii).ipKca) = obs(i).tp(ii).pKca;
            else
                yobs(i,sys.tp(ii).ipKca) = pKca;
                obs(i).tp(ii).pKca = nan;
            end
            if (isgood(obs(i).tp(ii).OmegaCa))
                yobs(i,sys.tp(ii).ipOmegaCa) = p(obs(i).tp(ii).OmegaCa); % Omega is dimensionless
            else
                yobs(i,sys.tp(ii).ipOmegaCa) = nan;
                obs(i).tp(ii).OmegaCa = nan;
            end
            if (isgood(obs(i).tp(ii).eOmegaCa))
                wobs(i,sys.tp(ii).ipOmegaCa) = w(obs(i).tp(ii).OmegaCa, obs(i).tp(ii).eOmegaCa);
            else
                wobs(i,sys.tp(ii).ipOmegaCa) = nan;
                obs(i).tp(ii).eOmegaCa = nan;
            end
            
            if (isgood(obs(i).tp(ii).ph_free)) % ph_free (ph on free scale)
                yobs(i,sys.tp(ii).iph_free) = obs(i).tp(ii).ph_free ;
            else
                yobs(i,sys.tp(ii).iph_free) = nan;
                obs(i).tp(ii).ph_free = nan;
            end
            if (isgood(obs(i).tp(ii).eph_free)) 
                wobs(i,sys.tp(ii).iph_free) = obs(i).tp(ii).eph_free.^(-2);
            else
                wobs(i,sys.tp(ii).iph_free) = nan;
                obs(i).tp(ii).eph_free = nan;
            end
            if (isgood(obs(i).tp(ii).pfh)) % pfh activity coefficient
                yobs(i,sys.tp(ii).ipfh) = obs(i).tp(ii).pfh ;
            else
                yobs(i,sys.tp(ii).ipfh) = pfh;
                obs(i).tp(ii).pfh = pfh;
            end
            if (isgood(obs(i).tp(ii).epfh)) % pfh activity coefficient
                wobs(i,sys.tp(ii).ipfh) = (obs(i).tp(ii).epfh).^(-2) ;
            else
                wobs(i,sys.tp(ii).ipfh) = (epfh).^(-2);
                obs(i).tp(ii).epfh = epfh;
            end
        end
    end
end

function [y,gy,ggy] = update_y(y,x,obs,sys,opt)
    nTP         = length(sys.tp);
    M           = sys.M;
    K           = sys.K;
    nv          = size(M,2);
    gy          = zeros(nv,length(y));
    ggy         = zeros(nv,length(y),length(y));
    sal = x(sys.isal);
    y = y;
    e3  = eps^3;
    ie3 = sqrt(-1)*e3;
    for i = 1:nTP
        % use complex step method to get ∂T, ∂S, ∂P
        [pK,gpK] = calc_pK(opt,x(sys.tp(i).iT), x(sys.isal), x(sys.tp(i).iP));

        [pK_T,gpK_T] = calc_pK(opt, x(sys.tp(i).iT) + ie3  , x(sys.isal)       , x(sys.tp(i).iP)       );
        [pK_S,gpK_S] = calc_pK(opt, x(sys.tp(i).iT)        , x(sys.isal) + ie3 , x(sys.tp(i).iP)       );
        [pK_P,gpK_P] = calc_pK(opt, x(sys.tp(i).iT)        , x(sys.isal)       , x(sys.tp(i).iP) + ie3 );
        %
        ggpK = zeros(length(pK),3,3);
        ggpK(:,1,:) = imag(gpK_T)/e3;
        ggpK(:,2,:) = imag(gpK_S)/e3;
        ggpK(:,3,:) = imag(gpK_P)/e3;
        
        iTSP        = [ sys.tp(i).iT, sys.isal, sys.tp(i).iP];
        isal = sys.isal;

        % update y
        if (isnan(obs.tp(i).pK0))
            y(sys.tp(i).ipK0)             = pK(1); 
            gy(sys.tp(i).ipK0,iTSP)       = gpK(1,:);
            ggy(sys.tp(i).ipK0,iTSP,iTSP) = ggpK(1,:,:);

        end
        if (isnan(obs.tp(i).pK1))
            y(sys.tp(i).ipK1)             = pK(2);
            gy(sys.tp(i).ipK1,iTSP)       = gpK(2,:);
            ggy(sys.tp(i).ipK1,iTSP,iTSP) = ggpK(2,:,:);
        end
        if (isnan(obs.tp(i).pK2))
            y(sys.tp(i).ipK2)             = pK(3);
            gy(sys.tp(i).ipK2,iTSP)       = gpK(3,:);
            ggy(sys.tp(i).ipK2,iTSP,iTSP) = ggpK(3,:,:);
        end
        
        if (isnan(obs.tp(i).pKb))
            y(sys.tp(i).ipKb)             = pK(4);
            gy(sys.tp(i).ipKb,iTSP)       = gpK(4,:);
            ggy(sys.tp(i).ipKb,iTSP,iTSP) = ggpK(4,:,:);
        end
        if (isnan(obs.TB))
            if (opt.TB == 1)
                % Uppstrom, L., Deep-Sea Research 21:161-162, 1974
                % ( copied from Orr's code )
                % TB = ( 0.000232/ 10.811) * (sal/1.80655)
                TB = 0.0004157 * sal / 35; 
                y(sys.ipTB)             = sys.p( TB );
                gy(sys.ipTB,isal)       = sys.dpdx( TB )   *  0.0004157 / 35 ;
                ggy(sys.ipTB,isal,isal) = sys.d2pdx2( TB ) * (0.0004157 / 35 ) ^ 2;
                
            elseif (opt.TB == 2)
                % Lee, Kim, Myrne, Millero, Feely, Yong-Ming Liu. 2010.
                % Geochemica Et Cosmochimica Acta 74 (6): 1801-1811.
                % ( copied from Sharp's code )
                % TB = (0.0002414/ 10.811) * (sal/1.80655)
                TB = 0.0004326 * sal / 35;
                y(sys.ipTB)             = sys.p(TB);
                gy(sys.ipTB,isal)       =  sys.dpdx( TB )  *   0.0004326 / 35;
                ggy(sys.ipTB,isal,isal) = sys.d2pdx2( TB ) * ( 0.0004326 / 35 ) ^ 2;
            elseif (opt.TB == 3)
                % Culkin, F., in Chemical Oceanography,
                % ed. Riley and Skirrow, 1965: GEOSECS references this
                % (copied from Orr's code)
                TB = 0.0004106 * sal / 35;
                y(sys.ipTB)             = p( TB );
                gy(sys.ipTB,isal)       = sys.dpdx( TB )   *   0.0004106 / 35 ;
                ggy(sys.ipTB,isal,isal) = sys.d2pdx2( TB ) * ( 0.0004106 / 35 ) ^ 2;
            end

        end
        if (isnan(obs.tp(i).pKw))
            y(sys.tp(i).ipKw)             = pK(5);
            gy(sys.tp(i).ipKw,iTSP)       = gpK(5,:);
            ggy(sys.tp(i).ipKw,iTSP,iTSP) = ggpK(5,:,:);
        end
        if (isnan(obs.tp(i).pKs))
            y(sys.tp(i).ipKs)             = pK(6);
            gy(sys.tp(i).ipKs,iTSP)       = gpK(6,:);
            ggy(sys.tp(i).ipKs,iTSP,iTSP) = ggpK(6,:,:);
        end
        if (isnan(obs.TS))
            % Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
            % copied from Orr's code
            TS     = ( 0.14 / 96.062 ) * ( sal / 1.80655 );
            y(sys.ipTS)             = sys.p( TS );
            gy(sys.ipTS,isal)       = sys.dpdx( TS )   *   (0.14 / 96.062 ) / 1.80655 ;
            ggy(sys.ipTS,isal,isal) = sys.d2pdx2( TS ) * ( (0.14 / 96.062 ) / 1.80655 ) ^ 2 ;
        end
        if (isnan(obs.tp(i).pKf))
            y(sys.tp(i).ipKf)             = pK(7);
            gy(sys.tp(i).ipKf,iTSP)       = gpK(7,:);
            ggy(sys.tp(i).ipKf,iTSP,iTSP) = ggpK(7,:,:);
        end
        if (isnan(obs.TF))
            % Riley, J. P., Deep-Sea Research 12:219-220, 1965:
            % this is .000068.*Sali./35. = .00000195.*Sali      % why don't we use this expression?
            TF  = ( 0.000067 / 18.998 ) * ( sal / 1.80655 ); 
            y(sys.ipTF)              = sys.p( TF );
            gy(sys.ipTF,isal)        = sys.dpdx( TF )   *   (0.000067 / 18.998 ) / 1.80655 ;
            ggy(sys.ipTF,isal,isal)  = sys.d2pdx2( TF ) * ( (0.000067 / 18.998 ) / 1.80655 ) ^ 2 ;
        end
        if (isnan(obs.tp(i).pK1p))
            y(sys.tp(i).ipK1p)             = pK(8); 
            gy(sys.tp(i).ipK1p,iTSP)       = gpK(8,:);
            ggy(sys.tp(i).ipK1p,iTSP,iTSP) = ggpK(8,:,:);
        end
        if (isnan(obs.tp(i).pK2p))
            y(sys.tp(i).ipK2p)             = pK(9);
            gy(sys.tp(i).ipK2p,iTSP)       = gpK(9,:);
            ggy(sys.tp(i).ipK2p,iTSP,iTSP) = ggpK(9,:,:);
        end
        if (isnan(obs.tp(i).pK3p))
            y(sys.tp(i).ipK3p)             = pK(10); 
            gy(sys.tp(i).ipK3p,iTSP)       = gpK(10,:);
            ggy(sys.tp(i).ipK3p,iTSP,iTSP) = ggpK(10,:,:);
        end
        
        if (isnan(obs.tp(i).pKsi))
            y(sys.tp(i).ipKsi)             = pK(11);  
            gy(sys.tp(i).ipKsi,iTSP)       = gpK(11,:);
            ggy(sys.tp(i).ipKsi,iTSP,iTSP) = ggpK(11,:,:);
        end        
        if (isnan(obs.tp(i).pKnh4))
            y(sys.tp(i).ipKnh4)             = pK(12); 
            gy(sys.tp(i).ipKnh4,iTSP)       = gpK(12,:);
            ggy(sys.tp(i).ipKnh4,iTSP,iTSP) = ggpK(12,:,:);
        end
        
        if (isnan(obs.tp(i).pKh2s))
            y(sys.tp(i).ipKh2s)             = pK(13); 
            gy(sys.tp(i).ipKh2s,iTSP)       = gpK(13,:); 
            ggy(sys.tp(i).ipKh2s,iTSP,iTSP) = ggpK(13,:,:); 
        end
        if (isnan(obs.tp(i).pp2f))
            y(sys.tp(i).ipp2f)             = pK(14);
            gy(sys.tp(i).ipp2f,iTSP)       = gpK(14,:);
            ggy(sys.tp(i).ipp2f,iTSP,iTSP) = ggpK(14,:,:);
        end
        if (isnan(obs.tp(i).pKar))
            y(sys.tp(i).ipKar)             = pK(15);
            gy(sys.tp(i).ipKar,iTSP)       = gpK(15,:);
            ggy(sys.tp(i).ipKar,iTSP,iTSP) = ggpK(15,:,:);
            y(sys.tp(i).ipKca)             = pK(16);
            gy(sys.tp(i).ipKca,iTSP)       = gpK(16,:);
            ggy(sys.tp(i).ipKca,iTSP,iTSP) = ggpK(16,:,:);
        end
        if (isnan(obs.TCa))
            if opt.K1K2 == 6 || opt.K1K2 == 7
                % Calculate Ca for GEOSECS, Riley and Skirrow 1965
                TCa    =  ( 0.01026 * sal / 35) ;
                y(sys.ipTCa)             = sys.p( TCa );
                gy(sys.ipTCa,isal)       = sys.dpdx( TCa )   *   0.01026 / 35 ;
                ggy(sys.iTCal,isal,isal) = sys.d2pdx2( TCa ) * ( 0.01026 / 35 ) ^ 2 ;
            else
                % Calculate Ca, Riley and Tongdui 1967
                % this is 0.010285.*obs.sal./35;
                TCa = ( 0.02128 / 40.087 * ( sal / 1.80655 ) );
                y(sys.ipTCa)             = sys.p( TCa ); 
                gy(sys.ipTCa,isal)       = sys.dpdx( TCa )   *   ( 0.02128 / 40.087 ) / 1.80655 ; 
                ggy(sys.ipTCa,isal,isal) = sys.d2pdx2( TCa ) * ( ( 0.02128 / 40.087 ) / 1.80655 ) ^ 2 ; 
            end
        end
    end
end
% ---------------------------------------------------------------------------------

function z0 = init(yobs,sys,opt)
    q   = sys.q;
    p = sys.p;
    
    y0  = yobs; 
    dic = q(yobs(sys.ipTC));
    alk = q(yobs(sys.ipTA));
    if (isnan(dic))
        dic         = 2200e-6;
        y0(sys.ipTC) = p(dic);
    end
    if (isnan(alk))
        alk         = 2200e-6;
        y0(sys.ipTA) = p(alk);
    end
    
    nTP = length(sys.tp);
    for i = 1:nTP
        % solve for the [H+] using only the carbonate alkalinity
        gam     = dic/alk;
        K0      = q(y0(sys.tp(i).ipK0));
        K1      = q(y0(sys.tp(i).ipK1));
        K2      = q(y0(sys.tp(i).ipK2));
        h       = 0.5*( ( gam - 1 ) * K1 + ( ( 1 - gam )^2 * ...
                                             K1^2 - 4 * K1 * K2 * ( 1 - 2 * gam ) ).^0.5 ) ;
        hco3    =  h * alk / (h + 2 * K2 );
        co2st   = h * hco3 / K1 ;
        co3     = dic*K1*K2/(K1*h + h*h + K1*K2) ;
        fco2    = co2st/K0;
        
        y0(sys.tp(i).iph)    = p(h);
        y0(sys.tp(i).iphco3)  = p(hco3);
        y0(sys.tp(i).ipco2st) = p(co2st);
        y0(sys.tp(i).ipco3)   = p(co3);
        y0(sys.tp(i).ipfco2)  = p(fco2);
        
        Kb      = q(y0(sys.tp(i).ipKb));
        TB      = q(yobs(sys.ipTB)); 
        boh4    = TB * Kb / (Kb + h) ;
        boh3    = TB - boh4;
        y0(sys.ipTB)         = p(TB);
        y0(sys.tp(i).ipboh3) = p(boh3);
        y0(sys.tp(i).ipboh4) = p(boh4);
        
        Kw      = q(y0(sys.tp(i).ipKw));
        oh      = Kw / h;
        y0(sys.tp(i).ipoh)   = p(oh);

        
        Ks      = q(y0(sys.tp(i).ipKs));
        TS      = q(yobs(sys.ipTS));
        
        h_tot = h * ( 1 + TS / Ks );
        h_free = h;
        hso4 = h_tot - h_free;
        so4     = Ks * hso4 / h;
        y0(sys.ipTS)            = p(TS);
        y0(sys.tp(i).iphso4)    = p(hso4);
        y0(sys.tp(i).ipso4)     = p(so4);
        
        Kf      = q(y0(sys.tp(i).ipKf));
        TF      = q(yobs(sys.ipTF));
        HF      = TF / ( 1 + Kf / h_free );
        F       = Kf * HF / h_free;
        h_sws = h_tot + HF;
        y0(sys.ipTF)         = p(TF);
        y0(sys.tp(i).ipF)    = p(F);
        y0(sys.tp(i).ipHF)   = p(HF);
        
        K1p     = q(y0(sys.tp(i).ipK1p));
        K2p     = q(y0(sys.tp(i).ipK2p));
        K3p     = q(y0(sys.tp(i).ipK3p));
        TP      = q(yobs(sys.ipTP));
        d = ( h^3 + K1p * h^2 + K1p * K2p * h + K1p * K2p * K3p);
        h3po4   = TP * h^3 / d;
        h2po4   = TP * K1p * h^2 / d;
        hpo4    = TP * K1p * K2p * h / d;
        po4     = TP * K1p * K2p * K3p / d;
        y0(sys.ipTP)             = p(TP);
        y0(sys.tp(i).iph3po4)    = p(h3po4);
        y0(sys.tp(i).iph2po4)    = p(h2po4);
        y0(sys.tp(i).iphpo4)     = p(hpo4);
        y0(sys.tp(i).ippo4)      = p(po4);
        
        Ksi     = q(y0(sys.tp(i).ipKsi));
        TSi     = q(yobs(sys.ipTSi));
        siooh3  = TSi / ( 1 + h / Ksi );
        sioh4   = TSi - siooh3;
        y0(sys.ipTSi)            = p(TSi);
        y0(sys.tp(i).ipsiooh3)   = p(siooh3);
        y0(sys.tp(i).ipsioh4)    = p(sioh4);
        
        
        Knh4    = q(y0(sys.tp(i).ipKnh4));
        TNH3    = q(yobs(sys.ipTNH3));
        nh3     = TNH3 / ( 1 + h / Knh4 );
        nh4     = TNH3 - nh3 ;
        y0(sys.ipTNH3)       = p(TNH3);
        y0(sys.tp(i).ipnh3)  = p(nh3);
        y0(sys.tp(i).ipnh4)  = p(nh4);
        
        Kh2s    = q(y0(sys.tp(i).ipKh2s));
        TH2S    = q(yobs(sys.ipTH2S));
        hs      = TH2S / ( 1 + h / Kh2s );
        h2s     = TH2S - hs ;
        y0(sys.ipTH2S)       = p(TH2S);
        y0(sys.tp(i).ipHS)   = p(hs);
        y0(sys.tp(i).ipH2S)  = p(h2s);
        
        p2f     = q(y0(sys.tp(i).ipp2f));
        pco2    = fco2/p2f;
        y0(sys.tp(i).ippco2)  = p(pco2);
        
        Kar     = q(y0(sys.tp(i).ipKar));
        TCa    = q(yobs(sys.ipTCa));
        OmegaAr = co3 * TCa / Kar;
        Kca     = q(y0(sys.tp(i).ipKca));
        OmegaCa = co3 * TCa / Kca ;
        y0(sys.ipTCa)           = p(TCa);
        y0(sys.tp(i).ipca)       = p(TCa);
        y0(sys.tp(i).ipOmegaAr)  = p(OmegaAr);
        y0(sys.tp(i).ipOmegaCa)  = p(OmegaCa);
        
        
        y0(sys.tp(i).iph_tot) = p(h_tot);
        y0(sys.tp(i).iph_sws) = p(h_sws);
        y0(sys.tp(i).iph_free) = p(h_free);
        fh = q(y0(sys.tp(i).ipfh));
        h_nbs = h_free*fh;
        y0(sys.tp(i).iph_nbs) = p(h_nbs);
    end
    nlam    = size(sys.M,1) + size(sys.K,1);
    lam     = zeros(nlam,1);
    z0      = [y0(:);lam(:)];
end


function [est] = parse_output(z,sigx,opt,sys)
    % populate est, output structure with best estimates
    %
    % INPUT:
    %   z  := system state including equilibrium constants and lagrange multipliers
    
    p       = sys.p;
    q       = sys.q;

    ebar    = @(j) (0.5 * ( q( z(j) - sigx(j) ) - q( z(j) + sigx(j) ) ) );
    ebar_l  = @(j) ( q( -sigx(j) ) ); % lower sigma
    ebar_u  = @(j) ( q( sigx(j) ) ); % upper sigma
        
    % populate 'est' structure with best estimate:
    %   1. p(value) and p(error) where p(x) = -log10(x)
    %   2. value and average error about the value in 'q' 
    %           where q(x) = x^(-10)
    %   3. upper and lower bounds in 'q' space, not symmetric
    %           about the value in 'q' space

    est.sal     = z(sys.isal);
    est.esal    = sigx(sys.isal);
    
    % TC
    est.pTC     = z(sys.ipTC);               
    est.epTC    = sigx(sys.ipTC);    
    est.TC      = q(z(sys.ipTC))*1e6; % 1e6  converts mol/kg to µmol/kg
    est.eTC     = ebar(sys.ipTC)*1e6;      
    est.eTC_l   = ebar_l(sys.ipTC)*1e6;    
    est.eTC_u   = ebar_u(sys.ipTC)*1e6;
    
    % TA
    est.pTA     = z(sys.ipTA);               
    est.epTA    = sigx(sys.ipTA);
    est.TA      = q(z(sys.ipTA))*1e6; % convt 
    est.eTA     = ebar(sys.ipTA)*1e6;      
    est.eTA_l   = ebar_l(sys.ipTA)*1e6;    
    est.eTA_u   = ebar_u(sys.ipTA)*1e6;

    % TB borate
    est.pTB     = z(sys.ipTB);
    est.epTB    = sigx(sys.ipTB);
    est.TB      = q(z(sys.ipTB))*1e6; % convt mol/kg to µmol/kg
    est.eTB     = ebar(sys.ipTB)*1e6;
    est.eTB_l   = ebar_l(sys.ipTB)*1e6;
    est.eTB_u   = ebar_u(sys.ipTB)*1e6;
    
    % TS sulfate
    est.pTS     = z(sys.ipTS);
    est.epTS    = sigx(sys.ipTS);
    est.TS      = q(z(sys.ipTS))*1e6;  % convt mol/kg to µmol/kg
    est.eTS     = ebar(sys.ipTS)*1e6;
    est.eTS_l   = ebar_l(sys.ipTS)*1e6;
    est.eTS_u   = ebar_u(sys.ipTS)*1e6;
    
    % TF fluoride
    est.pTF     = z(sys.ipTF);
    est.epTF        = sigx(sys.ipTF);
    est.TF      = q(z(sys.ipTF))*1e6; % convt mol/kg to µmol/kg
    est.eTF     = ebar(sys.ipTF)*1e6;
    est.eTF_l   = ebar_l(sys.ipTF)*1e6;
    est.eTF_u   = ebar_u(sys.ipTF)*1e6;

    % TP
    est.pTP     = z(sys.ipTP);           
    est.epTP    = sigx(sys.ipTP);
    est.TP      = q(z(sys.ipTP))*1e6;   % convt mol/kg to µmol/kg   
    est.eTP     = ebar(sys.ipTP)*1e6;
    est.eTP_l   = ebar_l(sys.ipTP)*1e6; 
    est.eTP_u   = ebar_u(sys.ipTP)*1e6;

    % TSi
    est.pTSi    = z(sys.ipTSi);         
    est.epTSi   = sigx(sys.ipTSi);
    est.TSi     = q(z(sys.ipTSi))*1e6;  % convt mol/kg to µmol/kg 
    est.eTSi    = ebar(sys.ipTSi)*1e6;
    est.eTSi_l  = ebar_l(sys.ipTSi)*1e6; 
    est.eTSi_u  = ebar_u(sys.ipTSi)*1e6;
    
    % TNH3
    est.pTNH3       = z(sys.ipTNH3);       
    est.epTNH3      = sigx(sys.ipTNH3);
    est.TNH3        = q(z(sys.ipTNH3))*1e6;   % convt mol/kg to µmol/kg
    est.eTNH3       = ebar(sys.ipTNH3)*1e6;
    est.eTNH3_l     = ebar_l(sys.ipTNH3)*1e6; 
    est.eTNH3_u     = ebar_u(sys.ipTNH3)*1e6;
    
    % TH2S
    est.pTH2S       = z(sys.ipTH2S);       
    est.epTH2S      = sigx(sys.ipTH2S);
    est.TH2S        = q(z(sys.ipTH2S))*1e6; % convt mol/kg to µmol/kg
    est.eTH2S       = ebar(sys.ipTH2S)*1e6;
    est.eTH2S_l     = ebar_l(sys.ipTH2S)*1e6; 
    est.eTH2S_u     = ebar_u(sys.ipTH2S)*1e6;

    % TCa
    est.pTCa       = z(sys.ipTCa);       
    est.epTCa      = sigx(sys.ipTCa);
    est.TCa        = q(z(sys.ipTCa)); % no conversion, mol/kg
    est.eTCa       = ebar(sys.ipTCa);
    est.eTCa_l     = ebar_l(sys.ipTCa); 
    est.eTCa_u     = ebar_u(sys.ipTCa);
    
    nTP = length(sys.tp);
    for i = 1:nTP
        % temp (deg C)
        est.tp(i).T     = z(sys.tp(i).iT);
        est.tp(i).eT    = sigx(sys.tp(i).iT);
        est.tp(i).eT_l  = z(sys.tp(i).iT)-sigx(sys.tp(i).iT);
        est.tp(i).eT_u  = z(sys.tp(i).iT)+sigx(sys.tp(i).iT);
        % pressure (dbar)
        est.tp(i).P     = z(sys.tp(i).iP);        
        est.tp(i).eP    = sigx(sys.tp(i).iP);
        est.tp(i).eP_l  = z(sys.tp(i).iP)-sigx(sys.tp(i).iP);
        est.tp(i).eP_u  = z(sys.tp(i).iP)+sigx(sys.tp(i).iP);
        
        % pH on the scale used to compute the pK values
        est.tp(i).phscale   = opt.phscale;
        est.tp(i).ph   = z(sys.tp(i).iph);
        est.tp(i).eph  = sigx(sys.tp(i).iph);
        est.tp(i).h    = q(z(sys.tp(i).iph)) * 1e6;
        est.tp(i).eh   = ebar(sys.tp(i).iph) * 1e6;
        est.tp(i).eh  = ebar_l(sys.tp(i).iph) * 1e6;
        est.tp(i).eh  = ebar_u(sys.tp(i).iph) * 1e6;

        % pH_free
        est.tp(i).ph_free    = z(sys.tp(i).iph_free);
        est.tp(i).eph_free   = sigx(sys.tp(i).iph_free);
        est.tp(i).h_free    = q(z(sys.tp(i).iph_free)) * 1e6;
        est.tp(i).eh_free   = ebar(sys.tp(i).iph_free) * 1e6;
        est.tp(i).eh_freel  = ebar_l(sys.tp(i).iph_free) * 1e6;
        est.tp(i).eh_freeu  = ebar_u(sys.tp(i).iph_free) * 1e6;

        % pH_tot
        est.tp(i).ph_tot    = z(sys.tp(i).iph_tot);
        est.tp(i).eph_tot   = sigx(sys.tp(i).iph_tot);

        % pH_sws
        est.tp(i).ph_sws    = z(sys.tp(i).iph_sws);
        est.tp(i).eph_sws   = sigx(sys.tp(i).iph_sws);

        % pH_nbs
        est.tp(i).ph_nbs    = z(sys.tp(i).iph_nbs);
        est.tp(i).eph_nbs   = sigx(sys.tp(i).iph_nbs);
        
        % fH = activity coefficient
        est.tp(i).pfh       = z(sys.tp(i).ipfh);
        est.tp(i).epfh      = sigx(sys.tp(i).ipfh);
        est.tp(i).fh        = q(z(sys.tp(i).ipfh))*1e6;
        est.tp(i).efh       = ebar(sys.tp(i).ipfh)*1e6;
        est.tp(i).efh_l     = ebar_l(sys.tp(i).ipfh)*1e6;
        est.tp(i).efh_u     = ebar_u(sys.tp(i).ipfh)*1e6;

        % fCO2
        est.tp(i).fco2      = q(z(sys.tp(i).ipfco2))*1e6; % convt atm to µatm
        est.tp(i).efco2     = ebar(sys.tp(i).ipfco2)*1e6;
        est.tp(i).efco2_l   = ebar_l(sys.tp(i).ipfco2)*1e6;
        est.tp(i).efco2_u   = ebar_u(sys.tp(i).ipfco2)*1e6;
        est.tp(i).pfco2     = z(sys.tp(i).ipfco2);
        est.tp(i).epfco2    = sigx(sys.tp(i).ipfco2);

        % pCO2
        est.tp(i).pco2      = q(z(sys.tp(i).ippco2))*1e6; % convt atm to µatm
        est.tp(i).epco2     = ebar(sys.tp(i).ippco2)*1e6;
        est.tp(i).epco2_l   = ebar_l(sys.tp(i).ippco2)*1e6;
        est.tp(i).epco2_u   = ebar_u(sys.tp(i).ippco2)*1e6;
        est.tp(i).ppco2     = z(sys.tp(i).ippco2);
        est.tp(i).eppco2    = sigx(sys.tp(i).ippco2);
        est.tp(i).pp2f      = z(sys.tp(i).ipp2f);
        est.tp(i).epp2f     = sigx(sys.tp(i).ipp2f);

        % HCO3
        est.tp(i).hco3      = q(z(sys.tp(i).iphco3))*1e6; % convt mol/kg to µmol/kg
        est.tp(i).ehco3     = ebar(sys.tp(i).iphco3)*1e6;
        est.tp(i).ehco3_l   = ebar_l(sys.tp(i).iphco3)*1e6;
        est.tp(i).ehco3_u   = ebar_u(sys.tp(i).iphco3)*1e6;
        est.tp(i).phco3     = z(sys.tp(i).iphco3);
        est.tp(i).ephco3    = sigx(sys.tp(i).iphco3);

        % CO3
        est.tp(i).co3       = q(z(sys.tp(i).ipco3))*1e6;
        est.tp(i).eco3      = ebar(sys.tp(i).ipco3)*1e6;
        est.tp(i).eco3_l    = ebar_l(sys.tp(i).ipco3)*1e6;
        est.tp(i).eco3_u    = ebar_u(sys.tp(i).ipco3)*1e6;
        est.tp(i).pco3      = z(sys.tp(i).ipco3);
        est.tp(i).epco3     = sigx(sys.tp(i).ipco3);

        % CO2*
        est.tp(i).pco2st    = z(sys.tp(i).ipco2st);
        est.tp(i).epco2st   = sigx(sys.tp(i).ipco2st);
        est.tp(i).co2st     = q(z(sys.tp(i).ipco2st));
        est.tp(i).eco2st    = ebar(sys.tp(i).ipco2st);
        est.tp(i).eco2st_l  = ebar_l(sys.tp(i).ipco2st);
        est.tp(i).eco2st_u  = ebar_u(sys.tp(i).ipco2st);

        % pK0 
        est.tp(i).pK0       = z(sys.tp(i).ipK0);
        est.tp(i).epK0      = sigx(sys.tp(i).ipK0);
        est.tp(i).K0        = q(z(sys.tp(i).ipK0));
        est.tp(i).eK0       = ebar(sys.tp(i).ipK0);
        est.tp(i).eK0_l     = ebar_l(sys.tp(i).ipK0);
        est.tp(i).eK0_u     = ebar_u(sys.tp(i).ipK0);

        % pK1
        est.tp(i).pK1   = z(sys.tp(i).ipK1);
        est.tp(i).epK1  = sigx(sys.tp(i).ipK1);
        est.tp(i).K1    = q(z(sys.tp(i).ipK1));
        est.tp(i).eK1   = ebar(sys.tp(i).ipK1);
        est.tp(i).eK1_l = ebar_l(sys.tp(i).ipK1); 
        est.tp(i).eK1_u = ebar_u(sys.tp(i).ipK1);

        % pK2
        est.tp(i).pK2   = z(sys.tp(i).ipK2);
        est.tp(i).epK2  = sigx(sys.tp(i).ipK2);
        est.tp(i).K2    = q(z(sys.tp(i).ipK2));
        est.tp(i).eK2   = ebar(sys.tp(i).ipK2);
        est.tp(i).eK2_l = ebar_l(sys.tp(i).ipK2);
        est.tp(i).eK2_u = ebar_u(sys.tp(i).ipK2);
        
        % OH 
        est.tp(i).oh    = q(z(sys.tp(i).ipoh))*1e6; % convt
        est.tp(i).eoh   = ebar(sys.tp(i).ipoh)*1e6;
        est.tp(i).eoh_l = ebar_l(sys.tp(i).ipoh)*1e6;
        est.tp(i).eoh_u = ebar_u(sys.tp(i).ipoh)*1e6;
        est.tp(i).poh   = z(sys.tp(i).ipoh);
        est.tp(i).epoh  = sigx(sys.tp(i).ipoh);

        % pKw 
        est.tp(i).pKw   = z(sys.tp(i).ipKw);
        est.tp(i).epKw  = sigx(sys.tp(i).ipKw);
        est.tp(i).Kw    = q(z(sys.tp(i).ipKw));
        est.tp(i).eKw   = ebar(sys.tp(i).ipKw);
        est.tp(i).eKw_l = ebar_l(sys.tp(i).ipKw);
        est.tp(i).eKw_u = ebar_u(sys.tp(i).ipKw);

        % BOH4 borate
        est.tp(i).boh4      = q(z(sys.tp(i).ipboh4))*1e6; % convt mol/kg to µmol/kg
        est.tp(i).eboh4     = ebar(sys.tp(i).ipboh4)*1e6;
        est.tp(i).eboh4_l   = ebar_l(sys.tp(i).ipboh4)*1e6;
        est.tp(i).eboh4_u   = ebar_u(sys.tp(i).ipboh4)*1e6;
        est.tp(i).pboh4     = z(sys.tp(i).ipboh4);
        est.tp(i).epboh4    = sigx(sys.tp(i).ipboh4);

        % BOH3
        est.tp(i).boh3      = q(z(sys.tp(i).ipboh3))*1e6;
        est.tp(i).eboh3     = ebar(sys.tp(i).ipboh3)*1e6;
        est.tp(i).eboh3_l   = ebar_l(sys.tp(i).ipboh3)*1e6;
        est.tp(i).eboh3_u   = ebar_u(sys.tp(i).ipboh3)*1e6;
        est.tp(i).pboh3     = z(sys.tp(i).ipboh3);
        est.tp(i).epboh3    = sigx(sys.tp(i).ipboh3);
            
        % pKb
        est.tp(i).pKb       = z(sys.tp(i).ipKb);
        est.tp(i).epKb      = sigx(sys.tp(i).ipKb);
        est.tp(i).Kb        = q(z(sys.tp(i).ipKb));
        est.tp(i).eKb       = ebar(sys.tp(i).ipKb);
        est.tp(i).eKb_l     = ebar_l(sys.tp(i).ipKb);
        est.tp(i).eKb_u     = ebar_u(sys.tp(i).ipKb);

        % SO4 sulfate
        est.tp(i).so4       = q(z(sys.tp(i).ipso4)); % mol/kg
        est.tp(i).eso4      = ebar(sys.tp(i).ipso4);
        est.tp(i).eso4_l    = ebar_l(sys.tp(i).ipso4);
        est.tp(i).eso4_u    = ebar_u(sys.tp(i).ipso4);
        est.tp(i).pso4      = z(sys.tp(i).ipso4);
        est.tp(i).epso4     = sigx(sys.tp(i).ipso4);

        % HSO4
        est.tp(i).hso4      = q(z(sys.tp(i).iphso4))*1e6;
        est.tp(i).ehso4     = ebar(sys.tp(i).iphso4)*1e6;
        est.tp(i).ehso4_l   = ebar_l(sys.tp(i).iphso4)*1e6;
        est.tp(i).ehso4_u   = ebar_u(sys.tp(i).iphso4)*1e6;
        est.tp(i).phso4     = z(sys.tp(i).iphso4);
        est.tp(i).ephso4    = sigx(sys.tp(i).iphso4);

        % pKs
        est.tp(i).pKs       = z(sys.tp(i).ipKs);
        est.tp(i).epKs      = sigx(sys.tp(i).ipKs);
        est.tp(i).Ks        = q(z(sys.tp(i).ipKs));
        est.tp(i).eKs       = ebar(sys.tp(i).ipKs);
        est.tp(i).eKs_l     = ebar_l(sys.tp(i).ipKs);
        est.tp(i).eKs_u     = ebar_u(sys.tp(i).ipKs);

        % F fluoride 
        est.tp(i).F         = q(z(sys.tp(i).ipF))*1e6; % convt
        est.tp(i).eF        = ebar(sys.tp(i).ipF)*1e6;
        est.tp(i).eF_l      = ebar_l(sys.tp(i).ipF)*1e6;
        est.tp(i).ef_u      = ebar_u(sys.tp(i).ipF)*1e6;
        est.tp(i).pF        = z(sys.tp(i).ipF);
        est.tp(i).epF       = sigx(sys.tp(i).ipF);

        % HF 
        est.tp(i).HF        = q(z(sys.tp(i).ipHF))*1e6;
        est.tp(i).eHF       = ebar(sys.tp(i).ipHF)*1e6;
        est.tp(i).eHF_l     = ebar_l(sys.tp(i).ipHF)*1e6;
        est.tp(i).eHF_u     = ebar_u(sys.tp(i).ipHF)*1e6;
        est.tp(i).pHF       = z(sys.tp(i).ipHF);
        est.tp(i).epHF      = sigx(sys.tp(i).ipHF);

        % pKf
        est.tp(i).pKf       = z(sys.tp(i).ipKf);
        est.tp(i).epKf      = sigx(sys.tp(i).ipKf);
        est.tp(i).Kf        = q(z(sys.tp(i).ipKf));
        est.tp(i).eKf       = ebar(sys.tp(i).ipKf);
        est.tp(i).eKf_l     = ebar_l(sys.tp(i).ipKf);
        est.tp(i).eKf_u     = ebar_u(sys.tp(i).ipKf);

        % PO4
        est.tp(i).po4       = q(z(sys.tp(i).ippo4))*1e6; % convt
        est.tp(i).epo4      = ebar(sys.tp(i).ippo4)*1e6;
        est.tp(i).epo4_l    = ebar_l(sys.tp(i).ippo4)*1e6;
        est.tp(i).epo4_u    = ebar_u(sys.tp(i).ippo4)*1e6;
        est.tp(i).ppo4      = z(sys.tp(i).ippo4);
        est.tp(i).eppo4     = sigx(sys.tp(i).ippo4);

        % HPO4
        est.tp(i).hpo4      = q(z(sys.tp(i).iphpo4))*1e6;
        est.tp(i).ehpo4     = ebar(sys.tp(i).iphpo4)*1e6;
        est.tp(i).ehpo4_l   = ebar_l(sys.tp(i).iphpo4)*1e6;
        est.tp(i).ehpo4_u   = ebar_u(sys.tp(i).iphpo4)*1e6;
        est.tp(i).phpo4     = z(sys.tp(i).iphpo4);
        est.tp(i).ephpo4    = sigx(sys.tp(i).iphpo4);

        % H2PO4
        est.tp(i).h2po4     = q(z(sys.tp(i).iph2po4))*1e6;
        est.tp(i).eh2po4    = ebar(sys.tp(i).iph2po4)*1e6;
        est.tp(i).eh2po4_l  = ebar_l(sys.tp(i).iph2po4)*1e6;
        est.tp(i).eh2po4_u  = ebar_u(sys.tp(i).iph2po4)*1e6;
        est.tp(i).ph2po4    = z(sys.tp(i).iph2po4);
        est.tp(i).eph2po4   = sigx(sys.tp(i).iph2po4);            

        % H3PO4
        est.tp(i).h3po4     = q(z(sys.tp(i).iph3po4))*1e6;
        est.tp(i).eh3po4    = ebar(sys.tp(i).iph3po4)*1e6;
        est.tp(i).eh3po4_l  = ebar_l(sys.tp(i).iph3po4)*1e6;
        est.tp(i).eh3po4_u  = ebar_u(sys.tp(i).iph3po4)*1e6;
        est.tp(i).ph3po4    = z(sys.tp(i).iph3po4);
        est.tp(i).eph3po4   = sigx(sys.tp(i).iph3po4);

        % pK1p
        est.tp(i).pK1p      = z(sys.tp(i).ipK1p);
        est.tp(i).epK1p     = sigx(sys.tp(i).ipK1p);
        est.tp(i).K1p       = q(z(sys.tp(i).ipK1p));
        est.tp(i).eK1p      = ebar(sys.tp(i).ipK1p);
        est.tp(i).eK1p_l    = ebar_l(sys.tp(i).ipK1p);
        est.tp(i).eK1p_u    = ebar_u(sys.tp(i).ipK1p);

        % pK2p
        est.tp(i).pK2p      = z(sys.tp(i).ipK2p);
        est.tp(i).epK2p     = sigx(sys.tp(i).ipK2p);
        est.tp(i).K2p       = q(z(sys.tp(i).ipK2p));
        est.tp(i).eK2p      = ebar(sys.tp(i).ipK2p);
        est.tp(i).pK2p_l    = ebar_l(sys.tp(i).ipK2p);
        est.tp(i).eK2p_u    = ebar_u(sys.tp(i).ipK2p);

        % pK3p
        est.tp(i).pK3p      = z(sys.tp(i).ipK3p);
        est.tp(i).epK3p     = sigx(sys.tp(i).ipK3p);
        est.tp(i).K3p       = q(z(sys.tp(i).ipK3p));
        est.tp(i).eK3p      = ebar(sys.tp(i).ipK3p);
        est.tp(i).eK3p_l    = ebar_l(sys.tp(i).ipK3p);
        est.tp(i).eK3p_u    = ebar_u(sys.tp(i).ipK3p);
        
        % SiOH4
        est.tp(i).sioh4     = q(z(sys.tp(i).ipsioh4))*1e6; % convt
        est.tp(i).esioh4    = ebar(sys.tp(i).ipsioh4)*1e6;
        est.tp(i).esioh4_l  = ebar_l(sys.tp(i).ipsioh4)*1e6;
        est.tp(i).esioh4_u  = ebar_u(sys.tp(i).ipsioh4)*1e6;
        est.tp(i).psioh4    = z(sys.tp(i).ipsioh4);
        est.tp(i).epsioh4   = sigx(sys.tp(i).ipsioh4);

        % SiOH3
        est.tp(i).siooh3    = q(z(sys.tp(i).ipsiooh3))*1e6;
        est.tp(i).esiooh3   = ebar(sys.tp(i).ipsiooh3)*1e6;
        est.tp(i).esiooh3_l = ebar_l(sys.tp(i).ipsiooh3)*1e6;
        est.tp(i).esiooh3_u = ebar_u(sys.tp(i).ipsiooh3)*1e6;
        est.tp(i).psiooh3   = z(sys.tp(i).ipsiooh3);
        est.tp(i).epsiooh3  = sigx(sys.tp(i).ipsiooh3);

        % pKsi
        est.tp(i).pKsi      = z(sys.tp(i).ipKsi);
        est.tp(i).epKsi     = sigx(sys.tp(i).ipKsi);
        est.tp(i).Ksi       = q(z(sys.tp(i).ipKsi));
        est.tp(i).eKsi      = ebar(sys.tp(i).ipKsi);
        est.tp(i).eKsi_l    = ebar_l(sys.tp(i).ipKsi);
        est.tp(i).eKsi_u    = ebar_u(sys.tp(i).ipKsi);


        % NH3
        est.tp(i).nh3       = q(z(sys.tp(i).ipnh3))*1e6; % convt
        est.tp(i).enh3      = ebar(sys.tp(i).ipnh3)*1e6;
        est.tp(i).enh3_l    = ebar_l(sys.tp(i).ipnh3)*1e6;
        est.tp(i).enh3_u    = ebar_u(sys.tp(i).ipnh3)*1e6;
        est.tp(i).pnh3      = z(sys.tp(i).ipnh3);
        est.tp(i).epnh3     = sigx(sys.tp(i).ipnh3);

        % NH4
        est.tp(i).nh4       = q(z(sys.tp(i).ipnh4))*1e6;
        est.tp(i).enh4      = ebar(sys.tp(i).ipnh4)*1e6;
        est.tp(i).enh4_l    = ebar_l(sys.tp(i).ipnh4)*1e6;
        est.tp(i).enh4_u    = ebar_u(sys.tp(i).ipnh4)*1e6;
        est.tp(i).pnh4      = z(sys.tp(i).ipnh4);
        est.tp(i).epnh4     = sigx(sys.tp(i).ipnh4);

        % pKNH4
        est.tp(i).pKnh4     = z(sys.tp(i).ipKnh4);
        est.tp(i).epKnh4    = sigx(sys.tp(i).ipKnh4);
        est.tp(i).Knh4      = q(z(sys.tp(i).ipKnh4));
        est.tp(i).eKnh4     = ebar(sys.tp(i).ipKnh4);
        est.tp(i).eKnh4_l   = ebar_l(sys.tp(i).ipKnh4);
        est.tp(i).eKnh4_u   = ebar_u(sys.tp(i).ipKnh4);

        % HS
        est.tp(i).HS        = q(z(sys.tp(i).ipHS))*1e6; % convt
        est.tp(i).eHS       = ebar(sys.tp(i).ipHS)*1e6;
        est.tp(i).eHS_l     = ebar_l(sys.tp(i).ipHS)*1e6;
        est.tp(i).eHS_u     = ebar_u(sys.tp(i).ipHS)*1e6;
        est.tp(i).pHS       = z(sys.tp(i).ipHS);
        est.tp(i).epHS      = sigx(sys.tp(i).ipHS);

        % H2S
        est.tp(i).H2S       = q(z(sys.tp(i).ipH2S))*1e6;
        est.tp(i).eH2S      = ebar(sys.tp(i).ipH2S)*1e6;
        est.tp(i).eH2S_l    = ebar_l(sys.tp(i).ipH2S)*1e6;
        est.tp(i).eHS2_u    = ebar_u(sys.tp(i).ipH2S)*1e6;
        est.tp(i).pH2S      = z(sys.tp(i).ipH2S);
        est.tp(i).epH2S     = sigx(sys.tp(i).ipH2S);

        % pKh2s
        est.tp(i).pKh2s     = z(sys.tp(i).ipKh2s);
        est.tp(i).epKh2s    = sigx(sys.tp(i).ipKh2s);
        est.tp(i).Kh2s      = q(z(sys.tp(i).ipKh2s));
        est.tp(i).eKh2s     = ebar(sys.tp(i).ipKh2s);
        est.tp(i).eKh2s_l   = ebar_l(sys.tp(i).ipKh2s);
        est.tp(i).eKh2s_u   = ebar_u(sys.tp(i).ipKh2s);

        % Ca
        est.tp(i).ca        = q(z(sys.tp(i).ipca))*1e6;
        est.tp(i).eca       = ebar(sys.tp(i).ipca)*1e6;
        est.tp(i).eca_l     = ebar_l(sys.tp(i).ipca)*1e6;
        est.tp(i).eca_u     = ebar_u(sys.tp(i).ipca)*1e6;
        est.tp(i).pca       = z(sys.tp(i).ipca);
        est.tp(i).epca      = sigx(sys.tp(i).ipca);

        % Omega_Ar
        est.tp(i).OmegaAr    = q(z(sys.tp(i).ipOmegaAr)); % unitless
        est.tp(i).eOmegaAr   = ebar(sys.tp(i).ipOmegaAr);
        est.tp(i).eOmegaAr_l = ebar_l(sys.tp(i).ipOmegaAr);
        est.tp(i).eOmegaAr_u = ebar_u(sys.tp(i).ipOmegaAr);
        est.tp(i).pOmegaAr   = z(sys.tp(i).ipOmegaAr);
        est.tp(i).epOmegaAr  = sigx(sys.tp(i).ipOmegaAr);

        % Omega_Ca
        est.tp(i).OmegaCa    = q(z(sys.tp(i).ipOmegaCa));
        est.tp(i).eOmegaCa   = ebar(sys.tp(i).ipOmegaCa);
        est.tp(i).eOmegaCa_l = ebar_l(sys.tp(i).ipOmegaCa);
        est.tp(i).eOmegaCa_u = ebar_u(sys.tp(i).ipOmegaCa);
        est.tp(i).pOmegaCa   = z(sys.tp(i).ipOmegaCa);
        est.tp(i).epOmegaCa  = sigx(sys.tp(i).ipOmegaCa);

        % pKar
        est.tp(i).pKar      = z(sys.tp(i).ipKar);
        est.tp(i).epKar     = sigx(sys.tp(i).ipKar);
        est.tp(i).Kar       = q(z(sys.tp(i).ipKar));
        est.tp(i).eKar      = ebar(sys.tp(i).ipKar);
        est.tp(i).eKar_l    = ebar_l(sys.tp(i).ipKar);
        est.tp(i).eKar_u    = ebar_u(sys.tp(i).ipKar);

        % pKca
        est.tp(i).pKca      = z(sys.tp(i).ipKca);
        est.tp(i).epKca     = sigx(sys.tp(i).ipKca);
        est.tp(i).Kca       = q(z(sys.tp(i).ipKca));
        est.tp(i).eKca      = ebar(sys.tp(i).ipKca);
        est.tp(i).eKca_l    = ebar_l(sys.tp(i).ipKca);
        est.tp(i).eKca_u    = ebar_u(sys.tp(i).ipKca);
    end
end





