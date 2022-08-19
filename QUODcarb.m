function [y,sigy,yobs,wobs,iflag] = QUODcarb(yobs,wobs,temp,sal,pres,sys)
% [y,sigy,yobs,wobs] = QUODcarb(yobs,wobs,temp,sal,pres,sys);
%
% OUTPUT:
%    y := posterior estimates of co2-system variables and equilibrium constants
% sigy := standard-deviation of marginalized posterior probability for the elements of y
% yobs := same as input except that the pK's have been added
% wobs := same as input except that the precisions of the pK's have been added
%
% INPUT:
%   yobs  := co2-system measured quantities 
%   wobs  := precisions of co2-system measured quantities
%   temp  := temperature in deg C
%   sal   := Salinity in PSU
%   pres  := Pressure in dbar
%   sys   := struct with stuff for co2-system solver (initialized using mksys.m)

    p = sys.p;
    q = sys.q;

    nv = length(sys.variables);
    
    [pK,gpK] = local_pK(temp,sal,pres);
    pK0   = pK(1);  pK1  = pK(2);  pK2   = pK(3);  pKb   = pK(4);  
    pKw   = pK(5);  pKs  = pK(6);  pKf   = pK(7);  pK1p  = pK(8);  
    pK2p  = pK(9);  pK3p = pK(10); pKsi  = pK(11); pKnh4 = pK(12); 
    pKh2s = pK(13); pp2f = pK(14);
    
    %
    % add "observations" for the equilibrium constants 
    %
    
    if (ismember('K0',sys.variables))
        %wK0 = 1/(1 + (0.002/pKsys(1)))^2 ; % 0.002 error on pK0, taken from literature
        wobs(sys.iK0) = (0.002).^(-2);
        yobs(sys.iK0) = pK0;        
    end
    if (ismember('K1',sys.variables))
        wobs(sys.iK1) = (0.01).^(-2);  % wK1 = 1/(1 + (0.01/pKsys(2)))^2 ;
        yobs(sys.iK1) = pK1;
    end
    if (ismember('K2',sys.variables))
        wobs(sys.iK2) = (0.02).^(-2);  % wK2 = 1/(1 + (0.02/pKsys(3)))^2 ;
        yobs(sys.iK2) = pK2;
    end
    if (ismember('p2f',sys.variables))
        wobs(sys.ip2f) = (0.001).^(-2);
        yobs(sys.ip2f) = pp2f;
    end
    if (ismember('Kb',sys.variables))
        wobs(sys.iKb) = (0.01).^(-2);  % wKb = 1/(1 + (0.01/pKsys(4)))^2 ;
        yobs(sys.iKb) = pKb;
    end
    if (ismember('Kw',sys.variables))
        wobs(sys.iKw) = (0.01).^(-2);  % wKw = 1/(1 + (0.01/pKsys(5)))^2 ;
        yobs(sys.iKw) = pKw;
    end
    if (ismember('Ks',sys.variables))
        wobs(sys.iKs) = (0.0021).^(-2); % wKs = 1/(1 + (0.0021/pKsys(6)))^2 ;
        yobs(sys.iKs) = pKs;
    end
    if (ismember('KF',sys.variables))
        wobs(sys.iKF) = (0.02).^(-2);   % wKF = 1/(p(1 + 0.02/KF))^2 ; % 2% relative error
        yobs(sys.iKF) = pKf;
    end
    if (ismember('K1p',sys.variables))
        wobs(sys.iK1p) = (0.09).^(-2);  % wK1p = 1/(1 + (0.09/pKsys(8)))^2 ;
        yobs(sys.iK1p) = pK1p;
    end
    if (ismember('K2p',sys.variables))
        wobs(sys.iK2p) = (0.03).^(-2);  % wK2p = 1/(1 + (0.03/pKsys(9)))^2 ;
        yobs(sys.iK2p) = pK2p;
    end
    if (ismember('K3p',sys.variables))
        wobs(sys.iK3p) = (0.02).^(-2);  % wK3p = 1/(1 + (0.02/pKsys(10)))^2 ;
        yobs(sys.iK3p) = pK3p;
    end
    if (ismember('Ksi',sys.variables))
        wobs(sys.iKsi) = (0.02).^(-2);  % wKSi = 1/(1 + (0.02/pKsys(11)))^2 ;
        yobs(sys.iKsi) = pKsi;
    end
    % add NH3 and H2S
    if (ismember('Knh4',sys.variables))
        wobs(sys.iKnh4) = (0.00017).^(-2);  % wKnh4 = 1/(1 + (0.00017/pKsys(11)))^2 ;
        yobs(sys.iKnh4) = pKnh4;
    end
    if (ismember('Kh2s',sys.variables))
        wobs(sys.iKh2s) = (0.033).^(-2);  % wKh2s = 1/(1 + (0.033/pKsys(11)))^2 ;
        yobs(sys.iKh2s) = pKh2s;
    end
    
    gun = @(z) grad_limpco2(z,yobs,wobs,pK,gpK,sys);
    % test limpco2 gradient and hessian using complex step method
    %{
      
      [f0,g0,H0] = limpco2(z0,yobs,wobs,sys);
      I = eye(length(z0));
      for k = 1:length(z0)
      [f,g,H] = limpco2(z0+sqrt(-1)*eps^3*I(:,k),yobs,wobs,sys);
      fprintf('%i %e %e \n',k,imag(f)/eps^3,real(g(k)))
      for kk = 1:length(z0)
      fprintf('(%i,%i) %e %e \n',k,kk,imag(g(kk))/eps^3,real(H(k,kk)));
      end
      keyboard
      end
    %}
    z0 = init(yobs,pK,sys);
    tol = 1e-9;
    [z,J,iflag] = newtn(z0,gun,tol);
    
    if (iflag ~=0)
        fprintf('Newton''s method iflag = %i\n',iflag);
    end
    % posterior uncertainty
    warning('off');
    C = inv(J);
    warning('on');
    C = C(1:nv,1:nv);
    y = z(1:nv);
    sigy = sqrt(diag(C));
end

%-----------------------------------------------------------------

function [g,H] = grad_limpco2(z,y,w,pK,gpK,sys)
    I = eye(length(z));
    H = zeros(length(z),length(z));
    for k = 1:length(z)
        [~,g] = limpco2(z+sqrt(-1)*(eps^3)*I(:,k),y,w,pK,gpK,sys);
        H(k,:) = imag(g(:))/(eps^3);
    end
    g = real(g(:));
    %H = imag(g(:)/eps^3); % complex step 
end

function [f,g] = limpco2(z,y,w,pK,gpK,sys)
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
    nv = size(M,2);
    nlam = size(M,1)+size(K,1);
    if (ismember('sulfate',sys.abr)) % MF doesn't understand why this is here
        nlam = nlam+1;
    end
    x   =  z(1:nv);      % measureable variables
    lam =  z(nv+1:end);  % Lagrange multipliers 
    
    % Make a vector of measured quantities
    
    i = find(~isnan(y));
    y = y(i);
    
    % Make a precision matrix
    W = diag(w(i));
    
    % Build a matrix that Picks out the measured components of x
    I = eye(nv); % for chain rule
    
    PP = I(i,:); % picking/pick out the measured ones
    e = PP*x - y;

    % column vector with zeros and pK's at end
    nrk = size(K,1) ;
    zpK = zeros(nrk,1);
    zgpK = zeros(nrk,length(sys.variables));
    %zggpK = zeros(nrk,length(sys.variables),length(sys.variables));

    % pK  = [pK0;pK1;pK2;pKb;pKw;pKs;pKf;pK1p;pK2p;pK3p;pKsi;pKnh4;pKh2s;pp2f];
    %       (1)  (2) (3) (4) (5) (6) (7) (8)  (9)  (10) (11)  (12) (13) (14)
    if (ismember('K0',sys.variables))
        zpK(sys.jK0)                               = pK(1); % pK0
        zgpK([sys.jK0],[sys.iT,sys.iS,sys.iP])     = gpK(1,:); % ∂T, ∂S, ∂P
        
        zpK(sys.jK1)                               = pK(2); % pK1
        zgpK([sys.jK1],[sys.iT,sys.iS,sys.iP])     = gpK(2,:);       
        
        zpK(sys.jK2)                               = pK(3); % pK2
        zgpK([sys.jK2],[sys.iT,sys.iS,sys.iP])     = gpK(3,:);  
    end

    if (ismember('Kb',sys.variables))
        zpK(sys.jKb)                               = pK(4); % pKb
        zgpK([sys.jKb],[sys.iT,sys.iS,sys.iP])     = gpK(4,:);
    end

    if (ismember('Kw',sys.variables))
        zpK(sys.jKw)                               = pK(5); % pKw
        zgpK([sys.jKw],[sys.iT,sys.iS,sys.iP])     = gpK(5,:);  
    end

    if (ismember('Ks',sys.variables))
        zpK(sys.jKs)                               = pK(6); % pKs
        zgpK([sys.jKs],[sys.iT,sys.iS,sys.iP])     = gpK(6,:);  
    end

    if (ismember('Kf',sys.variables))
        zpK(sys.jKf)                               = pK(7); % pKf
        zgpK([sys.jKf],[sys.iT,sys.iS,sys.iP])     = gpK(7,:);
    end

    if (ismember('K1p',sys.variables))
        zpK(sys.jK1p)                              = pK(8); % pK1p
        zgpK([sys.jK1p],[sys.iT,sys.iS,sys.iP])    = gpK(8,:); 

        zpK(sys.jK2p)                              = pK(9); % pK2p
        zgpK([sys.jK2p],[sys.iT,sys.iS,sys.iP])    = gpK(9,:); 

        zpK(sys.jK3p)                              = pK(10); % pK3p
        zgpK([sys.jK3p],[sys.iT,sys.iS,sys.iP])    = gpK(10,:); 
    end

    if (ismember('Ksi',sys.variables))
        zpK(sys.jKsi)                              = pK(11); % pKsi
        zgpK([sys.jKsi],[sys.iT,sys.iS,sys.iP])    = gpK(11,:); 
    end

    if (ismember('Knh4',sys.variables))
        zpK(sys.jKnh4)                             = pK(12); % pKnh4
        zgpK([sys.jKnh4],[sys.iT,sys.iS,sys.iP])   = gpK(12,:);  
    end

    if (ismember('Kh2s',sys.variables))
        zpK(sys.Kh2s)                              = pK(13); % pKh2s
        zgpK([sys.jKh2s],[sys.iT,sys.iS,sys.iP])   = gpK(13,:); 
    end

    if (ismember('p2f',sys.variables))
        zpK(sys.jp2f)                              = pK(14); % pp2f
        zgpK([sys.jp2f],[sys.iT,sys.iS,sys.iP])    = gpK(14,:); 
    end
    
    % constraint equations
    if (ismember('hf',sys.variables))
        c = [  M * q( x ); ... 
               (-K * x) - zpK;...
               sys.f2t(x) ] ; 
    else
        c = [  M * q( x ); ...
               (-K * x) - zpK  ] ; 
    end
    
    f = 0.5 *  e.' * W * e  + lam.' * c ;  % limp, method of lagrange multipliers
    
    % -(-1/2 sum of squares) + constraint eqns, minimize f => grad(f) = 0
    
    if ( nargout > 1 ) % compute the gradient
        if (ismember('hf',sys.variables))
            gf2t = zeros(1,nv);
            gf2t(1,[sys.ih, sys.iKs, sys.iTS, sys.ihf]) = sys.gf2t(x);
            dcdx = [ M * diag( sys.dqdx( x ) ); ...
                     (-K - zgpK) ;...
                     gf2t ]; % constraint eqns wrt -log10(concentrations)
        else
            dcdx = [ M * diag( sys.dqdx( x ) ); ...
                     (-K - zgpK) ]; % c'
        end    
        g = [ e.' * W * PP +  lam.' * dcdx ,  c.' ];
        
    end
%     
%     if ( nargout > 2 ) % compute the Hessian
%         ddq =  diag( sys.d2qdx2( x ) ); % q"
%         [nr,nc] = size(M);
%         gg = zeros(nc,1);
%         for row = (1:nr)
%             gg = gg + lam(row)*diag(M(row,:))*ddq;
%         end
%         for row = (nr+1):(nr+nrk)
%             gg = gg + lam(row)*(-zggpK((row-nr),:,:)); % ggpK
%         end
%         keyboard
%         % gg is 24x42x42 right now, idk how to get it down to 42x42 (-MF)
%         if (ismember('hf',sys.variables))
%             dhfdx2 = zeros(nc,nc);
%             ii = [sys.iKs,sys.iTS];
%             dhfdx2(ii,ii) = sys.ggf2t(x);
%             gg = gg + lam(nr+1)*dhfdx2;
%         end
%         H = [  PP.'*W*PP + gg , dcdx.'  ; ...
%                dcdx         , zeros(nlam)  ];
%     end
end




function z0 = init(yobs,pk,sys);
    q = sys.q;
    p = sys.p;
    k = q(pk);
    
    pK0 = pk(1);  pK1  = pk(2);  pK2  = pk(3);  pKb  = pk(4);   pKw  = pk(5);   pKs  = pk(6);   
    pKF = pk(7);  pK1p = pk(8);  pK2p = pk(9);  pK3p = pk(10);  pKsi = pk(11);  pp2f = pk(12);
    % add pKnh3, pKh2s
    
    K0 = k(1);  K1  = k(2);  K2  = k(3);  Kb  = k(4);   Kw  = k(5);   Ks  = k(6);   
    KF = k(7);  K1p = k(8);  K2p = k(9);  K3p = k(10);  Ksi = k(11);  p2f = k(12);

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

    % solve for the [H+] ion concentration using only the carbonate alkalinity
    gam = dic/alk;
    h = 0.5*( ( gam - 1 ) * K1 + ( ( 1 - gam )^2 * K1^2 - 4 * K1 * K2 * ( 1 - 2 * gam ) ).^0.5 ) ;    
    hco3 =  h * alk / (h + 2 * K2 );
    co2st = h * hco3 / K1 ;
    co3 = 0.5 * ( alk - hco3 ) ;
    fco2 = co2st/K0;
    pco2 = fco2/p2f;
    
    y0(sys.ih)     = p(h);
    y0(sys.ihco3)  = p(hco3);
    y0(sys.ico2st) = p(co2st);
    y0(sys.ico3)   = p(co3);
    y0(sys.ifco2)  = p(fco2);
    y0(sys.ipco2)  = p(pco2);

    if (ismember('Kw',sys.variables))
        oh = Kw / h;
        y0(sys.iKw) = pKw;
        y0(sys.ioh) = p(oh);
    end
    
    if (ismember('Kb',sys.variables));
        TB = q(yobs(sys.iTB));
        boh4 = TB / ( 1 + h / Kb );
        boh3 = TB - boh4;        
        y0(sys.iKb)   = pKb;
        y0(sys.iboh3) = p(boh3);
        y0(sys.iboh4) = p(boh4);
        y0(sys.iTB)   = p(TB);
    end
    
    if (ismember('Ks',sys.variables));
        TS = q(yobs(sys.iTS));
        hf = h / ( 1 + TS / Ks );
        hso4 = TS / ( 1 + Ks / hf);
        so4  = Ks * hso4 / hf;
        y0(sys.iKs)   = pKs;
        y0(sys.iTS)   = p(TS);
        y0(sys.ihf)   = p(hf);
        y0(sys.ihso4) = p(hso4);
        y0(sys.iso4)  = p(so4);
    end
    
    if (ismember('KF',sys.variables));
        TF = q(yobs(sys.iTF));
        HF = TF / ( 1 + KF / h );
        F  = KF * HF / h;
        y0(sys.iKF) = pKF;
        y0(sys.iTF) = p(TF);
        y0(sys.iF)  = p(F);
        y0(sys.iHF) = p(HF);
    end
    
    if (ismember('K1p',sys.variables));
        TP = q(yobs(sys.iTP));
        d = ( h^3 + K1p * h^2 + K1p * K2p * h + K1p * K2p * K3p);
        h3po4 = TP * h^3 / d;
        h2po4 = TP * K1p * h^2 / d;
        hpo4  = TP * K1p * K2p * h / d;
        po4   = TP * K1p * K2p * K3p / d;
        y0(sys.iK1p)   = pK1p;
        y0(sys.iK2p)   = pK2p;
        y0(sys.iK3p)   = pK3p;
        y0(sys.iTP)    = p(TP);
        y0(sys.ih3po4) = p(h3po4);
        y0(sys.ih2po4) = p(h2po4);
        y0(sys.ihpo4)  = p(hpo4);
        y0(sys.ipo4)   = p(po4);
    end
    
    if (ismember('Ksi',sys.variables));
        TSi = q(yobs(sys.iTSi));
        siooh3 = TSi / ( 1 + h / Ksi );
        sioh4  = TSi - siooh3;
        y0(sys.iKsi)    = pKsi;
        y0(sys.isiooh3) = p(siooh3);
        y0(sys.isioh4)  = p(sioh4);
        y0(sys.iTSi)    = p(TSi);
    end

    % add the Lagrange multipliers
    nlam = size(sys.M,1)+size(sys.K,1);
    if(ismember('sulfate',sys.abr))
        nlam = nlam+1;
    end
    lam = zeros(nlam,1);
    z0 = [y0(:);lam(:)];

end

