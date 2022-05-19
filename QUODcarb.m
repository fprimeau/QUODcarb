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
    
    [K0,K1,K2,Kb,Kw,Ks,KF,K1p,K2p,K3p,KSi,p2f] = local_K(temp,sal,pres);    
    pk = p([K0,K1,K2,Kb,Kw,Ks,KF,K1p,K2p,K3p,KSi,p2f]);


    %
    % add "observations" for the equilibrium constants 
    %
    if (ismember('K0',sys.variables))
        %wK0 = 1/(1 + (0.002/pKsys(1)))^2 ; % 0.002 error on pK0, taken from literature
        wobs(sys.iK0) = (0.002).^(-2);
        yobs(sys.iK0) = p(K0);        
    end
    if (ismember('K1',sys.variables))
        wobs(sys.iK1) = (0.01).^(-2);  %wK1 = 1/(1 + (0.01/pKsys(2)))^2 ;
        yobs(sys.iK1) = p(K1);
    end
    if (ismember('K2',sys.variables))
        wobs(sys.iK2) = (0.02).^(-2);  %wK2 = 1/(1 + (0.02/pKsys(3)))^2 ;
        yobs(sys.iK2) = p(K2);
    end
    if (ismember('Kb',sys.variables))
        wobs(sys.iKb) = (0.01).^(-2);  %wKb = 1/(1 + (0.01/pKsys(4)))^2 ;
        yobs(sys.iKb) = p(Kb);
    end
    if (ismember('Kw',sys.variables))
        wobs(sys.iKw) = (0.01).^(-2);  %wKw = 1/(1 + (0.01/pKsys(5)))^2 ;
        yobs(sys.iKw) = p(Kw);
    end
    if (ismember('Ks',sys.variables))
        wobs(sys.iKs) = (0.0021).^(-2); %wKs = 1/(1 + (0.0021/pKsys(6)))^2 ;
        yobs(sys.iKs) = p(Ks);
    end
    if (ismember('KF',sys.variables))
        wobs(sys.iKF) = (0.02).^(-2);   %wKF = 1/(p(1 + 0.02/KF))^2 ; % 2% relative error
        yobs(sys.iKF) = p(KF);
    end
    if (ismember('K1p',sys.variables))
        wobs(sys.iK1p) = (0.09).^(-2);  %wK1p = 1/(1 + (0.09/pKsys(8)))^2 ;
        yobs(sys.iK1p) = p(K1p);
    end
    if (ismember('K2p',sys.variables))
        wobs(sys.iK2p) = (0.03).^(-2);  %wK2p = 1/(1 + (0.03/pKsys(9)))^2 ;
        yobs(sys.iK2p) = p(K2p);
    end
    if (ismember('K3p',sys.variables))
        wobs(sys.iK3p) = (0.02).^(-2);  %wK3p = 1/(1 + (0.02/pKsys(10)))^2 ;
        yobs(sys.iK3p) = p(K3p);
    end
    if (ismember('KSi',sys.variables))
        wobs(sys.iKSi) = (0.02).^(-2);  %wKSi = 1/(1 + (0.02/pKsys(11)))^2 ;
        yobs(sys.iKSi) = p(KSi);
    end
    
    
    gun = @(z) grad_limpco2(z,yobs,wobs,sys);
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
    z0 = init(yobs,pk,sys);
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

function [g,H] = grad_limpco2(z,y,w,sys)
    [~,g,H] = limpco2(z,y,w,sys);
    g = g(:);
end

function [f,g,H] = limpco2(z,y,w,sys)
% [f,g,H] = limpco2(z,y,w,tc,s,P)
%
% Negative log probability for the co2 system  a.k.a. log improbability i.e., limp!
%
% INPUT:
%
%   z  := system state including equilibrium constants and lagrange multipliers
%   y  := measured components, non-measured components set to NaN
%   w  := measurement precisions, same size and order as y, non-measured components set to anything, they are ignored
%
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
    if (ismember('sulfate',sys.abr))
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
    
    
    P = I(i,:); % picking/pick out the measured ones
    e = P*x - y;
    
    % constraint equations
    
    if (ismember('hf',sys.variables));
        c = [  M * q( x ); ...
               -K * x;...
               sys.f2t(x) ] ; 
    else
        c = [  M * q( x ); ...
               -K * x  ] ; 
    end
    
    f = 0.5 *  e.' * W * e  + lam.' * c ;  % limp, method of lagrange multipliers
    
    % -(-1/2 sum of squares) + constraint eqns, minimize f => grad(f) = 0
    
    if ( nargout > 1 ) % compute the gradient
        if (ismember('hf',sys.variables))
            gf2t = zeros(1,nv);
            gf2t(1,[sys.ih, sys.iKs, sys.iTS, sys.ihf]) = sys.gf2t(x);
            dcdx = [ M * diag( sys.dqdx( x ) ); ...
                     -K;...
                     gf2t ]; % constraint eqns wrt -log10(concentrations)
        else
            dcdx = [ M * diag( sys.dqdx( x ) ); ...
                     -K ]; % constraint eqns wrt -log10(concentrations)
        end    
        g = [ e.' * W * P +  lam.' * dcdx ,  c.' ];
        
    end
    if ( nargout > 2 ) % compute the Hessian
        
        ddq =  diag( sys.d2qdx2( x ) ); % q"
        [nr,nc] = size(M);
        gg = zeros(nc,1);
        for row = 1:nr
            gg = gg + lam(row)*diag(M(row,:))*ddq;
        end
        if (ismember('hf',sys.variables))
            dhfdx2 = zeros(nc,nc);
            ii = [sys.iKs,sys.iTS];
            dhfdx2(ii,ii) = sys.ggf2t(x);
            gg = gg + lam(nr+1)*dhfdx2;
        end
        H = [  P.'*W*P + gg , dcdx.'  ; ...
               dcdx         , zeros(nlam)  ];
    end
end




function z0 = init(yobs,pk,sys);
    q = sys.q;
    p = sys.p;
    k = q(pk);
    
    pK0 = pk(1);  pK1  = pk(2);  pK2  = pk(3);  pKb  = pk(4);   pKw  = pk(5);   pKs  = pk(6);   
    pKF = pk(7);  pK1p = pk(8);  pK2p = pk(9);  pK3p = pk(10);  pKSi = pk(11);  pp2f = pk(12);

    K0 = k(1);  K1  = k(2);  K2  = k(3);  Kb  = k(4);   Kw  = k(5);   Ks  = k(6);   
    KF = k(7);  K1p = k(8);  K2p = k(9);  K3p = k(10);  KSi = k(11);  p2f = k(12);

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
    
    y0(sys.ih)     = p(h);
    y0(sys.ihco3)  = p(hco3);
    y0(sys.ico2st) = p(co2st);
    y0(sys.ico3)   = p(co3);
    y0(sys.ifco2)  = p(co2st/K0);
    
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
    
    if (ismember('KSi',sys.variables));
        TSi = q(yobs(sys.iTSi));
        siooh3 = TSi / ( 1 + h / KSi );
        sioh4  = TSi - siooh3;
        y0(sys.iKSi)    = pKSi;
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

% calculate equilibrium constants-------------------------------------
function [Kh,K1,K2,Kb,Kw,Ks,Kf,K1p,K2p,K3p,Ksi,p2f] = local_K(T,S,P)
% COPIED FROM co2sys.m Orr et al. (2018)  Github
% Originally from  van Heuven et al. (2011)
% Original co2sys is from Lewis and Wallace (1998)
% 
% T = Temp (deg C) input
% P = pressure (dbar)    
    
    TK = T + 273.15; % convert to Kelvin
    Rgas = 83.1451; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
    RT = Rgas.*TK;
    Pbar = P./10;
    IonS = 19.924 .* S ./ (1000-1.005 .* S); % from DOE handbook

    % pCO2 to fCO2 conversion (Weiss 1974) valid to within 0.1%
    Pstd = 1.01325; 
    delC = (57.7 - 0.118.*TK); 
    B = -1636.75 + 12.0408.*TK - 0.0327957.*TK.^2 + 3.16528.*0.00001.*TK.^3;
    p2f = exp((B + 2.*delC).*Pstd./(RT));

    % calculate TF (Riley 1965)--------------------------------------
    TF = (0.000067./18.998).*(S./1.80655); % mol/kg-SW

    % calculate TS (Morris & Riley 1966)-----------------------------
    TS = (0.14./96.062).*(S./1.80655); % mol/kg-SW
    
    % calculate Kh (Weiss 1974)--------------------------------------
    TK100 = TK./100;
    lnKh = -60.2409 + 93.4517 ./ TK100 + 23.3585 .* log(TK100) + S .* ...
        (0.023517 - 0.023656 .* TK100 + 0.0047036 .* TK100 .^2);
    Kh = exp(lnKh);

    % calculate Ks (Dickson 1990a)----------------------------------
    lnKs = -4276.1./TK + 141.328 - 23.093 .* log(TK) + ...
        (-13856./TK + 324.57 - 47.986 .* log(TK)) .* sqrt(IonS) + ...
        (35474./TK - 771.54 + 114.723 .* log(TK)) .* IonS + ...
        (-2698./TK) .* sqrt(IonS) .* IonS + (1776./TK) .* IonS.^2;
    Ks = exp(lnKs) .* (1 - 0.001005 .* S); % converted to mol/kg-SW

    % calculate Kf (Dickson 1979)----------------------------------
    lnKf = 1590.2/TK - 12.641 + 1.525 .* IonS.^0.5;
    Kf = exp(lnKf) .* (1 - 0.001005 .* S); % converted to mol/kg-SW
   
    % pH scale conversion factors (not pressure corrected)-----------
    SWS2tot = (1 + TS ./Ks)./(1 + TS ./Ks + TF./Kf);
    FREE2tot = 1 + TS./Ks;
    
    % calculate fH (Takahashi et al 1982)--------------------------
    fH = 1.2948 - 0.002036 .* TK + (0.0004607 - ...
        0.000001475 .* TK) .* S.^2 ;

    % calculate Kb (Dickson 1990)----------------------------------
    lnKbt = -8966.9 - 2890.53 .* sqrt(S) - 77.942 .* S + ...
        1.728 .* sqrt(S) .* S - 0.0996 .* S.^2;
    lnKb = lnKbt./ TK + 148.0248 + 137.1942 .* sqrt(S) + ...
        1.62142 .* S + (-24.4344 - 25.085 .* sqrt(S) - 0.2474 .* ...
        S) .* log(TK) + 0.053105 .* sqrt(S) .* TK;
    Kb = exp(lnKb)./SWS2tot; % SWS pH scale, mol/kg-SW

    % calculate Kw (Millero 1995)--------------------------------
    lnKw = 148.9802 - 13847.26 ./ TK - 23.6521 .* log(TK) + ...
        (-5.977 + 118.67 ./ TK + 1.0495 .* log(TK)) .* ...
        sqrt(S) - 0.01615 .* S;
    Kw = exp(lnKw);

    % calculate K1p, K2p, K3p, Ksi (Yao and Millero 1995)-------
    lnK1p = -4576.752 ./TK + 115.54 - 18.453 .* log(TK) + ...
        (-106.736./TK + 0.69171) .* sqrt(S) + (-0.65643./TK - 0.01844).*S;
    K1p = exp(lnK1p);

    lnK2p = -8814.715./TK + 172.1033 - 27.927.*log(TK) + ...
        (-160.34./TK + 1.3566).*sqrt(S) + (0.37335./TK - 0.05778).*S;
    K2p = exp(lnK2p);

    lnK3p = -3070.75./TK - 18.126 + (17.27039./TK + 2.81197).*sqrt(S) + ...
        (-44.99486./TK - 0.09984).*S;
    K3p = exp(lnK3p);

    lnKsi = -8904.2./TK + 117.4 - 19.334.*log(TK) + (-458.79./TK + ...
        3.5913).*sqrt(IonS) + (188.74/TK - 1.5998).*IonS + ...
        (-12.1652./TK + 0.07871).*IonS.^2;
    Ksi = exp(lnKsi).*(1 - 0.001005.*S); % convert to mol/kg-SW

    % calculate K1 & K2 (Mehrbach refit by Dickson and Millero 1987)---
    pK1 = 3670.7 ./TK - 62.008 + 9.7944 .* log(TK) - 0.0118.*S + ...
        0.000116.*S.^2;
    K1 = 10.^(-pK1); % SWS pH scale in mol/kg-SW

    pK2 = 1394.7./TK + 4.777 - 0.0184.*S + 0.000118.*S.^2;
    K2 = 10.^(-pK2); % SWS pH scale in mol/kg-SW

    % corrections for pressure---------------------------------------
    % sources: Millero 1995, 1992, 1982, 1979; Takahashi et al. 1982;
    %   Culberson & Pytkowicz 1968; Edmond & Gieskes 1970.

    dV = -25.5 + 0.1271.*T;
    Ka = (-3.08 + 0.0877 .* T) ./1000;
    lnK1fac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on K1

    dV = -15.82 - 0.0219 .* T;
    Ka = (1.13 - 0.1475 .*T)./1000;
    lnK2fac = (-dV + 0.5.*Ka .*Pbar).*Pbar./RT; % pressure effect on K2

    dV = -29.48 + 0.1622.*T - 0.002608.*T.^2;
    Ka = -2.84./1000;
    lnKbfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on Kb

    dV = -20.02 + 0.1119.*T - 0.001409 .*T.^2;
    Ka = (-5.13 + 0.0794.*T)./1000;
    lnKwfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on Kw

    dV = -9.78 - 0.009.*T - 0.000942.*T.^2;
    Ka = (-3.91 + 0.054.*T)./1000;
    lnKffac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on Kf

    dV = -18.03 + 0.0466.*T + 0.000316.*T.^2;
    Ka = (-4.53 + 0.09.*T)./1000;
    lnKsfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on Ks

    dV = -14.51 + 0.1211.*T - 0.000321.*T.^2;
    Ka  = (-2.67 + 0.0427.*T)./1000;
    lnK1pfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on K1p

    dV = -23.12 + 0.1758.*T - 0.002647.*T.^2;
    Ka  = (-5.15 + 0.09  .*T)./1000;
    lnK2pfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on K2p

    dV = -26.57 + 0.202 .*T - 0.003042.*T.^2;
    K  = (-4.08 + 0.0714.*T)./1000;
    lnK3pfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on K3p

    dV = -29.48 + 0.1622.*T - 0.002608.*T.^2;
    K  = -2.84./1000;
    lnKsifac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on Ksi

    % correct all Ks for pressure effects
    K1fac  = exp(lnK1fac);  K1  = K1 .*K1fac;
    K2fac  = exp(lnK2fac);  K2  = K2 .*K2fac;
    Kwfac  = exp(lnKwfac);  Kw  = Kw .*Kwfac;
    Kbfac  = exp(lnKbfac);  Kb  = Kb .*Kbfac;
    Kffac  = exp(lnKffac);  Kf  = Kf .*Kffac;
    Ksfac  = exp(lnKsfac);  Ks  = Ks .*Ksfac;
    K1pfac = exp(lnK1pfac); K1p = K1p.*K1pfac;
    K2pfac = exp(lnK2pfac); K2p = K2p.*K2pfac;
    K3pfac = exp(lnK3pfac); K3p = K3p.*K3pfac;
    Ksifac = exp(lnKsifac); Ksi = Ksi.*Ksifac;

    % CorrectpHScaleConversionsForPressure:
    % fH has been assumed to be independent of pressure.
    SWS2tot  = (1 + TS./Ks)./(1 + TS./Ks + TF./Kf);
    FREE2tot =  1 + TS./Ks;

end