% new version with subfunctions and derivatives -MF

% calculate equilibrium constants-------------------------------------
function [pK,gpK] = local_pK(TC,S,P)
% base equations COPIED FROM co2sys.m Orr et al. (2018)  Github
% Originally from  van Heuven et al. (2011)
% Original co2sys is from Lewis and Wallace (1998)
% 
% INPUT:
%   T  = Temp (deg C)
%   S  = Salinity
%   P  = pressure (dbar) 
%
% OUTPUT:
%    pK  = [pK0;pK1;pK2;pKb;pKw;pKs;pKf;pK1p;pK2p;pK3p;pKsi;pKnh4;pKh2s;pp2f]; % column vector of pK's
%   gpK  = [pK_T, pK_S, pK_P]; first derivatives (gradient of pK)
    
    TK = TC + 273.15; % convert to Kelvin
    Rgas = 83.1451; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
    RT = Rgas.*TK;
    Pbar = P./10;
    IonS = 19.924 .* S ./ (1000-1.005 .* S); % from DOE handbook
    % derivatives
    RT_T = Rgas;
    IonS_S = (19.924.*(1000-1.005.*S) + 1.005*(19.924.*S) )./...
            ((1000-1.005*S).^2); % NEW as of 7/25 (bug before)

    LOG10 = log(10);
    p = @(x) -log10(x); 
    q = @(x) 10.^(-x);  % use q, i.e., a backward p
    dpdx = @(x) -1 / (x .* LOG10); % p'
    dqdx = @(x) - LOG10 * 10.^( -x );  % q'

    % pCO2 to fCO2 conversion (Weiss 1974) valid to within 0.1% ----------
    function [pp2f,gpp2f] = calc_p2f(TK,RT,RT_T)
        Pstd = 1.01325; 
        delC = (57.7 - 0.118.*TK); 
        a = -1636.75; b = 12.0408; c = -0.0327957; 
        d = 3.16528*0.00001; h = -0.118;
        %B = -1636.75 + 12.0408.*TK - 0.0327957.*TK.^2 + 3.16528.*0.00001.*TK.^3;
        B = a + b.*TK + c.*TK.^2 + d.*TK.^3;
        B_T = b + 2*c*TK + 3*d*TK^2;
        
        pp2f = -((B + 2.*delC).*Pstd./(RT))/log(10);
        pp2f_T =  (-(B + 2.*delC).* (-Pstd.*(RT.^-2)*RT_T) - ...
             (B_T + 2*h) * Pstd/RT ) / log(10); % derivative wrt T
        
        gpp2f = [pp2f_T, 0, 0]; % gradient of p(p2f);
    end
    [pp2f,gpp2f] = calc_p2f(TK,RT,RT_T);

    % calculate K0 (Weiss 1974)-------------------------------------------
    function [pK0,gpK0] = calc_pK0(TK,S)
        TK100 = TK./100;
        TK100_T = 1/100;
        a = -60.2409; b = 93.4517; c = 23.3585; d = 0.023517; 
        g = -0.023656; h = 0.0047036;
        %lnK0 = -60.2409 + 93.4517 ./ TK100 + 23.3585 .* log(TK100) + S .* ...
        %(0.023517 - 0.023656 .* TK100 + 0.0047036 .* TK100 .^2);
        pK0 = -( a + b./TK100 + c.*log(TK100) + S.*( d + g.*TK100 + ...
                h.*(TK100.^2) ) ) ./ log(10);
        pK0_T =  -(-b.*TK100_T./(TK100.^2) + c.*TK100_T./TK100 + ...
                S.* (g.*TK100_T + 2.*h.*TK100.*TK100_T) ) ./log(10); % derivative wrt T
        pK0_S = -( d + g.*TK100 + h.*(TK100.^2)  ) ./log(10); % derivative wrt S

        gpK0 = [pK0_T, pK0_S, 0];
    end
    [pK0,gpK0] = calc_pK0(TK,S);

    % calculate Ks (Dickson 1990a)----------------------------------------
    function [pKs,gpKs] = calc_pKs(TC,TK,S,IonS,IonS_S,Pbar,RT,RT_T)
        a = -4276.1; b = 141.328; c = -23.093; d = -13856;
        g = 324.57; h = -47.986; l = 35474; m = -771.54;
        n = 114.723; o = -2698; r = 1776; 
        %lnKs = -4276.1./TK + 141.328 - 23.093 .* log(TK) + ...
        %    (-13856./TK + 324.57 - 47.986 .* log(TK)) .* sqrt(IonS) + ...
        %    (35474./TK - 771.54 + 114.723 .* log(TK)) .* IonS + ...
        %    (-2698./TK) .* sqrt(IonS) .* IonS + (1776./TK) .* IonS.^2;
        sqrtIonS_S = 0.5*IonS_S/sqrt(IonS); 
        pKs = - ( a./TK + b + c.*log(TK) + ...
                ( d./TK + g + h.*log(TK) ).*sqrt(IonS) + ...
                ( l./TK + m + n.*log(TK) ).*IonS + ...
                ( o./TK ).*sqrt(IonS).*IonS + ...
                ( r./TK ).*(IonS.^2) ) ./ log(10) + p(1 - 0.001005.*S);
        pKs_T = - ( -a./(TK.^2) + c./TK + ...
                ( -d./(TK.^2) + h./TK ).*sqrt(IonS) + ...
                ( -l./(TK.^2) + n./TK ).*IonS + ...
                ( -o./(TK.^2) ).* sqrt(IonS).*IonS + ...
                ( -r./(TK.^2) ).*(IonS.^2) ) ./ log(10); 
        pKs_S = -( (d./TK + g + h.*log(TK) ).*sqrtIonS_S + ...
                ( l./TK + m + n.*log(TK) ) .* IonS_S + ...
                ( o./TK ).*sqrtIonS_S.*IonS + ( o/TK ).*sqrt(IonS).*IonS_S + ...
                ( r./TK).*(2.*IonS.*IonS_S) )./ log(10) + ...
                dpdx(1 - 0.001005 .* S).*(-0.001005);
        
        % pressure correction
        aa = -18.03; bb = 0.0466; cc = 0.000316; 
        dd = -4.53; gg = 0.09; hh = 1000;
        %dV = -18.03 + 0.0466.*T + 0.000316.*T.^2;
        %Ka = (-4.53 + 0.09.*T)./1000;
        %lnKsfac = (-dV + 0.5.*Ka.*Pbar).*Pbar./RT; % pressure effect on Ks
        pKsfac = -( (-(aa + bb.*TC + cc.*TC.^2) + 0.5.*...
                    ((dd+gg.*TC)./hh).*Pbar ).* Pbar./RT ) ./ log(10);        
        pKsfac_P = -( (-(aa + bb.*TC + cc.*TC.^2)./RT ) + ...
                    (Pbar.*((dd+gg.*TC)/hh)./RT) ) ./ log(10);
        pKsfac_T = -( (-(bb + 2.*cc.*TC) + 0.5.*Pbar.*(gg./hh) ) .* ...
                    Pbar./RT + ( -(aa + bb.*TC + cc.*TC.^2) + 0.5.* ...
                    Pbar.*((dd+gg.*TC)./hh)).*(-Pbar.*RT_T./RT.^2) ) ...
                    ./ log(10);
        
        pKs = pKs + pKsfac;
        pKs_T = pKs_T + pKsfac_T;
        pKs_P = pKsfac_P;

        gpKs = [pKs_T, pKs_S, pKs_P];
    end
    [pKs,gpKs] = calc_pKs(TC,TK,S,IonS,IonS_S,Pbar,RT,RT_T);

    % calculate Kf (Dickson 1979)----------------------------------
    function [pKf,gpKf] = calc_pKf(TC,TK,S,IonS,IonS_S,Pbar,RT,RT_T)
        sqrtIonS_S = 0.5*IonS_S/sqrt(IonS); 
        a = 1590.2; b = -12.641; c = 1.525;
        %lnKf = 1590.2/TK - 12.641 + 1.525 .* IonS.^0.5;
        pKf = -( a/TK + b + c .* sqrt(IonS) ) ./log(10) ...
                    + p(1 - 0.001005 .* S) ; % converted to mol/kg-SW
        pKf_T = -(-a/(TK.^2) )./log(10); % derivative wrt T
        pKf_S =  -c.*sqrtIonS_S ./ log(10) + ...
                    dpdx(1 - 0.001005 .*S).*(-0.001005) ;  % derivative wrt S

        % pressure correction
        aa = -9.78; bb = -0.009; cc = -0.000942; 
        dd = -3.91; gg = 0.054; hh = 1000;
        %dV = -9.78 - 0.009.*T - 0.000942.*T.^2; %Ka = (-3.91 + 0.054.*T)./1000;
        pKffac = -( (-(aa + bb.*TC + cc.*TC.^2) + 0.5.* ...
                    ((dd + gg*TC)./hh).*Pbar ) .* Pbar./RT) ./ log(10); 
        pKffac_P =  -( (-(aa + bb.*TC + cc.*TC.^2)./RT ) + ...
                    (Pbar.*((dd+gg.*TC)/hh)./RT) ) ./ log(10);
        pKffac_T = -( (-(bb + 2.*cc.*TC) + 0.5.*Pbar.*(gg./hh) ).* ...
                    Pbar./RT + (-(aa + bb.*TC + cc.*TC.^2) + 0.5.*Pbar.* ...
                    ((dd+gg.*TC)./hh)).*(-Pbar .* RT_T ./ RT.^2)) ./log(10);

        pKf = pKf + pKffac;
        pKf_T = pKf_T + pKffac_T;
        pKf_P = pKffac_P;

        gpKf = [pKf_T, pKf_S, pKf_P];
    end
    [pKf,gpKf] = calc_pKf(TC,TK,S,IonS,IonS_S,Pbar,RT,RT_T);

    % calculate TF (Riley 1965)--------------------------------------
    TF = (0.000067./18.998).*(S./1.80655); % mol/kg-SW
    TF_S = 1.95217e-6; % derivative wrt S
    pTF = p(TF);
    pTF_S = dpdx(TF).*TF_S;

    % calculate TS (Morris & Riley 1966)-----------------------------
    TS = (0.14./96.062).*(S./1.80655); % mol/kg-SW
    TS_S = 0.000806727; % derivative wrt S
    pTS = p(TS);
    pTS_S = dpdx(TS).*TS_S;
    
    % calculate fH (Takahashi et al 1982)--------------------------
    %fH = 1.2948 - 0.002036 .* TK + (0.0004607 - ...
     %   0.000001475 .* TK) .* S.^2 ;

    % calculate Kb (Dickson 1990)----------------------------------
    function [pKb,gpKb] = calc_pKb(TC, TK,S,Pbar,RT, RT_T)
        a = -8966.9; b = -2890.53; c = -77.942; d = 1.728; g = -0.0996;
        h = 148.0248; l = 137.1942; m = 1.62142; n = -24.4344;
        o = -25.085; r = -0.2474; x = 0.053105;
        %lnKbt = -8966.9 - 2890.53 .* sqrt(S) - 77.942 .* S + ...
        %    1.728 .* sqrt(S) .* S - 0.0996 .* S.^2;
        %lnKb = lnKbt./ TK + 148.0248 + 137.1942 .* sqrt(S) + ...
        %    1.62142 .* S + (-24.4344 - 25.085 .* sqrt(S) - 0.2474 .* ...
        %    S) .* log(TK) + 0.053105 .* sqrt(S) .* TK;
        lnKbtop = a + b.*sqrt(S) + c.*S + d.*sqrt(S).*S + g.*S.^2;
        pKb = - (lnKbtop./TK + h + l.*sqrt(S) + m.*S + log(TK).* ...
                (n + o.*sqrt(S) + r.*S ) + x.*sqrt(S).*TK ) ./ log(10);
        pKb_T = -(-lnKbtop./(TK.^2) + (1./TK).*(n + o.*sqrt(S) + r.*S) + ...
                x.*sqrt(S) ) ./ log(10); % derivative wrt T
        pKb_S = - ((0.5.*b./(sqrt(S).*TK)) + (c./TK) + 1.5.*d.*sqrt(S)./TK + ...
                2.*g.*S./TK + (0.5.*l./sqrt(S)) + m + ...
                log(TK).*((0.5.*o./sqrt(S)) + r) + (0.5.*x.*TK./sqrt(S)) )...
                ./ log(10); % derivative wrt S

        % pressure correction
        aa = -29.48; bb = 0.1622; cc = -0.002608; dd = -2.84; hh = 1000;
        %dV = -29.48 + 0.1622.*T - 0.002608.*T.^2; % Ka = -2.84./1000;
        pKbfac = -( (-(aa + bb.*TC + cc.*TC.^2) + 0.5.*(dd./hh).*Pbar)...
                    .*Pbar./RT)./log(10); 
        pKbfac_P = -( (-(aa + bb.*TC + cc.*TC.^2)./RT ) + ...
                      (Pbar.*(dd./hh)./RT) ) ./ log(10);
        pKbfac_T = -( (-(bb + 2.*cc.*TC) ) .* Pbar  ./ RT + ...
                      (-(aa + bb.*TC + cc.*TC.^2) + 0.5.*Pbar.*(dd./hh)).* ... 
                      (-Pbar .* RT_T ./ RT.^2)) ./log(10);
        
        pKb = pKb + pKbfac;
        pKb_T = pKb_T + pKbfac_T;
        pKb_P = pKbfac_P;

        gpKb = [pKb_T, pKb_S, pKb_P];
    end
    [pKb,gpKb] = calc_pKb(TC,TK,S,Pbar,RT,RT_T);

    % calculate Kw (Millero 1995)--------------------------------
    function [pKw,gpKw] = calc_pKw(TC,TK,S,Pbar,RT,RT_T)
        a = 148.9802; b = -13847.26; c = -23.6521; d = -5.977;
        g = 118.67; h = 1.0495; l = -0.01615;
        %lnKw = 148.9802 - 13847.26 ./ TK - 23.6521 .* log(TK) + ...
        %    (-5.977 + 118.67 ./ TK + 1.0495 .* log(TK)) .* ...
        %    sqrt(S) - 0.01615 .* S;
        pKw = - (a + b./TK + c.*log(TK) + ...
                ( d + g./TK + h.*log(TK) ).*sqrt(S) + l.*S ) ./ log(10);
        pKw_T =  -( -b./(TK.^2) + c./TK + (-g./(TK.^2) + h./TK  ) .* ...
            sqrt(S) )./ log(10); % derivative wrt T
        pKw_S =  -( (d + g./TK + h.*log(TK)) .* (0.5./sqrt(S)) + l  ) ...
                ./ log(10); % derivative wrt S

        % pressure correction
        aa = -20.02; bb = 0.1119; cc = -0.001409; 
        dd = -5.13; gg = 0.0794; hh = 1000;
        % dV = -20.02 + 0.1119.*TC - 0.001409 .*TC.^2; % Ka = (-5.13 + 0.0794.*TC)./1000;
        pKwfac = - ( (-(aa + bb.*TC + cc.*TC.^2) + 0.5.*((dd+gg.*TC)./hh)...
                    .*Pbar).*Pbar./RT ) ./ log(10) ; 
        pKwfac_P =  -( (-(aa + bb.*TC + cc.*TC.^2)./RT ) + ...
                    (Pbar.*((dd+gg.*TC)/hh)./RT) ) ./ log(10);
        pKwfac_T = -( (-(bb + 2.*cc.*TC) + 0.5.*Pbar.*(gg./hh) ) ...
                     .* Pbar ./ RT + (-(aa + bb.*TC + cc.*TC.^2) + ...
                     0.5.*Pbar.*((dd+gg.*TC)./hh)).* ... 
                     (-Pbar .* RT_T ./ RT.^2)) ./log(10);
        pKw = pKw + pKwfac;
        pKw_T = pKw_T + pKwfac_T;
        pKw_P = pKwfac_P;

        gpKw = [pKw_T,pKw_S,pKw_P];
    end
    [pKw,gpKw] = calc_pKw(TC,TK,S,Pbar,RT,RT_T);

    % calculate K1p (Yao and Millero 1995)--------------------------------
    function [pK1p,gpK1p] = calc_pK1p(TC,TK,S,Pbar,RT,RT_T)
        a = -4576.752; b = 115.54; c = -18.453; d = -106.736;
        g = 0.69171; h = -0.65643; l = -0.01844;
        %lnK1p = -4576.752 ./TK + 115.54 - 18.453 .* log(TK) + ...
        %    (-106.736./TK + 0.69171) .* sqrt(S) + (-0.65643./TK - 0.01844).*S;
        
        pK1p = - (a./TK + b + c.*log(TK) + ...
                (d./TK + g ).*sqrt(S) + (h./TK + l).*S ) ./ log(10);
        pK1p_T = - ( -a./(TK.^2) + c./TK + sqrt(S) .* (-d./(TK.^2)) + ...
                (-h./(TK.^2)).*S ) ./ log(10); % derivative wrt T
        pK1p_S = - ( (d./TK + g).*(0.5./sqrt(S)) + ...
                (h./TK + l) ) ./ log(10);

        % pressure correction
        aa = -14.51; bb = 0.1211; cc = -0.000321; 
        dd = -2.67; gg = 0.0427; hh = 1000;
        % dV = -14.51 + 0.1211.*TC - 0.000321.*TC.^2; % Ka  = (-2.67 + 0.0427.*TC)./1000;
        pK1pfac = - ( ( -(aa + bb.*TC + cc.*TC.^2) + 0.5.* ...
                    ((dd + gg.*TC)./hh).*Pbar ) .*Pbar./RT ) ./ log(10); 
        pK1pfac_P =  -( (-(aa + bb.*TC + cc.*TC.^2)./RT ) + ...
                    (Pbar.*((dd+gg.*TC)/hh)./RT) ) ./ log(10);
        pK1pfac_T = -( (-(bb + 2.*cc.*TC) + 0.5.*Pbar.*(gg./hh) ) .* ...
                    Pbar./RT + (-(aa + bb.*TC + cc.*TC.^2) + 0.5.* ...
                    Pbar.*((dd+gg.*TC)./hh)).* ... 
                    (-Pbar .* RT_T ./ RT.^2)) ./log(10);
        pK1p = pK1p + pK1pfac;
        pK1p_T = pK1p_T + pK1pfac_T;
        pK1p_P = pK1pfac_P;

        gpK1p = [pK1p_T,pK1p_S,pK1p_P];
    end
    [pK1p,gpK1p] = calc_pK1p(TC,TK,S,Pbar,RT,RT_T);

    % calculate K2p (Yao and Millero 1995)--------------------------------
    function [pK2p,gpK2p] = calc_pK2p(TC,TK,S,Pbar,RT,RT_T)
        a = -8814.715; b = 172.1033; c = -27.927; d = -160.34;
        g = 1.3566; h = 0.37335; l = -0.05778;
        %lnK2p = -8814.715./TK + 172.1033 - 27.927.*log(TK) + ...
        %    (-160.34./TK + 1.3566).*sqrt(S) + (0.37335./TK - 0.05778).*S;
        pK2p = -( a./TK + b + c.*log(TK) + (d./TK + g).*sqrt(S) + ...
                (h./TK + l).*S ) ./log(10);
        pK2p_T = -( -a./(TK.^2) + c./TK + sqrt(S).* (-d./(TK.^2)) + ...
                (-h./(TK.^2)).*S ) ./ log(10);
        pK2p_S = -( (d./TK + g) .* ( 0.5./sqrt(S) ) + ...
                (h./TK + l) ) ./ log(10);

        % pressure correction
        aa = -23.12; bb = 0.1758; cc = -0.002647;
        dd = -5.15; gg = 0.09; hh = 1000;
        % dV = -23.12 + 0.1758.*TC - 0.002647.*TC.^2; % Ka  = (-5.15 + 0.09  .*TC)./1000;
        pK2pfac = -( (-(aa + bb.*TC + cc.*TC.^2) + 0.5.*((dd + gg.*TC)/hh) ...
                .*Pbar) .*Pbar./RT ) ./ log(10); 
        pK2pfac_P =  -( (-(aa + bb.*TC + cc.*TC.^2)./RT ) + ...
                      (Pbar.*((dd+gg.*TC)/hh)./RT) ) ./ log(10);
        pK2pfac_T = -( (-(bb + 2.*cc.*TC) + 0.5.*Pbar.*(gg./hh) ) .* Pbar  ./ RT + ...
                      (-(aa + bb.*TC + cc.*TC.^2) + 0.5.*Pbar.*((dd+gg.*TC)./hh)).* ... 
                      (-Pbar .* RT_T ./ RT.^2)) ./log(10);
        pK2p = pK2p + pK2pfac;
        pK2p_T = pK2p_T + pK2pfac_T;
        pK2p_P = pK2pfac_P;

        gpK2p = [pK2p_T, pK2p_S, pK2p_P];
    end
    [pK2p,gpK2p] = calc_pK2p(TC,TK,S,Pbar,RT,RT_T);

    % calculate K3p (Yao and Millero 1995)--------------------------------
    function [pK3p,gpK3p] = calc_pK3p(TC,TK,S,Pbar,RT,RT_T)
        a = -3070.75; b = -18.126; c = 17.27039; d = 2.81197;
        g = -44.99486; h = -0.09984;
        %lnK3p = -3070.75./TK - 18.126 + (17.27039./TK + 2.81197).*sqrt(S) + ...
        %    (-44.99486./TK - 0.09984).*S;
        pK3p = -( a./TK + b + (c./TK + d).*sqrt(S) +...
                (g./TK + h) .*S ) ./ log(10);
        pK3p_T = - (-a./(TK.^2) + (-c./(TK.^2).*sqrt(S) + ...
                (-g./(TK.^2)).*S )  )./ log(10);
        pK3p_S = - ( (c./TK + d) .* (0.5./sqrt(S)) + ...
                (g./TK + h) ) ./ log(10);
        
        % pressure correction
        aa = -26.57; bb = 0.202; cc = -0.003042; 
        dd = -4.08; gg = 0.0714; hh = 1000;
        % dV = -26.57 + 0.202 .*TC - 0.003042.*TC.^2; %Ka  = (-4.08 + 0.0714.*TC)./1000;
        pK3pfac = -( (-(aa + bb.*TC + cc.*TC.^2) + 0.5.* ...
                    ((dd + gg.*TC)./hh) .*Pbar) .* Pbar./RT ) ./log(10);
        pK3pfac_P =  -( (-(aa + bb.*TC + cc.*TC.^2)./RT ) + ...
                    (Pbar.*((dd+gg.*TC)/hh)./RT) ) ./ log(10);
        pK3pfac_T = -( (-(bb + 2.*cc.*TC) + 0.5.*Pbar.*(gg./hh) ) .* ...
                    Pbar  ./ RT + (-(aa + bb.*TC + cc.*TC.^2) + 0.5.* ...
                    Pbar.*((dd+gg.*TC)./hh)).* ... 
                    (-Pbar .* RT_T ./ RT.^2)) ./log(10);
        pK3p = pK3p + pK3pfac;
        pK3p_T = pK3p_T + pK3pfac_T;
        pK3p_P = pK3pfac_P;

        gpK3p = [pK3p_T, pK3p_S, pK3p_P]; % gradient pK3p
    end
   [pK3p,gpK3p] = calc_pK3p(TC,TK,S,Pbar,RT,RT_T);
   % S and P are fine but T is NOT
    
    % calculate Ksi (Yao and Millero 1995)--------------------------------
    function [pKsi,gpKsi] = calc_pKsi(TC,TK,S,IonS,IonS_S,Pbar,RT,RT_T)
        sqrtIonS_S = 0.5*IonS_S/sqrt(IonS);
        a = -8904.2; b = 117.4; c = -19.334; d = -458.79;
        g = 3.5913; h = 188.74; l = -1.5998; m = -12.1652; n = 0.07871;
        %lnKsi = -8904.2./TK + 117.4 - 19.334.*log(TK) + (-458.79./TK + ...
        %    3.5913).*sqrt(IonS) + (188.74/TK - 1.5998).*IonS + ...
        %    (-12.1652./TK + 0.07871).*IonS.^2;
        pKsi = - ( a./TK + b + c.*log(TK) + (d./TK + g).*sqrt(IonS) + ...
            (h./TK + l).*IonS + (m./TK + n).*IonS.^2 ) ./ log(10) + ...
            p(1 - 0.001005.*S) ; % convert to mol/kg-SW
        pKsi_T = - ( -a./(TK.^2) + c./TK + sqrt(IonS).*(-d./(TK.^2)) + ...
            IonS.*(-h./(TK.^2)) + IonS.^2.*(-m./(TK.^2)) ) ./ log(10) ;
       pKsi_S = - ( (d./TK + g).*(sqrtIonS_S) + (h./TK + l).*IonS_S + ...
            2.*IonS.*IonS_S.*(m./TK + n) ) ...
            ./ log(10) + dpdx(1 - 0.001005.*S).*(-0.001005) ;

        % pressure correction
        aa = -29.48; bb = 0.1622; cc = -0.002608;
        dd = -2.84; hh = 1000;
        % dV = -29.48 + 0.1622.*TC - 0.002608.*TC.^2; % Ka  = -2.84./1000;
        pKsifac = -( (-(aa + bb.*TC + cc.*TC.^2) + 0.5.* ...
                (dd/hh).*Pbar).*Pbar./RT ) ./ log(10);
        pKsifac_P = -( (-(aa + bb.*TC + cc.*TC.^2)./RT ) + ...
                      (Pbar.*(dd./hh)./RT) ) ./ log(10);
        pKsifac_T = -( (-(bb + 2.*cc.*TC) ) .* Pbar  ./ RT + ...
                      (-(aa + bb.*TC + cc.*TC.^2) + 0.5.*Pbar.*(dd./hh)).* ... 
                      (-Pbar .* RT_T ./ RT.^2)) ./log(10);
        
        pKsi = pKsi + pKsifac;
        pKsi_T = pKsi_T + pKsifac_T;
        pKsi_P = pKsifac_P;

        gpKsi = [pKsi_T, pKsi_S, pKsi_P]; % gradient pKsi
    end
    [pKsi,gpKsi] = calc_pKsi(TC,TK,S,IonS,IonS_S,Pbar,RT,RT_T);

    % calculate pK1 (Mehrbach refit by Dickson and Millero 1987)---
    function [pK1,gpK1] = calc_pK1(TC,TK,S,Pbar,RT,RT_T)
        a = 3670.7; b = -62.008; c = 9.7944; d = -0.0118; g = 0.000116;
        %pK1 = 3670.7 ./TK - 62.008 + 9.7944 .* log(TK) - 0.0118.*S + ...
        %    0.000116.*S.^2;
        pK1 = a./TK + b + c.*log(TK) + d.*S + g.*S.^2;
        pK1_T = -a./(TK.^2) + c./TK;
        pK1_S = d + 2.*g.*S;

        % pressure correction
        aa = -25.5; bb = 0.1271; dd = -3.08; gg = 0.0877; hh = 1000;
        % dV = -25.5 + 0.1271.*TC; % Ka = (-3.08 + 0.0877 .* TC) ./1000;
        pK1fac = - ( (-(aa + bb.*TC) + 0.5.*((dd + gg.*TC)./hh) ...
                .*Pbar) .*Pbar./RT ) ./ log(10);
        pK1fac_T = - ( (-(bb) + 0.5.*(gg./hh).*Pbar).*Pbar./RT + ...
                (-(aa + bb.*TC) + 0.5.*((dd + gg.*TC)./hh) ...
                .*Pbar) .*(-Pbar.*RT_T./(RT.^2))  ) ./ log(10);
        pK1fac_P = - ( 0.5.*((dd + gg.*TC)./hh) .* Pbar./RT + ...
                (-(aa + bb.*TC) + 0.5.*((dd + gg.*TC)./hh) ...
                .*Pbar) .*(1./RT) ) ./ log(10);
        pK1 = pK1 + pK1fac;
        pK1_T = pK1_T + pK1fac_T;
        pK1_P = pK1fac_P;

        gpK1 = [pK1_T, pK1_S, pK1_P];
    end
    [pK1,gpK1] = calc_pK1(TC,TK,S,Pbar,RT,RT_T);

    % calculate pK2 (Mehrbach refit by Dickson and Millero 1987)---
    function [pK2,gpK2] = calc_pK2(TC,TK,S,Pbar,RT,RT_T)
        a = 1394.7; b = 4.777; c = -0.0184; d = 0.000118;
        %pK2 = 1394.7./TK + 4.777 - 0.0184.*S + 0.000118.*S.^2;
        pK2 = a./TK + b + c.*S + d.*S.^2;
        pK2_T = -a./(TK.^2);
        pK2_S = c + 2*d*S;

        % pressure correction
        aa = -15.82; bb = -0.0219; dd = 1.13; gg = -0.1475; hh = 1000;
        % dV = -15.82 - 0.0219 .* TC; % Ka = (1.13 - 0.1475 .*TC)./1000;
        pK2fac = - ( (-(aa + bb.*TC) + 0.5.*((dd + gg.*TC)./hh) .*Pbar) ...
            .*Pbar./RT ) ./ log(10) ; 
        pK2fac_T = - ( (-(bb) + 0.5.*(gg./hh).*Pbar).*Pbar./RT + ...
                (-(aa + bb.*TC) + 0.5.*((dd + gg.*TC)./hh) ...
                .*Pbar) .*(-Pbar.*RT_T./(RT.^2))  ) ./ log(10);
        pK2fac_P = - ( 0.5.*((dd + gg.*TC)./hh) .* Pbar./RT + ...
                (-(aa + bb.*TC) + 0.5.*((dd + gg.*TC)./hh) ...
                .*Pbar) .*(1./RT) ) ./ log(10);
        
        pK2 = pK2 + pK2fac;
        pK2_T = pK2_T + pK2fac_T;
        pK2_P = pK2fac_P;

        gpK2 = [pK2_T, pK2_S, pK2_P]; % gradient pK2 with pressure correction
    end
    [pK2,gpK2] = calc_pK2(TC,TK,S,Pbar,RT,RT_T);

    % Ammonia, added by Sharp et al 2021, from Clegg and Whitfield (1995)
    function [pKnh4, gpKnh4] = calc_pKnh4(TC,TK,S,Pbar,RT,RT_T)
        a = 9.44605; b = -2729.33; c = 1/298.15; d = 0.04203362;
        g = -11.24742; h = -13.6416; l = 1.176949; m = -0.02860785;
        n = 545.4834; o = -0.1462507; r = 0.0090226468;
        t = -0.0001471361; v = 10.5425; w = 0.004669309;
        x = -0.0001691742; y = -0.5677934; z = -2.354039e-05;
        az = 0.009698623;
        % pKnh4 = 9.44605 - 2729.33 * (1/298.15-1/TK)) +...
                % (0.04203362 - 11.24742/TK) * S^0.25 +...
                % (-13.6416 + 1.176949 * sqrt(TK) - ...
                % 0.02860785*TK + 545.4834/TK)*sqrt(S) + ...
                % (-0.1462507 + 0.0090226468*sqrt(TK) - ...
                % 0.0001471361*TK + 10.5425/TK)*sqrt(S)*S + ...
                % (0.004669309 - 0.0001691742*sqrt(TK) - ...
                % 0.5677934/TK) * (S^2) + ...
                % (-2.354039e-05 + 0.009698623/TK)*(S^2.5); % total pH scale
         % pKnh4 = pKnh4 + p(1-0.001005*S); % convert to mol/kg-SW
         pKnh4 = a + b.* (c - 1./TK) + (d + g./TK).*S.^(0.25) + ...
                (h + l.*sqrt(TK) + m.*TK + n./TK ).*sqrt(S) + ...
                (o + r.*sqrt(TK) + t.*TK + v./TK ).*sqrt(S).*S + ...
                (w + x.*sqrt(TK) + y./TK ).*(S.^2) + ...
                (z + az./TK) .*(S.^(2.5)) ...
                + p(1 - 0.001005.*S); % convert to mol/kg-SW, on the total scale
         pKnh4_T = b./(TK.^2) - g./(TK.^2).*S.^(0.25) + ...
                (0.5*l./sqrt(TK) + m - n./(TK.^2) ).*sqrt(S) + ...
                (0.5*r./sqrt(TK) + t - v./(TK.^2) ).*sqrt(S).*S + ...
                (0.5*x./sqrt(TK) - y./(TK.^2) ).*(S.^2) + ...
                (-az./(TK.^2)).*(S.^(2.5));
         pKnh4_S = (d + g./TK).*( 0.25./(S.^(0.75)) ) + ...
                (h + l.*sqrt(TK) + m.*TK + n./TK ).*(0.5./sqrt(S)) + ...
                (o + r.*sqrt(TK) + t.*TK + v./TK ).*(1.5.*sqrt(S)) + ...
                (w + x.*sqrt(TK) + y./TK ).*(2.*S) + ...
                (z + az./TK) .*(2.5.*S.^(1.5)) ...
                + dpdx(1 - 0.001005.*S).*(-0.001005);

         % pressure correction
         aa = -26.43; bb = 0.0889; cc = -0.000905;
         dd = -5.03; gg = 0.0814; hh = 1000;
         % dV = -26.43 + 0.0889*TC - 0.000905*TC^2;
         % Ka = (-5.03 + 0.0814*TC)/1000;
         pKnh4fac = - ( ( -(aa + bb.*TC + cc.*TC.^2) + 0.5.*...
                    ((dd + gg.*TC)./hh).*Pbar ) .*Pbar./RT ) ./ log(10);
         pKnh4fac_P = -( (-(aa + bb.*TC + cc.*TC.^2)./RT ) + ...
                    (Pbar.*((dd+gg.*TC)/hh)./RT) ) ./ log(10);
         pKnh4fac_T = -( (-(bb + 2.*cc.*TC) + 0.5.*Pbar.*(gg./hh) ) ...
                    .* Pbar  ./ RT + (-(aa + bb.*TC + cc.*TC.^2) + ...
                    0.5.*Pbar.*((dd+gg.*TC)./hh)).* ... 
                    (-Pbar .* RT_T ./ RT.^2)) ./log(10);
        pKnh4 = pKnh4 + pKnh4fac;
        pKnh4_T = pKnh4_T + pKnh4fac_T;
        pKnh4_P = pKnh4fac_P;

        gpKnh4 = [pKnh4_T, pKnh4_S, pKnh4_P];
    end
    [pKnh4,gpKnh4] = calc_pKnh4(TC,TK,S,Pbar,RT,RT_T);

    % hydrogen sulfide, added by Sharp et al 2021, from Millero et al (1988)
    function [pKh2s,gpKh2s] = calc_pKh2s(TC,TK,S,Pbar,RT,RT_T)
        a = 225.838; b = -13275.3; c = -34.6435;
        d = 0.3449; h = -0.0274;
        % lnKh2s = (225.838 - 13275.3/TK - 34.6435*log(TK) + ...
        %           0.3449*sqrt(S) - 0.0274*S); % total scale
        pKh2s = - (a + b./TK + c.*log(TK) + d.*sqrt(S) + ...
                    h.*S) ./ log(10);
        pKh2s_T = - (-b./(TK.^2) + c./TK ) ./ log(10);
        pKh2s_S = - (0.5.*d./sqrt(S) + h) ./ log(10);

        % pressure correction
        aa = -11.07; bb = -0.009; cc = -0.000942;
        dd = -2.89; gg = 0.054; hh = 1000;
        % dV = -11.07 - 0.009*TC - 0.000942*TC^2; Ka = (-2.89 + 0.054*TC)/1000;
        pKh2sfac = - ( ( -(aa + bb.*TC + cc.*TC.^2) + 0.5.*...
                    ((dd + gg.*TC)./hh).*Pbar ) .* Pbar./RT ) ./ log(10);
        pKh2sfac_P = -( (-(aa + bb.*TC + cc.*TC.^2)./RT ) + ...
                    (Pbar.*((dd+gg.*TC)/hh)./RT) ) ./ log(10);
        pKh2sfac_T = -( (-(bb + 2.*cc.*TC) + 0.5.*Pbar.*(gg./hh) ) .* ...
                    Pbar  ./ RT + (-(aa + bb.*TC + cc.*TC.^2) + 0.5.* ...
                    Pbar.*((dd+gg.*TC)./hh)).* ... 
                    (-Pbar .* RT_T ./ RT.^2)) ./log(10);

        pKh2s = pKh2s + pKh2sfac;
        pKh2s_T = pKh2s_T + pKh2sfac_T;
        pKh2s_P = pKh2sfac_P;

        gpKh2s = [pKh2s_T, pKh2s_S, pKh2s_P];
    end
    [pKh2s,gpKh2s] = calc_pKh2s(TC,TK,S,Pbar,RT,RT_T);

    % corrections for pressure sources------------------------------------
    % Millero 1995, 1992, 1982, 1979; Takahashi et al. 1982;
    %   Culberson & Pytkowicz 1968; Edmond & Gieskes 1970.
    
    % pH scale conversion factors (now pressure corrected)-----------
    function [pSWS2tot,gpSWS2tot] = calc_pSWS2tot(pTS,pTS_S,pTF,pTF_S,pKs,gpKs,pKf,gpKf)
        % FREE2tot = 1 + TS./Ks; % not needed for right now
        top = 1 + q(pTS - pKs) ;
        top_T = dqdx(pTS - pKs) * (-gpKs(1)) ;
        top_S = dqdx(pTS - pKs) * (pTS_S - gpKs(2)) ;
        bot = top + q(pTF - pKf) ;
        bot_T = top_T + dqdx(pTF - pKf) * (-gpKf(1)) ;
        bot_S = top_S + dqdx(pTF - pKf) * (pTF_S - gpKf(2));

        pSWS2tot = p(top./bot);
        pSWS2tot_T = (-top_T / (top*log(10)) ) + (bot_T / (bot*log(10)) ); 
        pSWS2tot_S = (-top_S / (top*log(10)) ) + (bot_S / (bot*log(10)) );

        gpSWS2tot = [pSWS2tot_T, pSWS2tot_S,0];
    end
    [pSWS2tot,gpSWS2tot] = calc_pSWS2tot(pTS,pTS_S,pTF,pTF_S,pKs,gpKs,pKf,gpKf);     
    
    % convert from SWS to pH total scale ----------------------------------
    pK1 = pK1 + pSWS2tot;
    gpK1(1) = gpK1(1) + gpSWS2tot(1); gpK1(2) = gpK1(2) + gpSWS2tot(2);

    pK2 = pK2 + pSWS2tot;
    gpK2(1) = gpK2(1) + gpSWS2tot(1); gpK2(2) = gpK2(2) + gpSWS2tot(2);
    
    pKw = pKw + pSWS2tot;
    gpKw(1) = gpKw(1) + gpSWS2tot(1); gpKw(2) = gpKs(2) + gpSWS2tot(2);

    pK1p = pK1p + pSWS2tot;
    gpK1p(1) = gpK1p(1) + gpSWS2tot(1); gpK1p(2) = gpK1p(2) + gpSWS2tot(2);

    pK2p = pK2p + pSWS2tot;
    gpK2p(1) = gpK2p(1) + gpSWS2tot(1); gpK2p(2) = gpK2p(2) + gpSWS2tot(2);

    pK3p = pK3p + pSWS2tot;
    gpK3p(1) = gpK3p(1) + gpSWS2tot(1); gpK3p(2) = gpK3p(2) + gpSWS2tot(2);

    pKsi = pKsi + pSWS2tot;
    gpKsi(1) = gpKsi(1) + gpSWS2tot(1); gpKsi(2) = gpKsi(2) + gpSWS2tot(2);

    % p(K) = -log10(K); % K = equilibrium constants
    pK = [pK0;    pK1;  pK2;  pKb;  pKw;  pKs;  pKf;  pK1p;  pK2p; ...
             pK3p;  pKsi;  pKnh4;  pKh2s;  pp2f];
    % g(pK) = gradient(pK); % first derivatives of pK-> [d/dT, d/dS, d/dP];
    gpK = [gpK0; gpK1; gpK2; gpKb; gpKw; gpKs; gpKf; gpK1p; gpK2p; ...
            gpK3p; gpKsi; gpKnh4; gpKh2s; gpp2f];
 
end

function [ggpK] = local_ggpK(TC,S,P)
    % INPUT:
    %   T  = Temp (deg C)
    %   S  = Salinity
    %   P  = pressure (dbar) 

    % OUTPUT:
    % ggpK = second derivatives of pK using complex step method
    %        -> [d/dTdT, d/dTdS, d/dTdP, d/dSdS, d/dSdP, d/dPdP];

    % set up complex variables
    TCi = TC + sqrt(-1)*eps^3;
    Si = S + sqrt(-1)*eps^3;
    Pi = P + sqrt(-1)*eps^3;

    [~,gpK_T] = local_pK(TCi,S,P); % imaginary in T
    [~,gpK_S] = local_pK(TC,Si,P); % imaginary in S
    [~,gpK_P] = local_pK(TC,S,Pi); % imaginary in P

    ggpK = zeros(14,6);

    % ggpK0
    ggpK(:,1:3) = imag(gpK_T(:,1:3))./eps^3; % K_TT, K_ST, K_PT
    ggpK(:,4) = imag(gpK_S(:,2))./eps^3; % K_SS 
    ggpK(:,5:6) = imag(gpK_P(:,2:3))./eps^3; % K_SP, K_PP

end


