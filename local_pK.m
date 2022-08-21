function [pK,gpK] = local_pK(T,S,P)
% base equations COPIED FROM co2sys.m Orr et al. (2018)  Github
% Originally from  van Heuven et al. (2011)
% Original co2sys is from Lewis and Wallace (1998)

% INPUT:
%   T = Temp (deg C)
%   S  = Salinity
%   P  = pressure (dbar)

% OUTPUT:
%    pK  = [pK0;pK1;pK2;pKb;pKw;pKs;pKf;pK1p;pK2p;pK3p;pKsi;pKnh4;pKh2s;pp2f];
%   gpK  = [pK_T, pK_S, pK_P]; first derivatives (gradient of pK)

    TK   = T + 273.15; % convert to Kelvin
    Rgas = 83.1451; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
    RT   = Rgas * TK;
    RT_T = Rgas;
    
    Pbar = P / 10; % convert from dbar to bar
    ions   = @(S) 19.924 * S / ( 1000 - 1.005 * S); % from DOE handbook
    ions_S = @(S) 19.924 / ( 1000 - 1.005 * S ) + ...
             1.005 * 19.924 * S / ( 1000 -1.005 * S )^2 ; 
    
    LOG10 = log(10);
    p = @(x) -log10(x);
    q = @(x) 10.^(-x);  % inverse p, i.e., a backward p
    dpdx = @(x) -1 / (x * LOG10);     % p'
    dqdx = @(x) -LOG10 * 10.^( -x );  % q'
    
    
    % corrections for pressure------------------------------------
    % Millero 1995, 1992, 1982, 1979; Takahashi et al. 1982;
    %   Culberson & Pytkowicz 1968; Edmond & Gieskes 1970.
    dV   = @(T,a) a(1) + a(2) * T +     a(3) * T^2; 
    dV_T = @(T,a)        a(2)     + 2 * a(3) * T;
    Ka   = @(T,b) ( b(1) + b(2) * T ) / 1000;
    Ka_T = @(T,b) b(2) / 1000;
    ppfac   = @(T,Pbar,a,b)... 
              -( -dV(T,a)   * Pbar   + 0.5 * Ka(T,b)   * Pbar^2 )      / ( RT   * LOG10 );
    ppfac_T = @(T,Pbar,a,b)... 
              -( -dV_T(T,a) * Pbar   + 0.5 * Ka_T(T,b) * Pbar^2 )      / ( RT   * LOG10 ) ...
              -( -dV(T,a)   * Pbar   + 0.5 * Ka(T,b)   * Pbar^2 ) * (-RT_T / ( RT^2 * LOG10 ));
    ppfac_P = @(T,Pbar,a,b)... 
              -( -dV(T,a)         +       Ka(T,b)   * Pbar ) / (10 * RT * LOG10 ) ;
    

    

    
    % compute the pK's and their derivatives w.r.t. T,P,and S
    [pp2f    , gpp2f    ] = calc_p2f(T,RT,RT_T);
    [pK0     , gpK0     ] = calc_pK0(T,S);
    [pKs     , gpKs     ] = calc_pKs(T,S,Pbar,RT,RT_T);
    [pKf     , gpKf     ] = calc_pKf(T,S,Pbar,RT,RT_T);
    [pKb     , gpKb     ] = calc_pKb(T,S,Pbar,RT,RT_T);
    [pKw     , gpKw     ] = calc_pKw(T,S,Pbar,RT,RT_T);
    [pK1p    , gpK1p    ] = calc_pK1p(T,S,Pbar,RT,RT_T);
    [pK2p    , gpK2p    ] = calc_pK2p(T,S,Pbar,RT,RT_T);
    [pK3p    , gpK3p    ] = calc_pK3p(T,S,Pbar,RT,RT_T);
    [pKsi    , gpKsi    ] = calc_pKsi(T,S,Pbar,RT,RT_T);
    [pK1     , gpK1     ] = calc_pK1(T,S,Pbar,RT,RT_T);
    [pK2     , gpK2     ] = calc_pK2(T,S,Pbar,RT,RT_T);
    [pKnh4   , gpKnh4   ] = calc_pKnh4(T,S,Pbar,RT,RT_T);
    [pKh2s   , gpKh2s   ] = calc_pKh2s(T,S,Pbar,RT,RT_T);
    [pSWS2tot, gpSWS2tot] = calc_pSWS2tot(S,pKs,gpKs,pKf,gpKf);

    % convert from SWS to pH total scale ----------------------------------
    pK1  = pK1  + pSWS2tot;  gpK1  = gpK1  + gpSWS2tot;
    pK2  = pK2  + pSWS2tot;  gpK2  = gpK2  + gpSWS2tot;
    pKw  = pKw  + pSWS2tot;  gpKw  = gpKw  + gpSWS2tot;
    pK1p = pK1p + pSWS2tot;  gpK1p = gpK1p + gpSWS2tot;
    pK2p = pK2p + pSWS2tot;  gpK2p = gpK2p + gpSWS2tot;
    pK3p = pK3p + pSWS2tot;  gpK3p = gpK3p + gpSWS2tot;
    pKsi = pKsi + pSWS2tot;  gpKsi = gpKsi + gpSWS2tot;

    
    %
    % output
    %
    pK = [pK0;    pK1;  pK2;  pKb;  pKw;  pKs;  pKf;  pK1p;  pK2p; ...
          pK3p;  pKsi;  pKnh4;  pKh2s;  pp2f];
    gpK = [gpK0; gpK1; gpK2; gpKb; gpKw; gpKs; gpKf; gpK1p; gpK2p; ...
           gpK3p; gpKsi; gpKnh4; gpKh2s; gpp2f];

    
    
    %
    % subfunctions
    %
 
    function [pp2f,gpp2f] = calc_p2f(T,RT,RT_T)
    % pCO2 to fCO2 conversion (Weiss 1974) valid to within 0.1% ----------

        TK = T + 273.15; % convert to Kelvin
        Pstd = 1.01325;
        delC = (57.7 - 0.118.*TK);
        delC_T = -0.118;

        a = [ -1636.75; 12.0408; 0.0327957; 3.16528e-5 ]; 
        
        f  = @(T,a) a(1) + a(2) * T +     a(3) * T^2 +     a(4) * T^3;
        df = @(T,a)        a(2)     + 2 * a(3) * T   + 3 * a(4) * T^2;
  
        pp2f   = -( ( f(TK,a) + 2 * delC ) * Pstd / RT ) / LOG10;
        pp2f_T = -( ( df(TK,a) + 2 * delC_T ) * Pstd / RT - ...
                    ( f(TK,a) + 2 * delC ) * Pstd * RT_T / RT^2 ) / LOG10;
        pp2f_S = 0;
        pp2f_P = 0;
        
        gpp2f = [pp2f_T, pp2f_S, pp2f_P]; % gradient of p(p2f);
    end
    
    function [pK0,gpK0] = calc_pK0(T,S)
    % calculate K0 (Weiss 1974)-------------------------------------------
        TK = T + 273.15; % convert to Kelvin
        TK100 = TK./100;
        TK100_T = 1/100;

        a = [  -60.2409;   93.4517;   23.3585    ];
        b = [    0.023517; -0.023656;  0.0047036 ];

        f  = @(T,a) a(1)  + a(2) / T   + a(3) * log(T);
        df = @(T,a)       - a(2) / T^2 + a(3) / T;
        
        g  = @(T,b) b(1) + b(2) * T +     b(3) * T^2;
        dg = @(T,b)        b(2)     + 2 * b(3) * T; 
        
        pK0   = -(  f(TK100,a) +  g(TK100,b) * S ) / LOG10;
        pK0_T = -( df(TK100,a) + dg(TK100,b) * S ) * TK100_T / LOG10;
        pK0_S = -(                g(TK100,b)     ) / LOG10;
        
        gpK0 = [pK0_T, pK0_S, 0];
    end
    
    function [pKs,gpKs] = calc_pKs(T,S,Pbar,RT,RT_T)
    % calculate Ks (Dickson 1990a)----------------------------------------
        TK = T + 273.15; % convert to Kelvin
        IonS = ions(S);
        IonS_S = ions_S(S);
        sqrtIonS_S = 0.5*IonS_S/sqrt(IonS);

        a1 = [  -4276.1;  141.328; -23.093 ]; 
        a2 = [ -13856.0;  324.57;  -47.986 ];
        a3 = [  35474.0; -771.54;  114.723 ];
        a4 = [  -2698.0;    0.0;     0.0   ];
        a5 = [   1776.0;    0.0;     0.0   ];

        f  = @(T,a)  a(1) / T    + a(2) + a(3) * log(T);
        df = @(T,a) -a(1) / T^2         + a(3) / T;
        
        pKs = -( f(TK,a1) + ...
                 f(TK,a2) * sqrt(IonS) + ...
                 f(TK,a3) * IonS + ...
                 f(TK,a4) * IonS^(1.5) + ...
                 f(TK,a5) * IonS^2 ) / LOG10 + p(1 - 0.001005 * S );
        pKs_T = -( df(TK,a1) + ...
                   df(TK,a2) * sqrt(IonS) + ...
                   df(TK,a3) * IonS + ...
                   df(TK,a4) * IonS^(1.5) + ...
                   df(TK,a5) * IonS^2 ) / LOG10;
        pKs_S = -( f(TK,a2) * IonS_S * 0.5 / sqrt(IonS) + ...
                   f(TK,a3) * IonS_S + ...
                   f(TK,a4) * IonS_S * 1.5 * IonS^(0.5) + ...
                   f(TK,a5) * IonS_S * 2.0 * IonS ) / LOG10 - 0.001005 * dpdx(1 - 0.001005 * S );

        
        % pressure correction
        a = [ -18.03; 0.0466; 0.000316 ];
        b = [-4.53; 0.09 ];
        
        pKs = pKs + ppfac(T,Pbar,a,b);
        pKs_T = pKs_T + ppfac_T(T,Pbar,a,b);
        pKs_P = ppfac_P(T,Pbar,a,b);
     
        gpKs = [pKs_T, pKs_S, pKs_P];
    end
    
    function [pKf,gpKf] = calc_pKf(T,S,Pbar,RT,RT_T)
    % calculate Kf (Dickson 1979)----------------------------------
        
        TK = T + 273.15; % convert to Kelvin
        IonS = ions(S);
        IonS_S = ions_S(S);
        sqrtIonS_S = 0.5*IonS_S/sqrt(IonS);
        a = 1590.2; b = -12.641; c = 1.525;

        pKf = -( a/TK + b + c * sqrt(IonS) ) ./LOG10 + p(1 - 0.001005 * S) ;
        pKf_T = -(-a/(TK.^2) ) / LOG10; 
        pKf_S =  -c * sqrtIonS_S / LOG10 + dpdx(1 - 0.001005  *S) * (-0.001005) ;  
        
        % pressure correction
        a = [ -9.78; -0.009; -0.000942 ];
        b = [ -3.91; 0.054 ];
        
        pKf = pKf + ppfac(T,Pbar,a,b);
        pKf_T = pKf_T + ppfac_T(T,Pbar,a,b);
        pKf_P = ppfac_P(T,Pbar,a,b);
        
        gpKf = [pKf_T, pKf_S, pKf_P];
    end
    
    
     
    function [pKb,gpKb] = calc_pKb(T,S,Pbar,RT,RT_T)
    % calculate Kb (Dickson 1990)----------------------------------
        TK = T + 273.15; % convert to Kelvin

        a1 = [ -8966.9;    -2890.53;     -77.942;   1.728; -0.0996 ];
        a2 = [   148.0248;   137.1942;     1.62142; 0.0;    0.0    ];
        a3 = [   -24.4344;   -25.085;     -0.2474;  0.0;    0.0    ];
        a4 = [     0.0;        0.053105;   0.0;     0.0;    0.0    ];

        f  = @(S,a) a(1)      + a(2) * S^(0.5)  + a(3) * S +       a(4) * S^(1.5) +     a(5) * S^2;
        df = @(S,a)       0.5 * a(2) * S^(-0.5) + a(3)     + 1.5 * a(4) * S^(0.5) + 2 * a(5) * S;

        pKb   = -(  f(S,a1) / TK +  f(S,a2)  +  f(S,a3) * log(TK)  +  f(S,a4) * TK ) / LOG10;
        pKb_T = -( -f(S,a1) / TK^2           +  f(S,a3) / TK       +  f(S,a4)      ) / LOG10;
        pKb_S = -( df(S,a1) / TK + df(S,a2)  + df(S,a3) * log(TK)  + df(S,a4) * TK ) / LOG10;
        
        % pressure correction
        a = [ -29.48; 0.1622; -0.002608 ];
        b = [ -2.84; 0];
        
        pKb = pKb + ppfac(T,Pbar,a,b);
        pKb_T = pKb_T + ppfac_T(T,Pbar,a,b);
        pKb_P = ppfac_P(T,Pbar,a,b);
        
        gpKb = [pKb_T, pKb_S, pKb_P];
    end
    
    function [pKw,gpKw] = calc_pKw(T,S,Pbar,RT,RT_T)
    % calculate Kw (Millero 1995)--------------------------------
        TK = T + 273.15; % convert to Kelvin

        a1 = [ 148.9802; -13847.26; -23.6521 ];
        a2 = [  -5.977;     118.67;   1.0495 ];
        a3 = [  -0.01615;     0.0;    0.0    ];

        f  = @(T,a)   a(1) + a(2) / T    + a(3) * log(T);
        df = @(T,a)        - a(2) / T^2  + a(3) / T;

        pKw   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
        pKw_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
        pKw_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;


        % pressure correction
        a = [ -20.02; 0.1119; -0.001409];
        b = [ -5.13; 0.0794 ];

        pKw   = pKw + ppfac(T,Pbar,a,b);
        pKw_T = pKw_T + ppfac_T(T,Pbar,a,b);
        pKw_P = ppfac_P(T,Pbar,a,b);
        
        gpKw = [pKw_T,pKw_S,pKw_P];
    end
    
    function [pK1p,gpK1p] = calc_pK1p(T,S,Pbar,RT,RT_T)
    % calculate K1p (Yao and Millero 1995)--------------------------------
        TK = T + 273.15; % convert to Kelvin
        
        a1 = [ -4576.752;    115.54;    -18.453 ];
        a2 = [  -106.736;      0.69171;   0.0   ];
        a3 = [    -0.65643;   -0.01844;   0.0   ];

        f  = @(T,a)   a(1) / T   + a(2) + a(3) * log(T);
        df = @(T,a) - a(1) / T^2        + a(3) / T;

        pK1p   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
        pK1p_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
        pK1p_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;

        
        % pressure correction       
        a = [ -14.51; 0.1211; -0.000321 ];
        b = [  -2.67; 0.0427 ];

        pK1p = pK1p + ppfac(T,Pbar,a,b);
        pK1p_T = pK1p_T + ppfac_T(T,Pbar,a,b);
        pK1p_P = ppfac_P(T,Pbar,a,b);

        gpK1p = [pK1p_T,pK1p_S,pK1p_P];

    end
    
    function [pK2p,gpK2p] = calc_pK2p(T,S,Pbar,RT,RT_T)
    % calculate K2p (Yao and Millero 1995)--------------------------------
        TK = T + 273.15; % convert to Kelvin

        a1 = [-8814.715;   172.1033;  -27.927 ]; 
        a2 = [ -160.34;      1.3566;    0.0   ];
        a3 = [    0.37335;  -0.05778;   0.0   ];
        
        f  = @(T,a)   a(1) / T   + a(2) + a(3) * log(T);
        df = @(T,a) - a(1) / T^2        + a(3) / T;

        pK2p   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
        pK2p_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
        pK2p_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;

        % pressure correction
        a = [ -23.12; 0.1758; -0.002647 ];
        b = [ -5.15; 0.09 ]; 

        pK2p = pK2p + ppfac(T,Pbar,a,b);
        pK2p_T = pK2p_T + ppfac_T(T,Pbar,a,b);
        pK2p_P = ppfac_P(T,Pbar,a,b);
        
        gpK2p = [pK2p_T, pK2p_S, pK2p_P];
    end
    
    function [pK3p,gpK3p] = calc_pK3p(T,S,Pbar,RT,RT_T)
    % calculate K3p (Yao and Millero 1995)--------------------------------
        TK = T + 273.15; % convert to Kelvin
        
        a1 = [ -3070.75;   -18.126   ];
        a2 = [    17.27039;  2.81197 ];
        a3 = [   -44.99486; -0.09984 ];
                
        f  = @(T,a)   a(1) / T   + a(2);
        df = @(T,a) - a(1) / T^2;
        
        pK3p   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S) +  f(TK,a3) * S ) / LOG10;
        pK3p_T = -( df(TK,a1) + df(TK,a2) * sqrt(S) + df(TK,a3) * S ) / LOG10;
        pK3p_S = -(              f(TK,a2) * 0.5 / sqrt(S) + f(TK,a3) ) / LOG10;
                
        % pressure correction
        a = [ -26.57; 0.202; -0.003042 ];
        b = [ -4.08; 0.0714 ];

        pK3p = pK3p + ppfac(T,Pbar,a,b);
        pK3p_T = pK3p_T + ppfac_T(T,Pbar,a,b);
        pK3p_P = ppfac_P(T,Pbar,a,b);
        
        gpK3p = [pK3p_T, pK3p_S, pK3p_P]; 
    end

    
    function [pKsi,gpKsi] = calc_pKsi(T,S,Pbar,RT,RT_T)
    % calculate Ksi (Yao and Millero 1995)--------------------------------
        TK = T + 273.15; % convert to Kelvin
        IonS = ions(S);
        IonS_S = ions_S(S);
        sqrtIonS_S = 0.5 * IonS_S / sqrt(IonS);
        
        f  = @(T,a)   a(1) / T   + a(2) + a(3) * log(T);
        df = @(T,a) - a(1) / T^2        + a(3) / T;
        
        a1 = [ -8904.2;    117.4;   -19.334 ];
        a2 = [  -458.79;     3.5913;  0.0   ]; 
        a3 = [   188.74;    -1.5998;  0.0   ];
        a4 = [   -12.1652;   0.07871; 0.0   ];
        
        pKsi = -( f(TK,a1) + ...
                  f(TK,a2) * sqrt(IonS) + ...
                  f(TK,a3) * IonS + ...
                  f(TK,a4) * IonS^2  ) / LOG10 + p(1 - 0.001005 * S);

        pKsi_T = -( df(TK,a1) + ...
                    df(TK,a2) * sqrt(IonS) + ...
                    df(TK,a3) * IonS + ...
                    df(TK,a4) * IonS^2 ) / LOG10;
        
        pKsi_S = -( f(TK,a2) * sqrtIonS_S + ...
                    f(TK,a3) * IonS_S + ...
                    f(TK,a4) * 2 * IonS_S * IonS) / LOG10 ...
                 -0.001005 * dpdx(1 - 0.001005 * S);
        
        % pressure correction
        a = [ -29.48; 0.1622; -0.002608 ];
        b =[ -2.84; 0];
        
        pKsi = pKsi + ppfac(T,Pbar,a,b);
        pKsi_T = pKsi_T + ppfac_T(T,Pbar,a,b);
        pKsi_P = ppfac_P(T,Pbar,a,b);
        
        gpKsi = [pKsi_T, pKsi_S, pKsi_P]; % gradient pKsi
    end
    
    function [pK1,gpK1] = calc_pK1(T,S,Pbar,RT,RT_T)
    % calculate pK1 (Mehrbach refit by Dickson and Millero 1987)---
        TK = T + 273.15; % convert to Kelvin
        a = 3670.7; b = -62.008; c = 9.7944; d = -0.0118; g = 0.000116;
        pK1   =  a / TK   + b + c * log(TK) + d * S + g * S^2;
        pK1_T = -a / TK^2     + c / TK;
        pK1_S = d + 2.*g.*S;
        
        % pressure correction
        a = [-25.5; 0.1271; 0];
        b = [ -3.08; 0.0877 ];

        pK1   = pK1 + ppfac(T,Pbar,a,b);
        pK1_T = pK1_T + ppfac_T(T,Pbar,a,b);
        pK1_P = ppfac_P(T,Pbar,a,b);

        gpK1 = [pK1_T, pK1_S, pK1_P];

    end
    
    function [pK2,gpK2] = calc_pK2(T,S,Pbar,RT,RT_T)
    % calculate pK2 (Mehrbach refit by Dickson and Millero 1987)---
        TK = T + 273.15; % convert to Kelvin
        a = 1394.7; b = 4.777; c = -0.0184; d = 0.000118;
        pK2 = a / TK + b + c * S + d * S^2;
        pK2_T = -a / (TK^2);
        pK2_S = c + 2 * d * S;
        
        % pressure correction
        a = [ -15.82; -0.0219; 0 ];
        b = [ 1.13; -0.1475 ];
        
        pK2 = pK2 + ppfac(T,Pbar,a,b);
        pK2_T = pK2_T + ppfac_T(T,Pbar,a,b);
        pK2_P = ppfac_P(T,Pbar,a,b);
        
        gpK2 = [pK2_T, pK2_S, pK2_P]; 
    end
    
    function [pKnh4, gpKnh4] = calc_pKnh4(T,S,Pbar,RT,RT_T)
    % Ammonia, added by Sharp et al 2021, from Clegg and Whitfield (1995)
        TK = T + 273.15; % convert to Kelvin
        a = 9.44605; b = -2729.33; c = 1/298.15; d = 0.04203362; g = -11.24742; 

        f  = @(T,a) a(1)  + a(2) * sqrt(T) + a(3) * T + a(4) / T;
        df = @(T,a)    0.5 *a(2) / sqrt(T) + a(3)     - a(4) / T^2;

        a1 = [ -13.6416;      1.176949;     -0.02860785;   545.4834       ];
        a2 = [  -0.1462507;    0.0090226468; -0.0001471361;  10.5425      ];
        a3 = [   0.004669309; -0.0001691742;  0;             -0.5677934   ];
        a4 = [  -2.354039e-05; 0;             0;              0.009698623 ];

        pKnh4 = a + b * (c - 1 / TK) + (d + g / TK) * S^(0.25) + ...
                f(TK,a1) * S^0.5 + ...
                f(TK,a2) * S^1.5 + ...
                f(TK,a3) * S^2   + ...
                f(TK,a4) * S^2.5 + ...
                + p(1 - 0.001005 * S); 
        
        pKnh4_T = b / (TK^2) - g / (TK^2) * S^(0.25) + ...
                  df(TK,a1) * S^0.5 + ...
                  df(TK,a2) * S^1.5 + ...
                  df(TK,a3) * S^2   + ...
                  df(TK,a4) * S^2.5;

        pKnh4_S = 0.25 * (d + g / TK) * S^(-0.75) + ...
                  f(TK,a1) * 0.5 * S^-0.5 + ...
                  f(TK,a2) * 1.5 * S^0.5 + ...
                  f(TK,a3) * 2.0 * S   + ...
                  f(TK,a4) * 2.5 * S^1.5 - ...
                  0.001005 * dpdx(1 - 0.001005 * S); 
        
        % pressure correction
        a = [ -26.43; 0.0889; -0.000905 ];
        b = [ -5.03; 0.0814 ];

        pKnh4 = pKnh4 + ppfac(T,Pbar,a,b);
        pKnh4_T = pKnh4_T + ppfac_T(T,Pbar,a,b);
        pKnh4_P = ppfac_P(T,Pbar,a,b);
        
        gpKnh4 = [pKnh4_T, pKnh4_S, pKnh4_P];
    end
    
    function [pKh2s,gpKh2s] = calc_pKh2s(T,S,Pbar,RT,RT_T)
    % hydrogen sulfide, added by Sharp et al 2021, from Millero et al (1988)
        TK = T + 273.15; % convert to Kelvin
        a = 225.838; b = -13275.3; c = -34.6435;
        d = 0.3449; h = -0.0274;
        pKh2s = - (a + b./TK + c.*log(TK) + d.*sqrt(S) + ...
                   h.*S) ./ LOG10;
        pKh2s_T = - (-b./(TK.^2) + c./TK ) ./ LOG10;
        pKh2s_S = - (0.5.*d./sqrt(S) + h) ./ LOG10;
        
        % pressure correction
        a = [ -11.07; -0.009; -0.000942 ];
        b = [ -2.89;  0.054 ];
        
        pKh2s = pKh2s + ppfac(T,Pbar,a,b);
        pKh2s_T = pKh2s_T + ppfac_T(T,Pbar,a,b);
        pKh2s_P = ppfac_P(T,Pbar,a,b);
        
        gpKh2s = [pKh2s_T, pKh2s_S, pKh2s_P];
    end
    
    
    function [pSWS2tot,gpSWS2tot] = calc_pSWS2tot(S,pKs,gpKs,pKf,gpKf)
    % pH scale conversion factors (not pressure corrected)-----------
    % FREE2tot = 1 + TS./Ks; 

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
        
        top = 1 + q(pTS - pKs) ;
        top_T = dqdx(pTS - pKs) * (-gpKs(1)) ;
        top_S = dqdx(pTS - pKs) * (pTS_S - gpKs(2)) ;
        top_P = dqdx(pTS - pKs) * (-gpKs(3)) ;
        bot = top + q(pTF - pKf) ;
        bot_T = top_T + dqdx(pTF - pKf) * (-gpKf(1)) ;
        bot_S = top_S + dqdx(pTF - pKf) * (pTF_S - gpKf(2));
        bot_P = top_P + dqdx(pTF - pKf) * ( -gpKf(3));
        pSWS2tot = p(top./bot);
        pSWS2tot_T = (-top_T / (top*LOG10) ) + (bot_T / (bot*LOG10) );
        pSWS2tot_S = (-top_S / (top*LOG10) ) + (bot_S / (bot*LOG10) );
        pSWS2tot_P = (-top_P / (top*LOG10) ) + (bot_P / (bot*LOG10) );
        gpSWS2tot = [pSWS2tot_T, pSWS2tot_S, pSWS2tot_P];
    end
    
end

