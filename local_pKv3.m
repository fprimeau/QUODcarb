function [pK,gpK] = local_pKv3(obs,T,S,P)
% base equations COPIED FROM co2sys.m Orr et al. (2018)  Github
% Originally from  van Heuven et al. (2011)
% Original co2sys is from Lewis and Wallace (1998)

% INPUT:
%   T = Temp (deg C)
%   S  = Salinity
%   P  = pressure (dbar)

% OUTPUT:
%    pK  = [pK0;pK1;pK2;pKb;pKw;pKs;pKf;pK1p;pK2p;pK3p;pKsi;pKnh4;pKh2s;pp2f;pKar;pKca];
%   gpK  = [pK_T, pK_S, pK_P]; first derivatives (gradient of pK)

    TK   = T + 273.15; % convert to Kelvin
    Rgas = 83.14462618; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
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
    
    
    % corrections for pressure---------------------------------------------
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
              -( -dV(T,a)   * Pbar   + 0.5 * Ka(T,b)   * Pbar^2 ) * (RT_T / ( RT^2 * LOG10 ));
    ppfac_P = @(T,Pbar,a,b)... 
              -( -dV(T,a)         +       Ka(T,b)   * Pbar ) / (10 * RT * LOG10 ) ;
    

    % compute the pK's and their derivatives w.r.t. T,P,and S -------------
    [pp2f    , gpp2f    ] = calc_p2f(T,RT,RT_T); 
    [pKs     , gpKs     ] = calc_pKs(obs,T,S,Pbar); 
    [pKf     , gpKf     ] = calc_pKf(obs,T,S,Pbar); 
    [pSWS2tot, gpSWS2tot, pFREE2tot, gpFREE2tot] = calc_pSWS2tot(S,pKs,gpKs,pKf,gpKf);
    [phf     , gphf     ] = calc_phf(obs,T,S);
    [pK0     , gpK0     ] = calc_pK0(T,S); 
    [pKb     , gpKb     ] = calc_pKb(obs,T,S,Pbar,pSWS2tot,gpSWS2tot,phf,gphf); 
    [pKw     , gpKw     ] = calc_pKw(obs,T,S,Pbar); 
    [pK1p    , gpK1p    ] = calc_pK1p(obs,T,S,Pbar); 
    [pK2p    , gpK2p    ] = calc_pK2p(obs,T,S,Pbar,phf,gphf); 
    [pK3p    , gpK3p    ] = calc_pK3p(obs,T,S,Pbar,phf,gphf); 
    [pKsi    , gpKsi    ] = calc_pKsi(obs,T,S,Pbar,phf,gphf); 
    [pK1     , gpK1     ] = calc_pK1(obs,T,S,Pbar,phf,gphf,pSWS2tot,gpSWS2tot); 
    [pK2     , gpK2     ] = calc_pK2(obs,T,S,Pbar,phf,gphf,pSWS2tot,gpSWS2tot);
    [pKnh4   , gpKnh4   ] = calc_pKnh4(T,S,Pbar);
    [pKh2s   , gpKh2s   ] = calc_pKh2s(T,S,Pbar);
    [pKar    , gpKar    ] = calc_pKar(obs,T,S,Pbar,phf,gphf);
    [pKca    , gpKca    ] = calc_pKca(obs,T,S,Pbar,phf,gphf); %??

    % convert from SWS to pH total scale ----------------------------------
    pK1  = pK1  + pSWS2tot;  gpK1  = gpK1  + gpSWS2tot;
    pK2  = pK2  + pSWS2tot;  gpK2  = gpK2  + gpSWS2tot;
    pKw  = pKw  + pSWS2tot;  gpKw  = gpKw  + gpSWS2tot;
    pK1p = pK1p + pSWS2tot;  gpK1p = gpK1p + gpSWS2tot;
    pK2p = pK2p + pSWS2tot;  gpK2p = gpK2p + gpSWS2tot;
    pK3p = pK3p + pSWS2tot;  gpK3p = gpK3p + gpSWS2tot;
    pKsi = pKsi + pSWS2tot;  gpKsi = gpKsi + gpSWS2tot;
    pKar = pKar + pSWS2tot;  gpKar = gpKar + gpSWS2tot;
    pKca = pKca + pSWS2tot;  gpKca = gpKca + gpSWS2tot;

    %
    % output
    %
    pK = [pK0;    pK1;  pK2;  pKb;  pKw;  pKs;  pKf;  pK1p;  pK2p; ...
          pK3p;  pKsi;  pKnh4;  pKh2s;  pp2f; pKar; pKca];
    gpK = [gpK0; gpK1; gpK2; gpKb; gpKw; gpKs; gpKf; gpK1p; gpK2p; ...
           gpK3p; gpKsi; gpKnh4; gpKh2s; gpp2f; gpKar; gpKca];
    

    
    %
    % subfunctions
    %
 
    function [pp2f,gpp2f] = calc_p2f(T,RT,RT_T)
    % pCO2 to fCO2 conversion (Weiss 1974) valid to within 0.1% -----------
        TK = T + 273.15; % convert to Kelvin
        Pstd = 1.01325;
        delC = (57.7 - 0.118.*TK);
        delC_T = -0.118;

        a = [ -1636.75; 12.0408; -0.0327957; 3.16528e-5 ]; 
        
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
    % calculate K0 (Weiss 1974)--------------------------------------------
        TK = T + 273.15; % convert to Kelvin
        TK100 = TK./100;
        TK100_T = 1/100;

        a = [  -60.2409;   93.4517;   23.3585  ];
        b = [  0.023517; -0.023656;  0.0047036 ];

        f  = @(T,a) a(1)  + a(2) / T   + a(3) * log(T);
        df = @(T,a)       - a(2) / T^2 + a(3) / T;
        
        g  = @(T,b) b(1) + b(2) * T +     b(3) * T^2;
        dg = @(T,b)        b(2)     + 2 * b(3) * T; 
        
        pK0   = -(  f(TK100,a) +  g(TK100,b) * S ) / LOG10;
        pK0_T = -( df(TK100,a) + dg(TK100,b) * S ) * TK100_T / LOG10;
        pK0_S = -(                g(TK100,b)     ) / LOG10;
        
        gpK0 = [pK0_T, pK0_S, 0];
    end
    
    function [pKs,gpKs] = calc_pKs(obs,T,S,Pbar)
        TK = T + 273.15; % convert to Kelvin
        IonS = ions(S);
        IonS_S = ions_S(S);
        sqrtIonS_S = 0.5*IonS_S/sqrt(IonS);

        % all calculated on free pH scale
        % stay on free pH scale (no conversion to SWS or total)
        if obs.cKSO4 == 1
            % calculate Ks (Dickson 1990a)---------------------------------
            a1 = [  -4276.1;  141.328; -23.093 ]; 
            a2 = [ -13856.0;   324.57; -47.986 ];
            a3 = [  35474.0;  -771.54; 114.723 ];
            a4 = [  -2698.0;      0.0;    0.0  ];
            a5 = [   1776.0;      0.0;    0.0  ];

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
        elseif obs.cKSO4 == 2
            % calculate Ks (Khoo et al 1977)-------------------------------
            a1 = [   647.59;  -6.3451;  0.0 ]; 
            a2 = [ 0.019085;  -0.5208;  0.0 ];

            f  = @(T,a)  a(1) / T    + a(2) + a(3) * log(T);
            df = @(T,a) -a(1) / T^2         + a(3) / T;
        
            pKs = ( f(TK,a1) + ...
                    a2(1) * TK + ...
                    ( a2(2) * sqrt(IonS)) ) + ...
                    p(1 - 0.001005 * S );
            pKs_T = -( df(TK,a1) + ...
                    a2(1) );
            pKs_S = ( a2(2) * IonS_S * 0.5 / sqrt(IonS) ) - ...
                   0.001005 * dpdx(1 - 0.001005 * S );
        elseif obs.cKSO4 == 3
            % calculate Ks (Waters and Millero, 2013) ---------------------
            % log(K^(.)_HSO4) = logKS0
            % log((K^(*)_HSO4)/(K^(.)_HSO4)) = logKSK0
            % KS = (10^logKSK0))*(10^logKS0)
            % (^ from Sharp's code)

            a1 = [  0.0         ;    562.69486;  -102.5154 ]; 
            a2 = [ -0.0001117033;    0.2477538;  -13273.76 ];
            a3 = [  4.24666     ;    -0.152671;  0.0267059 ];
            a4 = [ -0.000042128 ;          0.0;        0.0 ];
            a5 = [  0.2542181   ;  -0.00509534; 0.00071589 ];
            a6 = [ -0.00291179  ; 0.0000209968;        0.0 ];
            a7 =   -0.0000403724;

            f1  = @(T,a)  a(1) / T    + a(2) + a(3) * log(T);
            df1 = @(T,a) -a(1) / T^2         + a(3) / T;
            f2  = @(T,a)  a(1) * (T^2) + a(2) * T + a(3) / T ;
            df2 = @(T,a) 2*a(1)*T + a(2) - a(3) / T^2 ;
            f3  = @(T,a)  a(1) + a(2) * T + a(3) * T * log(T) ;
            df3 = @(T,a)  a(2) + a(3)*log(T) + a(3) ;

            logKS0  =   f1(TK,a1) + f2(TK,a2) ;
            logKSK0 = ( f3(TK,a3) + f2(TK,a4) ) * sqrt(S) + ...
                      ( f3(TK,a5) ) * S + ...
                      ( f3(TK,a6) ) * S^(1.5) + ...
                      ( a7 * S^2 ) ;
            
            pKs = (-logKS0 - logKSK0) + p(1 - 0.001005 * S ); % our 'p' form is -log10
         
            logKS0_T  = df1(TK,a1) + df2(TK,a2) ;
            logKSK0_T = ( df3(TK,a3) + df2(Tk,a4) ) + ...
                          df3(TK,a5) + ...
                          df3(TK,a6) ;
            logKSK0_S = ( f3(TK,a3) + f2(Tk,a4) ) * (-0.5/sqrt(S)) + ...
                        ( f3(TK,a6) * 0.5 * sqrt(S)) + ...
                        ( a7 * 2 * S );

            pKs_T = (-logKS0_T - logKSK0_T) ;
            pKs_S = (-logKSK0_S) - 0.001005 * dpdx(1-0.001005*S);
        end
        
        % pressure correction
        a = [ -18.03; 0.0466; 0.000316 ];
        b = [-4.53; 0.09 ];
        
        pKs = pKs + ppfac(T,Pbar,a,b);
        pKs_T = pKs_T + ppfac_T(T,Pbar,a,b);
        pKs_P = ppfac_P(T,Pbar,a,b);
     
        gpKs = [pKs_T, pKs_S, pKs_P];
    end
    
    function [pKf,gpKf] = calc_pKf(obs,T,S,Pbar)
        TK = T + 273.15; % convert to Kelvin
        IonS = ions(S);
        IonS_S = ions_S(S);
        sqrtIonS_S = 0.5*IonS_S/sqrt(IonS);

        % all calculated on free pH scale
        % stay on free pH scale (no conversion to SWS or total)
        if obs.cKF == 1
        % calculate Kf (Dickson 1979)--------------------------------------   
            a = 1590.2; b = -12.641; c = 1.525;

            pKf = -( a/TK + b + c * sqrt(IonS) ) ./LOG10 + p(1 - 0.001005 * S) ; % free pH scale
            pKf_T = -(-a/(TK.^2) ) / LOG10; 
            pKf_S =  -c * sqrtIonS_S / LOG10 + dpdx(1 - 0.001005  *S) * (-0.001005) ;  
        elseif obs.cKF == 2
        % calculate Perez and Fraga (1987)---------------------------------
        % (to be used for S: 10-40, T: 9-33) (from Sharp's code)
            a = 874;  b = -9.68;  c = 0.111;

            pKf = -( a/TK + b + c * sqrt(S) ) ./LOG10 ; % free pH scale
            pKf_T = -(-a/(TK.^2) ) / LOG10; 
            pKf_S =  -c * (-0.5/sqrt(S)) / LOG10 ;  
        end

        % pressure correction
        a = [ -9.78; -0.009; -0.000942 ];
        b = [ -3.91; 0.054 ];
        
        pKf = pKf + ppfac(T,Pbar,a,b);
        pKf_T = pKf_T + ppfac_T(T,Pbar,a,b);
        pKf_P = ppfac_P(T,Pbar,a,b);
        
        gpKf = [pKf_T, pKf_S, pKf_P];
    end
    
        
    function [pKb,gpKb] = calc_pKb(obs,T,S,Pbar,pSWS2tot,gpSWS2tot,phf,gphf)   
        TK = T + 273.15; % convert to Kelvin

        if obs.cK1K2 == 6 || obs.cK1K2 == 7
        % GEOSECS Peng et al.
        % Lyman, John, UCLA Thesis, 1957
        % fit by Li et al, JGR 74:5507-5525, 1969
        % (copied from Orr's code)
            pKb = -( -9.26 + 0.00886 * S + 0.01*T) ; % our 'p' is -log10
            pKb_T = -0.01;
            pKb_S = -0.00886;

        % convert pKb from NBS scale to SWS scale, then SWS to total
            pKb = pKb - phf ;
            pKb = pKb + pSWS2tot;
            pKb_T = pKb_T - gphf(1) + gpSWS2tot(1);
            pKb_S = pKb_S - gphf(2) + gpSWS2tot(2);

        else
        % calculate Kb (Dickson 1990)--------------------------------------
            a1 = [ -8966.9;    -2890.53;     -77.942;   1.728; -0.0996 ];
            a2 = [   148.0248;   137.1942;     1.62142; 0.0;    0.0    ];
            a3 = [   -24.4344;   -25.085;     -0.2474;  0.0;    0.0    ];
            a4 = [     0.0;        0.053105;   0.0;     0.0;    0.0    ];

            f  = @(S,a) a(1)      + a(2) * S^(0.5)  + a(3) * S +       a(4) * S^(1.5) +     a(5) * S^2;
            df = @(S,a)       0.5 * a(2) * S^(-0.5) + a(3)     + 1.5 * a(4) * S^(0.5) + 2 * a(5) * S;

            pKb   = -(  f(S,a1) / TK +  f(S,a2)  +  f(S,a3) * log(TK)  +  f(S,a4) * TK ) / LOG10;
            pKb_T = -( -f(S,a1) / TK^2           +  f(S,a3) / TK       +  f(S,a4)      ) / LOG10;
            pKb_S = -( df(S,a1) / TK + df(S,a2)  + df(S,a3) * log(TK)  + df(S,a4) * TK ) / LOG10;
        end
        % pressure correction
        a = [ -29.48; 0.1622; -0.002608 ];
        b = [ -2.84; 0];
        
        pKb = pKb + ppfac(T,Pbar,a,b);
        pKb_T = pKb_T + ppfac_T(T,Pbar,a,b);
        pKb_P = ppfac_P(T,Pbar,a,b);
        
        gpKb = [pKb_T, pKb_S, pKb_P];
    end
    
    function [pKw,gpKw] = calc_pKw(obs,T,S,Pbar)
    TK = T + 273.15; % convert to Kelvin

        if obs.cK1K2 == 7
        % calculate Kw (Millero, 1979)-------------------------------------
            a1 = [ 148.9802; -13847.26; -23.6521 ];
            a2 = [ -79.2447;   3298.72;  12.0408 ];
            a3 = [ -0.019813;      0.0;    0.0   ];

            f  = @(T,a)   a(1) + a(2) / T    + a(3) * log(T);
            df = @(T,a)        - a(2) / T^2  + a(3) / T;

            pKw   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
            pKw_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
            pKw_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;

        elseif obs.cK1K2 == 8
        % calculate Kw (Millero 1979)--------------------------------------
        % refit data of Harned and Owen, 1958
            a1 = [ 148.9802; -13847.26; -23.6521 ];

            f  = @(T,a)   a(1) + a(2) / T    + a(3) * log(T);
            df = @(T,a)        - a(2) / T^2  + a(3) / T;

            pKw   = -(  f(TK,a1) ) / LOG10;
            pKw_T = -( df(TK,a1) ) / LOG10;
            pKw_S = -(   0   ) / LOG10;
        else
        % calculate Kw (Millero, 1995)-------------------------------------
            a1 = [ 148.9802; -13847.26; -23.6521 ];
            a2 = [  -5.977;     118.67;   1.0495 ];
            a3 = [  -0.01615;     0.0;    0.0    ];

            f  = @(T,a)   a(1) + a(2) / T    + a(3) * log(T);
            df = @(T,a)        - a(2) / T^2  + a(3) / T;

            pKw   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
            pKw_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
            pKw_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;
        end 

        % pressure correction
        a = [ -20.02; 0.1119; -0.001409];
        b = [ -5.13; 0.0794 ];

        pKw   = pKw + ppfac(T,Pbar,a,b);
        pKw_T = pKw_T + ppfac_T(T,Pbar,a,b);
        pKw_P = ppfac_P(T,Pbar,a,b);
        
        gpKw = [pKw_T,pKw_S,pKw_P];
    end
    
    function [pK1p,gpK1p] = calc_pK1p(obs,T,S,Pbar)
    TK = T + 273.15; % convert to Kelvin
        if obs.cK1K2 == 7
        % calculate K1p----------------------------------------------------
        % Peng et al don't include the contribution from this term,
        % but it is so small it doesn't contribute. It needs to be kept
        % so that the routines work ok. (From Orr's code)
            pK1p = p(0.02);
            pK1p_T = 0; 
            pK1p_S = 0;        
        else
        % calculate K1p (Yao and Millero 1995)-----------------------------       
            a1 = [ -4576.752;    115.54;  -18.453 ];
            a2 = [  -106.736;   0.69171;     0.0  ];
            a3 = [  -0.65643;  -0.01844;     0.0  ];

            f  = @(T,a)   a(1) / T   + a(2) + a(3) * log(T);
            df = @(T,a) - a(1) / T^2        + a(3) / T;

            pK1p   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
            pK1p_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
            pK1p_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;
        end
        
        % pressure correction       
        a = [ -14.51; 0.1211; -0.000321 ];
        b = [  -2.67; 0.0427 ];

        pK1p = pK1p + ppfac(T,Pbar,a,b);
        pK1p_T = pK1p_T + ppfac_T(T,Pbar,a,b);
        pK1p_P = ppfac_P(T,Pbar,a,b);

        gpK1p = [pK1p_T,pK1p_S,pK1p_P];
    end
    
    function [pK2p,gpK2p] = calc_pK2p(obs,T,S,Pbar,phf,gphf)
        TK = T + 273.15; % convert to Kelvin

        if obs.cK1K2 == 7
        % calculate K2p (Kester and Pytkowicz, 1967)-----------------------
            pK2p = -(-9.039 - 1450/TK) / LOG10  ; 
            pK2p = pK2p - phf ; % convert from NBS to SWS pH scale

            pK2p_T = -( 1450 / TK^2 ) / LOG10 - gphf(1) ;
            pK2p_S = -gphf(2) ;
            pK2p_P = -gphf(3) ;

        else
        % calculate K2p (Yao and Millero 1995)-----------------------------
            a1 = [ -8814.715;  172.1033; -27.927 ]; 
            a2 = [   -160.34;    1.3566;    0.0  ];
            a3 = [   0.37335;  -0.05778;    0.0  ];
        
            f  = @(T,a)   a(1) / T   + a(2) + a(3) * log(T);
            df = @(T,a) - a(1) / T^2        + a(3) / T;

            pK2p   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
            pK2p_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
            pK2p_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;
            pK2p_P = 0 ;
        end

        % pressure correction
        a = [ -23.12; 0.1758; -0.002647 ];
        b = [ -5.15; 0.09 ]; 

        pK2p = pK2p + ppfac(T,Pbar,a,b);
        pK2p_T = pK2p_T + ppfac_T(T,Pbar,a,b);
        pK2p_P = pK2p_P + ppfac_P(T,Pbar,a,b);
        
        gpK2p = [pK2p_T, pK2p_S, pK2p_P];
    end
    
    function [pK3p,gpK3p] = calc_pK3p(obs,T,S,Pbar,phf,gphf)
        TK = T + 273.15; % convert to Kelvin

        if obs.cK1K2 == 7
        % calculate K3p (Kester and Pytkowicz, 1967)-----------------------
            pK3p = -(4.466 - 7276/TK) / LOG10;
            pK3p = pK3p - phf ; % convert from NBS to SWS pH scale

            pK3p_T = -( 7276 / TK^2 ) / LOG10 - gphf(1) ;
            pK3p_S = -gphf(2) ;
            pK3p_P = -gphf(3) ;    
        else
        % calculate K3p (Yao and Millero 1995)-----------------------------
            a1 = [    -3070.75; -18.126  ];
            a2 = [    17.27039;  2.81197 ];
            a3 = [   -44.99486; -0.09984 ];
                
            f  = @(T,a)   a(1) / T   + a(2);
            df = @(T,a) - a(1) / T^2;
        
            pK3p   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S) +  f(TK,a3) * S ) / LOG10;
            pK3p_T = -( df(TK,a1) + df(TK,a2) * sqrt(S) + df(TK,a3) * S ) / LOG10;
            pK3p_S = -(              f(TK,a2) * 0.5 / sqrt(S) + f(TK,a3) ) / LOG10;
            pK3p_P = 0 ;
          
        end
        % pressure correction
        a = [ -26.57; 0.202; -0.003042 ];
        b = [ -4.08; 0.0714 ];

        pK3p = pK3p + ppfac(T,Pbar,a,b);
        pK3p_T = pK3p_T + ppfac_T(T,Pbar,a,b);
        pK3p_P = pK3p_P + ppfac_P(T,Pbar,a,b);
        
        gpK3p = [pK3p_T, pK3p_S, pK3p_P]; 
    end

    
    function [pKsi,gpKsi] = calc_pKsi(obs,T,S,Pbar,phf,gphf)
    TK = T + 273.15; % convert to Kelvin
    IonS = ions(S);
    IonS_S = ions_S(S);
    sqrtIonS_S = 0.5 * IonS_S / sqrt(IonS);
    
        if obs.cK1K2 == 7
        % calculate Ksi (Sillen, Martell, and Bjerrum, 1964)---------------
            pKsi = p(4e-10) - phf ; % also convert NBS scale to SWS pH scale
            pKsi_T = -gphf(1);
            pKsi_S = -gphf(2);

        else
        % calculate Ksi (Yao and Millero 1995)-----------------------------         
            a1 = [  -8904.2;    117.4; -19.334 ];
            a2 = [  -458.79;   3.5913;   0.0   ]; 
            a3 = [   188.74;  -1.5998;   0.0   ];
            a4 = [ -12.1652;  0.07871;   0.0   ];
                        
            f  = @(T,a)   a(1) / T   + a(2) + a(3) * log(T);
            df = @(T,a) - a(1) / T^2        + a(3) / T;
        
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
        end
        % pressure correction
        a = [ -29.48; 0.1622; -0.002608 ];
        b =[ -2.84; 0];
        
        pKsi = pKsi + ppfac(T,Pbar,a,b);
        pKsi_T = pKsi_T + ppfac_T(T,Pbar,a,b);
        pKsi_P = ppfac_P(T,Pbar,a,b);
        
        gpKsi = [pKsi_T, pKsi_S, pKsi_P]; % gradient pKsi
    end
    
    function [pK1,gpK1] = calc_pK1(obs,T,S,Pbar,phf,gphf,pSWS2tot,gpSWS2tot)
        % calculate pK1 based on user choice input with obs.cK1K2
        TK = T + 273.15; % convert to Kelvin

        if obs.cK1K2 == 1
        % Roy et al, Marine Chemistry, 44:249-267, 1993 -------------------
            a1 = [ 2.83655     ; -2307.1266; -1.5529413];
            a2 = [ -0.020760841; -4.0484   ;  0.0      ];
            a3 = [  0.08468345 ;    0.0    ;  0.0      ];
            a4 = [ -0.00654208 ;    0.0    ;  0.0      ];

            f  = @(T,a) a(1) +  a(2) / T   + a(3) * log(T);
            df = @(T,a)       - a(2) / T^2 + a(3) / T;
            
            pK1 = - (f(TK,a1) + f(TK,a2) .* sqrt(S) + ...
                f(TK,a3) .* S +  f(TK,a4) .* (S.^(3/2)) )  ... % total scale
                / LOG10 + p(1 - 0.001005 * S) ; % convert to mol/kg-SW
            pK1 = pK1 - pSWS2tot; % convert from pH total to SWS scale

            pK1_T = -( df(TK,a1) + df(TK,a2) .* sqrt(S) + ...
                df(TK,a3) .* S + df(TK,a4) .* (S.^(3/2)) ) / LOG10;
            pK1_T = pK1_T - gpSWS2tot(1); % convert tot to SWS scale
            pK1_S = -(  f(TK,a2) .* 0.5 ./ sqrt(S) +  f(TK,a3) + ...
                (3/2) .* f(TK,a4) .* sqrt(S) ) / LOG10 ...
                -0.001005 * dpdx(1 - 0.001005 * S) ;
            pK1_S = pK1_S - gpSWS2tot(2); % convert tot to SWS scale
            % pass all pK1 out of this function as SWS scale

        elseif obs.cK1K2 == 2 
        % Goyet and Poisson, Deep-Sea Research, 36(11):1635-1654, 1989 ----
            a = 812.27; b = 3.356; c = -0.00171; d = 0.000091;
            pK1 = a ./ TK + b + c .* S .* log(TK) + d .* (S.^2) ; % SWS scale
            pK1_T = - a ./ (TK.^2) + c .* S ./ TK ;
            pK1_S = c .* log(TK) + 2 .* d .* S ;

        elseif obs.cK1K2 == 3
        % Hansson refit by Dickson and Millero, 1987 ----------------------
            a = 851.4; b = 3.237; c = -0.0106; d = 0.000105;
            pK1 = a ./ TK + b + c .* S + d .* (S.^2) ; % SWS scale
            pK1_T = - a ./ (TK.^2) ;
            pK1_S = c + 2 .* d .* S ;

        elseif obs.cK1K2 == 4
        % Mehrbach refit by Dickson and Millero, 1987 ---------------------
            a = 3670.7; b = -62.008; c = 9.7944; d = -0.0118; g = 0.000116;
            pK1   =  a ./ TK   + b + c .* log(TK) + d .* S + g .* S^2;
            pK1_T = -a ./ TK^2     + c ./ TK;
            pK1_S = d + 2.*g.*S;

        elseif obs.cK1K2 == 5
        % Hansson and Mehrbach refit by Dickson and Millero, 1987 ---------
            a = 845; b = 3.248; c = -0.0098; d = 0.000087;
            pK1   = a ./ TK + b + c .* S + d .* (S.^2) ;
            pK1_T = -a ./ (TK.^2) ;
            pK1_S = c + 2 .* d .* S ;
        
        elseif obs.cK1K2 == 6 || obs.cK1K2 == 7
        % GEOSECS and Peng et al use K1K2 from Mehrbach et al, 1973 -------
        % I.e., these are the original Mehrbach dissociation constants ----
            a = -13.7201; b = 0.031334; c = 3235.76; 
            d = 1.3e-5; g = -0.1032;

            pK1 = a + b .* TK + c ./ TK + d .* S .* TK + g .* sqrt(S) ; % NBS scale
            pK1 = pK1 - phf ; % convert from NBS to SWS pH scale

            pK1_T = ( b - c ./ (TK.^2) + d .* S ) - gphf(1) ; % SWS scale
            pK1_S = ( d .* TK + 0.5 .* g ./ sqrt(S) ) - gphf(2) ; % SWS scale

        elseif obs.cK1K2 == 8
        % Millero 1979, pure water case -----------------------------------
            a = 290.9097; b = -14554.21; c = -45.0575;
            pK1 = - (a + b ./ TK + c .* log(TK) ) / LOG10;
            pK1_T = - ( -b ./ (TK.^2) + c ./ TK ) / LOG10;
            pK1_S = 0;

        elseif obs.cK1K2 == 9
        % Cai and Wang 1998, for estuarine use ----------------------------
            a = 3404.71; b = 0.032786; c = -14.8435; d = -0.071692; 
            g = 0.0021487; h = 200.1; l = 0.3220;
            f1 = h ./ TK + l;

            pK1 = a ./ TK + b .* TK + c + d .* f1 .* sqrt(S) + g .* S ;
            pK1 = pK1 - phf; % convert from NBS scale to SWS scale

            pK1_T = ( -a ./ (TK.^2) + b + d .* sqrt(S) .* (-h ./ (TK.^2)) )...
                - gphf(1) ; % Temp derivative convert to sws scale
            pK1_S = ( 0.5 .* d .* f1 ./ sqrt(S) + g ) - gphf(2) ;

        elseif obs.cK1K2 == 10
        % Leuker, Dickson, Keeling, 2000 ----------------------------------
            a = 3633.86; b = -61.2172; c = 9.6777;
            d = -0.011555; g = 0.0001152;

            pK1 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2);
            pK1 = pK1 - pSWS2tot; % convert from total scale to SWS scale

            pK1_T = ( -a ./ (TK.^2) + c ./ TK ) - gpSWS2tot(1);
            pK1_S = ( d + 2 .* g .* S ) - gpSWS2tot(2);

        elseif obs.cK1K2 == 11
        % Mojica Prieto and Millero, 2002 ---------------------------------
            a = -43.6977; b = -0.0129037; c = 1.364e-4;
            d = 2885.378; g = 7.045159;

            pK1 = a + b .* S + c .* (S.^2) + d ./ TK + g .* log(TK) ;

            pK1_T = -d ./ (TK.^2) + g ./ TK ;
            pK1_S = b + 2 .* c .* S ;

        elseif obs.cK1K2 == 12
        % Millero et al, 2002 ---------------------------------------------
            a = 6.359; b = -0.00664; c = -0.01322; d = 4.989e-5;
            pK1 = a + b .* S + c .* T + d .* (T.^2) ; % tempC
            pK1_T = c + 2 .* d .* T ;
            pK1_S = b ;

        elseif obs.cK1K2 == 13
        % Millero 2006 work from titrations, S = 1 - 50, T = 0 - 50  ------
        % Millero, Graham, Huang, Bustos-Serrano, Pierrot, 2006
            a1 = [-126.3404;  13.4191; 0.0331; -5.33e-5] ;
            a2 = [ 6320.813; -530.123; -6.103;    0.0  ] ;
            a3 = [19.568224; -2.06950;   0.0 ;    0.0  ] ;
            
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK1 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK1_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK1_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);

        elseif obs.cK1K2 == 14
        % Millero 2010, for estuarine use ---------------------------------
            a1 = [-126.34048;  13.4038; 0.03206; -5.242e-5] ; 
            a2 = [  6320.813; -530.659; -5.8210;    0.0   ] ;
            a3 = [ 19.568224;  -2.0664;   0.0  ;    0.0   ] ;
                        
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK1 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK1_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK1_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);

        elseif obs.cK1K2 == 15
        % Waters, Millero, and Woosley, 2014 ------------------------------
            a1 = [-126.34048;  13.409160; 0.031646; -5.1895e-5] ;
            a2 = [  6320.813;  -531.3642;   -5.713;    0.0    ] ;
            a3 = [ 19.568224; -2.0669166;   0.0   ;    0.0    ] ;

            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK1 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK1_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK1_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);

        elseif obs.cK1K2 == 16
        % Sulpis et al, 2020 ----------------------------------------------
            a = 8510.63; b = -172.4493; c = 26.32996; 
            d = -0.011555; g = 0.0001152; 

            pK1 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2) ;
            pK1 = pK1 - pSWS2tot; % convert from tot to SWS scale

            pK1_T = ( -a ./ (TK.^2) + c ./ TK ) - gpSWS2tot(1);
            pK1_S = ( d + 2 .* g .* S ) - gpSWS2tot(2);

        elseif obs.cK1K2 == 17
        % Schockman and Byrne, 2021 ---------------------------------------
            a1 = [-126.34048;  13.568513; 0.031645; -5.3834e-5] ;
            a2 = [  6320.813;  -539.2304;   -5.635;    0.0    ] ;
            a3 = [ 19.568224; -2.0901396;    0.0  ;    0.0    ] ;
                        
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK1 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK1 = pK1 - pSWS2tot ; % convert from pHtot to SWS
            pK1_T = ( -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ) ...
                - gpSWS2tot(1) ;
            pK1_S = ( df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK) ) ...
                - gpSWS2tot(2);

        end

        % pressure correction for obs.cK1K2 ~=6 & ~=7 & ~=8
        a = [-25.5; 0.1271; 0];
        b = [ -3.08; 0.0877 ];

        pK1   = pK1 + ppfac(T,Pbar,a,b);
        pK1_T = pK1_T + ppfac_T(T,Pbar,a,b);
        pK1_P = ppfac_P(T,Pbar,a,b);

        gpK1 = [pK1_T, pK1_S, pK1_P];
        
    end
    
    function [pK2,gpK2] = calc_pK2(obs,T,S,Pbar,phf,gphf,pSWS2tot,gpSWS2tot)
    % calculate pK2 based on user choice input with obs.cK1K2
        TK = T + 273.15; % convert to Kelvin

        if obs.cK1K2 == 1
        % Roy et al, Marine Chemistry, 44:249-267, 1993 ------------------
            a1 = [  -9.226508 ; -3351.6106 ; -0.2005743 ];
            a2 = [-0.106901773;  -23.9722  ;   0.0      ];
            a3 = [  0.1130822 ;    0.0     ;   0.0      ];
            a4 = [ -0.00846934;    0.0     ;   0.0      ];

            f  = @(T,a) a(1) +  a(2) / T   + a(3) * log(T);
            df = @(T,a)       - a(2) / T^2 + a(3) / T;
            
            pK2 = - (f(TK,a1) + f(TK,a2) .* sqrt(S) + ...
                f(TK,a3) .* S +  f(TK,a4) .* (S.^(3/2)) )  ... % total scale
                / LOG10 + p(1 - 0.001005 * S) ; % convert to mol/kg-SW
            pK2 = pK2 - pSWS2tot; % convert from pH total to SWS scale

            pK2_T = -( df(TK,a1) + df(TK,a2) .* sqrt(S) + ...
                df(TK,a3) .* S + df(TK,a4) .* (S.^(3/2)) ) / LOG10;
            pK2_T = pK2_T - gpSWS2tot(1); % convert tot to SWS scale
            pK2_S = -(  f(TK,a2) .* 0.5 ./ sqrt(S) +  f(TK,a3) + ...
                (3/2) .* f(TK,a4) .* sqrt(S) ) / LOG10 ...
                -0.001005 * dpdx(1 - 0.001005 * S) ;
            pK2_S = pK2_S - gpSWS2tot(2); % convert tot to SWS scale
            % pass all pK2 out of this function as SWS scale

        elseif obs.cK1K2 == 2 
        % Goyet and Poisson, Deep-Sea Research, 36(11):1635-1654, 1989 ----
            a = 1450.87; b = 4.604; c = -0.00385; d = 0.000182;
            pK2 = a ./ TK + b + c .* S .* log(TK) + d .* (S.^2) ; % SWS scale
            pK2_T = - a ./ (TK.^2) + c .* S ./ TK ;
            pK2_S = c .* log(TK) + 2 .* d .* S ;

        elseif obs.cK1K2 == 3
        % Hansson refit by Dickson and Millero, 1987 ----------------------
            a = -3885.4; b = 125.844; c = -18.141; 
            d = -0.0192; g = 0.000132;

            pK2 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2) ; % SWS scale

            pK2_T = - a ./ (TK.^2) + c ./ TK ;
            pK2_S = d + 2 .* g .* S ;

        elseif obs.cK1K2 == 4
        % Mehrbach refit by Dickson and Millero, 1987 ---------------------
            a = 1394.7; b = 4.777; c = -0.0184; d = 0.000118;
            pK2 = a ./ TK + b + c .* S + d .* (S.^2);
            pK2_T = -a ./ (TK.^2);
            pK2_S = c + 2 .* d .* S;

        elseif obs.cK1K2 == 5
        % Hansson and Mehrbach refit by Dickson and Millero, 1987 ---------
            a = 1377.3; b = 4.824; c = -0.0185; d = 0.000122;
            pK2   = a ./ TK + b + c .* S + d .* (S.^2) ;
            pK2_T = -a ./ (TK.^2) ;
            pK2_S = c + 2 .* d .* S ;
        
        elseif obs.cK1K2 == 6 || obs.cK1K2 == 7
        % GEOSECS and Peng et al use K1K2 from Mehrbach et al, 1973 -------
        % I.e., these are the original Mehrbach dissociation constants ----
            a1 = [ 5371.9654; -128375.28;  1.671221 ] ;
            a2 = [  0.22913 ;    2.136  ; -8.0944e-4] ;
            a3 = [  18.3802 ;  -5617.11 ;     0.0   ] ;
            b = -2194.3055;

            f  = @(T,a) a(1) +  a(2) ./ T   + a(3) .* T ;
            df = @(T,a)    - a(2) ./ (T.^2) + a(3) ./ T ;
            
            pK2 = ( f(TK,a1) + f(TK,a2) .* S + f(TK,a3) .* log10(S) + ...
                b .* log10(TK) ) - phf ; % convert from NBS to SWS scale

            pK2_T = ( df(TK,a1) + df(TK,a2) .* S + df(TK,a3) .* log10(S) + ...
                b ./ (TK .* log(10)) ) - gphf(1) ;
            pK2_S = ( f(TK,a2) + f(TK,a3) ./ (S .* log(10)) ) - gphf(2) ;

        elseif obs.cK1K2 == 8
        % Millero 1979, pure water case -----------------------------------
            a = 207.6548; b = -11843.79; c = -33.6485;
            pK2 = - (a + b ./ TK + c .* log(TK) ) / LOG10;
            pK2_T = - ( -b ./ (TK.^2) + c ./ TK ) / LOG10;
            pK2_S = 0;

        elseif obs.cK1K2 == 9
        % Cai and Wang 1998, for estuarine use ----------------------------
            a = 2902.39; b = 0.02379; c = -6.4980; d = -0.3191; 
            g = 0.0198; h = -129.24; l = 1.4381;
            f1 = h ./ TK + l;

            pK2 = a ./ TK + b .* TK + c + d .* f1 .* sqrt(S) + g .* S ;
            pK2 = pK2 - phf; % convert from NBS scale to SWS scale

            pK2_T = ( -a ./ (TK.^2) + b + d .* sqrt(S) .* (-h ./ (TK.^2)) )...
                - gphf(1) ; % Temp derivative convert to sws scale
            pK2_S = ( 0.5 .* d .* f1 ./ sqrt(S) + g ) - gphf(2) ;

        elseif obs.cK1K2 == 10
        % Leuker, Dickson, Keeling, 2000 ----------------------------------
            a = 471.78; b = 25.929; c = -3.16967;
            d = -0.01781; g = 0.0001122;

            pK2 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2);
            pK2 = pK2 - pSWS2tot; % convert from total scale to SWS scale

            pK2_T = ( -a ./ (TK.^2) + c ./ TK ) - gpSWS2tot(1);
            pK2_S = ( d + 2 .* g .* S ) - gpSWS2tot(2);

        elseif obs.cK1K2 == 11
        % Mojica Prieto and Millero, 2002 ---------------------------------
            a1 = [ -452.0940;  21263.61; 68.483143 ] ;
            a2 = [ 13.142162; -581.4428; -1.967035 ] ;
            a3 = [ -8.101e-4;  0.259601;    0.0    ] ;

            f = @(T,a) a(1) + a(2) ./ TK + a(3) .* log(TK) ;
            df = @(T,a) -a(2) ./ (TK.^2) + a(3) ./ TK ;

            pK2 = f(TK,a1) + f(TK,a2) .* S + f(TK,a3) .* (S.^2) ;

            pK2_T = df(TK,a1) + df(TK,a2) .* S + df(TK,a3) .* (S.^2) ;
            pK2_S = f(TK,a2) + 2 .* f(TK,a3) .* S ;

        elseif obs.cK1K2 == 12
        % Millero et al, 2002 ---------------------------------------------
            a = 9.867; b = -0.01314; c = -0.01904; d = 2.448e-5;
            pK2 = a + b .* S + c .* T + d .* (T.^2) ; % tempC
            pK2_T = c + 2 .* d .* T ;
            pK2_S = b ;

        elseif obs.cK1K2 == 13
        % Millero 2006 work from titrations, S = 1 - 50, T = 0 - 50  ------
        % Millero, Graham, Huang, Bustos-Serrano, Pierrot, 2006
            a1 = [ -90.1833;  21.0894;  0.1248; -3.687e-4 ] ;
            a2 = [ 5143.692; -772.483; -20.051;    0.0    ] ;
            a3 = [14.613358;  -3.3336;    0.0 ;    0.0    ] ;
            
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK2 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;

            pK2_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK2_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);

        elseif obs.cK1K2 == 14
        % Millero 2010, for estuarine use ---------------------------------
            a1 = [ -90.18333;  21.3728;  0.1218; -3.688e-4] ; 
            a2 = [  5143.692; -788.289; -19.189;    0.0   ] ;
            a3 = [ 14.613358;   -3.374;   0.0  ;    0.0   ] ;
                        
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK2 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK2_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK2_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);

        elseif obs.cK1K2 == 15
        % Waters, Millero, and Woosley, 2014 ------------------------------
            a1 = [  -90.1833;  21.225890; 0.12450870; -3.7243e-4] ;
            a2 = [  5143.692;  -779.3444;  -19.91739;    0.0    ] ;
            a3 = [ 14.613358; -3.3534679;     0.0   ;    0.0    ] ;

            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK2 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK2_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK2_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);

        elseif obs.cK1K2 == 16
        % Sulpis et al, 2020 ----------------------------------------------
            a = 4226.23; b = -59.4636; c = 9.60817; 
            d = -0.01781; g = 0.0001122; 

            pK2 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2) ;
            pK2 = pK2 - pSWS2tot; % convert from tot to SWS scale

            pK2_T = ( -a ./ (TK.^2) + c ./ TK ) - gpSWS2tot(1);
            pK2_S = ( d + 2 .* g .* S ) - gpSWS2tot(2);

        elseif obs.cK1K2 == 17
        % Schockman and Byrne, 2021 ---------------------------------------
            a = 116.8067; b = -3655.02; c = -16.45817; 
            d = 0.04523; g = -0.615; h = -0.0002799; l = 4.969;

            pK2 = a + b ./ TK + c .* log(TK) + d .* S + g .* sqrt(S) + ...
                h .* (S.^2) + l .* (S ./ TK) ;
            pK2 = pK2 - pSWS2tot ; % convert from pHtot to SWS

            pK2_T = ( -b ./ (TK.^2) + c ./ TK - l .* S ./ (TK.^2) ) ...
                - gpSWS2tot(1) ;
            pK2_S = ( d + 0.5 .* g ./ sqrt(S) + 2 .* h .* S + l ./ TK ) ...
                - gpSWS2tot(2);

        end
        
        % pressure correction for obs.cK1K2 ~=6 & ~=7 & ~=8
        a = [ -15.82; -0.0219; 0 ];
        b = [ 1.13; -0.1475 ];
        
        pK2 = pK2 + ppfac(T,Pbar,a,b);
        pK2_T = pK2_T + ppfac_T(T,Pbar,a,b);
        pK2_P = ppfac_P(T,Pbar,a,b);
        
        gpK2 = [pK2_T, pK2_S, pK2_P]; 
    end
    
    function [pKnh4, gpKnh4] = calc_pKnh4(T,S,Pbar)
    % Ammonia, added by Sharp et al 2021, from Clegg and Whitfield (1995)
        TK = T + 273.15; % convert to Kelvin
        a = 9.244605; b = -2729.33; c = 1/298.15; d = 0.04203362; g = -11.24742; 

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
    
    function [pKh2s,gpKh2s] = calc_pKh2s(T,S,Pbar)
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
    
    function [pKar, gpKar] = calc_pKar(obs,T,S,Pbar,phf,gphf)
    % Aragonite solubility
        TK = T + 273.15;
        Rgas = 83.14462618; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
        RT   = Rgas * TK;
        RT_T = Rgas;

        if obs.cK1K2 == 6 || obs.cK1K2 == 7
        % calculate Kar-Aragonite (Berner, 1976) --------------------------
            a = 0.0000001; b = -34.452;      c = -39.866; 
            d = 110.21;    g = 0.0000075752;

            Kar = 1.45 .* a .* ( b + c .* S^(1/3) + d .* log10(S) ...
                    + g .* (TK.^2) ) ;
            Kar_T = 1.45 .* a .* 2 .* g .* TK ; 
            Kar_S = 1.45 .* a .* ( (1/3) .* c .* S ^ (-2/3) + ...
                    d ./ (S .* log(10)) ) ;

            pKar = p(Kar); % - phf ; % 'p' it and convert to SWS pH scale            
            pKar_T = dpdx(Kar) .* Kar_T ; %- gphf(1) ; 
            pKar_S = dpdx(Kar) .* Kar_S ; %- gphf(2) ;

            % pressure correction
            pKar = pKar - ((33.3 - 0.22 .* T) .* Pbar ./ RT ) ./ LOG10 ; % T = tempC
            pKar_T = pKar_T - (Pbar .* Rgas .* ...
                    (59.8730 .* T - 9155.988) ./ (RT)^2 ) ./ LOG10;
            pKar_P = - ((33.3 - 0.22 .* T) ./ RT ) ./ LOG10;

        else
        % calculate Kar-Aragonite (Mucci, 1983) ---------------------------
            a1 = [ -171.945; -0.077993; 2903.293; 71.595] ;
            a2 = [-0.068393; 0.0017276;   88.135;   0.0 ] ;
            b = -0.10018; c = 0.0059415;

            f  = @(T,a) a(1) + a(2) .* T + a(3) ./ T + a(4) .* log10(T) ;
            df = @(T,a)      a(2) - a(3)./ (T.^2) + a(4) ./ (T .* log(10)); 

            log10Kar = f(TK,a1) + f(TK,a2) .* sqrt(S) + ...
                        b .* S + c .* S ^(3/2) ;
            pKar = -log10Kar ; % pK = -log10(K);
            %pKar = pKar - phf ; % convert from NBS to SWS pH scale

            pKar_T = - ( df(TK,a1) + df(TK,a2) .* sqrt(S) ) ;%- gphf(1);
            pKar_S = - ( 0.5 .* f(TK,a2) ./ sqrt(S) + ...
                        b + (3/2) .* c .* sqrt(S) ) ;%- gphf(2);

            % pressure correction
            d1 = [(-48.76 + 2.8); 0.5304; 0.0 ] ;
            d2 = [-11.76; 0.3692] ;

            pKar = pKar + ppfac(T,Pbar,d1,d2);
            pKar_T = pKar_T + ppfac_T(T,Pbar,d1,d2);
            pKar_P = ppfac_P(T,Pbar,d1,d2);
        end
        gpKar = [pKar_T; pKar_S; pKar_P];
    end

    function [pKca, gpKca] = calc_pKca(obs,T,S,Pbar,phf,gphf)
    % Calcite solubility
        TK = T + 273.15;
        Rgas = 83.14462618; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
        RT   = Rgas * TK;
        RT_T = Rgas;

        if obs.cK1K2 == 6 || obs.cK1K2 == 7
        % calculate Kca-Calcite (Berner, 1976) ----------------------------
            a = 0.0000001; b = -34.452;      c = -39.866; 
            d = 110.21;    g = 0.0000075752;

            Kca = a .* ( b + c .* S^(1/3) + d .* log10(S) ...
                    + g .* (TK.^2) ) ;
            Kca_T = a .* 2 .* g .* TK ; 
            Kca_S = a .* ( (1/3) .* c .* S ^ (-2/3) + ...
                    d ./ (S .* log(10)) ) ;

            pKca = p(Kca) ;%- phf ; % 'p' it and convert to SWS pH scale            
            pKca_T = dpdx(Kca) .* Kca_T ;% - gphf(1) ; 
            pKca_S = dpdx(Kca) .* Kca_S ; %- gphf(2) ;

            % pressure correction
            pKca = pKca - ((36 - 0.2 .* T) .* Pbar ./ RT ) ./ LOG10 ; % T = tempC
            pKca_T = pKca_T - (Pbar .* Rgas .* ...
                    (54.43 .* T - 9888.03) ./ (RT)^2 ) ./ LOG10;
            pKca_P = - ((36 - 0.2 .* T) ./ RT ) ./ LOG10;

        else
        % calculate Kca-Calcite (Mucci, 1983) -----------------------------
            a1 = [-171.9065; -0.077993; 2839.319; 71.595] ;
            a2 = [ -0.77712; 0.0028426;   178.34;   0.0 ] ;
            b = -0.07711; c = 0.0041249;

            f  = @(T,a) a(1) + a(2) .* T + a(3) ./ T + a(4) .* log10(T) ;
            df = @(T,a)      a(2) - a(3)./ (T.^2) + a(4) ./ (T .* log(10)); 

            log10Kca = f(TK,a1) + f(TK,a2) .* sqrt(S) + ...
                        b .* S + c .* S ^(3/2) ;
            pKca = -log10Kca ; % pK = -log10(K);
            %pKca = pKca - phf ; % convert from NBS to SWS pH scale

            pKca_T = - ( df(TK,a1) + df(TK,a2) .* sqrt(S) );% - gphf(1);
            pKca_S = - ( 0.5 .* f(TK,a2) ./ sqrt(S) + ...
                        b + (3/2) .* c .* sqrt(S) ) ;%- gphf(2);

            % pressure correction
            d1 = [-48.76; 0.5304; 0.0 ] ;
            d2 = [-11.76; 0.3692] ;

            pKca = pKca + ppfac(T,Pbar,d1,d2);
            pKca_T = pKca_T + ppfac_T(T,Pbar,d1,d2);
            pKca_P = ppfac_P(T,Pbar,d1,d2);
        end
        gpKca = [pKca_T; pKca_S; pKca_P];
    end
    
    function [pSWS2tot,gpSWS2tot,pFREE2tot,gpFREE2tot] = calc_pSWS2tot(S,pKs,gpKs,pKf,gpKf)
    % pH scale conversion factors (not pressure corrected)----------- 

    % calculate TF (Riley 1965)--------------------------------------------
        TF = (0.000067./18.998).*(S./1.80655); % mol/kg-SW
        TF_S = 1.95217e-6; % derivative wrt S
        pTF = p(TF);
        pTF_S = dpdx(TF).*TF_S;
        
    % calculate TS (Morris & Riley 1966)-----------------------------------
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

     % FREE2tot = 1 + TS./Ks; ---------------------------------------------
        pFREE2tot = p(top);
        pFREE2tot_T = (-top_T / (top*LOG10) ) ;
        pFREE2tot_S = (-top_S / (top*LOG10) ) ;
        pFREE2tot_P = (-top_P / (top*LOG10) ) ;
        gpFREE2tot = [pFREE2tot_T, pFREE2tot_S, pFREE2tot_P];
    end

    function [phf,gphf] = calc_phf(obs,T,S)
        TK = T + 273.15; % convert to Kelvin
        if obs.cK1K2 == 7
        % fH def'n: Peng et al, Tellus 39B: 439-458, 1987: ----------------
            fH   = ( 1.29 - 0.00204.*(TK) + ...
                (0.00046 - (0.00000148 * (TK))) * (S).^2 ) ;
            phf = p(fH);
            gphf_T = dpdx(phf) - 0.00204 - 0.00000148 * (S.^2)  ;
            gphf_S = dpdx(phf) * 2 * (0.00046 - (0.00000148 * (TK))) * S ;
        else
        % fH def'n: Takahashi et al, Ch 3 in GEOSECS v.3, 1982 ------------
            fH = 1.2948 - 0.002036*(TK) + ...
                (0.0004607 - (0.000001475 * (TK))) * (S).^2 ;
            phf = p(fH);
            gphf_T = dpdx(phf) - 0.002036 - 0.000001475 * (S.^2) ;
            gphf_S = dpdx(phf) * 2 * (0.0004607 - (0.000001475 * (TK))) * S ;
        end
        gphf = [gphf_T, gphf_S, 0];
    end  
%     LOG10 = log(10);
%     p = @(x) -log10(x);
%     q = @(x) 10.^(-x);  % inverse p, i.e., a backward p
%     dpdx = @(x) -1 / (x * LOG10);     % p'
%     dqdx = @(x) -LOG10 * 10.^( -x );  % q'
    
end

