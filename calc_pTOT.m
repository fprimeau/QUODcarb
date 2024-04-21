function [pT,gpT,ggpT,epT] = calc_pTOT(opt,S)
% base equations COPIED from co2sys.m Orr et al. (2018) Github
% Originally from van Heuven et al. (2011)
% Original co2sys is from Lewis and Wallace (1998)

% INPUT:
%   S   = Salinity

% OUTPUT:
%   pT  = [   pTB;   pTS;   pTF;   pTCa ]
%           -log10 values of totals (pT)
%  gpT  = [  gpTB;  gpTS;  gpTF;  gpTCa ]
%           first derivatives (gradient of pT wrt S) 
% ggpT  = [ ggpTB; ggpTS; ggpTF; ggpTCa ]
%           second derivatives (Hessian of pT)
%  epT  = [  epTB;  epTS;  epTF;  epTCa ]
%           precisions of pT (w of errors)

    % utility functions
    LOG10   = log(10);
    p       = @(x) -log10(x);
    q       = @(x) 10.^(-x);  % inverse p, i.e., a backward p
    dpdx    = @(x) -1 / (x * LOG10);        % p'
    d2pdx2  = @(x) 1 / (x^2 * LOG10);       % p''
    dqdx    = @(x) -LOG10 * 10.^( -x );     % q'
    d2qdx2  = @(x) LOG10^2 * 10.^( -x );    % q''
    my_abs  = @(x) sqrt(x*x);

    % compute the totals and their derivatives
    [pTB  , gpTB  , ggpTB  , epTB  ] = calc_pTB(opt,S);
    [pTS  , gpTS  , ggpTS  , epTS  ] = calc_pTS(opt,S);
    [pTF  , gpTF  , ggpTF  , epTF  ] = calc_pTF(opt,S);
    [pTCa , gpTCa , ggpTCa , epTCa ] = calc_pTCa(opt,S);

    % ---------------------------------------------------------------------
    % output
    % ---------------------------------------------------------------------
    pT   = [   pTB;   pTS;   pTF;   pTCa ];
    gpT  = [  gpTB;  gpTS;  gpTF;  gpTCa ];
    ggpT = [ ggpTB; ggpTS; ggpTF; ggpTCa ];
    epT  = [  epTB;  epTS;  epTF;  epTCa ];

    % ---------------------------------------------------------------------
    % subfunctions
    % ---------------------------------------------------------------------

    function [pTB,gpTB,ggpTB,epTB] = calc_pTB(opt,S)
        if (opt.TB == 1)
            % Uppstrom, L., Deep-Sea Research 21:161-162, 1974
            % ( copied from Orr's code )
            % TB = ( 0.000232/ 10.811) * (sal/1.80655)
            TB     = 0.0004157 * S / 35;
            pTB    =  p( TB );
            gpTB   = dpdx( TB )   *  0.0004157 / 35 ;
            ggpTB  = d2pdx2( TB ) * (0.0004157 / 35 ) ^ 2;

            % std 5e-6 on avg 2.32e-4 for (B mg kg^-1)/(Cl o/oo)
            TBu     = ( ( (2.32e-4 + 5e-6)/10.811) * S/1.80655 ); 
            TBl     = ( ( (2.32e-4 - 5e-6)/10.811) * S/1.80655 );
            eTB     = (TBu - TBl) /2 ;
            epTB    = my_abs( p(TB + eTB) - pTB ); % mol/kg
            
        elseif (opt.TB == 2)
            % Lee, Kim, Myrne, Millero, Feely, Yong-Ming Liu. 2010.
            % Geochemica Et Cosmochimica Acta 74 (6): 1801-1811.
            % ( copied from Sharp's code )
            % TB = (0.0002414/ 10.811) * (sal/1.80655)
            TB      = 0.0004326 * S / 35;
            pTB     = p(TB);
            gpTB    = dpdx( TB )   *   0.0004326 / 35;
            ggpTB   = d2pdx2( TB ) * ( 0.0004326 / 35 ) ^ 2;
            
            % std 9e-7 on avg 2.414e-4
            TBu     = ( ( (2.414e-4 + 9e-7)/10.811) * S/1.80655);
            TBl     = ( ( (2.414e-4 - 9e-7)/10.811) * S/1.80655);
            eTB     = (TBu - TBl) /2;
            epTB    = my_abs( p(TB + eTB) - pTB ); % mol/kg
            
        elseif (opt.K1K2 == 6) || (opt.K1K2 == 7)
            % this is about 1% lower than Uppstrom's value
            % Culkin, F., in Chemical Oceanography,
            % ed. Riley and Skirrow, 1965: GEOSECS references this
            % (copied from Orr's code)
            TB      = 0.0004106 * S / 35;
            pTB     = p( TB );
            gpTB    = dpdx( TB )      *    0.0004106 / 35 ;
            ggpTB   = sys.d2pdx2( TB ) * ( 0.0004106 / 35 ) ^ 2;

            % can't find paper, assume same as Uppstrom
            % std 5e-6 on avg 2.32e-4
            TBu     = ( ( (2.32e-4 + 5e-6)/10.811) * S/1.80655 );
            TBl     = ( ( (2.32e-4 - 5e-6)/10.811) * S/1.80655 );
            eTB     = (TBu - TBl) /2 ;
            epTB    = my_abs( p(TB + eTB) - pTB ); % mol/kg
        end 
    end

    function [pTS,gpTS,ggpTS,epTS] = calc_pTS(opt,S)
        % Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
        % copied from Orr's code
        TS      = ( 0.14 / 96.062 ) * ( S / 1.80655 );
        pTS     = p( TS );
        gpTS    = dpdx( TS )   *   (0.14 / 96.062 ) / 1.80655 ;
        ggpTS   = d2pdx2( TS ) * ( (0.14 / 96.062 ) / 1.80655 ) ^ 2 ;
        
        % 0.14000 ± 0.00023
        TSu     = ( ( (0.14+0.00023)/96.062 ) * S/ 1.80655 );
        TSl     = ( ( (0.14-0.00023)/96.062 ) * S/ 1.80655 );
        eTS     = (TSu - TSl) / 2;
        my_abs  = @(x) sqrt(x*x);
        epTS    = my_abs( p(TS + eTS) - pTS );
    end

    function [pTF,gpTF,ggpTF,epTF] = calc_pTF(opt,S)
        % Riley, J. P., Deep-Sea Research 12:219-220, 1965:
        % this is .000068.*Sali./35. = .00000195.*Sali   
        TF      = ( 0.000067 / 18.998 ) * ( S / 1.80655 ); 
        pTF     = p( TF );
        gpTF    = dpdx( TF )   *   (0.000067 / 18.998 ) / 1.80655 ;
        ggpTF   = d2pdx2( TF ) * ( (0.000067 / 18.998 ) / 1.80655 ) ^ 2 ;

        % 6.7 ± 0.1 e-5
        TFu     = ( ( (6.7e-5 + 0.1e-5)/18.998) * S/1.80655 );
        TFl     = ( ( (6.7e-5 - 0.1e-5)/18.998) * S/1.80655 );
        eTF     = (TFu - TFl) / 2;
        my_abs  = @(x) sqrt(x*x);
        epTF    = my_abs( p(TF + eTF) - pTF );
    end

    function [pTCa,gpTCa,ggpTCa,epTCa] = calc_pTCa(opt,S)
        if opt.K1K2 == 6 || opt.K1K2 == 7
            % Calculate Ca for GEOSECS, Riley and Skirrow 1965
            TCa     =  ( 0.01026 * S / 35) ;
            pTCa    = p( TCa );
            gpTCa   = dpdx( TCa )   *   0.01026 / 35 ;
            ggpTCa  = d2pdx2( TCa ) * ( 0.01026 / 35 ) ^ 2 ;
        else
            % Calculate Ca, Riley and Tongudai 1967
            % this is 0.010285.*obs.sal./35;
            TCa     = ( 0.02128 / 40.087 * ( S / 1.80655 ) );
            pTCa    = p( TCa ); 
            gpTCa   = dpdx( TCa )   *   ( 0.02128 / 40.087 ) / 1.80655 ; 
            ggpTCa  = d2pdx2( TCa ) * ( ( 0.02128 / 40.087 ) / 1.80655 ) ^ 2 ; 
        end
        % mean 0.02128 ± 0.00006 Ca/Cl ratio (g/kg)/(o/oo)
        TCau    = ( (0.02128 + 6e-5)/ 40.087 * ( S / 1.80655 ) );
        TCal    = ( (0.02128 - 6e-5)/ 40.087 * ( S / 1.80655 ) );
        eTCa    = (TCau - TCal) / 2;
        my_abs  = @(x) sqrt(x*x);
        epTCa   = my_abs( p(TCa + eTCa) - pTCa );
    end

end



