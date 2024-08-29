function [out] = run_co2sys(obs,opt,tpopt,pair,nD)
    opt = check_opt(opt);
    tic
    for i = 1:nD
        par1type    = 1; % alkalinity
        par1        = obs(i).TA; % TA umol/kg
        epar1       = obs(i).eTA; % error alkalinity
        par2type    = 2; % DIC
        par2        = obs(i).TC; % DIC_in µmol/kg
        epar2       = obs(i).eTC; % error DIC
        par3type    = 3; % pH
        par3        = obs(i).tp(tpopt).ph;
        epar3       = obs(i).tp(tpopt).eph;
        par4        = obs(i).tp(tpopt).pco2;  % pCO2 converted to pH's temp
        epar4       = obs(i).tp(tpopt).epco2; %
        par4type    = 4;
        par5        = obs(i).tp(tpopt).co3; % CO3
        epar5       = obs(i).tp(tpopt).eco3;
        par5type    = 7;

        sal     = obs(i).sal; % salinity of sample
        tempin  = obs(i).tp(tpopt).T; % temp of sample
        presin  = obs(i).tp(tpopt).P; % pressure of sample (dbars)

        % tempout = obs(i).tp(tpopt).T;
        % presout = obs(i).tp(tpopt).P;
        % tempout = obs(i).tp(2).T; % FOR pH ZSCORE STUFF
        % presout = obs(i).tp(2).P; % FOR pH ZSCORE STUFF
        tempout = obs(i).tp(3).T; % FOR pCO2 ZSCORE STUFF
        presout = obs(i).tp(3).P; % FOR pCO2 ZSCORE STUFF

        sil     = obs(i).TSi; % total Si of sample
        po4     = obs(i).TP; % phosphate of sample
        pHscale = opt.phscale; % 1 = tot; 2 = sws; 3 = free; 4 = nbs
        k1k2c   = opt.K1K2;
        kso4c   = opt.KSO4; % 1 or 3 is Dickson 1990
        r       = 0; % correlation coefficient
        nh4     = 0;
        h2s     = 0;
        kfc     = opt.KF;
        tbc     = opt.TB;
        X       = opt.co2press; % calls for 'co2press,1' not just '1'

        % calculate some things
        p       = @(x) -log10(x);
        q       = @(x) 10.^(-x);  % inverse p, i.e., a backward p
        ebar    = @(px,pe) (0.5 * ( q( px - pe ) - q( px + pe ) ) );

        [pT,~,~,epT] = calc_pTOT(opt,sal);
        TB = q(pT(1)); 
        eTB = ebar(pT(1),epT(1));
        % Cal = q(pT(4));
        eTCa = ebar(pT(4),epT(4));

        [~,~,epKi] = calc_pK(opt,tempin,sal,presin);

        % errors
        esal    = obs(i).esal; % salinity error
        etemp   = obs(i).tp(tpopt).eT; % error temp
        epo4    = obs(i).eTP; % error PO4
        esi     = obs(i).eTSi; % error Si
        eBt     = eTB/TB; % wants it in frac form
        enh4    = 5e-4;
        eh2s    = 5e-4;
        ecal    = eTCa; % 60 umol/kg = 6e-5 mol/kg = 0.06 mmol/kg
        epK0    = epKi(1);
        epK1    = epKi(2);
        epK2    = epKi(3);
        epKb    = epKi(4);
        epKw    = epKi(5);
        epKar   = epKi(15);
        epKca   = epKi(16);
        epK     = [epK0, epK1, epK2, epKb, epKw, epKar, epKca];
        %    = [0.02,0.0075,0.015,0.01, 0.01, 0.02, 0.02];

        % par1 = alk, par2 = dic, par3 = ph, par4 = pco2; par5 = co3
        if pair == 1 % TA TC
            OUT = CO2SYS(par1,par2,par1type,par2type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,pHscale, ...
                k1k2c,kso4c,kfc,tbc,'co2_press',X);
            err = errors(par1,par2,par1type,par2type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,epar1,epar2, ...
                esal,etemp,esi,epo4,enh4,eh2s,epK,eBt,r,pHscale, ...
                k1k2c,kso4c,kfc,tbc,ecal);
        elseif pair == 2 % TC pH
            OUT = CO2SYS(par2,par3,par2type,par3type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,pHscale, ...
                k1k2c,kso4c,kfc,tbc,'co2_press',X);
            err = errors(par2,par3,par2type,par3type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,epar2,epar3, ...
                esal,etemp,esi,epo4,enh4,eh2s,epK,eBt,r,pHscale, ...
                k1k2c,kso4c,kfc,tbc,ecal);
        elseif pair == 3 % TA pH
            OUT = CO2SYS(par1,par3,par1type,par3type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,pHscale, ...
                k1k2c,kso4c,kfc,tbc,'co2_press',X);
            err = errors(par1,par3,par1type,par3type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,epar1,epar3, ...
                esal,etemp,esi,epo4,enh4,eh2s,epK,eBt,r,pHscale, ...
                k1k2c,kso4c,kfc,tbc,ecal);
        elseif pair == 6 % TC CO3
            OUT = CO2SYS(par2,par5,par2type,par5type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,pHscale, ...
                k1k2c,kso4c,kfc,tbc,'co2_press',X);
            err = errors(par2,par5,par2type,par5type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,epar2,epar5, ...
                esal,etemp,esi,epo4,enh4,eh2s,epK,eBt,r,pHscale, ...
                k1k2c,kso4c,kfc,tbc,ecal);
        elseif pair == 7 % TC pCO2
            OUT = CO2SYS(par2,par4,par2type,par4type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,pHscale, ...
                k1k2c,kso4c,kfc,tbc,'co2_press',X);
            err = errors(par2,par4,par2type,par4type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,epar2,epar4, ...
                esal,etemp,esi,epo4,enh4,eh2s,epK,eBt,r,pHscale, ...
                k1k2c,kso4c,kfc,tbc,ecal);
        elseif pair == 8 % TA pCO2
            OUT = CO2SYS(par1,par4,par1type,par4type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,pHscale, ...
                k1k2c,kso4c,kfc,tbc,'co2_press',X);
            err = errors(par1,par4,par1type,par4type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,epar1,epar4, ...
                esal,etemp,esi,epo4,enh4,eh2s,epK,eBt,r,pHscale, ...
                k1k2c,kso4c,kfc,tbc,ecal);
        elseif pair == 9 % TA CO3
            OUT = CO2SYS(par1,par5,par1type,par5type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,pHscale, ...
                k1k2c,kso4c,kfc,tbc,'co2_press',X);
            err = errors(par1,par5,par1type,par5type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,epar1,epar5, ...
                esal,etemp,esi,epo4,enh4,eh2s,epK,eBt,r,pHscale, ...
                k1k2c,kso4c,kfc,tbc,ecal);
        elseif pair == 10 % pH CO3
            OUT = CO2SYS(par3,par5,par3type,par5type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,pHscale, ...
                k1k2c,kso4c,kfc,tbc,'co2_press',X);
            err = errors(par3,par5,par3type,par5type,sal,tempin, ...
                tempout,presin,presout,sil,po4,nh4,h2s,epar3,epar5, ...
                esal,etemp,esi,epo4,enh4,eh2s,epK,eBt,r,pHscale, ...
                k1k2c,kso4c,kfc,tbc,ecal);
        end

        % exit struct
        out(i).TA      = OUT(1);
        out(i).eTA     = err(1); % sigma, not precision
        out(i).TC      = OUT(2);
        out(i).eTC     = err(2);
        out(i).TB      = OUT(93);
        out(i).TF      = OUT(94);
        out(i).TS      = OUT(95);
        out(i).TP      = OUT(96);
        out(i).TSi     = OUT(97);
        out(i).TNH4    = OUT(98);
        out(i).TH2S    = OUT(99);
        out(i).CAL     = OUT(100); % TCa
        out(i).ph      = OUT(21);
        out(i).eh      = (err(13)*1e-9); % was nano units
        out(i).eph     = (out(i).eh/(10^(-out(i).ph)))*(1/log(10)); % dph = (-1/ln(10))*(dH/H) from Orr's 'errors'
        out(i).ph_sws  = OUT(44);
        out(i).ph_free = OUT(45);
        out(i).ph_nbs  = OUT(46);
        out(i).ph_tot  = OUT(43);
        out(i).pfH     = -log10(OUT(101));
        out(i).pco2    = OUT(22);
        out(i).epco2   = err(14);
        out(i).co3     = OUT(25);
        out(i).eco3    = err(17);
        out(i).pK0     = p(OUT(78));
        out(i).epK0    = epK0; % not an output of 'errors'
        out(i).pK1     = p(OUT(79));
        out(i).epK1    = epK1; % not an output of 'errors'
        out(i).pK2     = p(OUT(80));
        out(i).epK2    = epK2; % not an output of 'errors'
        out(i).pKw     = p(OUT(83));
        out(i).epKw    = epKw; % not an output of 'errors'
        out(i).pKb     = p(OUT(84));
        out(i).epKb    = epKb; % not an output of 'errors'
        out(i).pKs     = p(OUT(86));
        out(i).pKf     = p(OUT(85));
        out(i).pKp1    = p(OUT(87));
        out(i).pKp2    = p(OUT(88));
        out(i).pKp3    = p(OUT(89));
        out(i).pKSi    = p(OUT(90));
        out(i).pKnh4   = p(OUT(91));
        out(i).pKh2s   = p(OUT(92));
        out(i).OmegaAr     = OUT(36);
        out(i).eOmegaAr    = err(21);
        out(i).OmegaCa     = OUT(35);
        out(i).eOmegaCa    = err(20);
        out(i).Rev         = OUT(34);
    end
    toc
end



% CO2SYSv3 ----------------------------------------------------------------
%  --> 01/30/23 MF edited it slightly to output CAL so output
%       vector is slightly different than listed in documentation
            % '100 - TCal            (umol/kgSW) ';
            % '101 - fH              (umol/kgSW) '};


function [DATA,HEADERS,NICEHEADERS]=CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE, ...
    SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,pHSCALEIN, ...
    K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON,varargin)
%**************************************************************************
%
% Current: CO2SYS.m v3.1.2   (Jan  2023: https://github.com/jonathansharp/CO2-System-Extd)
%          CO2SYS.m v2       (Dec  2016: https://github.com/jamesorr/CO2SYS-MATLAB)
%          CO2SYS.m v1       (Sept 2011: https://cdiac.ess-dive.lbl.gov/ftp/co2sys/CO2SYS_calc_MATLAB_v1.1/)
%
% CO2SYS is a MATLAB-version of the original CO2SYS for DOS. 
% CO2SYS calculates and returns the state of the carbonate system of 
%    oceanographic water samples, if supplied with enough input.
%
% Please note that, besides certain extended capabilities and minor
%    corrections, this software is intended to be nearly identical to the
%    DOS and Excel versions that have been released previously, meaning
%    that results obtained should be very nearly identical for identical
%    input.
% Additionally, several of the dissociation constants K1 and K2 that have 
%    been published since the original DOS version was written are implemented.
%    For a complete list of changes since version 1.0, see below.
%
% For much more info please have a look at:
%    Lewis, E., and D. W. R. Wallace. 1998. Program Developed for
%    CO2 System Calculations. ORNL/CDIAC-105. Carbon Dioxide Information
%    Analysis Center, Oak Ridge National Laboratory, U.S. Department of Energy,
%    Oak Ridge, Tennessee. http://cdiac.ornl.gov/oceans/co2rprt.html
%
%**************************************************************************
%
%  **** SYNTAX:
%  [RESULT,HEADERS,NICEHEADERS]=CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,...
%        ...SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,pHSCALEIN,...
%        ...K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON,varargin)
% 
%  **** SYNTAX EXAMPLES:
%  [Result]                     = CO2SYS(2400,2200,1,2,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [Result,Headers]             = CO2SYS(2400,   8,1,3,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [Result,Headers,Niceheaders] = CO2SYS( 500,   8,5,3,35,0,25,4200,0,15,1,0,0,1,4,1,1,1,'co2_press',1)
%  [A]                          = CO2SYS(2400,2000:10:2400,1,2,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [A]                          = CO2SYS(2400,2200,1,2,0:1:35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [A]                          = CO2SYS(2400,2200,1,2,35,0,25,0:100:4200,0,15,1,0,0,1,4,1,1,1,'co2_press',1)
%  
%  **** APPLICATION EXAMPLE (copy and paste this into command window):
%  tmps=0:40; sals=0:40; [X,Y]=meshgrid(tmps,sals);
%  A = CO2SYS(2300,2100,1,2,Y(:),X(:),nan,0,nan,1,1,0,0,1,9,1,1,1);
%  Z=nan(size(X)); Z(:)=A(:,4); figure; contourf(X,Y,Z,20); caxis([0 1200]); colorbar;
%  ylabel('Salinity [psu]'); xlabel('Temperature [degC]'); title('Dependence of pCO2 [uatm] on T and S')
% 
%**************************************************************************
%
% INPUT:
%
%   PAR1  (some unit) : scalar or vector of size n
%   PAR2  (some unit) : scalar or vector of size n
%   PAR1TYPE       () : scalar or vector of size n (Footnote *1 below)
%   PAR2TYPE       () : scalar or vector of size n (Footnote *1 below)
%   SAL            () : scalar or vector of size n
%   TEMPIN  (degr. C) : scalar or vector of size n 
%   TEMPOUT (degr. C) : scalar or vector of size n 
%   PRESIN     (dbar) : scalar or vector of size n 
%   PRESOUT    (dbar) : scalar or vector of size n
%   SI    (umol/kgSW) : scalar or vector of size n
%   PO4   (umol/kgSW) : scalar or vector of size n
%   NH4   (umol/kgSW) : scalar or vector of size n
%   H2S   (umol/kgSW) : scalar or vector of size n
%   pHSCALEIN         : scalar or vector of size n (Footnote *2 below)
%   K1K2CONSTANTS     : scalar or vector of size n (Footnote *3 below)
%   KSO4CONSTANT      : scalar or vector of size n (Footnote *4 below)
%   KFCONSTANT        : scalar or vector of size n (Footnote *5 below)
%   BORON             : scalar or vector of size n (Footnote *6 below)
%
% OPTIONAL INPUT:
%   'co2_press',X     : string,scalar pair
%               X = 0 : K0 and FugFac are not corrected for in situ pressure
%               X = 1 : K0 and FugFac are corrected for in situ pressure
%
%  (*1) Each element must be an integer, 
%      indicating that PAR1 (or PAR2) is of type: 
%  1 = TA
%  2 = DIC
%  3 = pH
%  4 = pCO2
%  5 = fCO2
%  6 = HCO3
%  7 = CO3
%  8 = CO2
% 
%  (*2) Each element must be an integer, 
%       indicating that the pH-input (PAR1 or PAR2, if any) is at:
%  1 = Total scale
%  2 = Seawater scale
%  3 = Free scale
%  4 = NBS scale
% 
%  (*3) Each element must be an integer, 
%       indicating the K1 and K2 dissociation constants that are to be used:
%   1 = Roy, 1993                                           T:    0-45  S:  5-45. Total scale. Artificial seawater.
%   2 = Goyet & Poisson                                     T:   -1-40  S: 10-50. Seaw. scale. Artificial seawater.
%   3 = HANSSON              refit BY DICKSON AND MILLERO   T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   4 = MEHRBACH             refit BY DICKSON AND MILLERO   T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   5 = HANSSON and MEHRBACH refit BY DICKSON AND MILLERO   T:    2-35  S: 20-40. Seaw. scale. Artificial seawater.
%   6 = GEOSECS (i.e., original Mehrbach)                   T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   7 = Peng	(i.e., originam Mehrbach but without XXX)   T:    2-35  S: 19-43. NBS scale.   Real seawater.
%   8 = Millero, 1979, FOR PURE WATER ONLY (i.e., Sal=0)    T:    0-50  S:     0. 
%   9 = Cai and Wang, 1998                                  T:    2-35  S:  0-49. NBS scale.   Real and artificial seawater.
%  10 = Lueker et al, 2000                                  T:    2-35  S: 19-43. Total scale. Real seawater.
%  11 = Mojica Prieto and Millero, 2002                     T:    0-45  S:  5-42. Seaw. scale. Real seawater
%  12 = Millero et al, 2002                                 T: -1.6-35  S: 34-37. Seaw. scale. Field measurements.
%  13 = Millero et al, 2006                                 T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  14 = Millero, 2010                                       T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  15 = Waters, Millero, & Woosley, 2014                    T:    0-50  S:  1-50. Seaw. scale. Real seawater.
%  16 = Sulpis et al, 2020                                  T: -1.7-32  S: 31-38. Total scale. Field measurements.
%  17 = Schockman & Byrne, 2021                             T:   15-35  S: 19-41. Total scale. Real seawater.
% 
%  (*4) Each element must be an integer that
%       indicates the KSO4 dissociation constant that is to be used:
%  1 = KSO4 of Dickson   (PREFERRED)
%  2 = KSO4 of Khoo
%  3 = KSO4 of Waters and Millero
%
%  (*5) Each element must be an integer that 
%       indicates the KHF dissociation constant that is to be used:
%  1 = KF of Dickson & Riley 1979  
%  2 = KF of Perez & Fraga, 1987  (PREFERRED)
%
%  (*6) Each element must be an integer that 
%       indicates the the formulation of the borate-to-salinity ratio to be used:
%  1 = TB of Uppstrom 1979
%  2 = TB of Lee 2010  (PREFERRED)
%
%**************************************************************************
%
% OUTPUT: *  an array containing the following parameter values (one row per sample):
%         *  a cell-array containing crudely formatted headers
%         *  a cell-array containing nicely formatted headers
%
%   POS   PARAMETER            UNIT
%
%    01 - TAlk                 (umol/kgSW)
%    02 - TCO2                 (umol/kgSW)
%    03 - pHin                 ()
%    04 - pCO2 input           (uatm)
%    05 - fCO2 input           (uatm)
%    06 - HCO3 input           (umol/kgSW)
%    07 - CO3 input            (umol/kgSW)
%    08 - CO2 input            (umol/kgSW)
%    09 - BAlk input           (umol/kgSW)
%    10 - OH input             (umol/kgSW)
%    11 - PAlk input           (umol/kgSW)
%    12 - SiAlk input          (umol/kgSW)
%    13 - AmmAlk input         (umol/kgSW)
%    14 - HSAlk input          (umol/kgSW)
%    15 - Hfree input          (umol/kgSW)
%    16 - RevelleFactor input  ()
%    17 - OmegaCa input        ()
%    18 - OmegaAr input        ()
%    19 - xCO2 input           (ppm)
%    20 - SIR input            ()
%    21 - pH output            ()
%    22 - pCO2 output          (uatm)
%    23 - fCO2 output          (uatm)
%    24 - HCO3 output          (umol/kgSW)
%    25 - CO3 output           (umol/kgSW)
%    26 - CO2 output           (umol/kgSW)
%    27 - BAlk output          (umol/kgSW)
%    28 - OH output            (umol/kgSW)
%    29 - PAlk output          (umol/kgSW)
%    30 - SiAlk output         (umol/kgSW)
%    31 - AmmAlk output        (umol/kgSW)
%    32 - HSAlk output         (umol/kgSW)
%    33 - Hfree output         (umol/kgSW)
%    34 - RevelleFactor output ()
%    35 - OmegaCa output       ()
%    36 - OmegaAr output       ()
%    37 - xCO2 output          (ppm)
%    38 - SIR output            ()
%    39 - pH input (Total)     ()          
%    40 - pH input (SWS)       ()          
%    41 - pH input (Free)      ()          
%    42 - pH input (NBS)       ()          
%    43 - pH output (Total)    ()          
%    44 - pH output (SWS)      ()          
%    45 - pH output (Free)     ()          
%    46 - pH output (NBS)      () 
%    47 - TEMP input           (deg C)     ***    
%    48 - TEMPOUT              (deg C)     ***
%    49 - PRES input           (dbar or m) ***
%    50 - PRESOUT              (dbar or m) ***
%    51 - PAR1TYPE             (integer)   ***
%    52 - PAR2TYPE             (integer)   ***
%    53 - K1K2CONSTANTS        (integer)   ***
%    54 - KSO4CONSTANT         (integer)   *** 
%    55 - KFCONSTANT           (integer)   *** 
%    56 - BORON                (integer)   *** 
%    57 - pHSCALE of input     (integer)   ***
%    58 - SAL                  (psu)       ***
%    59 - PO4                  (umol/kgSW) ***
%    60 - SI                   (umol/kgSW) ***
%    61 - NH4                  (umol/kgSW) ***
%    62 - H2S                  (umol/kgSW) ***
%    63 - K0  input            ()          
%    64 - K1  input            ()          
%    65 - K2  input            ()          
%    66 - pK1 input            ()          
%    67 - pK2 input            ()          
%    68 - KW  input            ()          
%    69 - KB  input            ()          
%    70 - KF  input            ()          
%    71 - KS  input            ()          
%    72 - KP1 input            ()          
%    73 - KP2 input            ()          
%    74 - KP3 input            ()          
%    75 - KSi input            ()              
%    76 - KNH4input            ()    
%    77 - KH2Sinput            ()    
%    78 - K0  output           ()          
%    79 - K1  output           ()          
%    80 - K2  output           ()          
%    81 - pK1 output           ()          
%    82 - pK2 output           ()          
%    83 - KW  output           ()          
%    84 - KB  output           ()          
%    85 - KF  output           ()          
%    86 - KS  output           ()          
%    87 - KP1 output           ()          
%    88 - KP2 output           ()          
%    89 - KP3 output           ()          
%    90 - KSi output           ()              
%    91 - KNH4output           ()    
%    92 - KH2Soutput           () 
%    93 - TB                   (umol/kgSW) 
%    94 - TF                   (umol/kgSW) 
%    95 - TS                   (umol/kgSW) 
%    96 - TP                   (umol/kgSW) 
%    97 - TSi                  (umol/kgSW)
%    98 - TNH4                 (umol/kgSW) 
%    99 - TH2S                 (umol/kgSW)
%
%    *** SIMPLY RESTATES THE INPUT BY USER 
%
% In all the above, the terms "input" and "output" may be understood
%    to refer to the 2 scenarios for which CO2SYS performs calculations, 
%    each defined by its own combination of temperature and pressure.
%    For instance, one may use CO2SYS to calculate, from measured DIC and
%    TAlk, the pH that that sample will have in the lab (e.g., T=25 degC, P=0
%    dbar), and what the in situ pH would have been (e.g., at T=1 degC, P=4500).
%    A = CO2SYS(2400,2200,1,2,35,25,1,0,4200,1,1,0,0,1,4,1,1,1)
%    pH_lab = A(3);  % 7.8429
%    pH_sea = A(20); % 8.0503
% 
%**************************************************************************
%
% **** Changes since 3.1 by JD Sharp.
%   - rigorous validation performed against PyCO2SYS
%     (https://github.com/mvdh7/PyCO2SYS)
%   - initial pH estimates obtained via the approach of Munhoven (2013)
%   - correction to solution for free scale pH within iterative pH solvers
%   - correction to uncertainty calculation for parameters at output conditions
%   - consitency implemented for [CO2(aq)] calculations
%   - substrate-inhibitor ratio (SIR; Bach, 2015) included as an output argument
%   - input uncertainty in [CO2], [HCO3], and [CO3] should now be in mol/kg
%   - option added for pressure corrections to K0 and fugacity factor
%
% **** Changes since 3.0 by JD Sharp.
%   - added KSO4 of Waters and Millero (2013)
%   - added K1 and K2 of Sulpis et al. (2020)
%   - added K1 and K2 of Schockman and Byrne (2021)
%
% **** Changes since 3.0 by JD Sharp based on code from D Pierrot.
%   - changed code to set pH values that don't converge to NaN. All	
%     subsequent calculated values also set to NaN.
%   - modified input function to separate KHSO4 and TB choices
%   - added KHF of Perez & Fraga as choice for HF dissociation constant
%   - modified output to reflect all changes mentioned above
%
% **** Changes since 3.0 by MP Humphreys.
%   - include Peng correction for Icase 16 and 17.
%   - fix Icase typo for CO2-HCO3 input pair.
%   - make corrections to (F) indexing in a few places.
%
% **** Changes since 2.1 (uploaded to GitHub Jul 2019) by JD Sharp
%	- now allows for input of NH4+ and H2S concentrations
%
% **** Additional changes since 2.0 by JD Sharp
%	- now allows for HCO3, CO3, and CO2 as input parameters for calculation and
%     for error propagation
%
% **** Changes since 2.0
%	- slight changes to allow error propagation
%	- new option to choose K1 & K2 from Waters et al. (2014): fixes inconsistencies with Millero (2010) identified by Orr et al. (2015)
%
% **** Changes since 1.01 (uploaded to CDIAC at June 11th, 2009):
% - Function cleans up its global variables when done (if you lose variables, this may be the cause -- see around line 570)
% - Added the outputting of K values
% - Implementation of constants of Cai and Wang, 1998
% - Implementation of constants of Lueker et al., 2000
% - Implementation of constants of Mojica-Prieto and Millero, 2002
% - Implementation of constants of Millero et al., 2002 (only their eqs. 19, 20, no TCO2 dependency)
% - Implementation of constants of Millero et al., 2006
% - Implementation of constants of Millero et al., 2010
% - Properly listed Sal and Temp limits for the available constants
% - added switch for using the new Lee et al., (2010) formulation of Total Borate
% - Minor corrections to the GEOSECS constants (gave NaN for some output in earlier version)
% - Fixed decimal point error on [H+] (did not get converted to umol/kgSW from mol/kgSW).
% - Changed 'Hfreein' to 'Hfreeout' in the 'NICEHEADERS'-output (typo)
%
% **** Changes since 1.00 (uploaded to CDIAC at May 29th, 2009):
% - added a note explaining that all known bugs were removed before release of 1.00
%
%**************************************************************************
%
% CO2SYS originally by Lewis and Wallace 1998
%
% Converted to MATLAB by Denis Pierrot at
% CIMAS, University of Miami, Miami, Florida
%
% Vectorization, internal refinements and speed improvements by
% Steven van Heuven, University of Groningen, The Netherlands.
% Questions, bug reports et cetera: svheuven@gmail.com
%
% Modifications for error propagation by JM Epitalon
%
% Extension to include input of CO2, HCO3, CO3, NH4, and H2S by
% Jonathan Sharp, University of South Florida.
%
% Modification to set pH values that do not converge to NaN, separate
% KHSO4 and TB, and to add the KHF of Perez & Fraga by Denis Pierrot,
% implemented in this version by Jonathan Sharp, University of Washington
%
% Bug fixes by Matthew Humphreys, NIOZ Texel, the Netherlands.
%
% Additional modifications for consistency with PyCO2SYS and other added
% options and input/output arguments by Jonathan Sharp, University of
% Washington
%
%**************************************************************************


%**************************************************************************
% NOTHING BELOW THIS SHOULD REQUIRE EDITING BY USER!
%**************************************************************************


% Declare global variables
global pHScale WhichKs WhoseKSO4 WhoseKF WhoseTB Pbar
global Sal sqrSal TempK logTempK TempCi TempCo Pdbari Pdbaro;
global FugFac VPFac PengCorrection ntps RGasConstant;
global fH RT;
global K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S CAL F;

% Added by JM Epitalon
% For computing derivative with respect to Ks, one has to call CO2sys with a perturbed K
% Requested perturbation is passed through the following global variables
global PertK    % Id of perturbed K
global Perturb  % perturbation

% Input conditioning

% set default for optional input argument
global p_opt
p_opt = 0;
% parse optional input argument
for i = 1:2:length(varargin)-1
    if strcmpi(varargin{i}, 'co2_press')
        p_opt = varargin{i+1};
    end
end

% Determine lengths of input vectors
veclengths=[length(PAR1) length(PAR2) length(PAR1TYPE)...
            length(PAR2TYPE) length(SAL) length(TEMPIN)...
            length(TEMPOUT) length(PRESIN) length(PRESOUT)...
            length(SI) length(PO4) length(NH4) length(H2S)...
            length(pHSCALEIN) length(K1K2CONSTANTS) length(KSO4CONSTANT)...
	        length(KFCONSTANT) length(BORON)];

if length(unique(veclengths))>2
	disp(' '); disp('*** INPUT ERROR: Input vectors must all be of same length, or of length 1. ***'); disp(' '); return
end

% Make column vectors of all input vectors
PAR1         =PAR1         (:);
PAR2         =PAR2         (:);
PAR1TYPE     =PAR1TYPE     (:);
PAR2TYPE     =PAR2TYPE     (:);
SAL          =SAL          (:);
TEMPIN       =TEMPIN       (:);
TEMPOUT      =TEMPOUT      (:);
PRESIN       =PRESIN       (:);
PRESOUT      =PRESOUT      (:);
SI           =SI           (:);
PO4          =PO4          (:);
NH4          =NH4          (:);
H2S          =H2S          (:);
pHSCALEIN    =pHSCALEIN    (:);
K1K2CONSTANTS=K1K2CONSTANTS(:);
KSO4CONSTANT =KSO4CONSTANT (:);
KFCONSTANT   =KFCONSTANT   (:);
BORON        =BORON        (:);

% Find the longest column vector:
ntps = max(veclengths);

% Populate column vectors
PAR1(1:ntps,1)          = PAR1(:)          ;
PAR2(1:ntps,1)          = PAR2(:)          ;
PAR1TYPE(1:ntps,1)      = PAR1TYPE(:)      ;
PAR2TYPE(1:ntps,1)      = PAR2TYPE(:)      ;
SAL(1:ntps,1)           = SAL(:)           ;
TEMPIN(1:ntps,1)        = TEMPIN(:)        ;
TEMPOUT(1:ntps,1)       = TEMPOUT(:)       ;
PRESIN(1:ntps,1)        = PRESIN(:)        ;
PRESOUT(1:ntps,1)       = PRESOUT(:)       ;
SI(1:ntps,1)            = SI(:)            ;
PO4(1:ntps,1)           = PO4(:)           ;
NH4(1:ntps,1)           = NH4(:)           ;
H2S(1:ntps,1)           = H2S(:)           ;
pHSCALEIN(1:ntps,1)     = pHSCALEIN(:)     ;
K1K2CONSTANTS(1:ntps,1) = K1K2CONSTANTS(:) ;
KSO4CONSTANT(1:ntps,1)  = KSO4CONSTANT(:)  ;
KFCONSTANT(1:ntps,1)    = KFCONSTANT(:)    ;
BORON(1:ntps,1)         = BORON(:)         ;

% Assign input to the 'historical' variable names.
pHScale      = pHSCALEIN;
WhichKs      = K1K2CONSTANTS;
WhoseKSO4    = KSO4CONSTANT;
WhoseKF      = KFCONSTANT;
WhoseTB      = BORON;
p1           = PAR1TYPE;
p2           = PAR2TYPE;
TempCi       = TEMPIN;
TempCo       = TEMPOUT;
Pdbari       = PRESIN;
Pdbaro       = PRESOUT;
Sal          = SAL;
sqrSal       = sqrt(SAL);
TP           = PO4;
TSi          = SI;
TNH4         = NH4;
TH2S         = H2S;

RGasConstant = 83.14462618; % ml bar-1 K-1 mol-1,
%                             recommended by NIST
%                             https://physics.nist.gov/cgi-bin/cuu/Value?r
%RGasConstant = 83.1451;  % ml bar-1 K-1 mol-1, DOEv2
%                           Compatible w/ CO2SYSv2.0.5
%RGasConstant = 83.14472; % ml bar-1 K-1 mol-1, DOEv3

% Generate empty vectors for...
TA   = nan(ntps,1); % Talk
TC   = nan(ntps,1); % DIC
PH   = nan(ntps,1); % pH
PC   = nan(ntps,1); % pCO2
FC   = nan(ntps,1); % fCO2
HCO3 = nan(ntps,1); % [HCO3]
CO3  = nan(ntps,1); % [CO3]
CO2  = nan(ntps,1); % [CO2*]

% Assign values to empty vectors.
F=(p1==1 & PAR1~=-999);   TA(F)=PAR1(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p1==2 & PAR1~=-999);   TC(F)=PAR1(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p1==3 & PAR1~=-999);   PH(F)=PAR1(F);
F=(p1==4 & PAR1~=-999);   PC(F)=PAR1(F)/1e6; % Convert from microatm. to atm.
F=(p1==5 & PAR1~=-999);   FC(F)=PAR1(F)/1e6; % Convert from microatm. to atm.
F=(p1==6 & PAR1~=-999); HCO3(F)=PAR1(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p1==7 & PAR1~=-999);  CO3(F)=PAR1(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p1==8 & PAR1~=-999);  CO2(F)=PAR1(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p2==1 & PAR2~=-999);   TA(F)=PAR2(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p2==2 & PAR2~=-999);   TC(F)=PAR2(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p2==3 & PAR2~=-999);   PH(F)=PAR2(F);
F=(p2==4 & PAR2~=-999);   PC(F)=PAR2(F)/1e6; % Convert from microatm. to atm.
F=(p2==5 & PAR2~=-999);   FC(F)=PAR2(F)/1e6; % Convert from microatm. to atm.
F=(p2==6 & PAR2~=-999); HCO3(F)=PAR2(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p2==7 & PAR2~=-999);  CO3(F)=PAR2(F)/1e6; % Convert from micromol/kg to mol/kg
F=(p2==8 & PAR2~=-999);  CO2(F)=PAR2(F)/1e6; % Convert from micromol/kg to mol/kg

% Generate the columns holding Si, Phos, Amm, H2S and Sal.
% Pure Water case:
F=(WhichKs==8);
Sal(F) = 0;
% GEOSECS and Pure Water:
F=(WhichKs==8 | WhichKs==6);  
TP(F)  = 0;
TSi(F) = 0;
TNH4(F)  = 0;
TH2S(F)  = 0;
% All other cases
F=~F;                         
TP(F)   = TP(F)./1e6;
TSi(F)  = TSi(F)./1e6;
TNH4(F) = TNH4(F)./1e6;
TH2S(F) = TH2S(F)./1e6;

% The vector 'PengCorrection' is used to modify the value of TA, for those
% cases where WhichKs==7, since PAlk(Peng) = PAlk(Dickson) + TP.
% Thus, PengCorrection is 0 for all cases where WhichKs is not 7
PengCorrection=zeros(ntps,1); F=WhichKs==7; PengCorrection(F)=TP(F);

% Calculate the constants for all samples at input conditions
% The constants calculated for each sample will be on the appropriate pH scale!
Constants(TempCi,Pdbari);

% Added by JM Epitalon
% For computing derivative with respect to Ks, one has to perturb the value of one K
% Requested perturbation is passed through global variables PertK and Perturb
if (~ isempty(PertK))
    switch PertK
        case {'K0'}
            K0 = K0 + Perturb;
        case {'K1'}
            K1 = K1 + Perturb;
        case {'K2'}
            K2 = K2 + Perturb;
        case {'KB'}
            KB = KB + Perturb;
        case {'KW'}
            KW = KW + Perturb;
        case {'BOR'}
            TB = TB + Perturb;
    end
end


% Make sure fCO2 is available for each sample that has pCO2 or CO2.
F = (~isnan(PC) & (p1==4 | p2==4));  FC(F) = PC(F).*FugFac(F);
F = (~isnan(CO2) & (p1==8 | p2==8)); FC(F) = CO2(F)./K0(F);

% Generate vectors for results, and copy the raw input values into them
TAc    = TA;
TCc    = TC;
PHic   = PH;
PCic   = PC;
FCic   = FC;
HCO3ic = HCO3;
CO3ic  = CO3;
CO2ic  = CO2;

% Generate vector describing the combination of input parameters
% So, the valid ones are:
% 12,13,15,16,17,18,23,25,26,27,28,35,36,37,38,56,57,67,68,78
Icase = 10*min(p1,p2) + max(p1,p2);

% Calculate missing values for AT,CT,PH,FC,HCO3,CO3,CO2:
% pCO2 will be calculated later on, routines work with fCO2.
F=Icase==12; % input TA, TC
if any(F)
F=(~isnan(TAc) & ~isnan(TCc) & F);
    PHic(F)                = CalculatepHfromTATC(TAc(F)-PengCorrection(F),TCc(F));
    F=(~isnan(PHic) & F);
    if any(F)
       FCic(F)              = CalculatefCO2fromTCpH(TCc(F), PHic(F));
       [CO3ic(F),HCO3ic(F)] = CalculateCO3HCO3fromTCpH(TCc(F),PHic(F));
    end
end
F=Icase==13; % input TA, pH
if any(F)
F=(~isnan(TAc) & ~isnan(PHic) & F);
    TCc(F)                 = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
    FCic(F)                = CalculatefCO2fromTCpH(TCc(F),PHic(F));
    [CO3ic(F),HCO3ic(F)]   = CalculateCO3HCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==14 | Icase==15 | Icase==18; % input TA, (pCO2 or fCO2 or CO2)
if any(F)
F=(~isnan(TAc) & ~isnan(FCic) & F);
    PHic(F)                = CalculatepHfromTAfCO2(TAc(F)-PengCorrection(F),FCic(F));
    F=(~isnan(PHic) & F);
    if any(F)
       TCc(F)              = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
       [CO3ic(F),HCO3ic(F)]= CalculateCO3HCO3fromTCpH(TCc(F),PHic(F));
    end
end
F=Icase==16; % input TA, HCO3
if any(F)
F=(~isnan(TAc) & ~isnan(HCO3ic) & F);
    PHic(F)                = CalculatepHfromTAHCO3(TAc(F)-PengCorrection(F),HCO3ic(F));  % added Peng correction // MPH
    F=(~isnan(PHic) & F);
    if any(F)
       TCc(F)              = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
       FCic(F)             = CalculatefCO2fromTCpH(TCc(F),PHic(F)); 
       CO3ic(F)            = CalculateCO3fromTCpH(TCc(F),PHic(F));
    end
end
F=Icase==17; % input TA, CO3
if any(F)
F=(~isnan(TAc) & ~isnan(CO3ic) & F);
    PHic(F)                 = CalculatepHfromTACO3(TAc(F)-PengCorrection(F),CO3ic(F));  % added Peng correction // MPH
    F=(~isnan(PHic) & F);
    if any(F)
       TCc(F)               = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
       FCic(F)              = CalculatefCO2fromTCpH(TCc(F),PHic(F)); 
       HCO3ic(F)            = CalculateHCO3fromTCpH(TCc(F),PHic(F));
    end
end
F=Icase==23; % input TC, pH
if any(F)
F=(~isnan(TCc) & ~isnan(PHic) & F);
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    FCic(F)                 = CalculatefCO2fromTCpH(TCc(F),PHic(F));
    [CO3ic(F),HCO3ic(F)]    = CalculateCO3HCO3fromTCpH(TCc(F), PHic(F));
end
F=Icase==24 | Icase==25 | Icase==28;  % input TC, (pCO2 or fCO2 or CO2)
if any(F)
F=(~isnan(TCc) & ~isnan(FCic) & F);
    PHic(F)                 = CalculatepHfromTCfCO2(TCc(F),FCic(F));
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    [CO3ic(F),HCO3ic(F)]    = CalculateCO3HCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==26; % input TC, HCO3
if any(F)
F=(~isnan(TCc) & ~isnan(HCO3ic) & F);
    [PHic(F),FCic(F)]       = CalculatepHfCO2fromTCHCO3(TCc(F),HCO3ic(F));
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    CO3ic(F)                = CalculateCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==27; % input TC, CO3
if any(F)
F=(~isnan(TCc) & ~isnan(CO3ic) & F);
    [PHic(F),FCic(F)]       = CalculatepHfCO2fromTCCO3(TCc(F),CO3ic(F));
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    HCO3ic(F)               = CalculateHCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==34 | Icase==35 | Icase==38; % input pH, (pCO2 or fCO2 or CO2)
if any(F)
F=(~isnan(PHic) & ~isnan(FCic) & F);
    TCc(F)                  = CalculateTCfrompHfCO2(PHic(F),FCic(F));
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    [CO3ic(F),HCO3ic(F)]    = CalculateCO3HCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==36; % input pH, HCO3
if any(F)
F=(~isnan(PHic) & ~isnan(HCO3ic) & F);
    TAc(F)                  = CalculateTAfrompHHCO3(PHic(F),HCO3ic(F)) + PengCorrection(F);
    TCc(F)                  = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
    FCic(F)                 = CalculatefCO2fromTCpH(TCc(F),PHic(F));
    CO3ic(F)                = CalculateCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==37; % input pH, CO3
if any(F)
F=(~isnan(PHic) & ~isnan(CO3ic) & F);
    TAc(F)                  = CalculateTAfrompHCO3(PHic(F),CO3ic(F)) + PengCorrection(F);
    TCc(F)                  = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
    FCic(F)                 = CalculatefCO2fromTCpH(TCc(F),PHic(F));
    HCO3ic(F)               = CalculateHCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==46 | Icase==56 | Icase==68; % input (pCO2 or fCO2 or CO2), HCO3
if any(F)
F=(~isnan(FCic) & ~isnan(HCO3ic) & F);
    PHic(F)                 = CalculatepHfromfCO2HCO3(FCic(F),HCO3ic(F));
    TCc(F)                  = CalculateTCfrompHfCO2(PHic(F),FCic(F));
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    CO3ic(F)                = CalculateCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==47 | Icase==57 | Icase==78; % input (pCO2 or fCO2 or CO2), CO3
if any(F)
F=(~isnan(FCic) & ~isnan(CO3ic) & F);
    PHic(F)                 = CalculatepHfromfCO2CO3(FCic(F),CO3ic(F));
    TCc(F)                  = CalculateTCfrompHfCO2 (PHic(F),FCic(F));
    TAc(F)                  = CalculateTAfromTCpH(TCc(F),PHic(F)) + PengCorrection(F);
    HCO3ic(F)               = CalculateHCO3fromTCpH(TCc(F),PHic(F));
end
F=Icase==67; % input HCO3, CO3
if any(F)
F=(~isnan(HCO3ic) & ~isnan(CO3ic) & F);
    PHic(F)                 = CalculatepHfromCO3HCO3(CO3ic(F),HCO3ic(F));
    TAc(F)                  = CalculateTAfrompHCO3(PHic(F),CO3ic(F)) + PengCorrection(F);
    TCc(F)                  = CalculateTCfromTApH(TAc(F)-PengCorrection(F),PHic(F));
    FCic(F)                 = CalculatefCO2fromTCpH(TCc(F),PHic(F));
    %CO2ic(F)                = CalculateCO2fromTCpH(TCc(F),PHic(F));
end

% By now, an fCO2 value is available for each sample.
% Generate the associated pCO2 value:
F = (isnan(PCic) & (p1~=4 | p2~=4)); PCic(F)  = FCic(F)./FugFac(F);
% Generate the associated CO2 value:
F = (isnan(CO2ic) & (p1~=8 | p2~=8)); CO2ic(F) = FCic(F).*K0(F);

% Calculate Other Params At Input Conditions:
BAlkinp    = nan(ntps,1); % Generate empty vectors
[OHinp,PAlkinp,SiAlkinp,AmmAlkinp,HSAlkinp,Hfreeinp,HSO4inp,HFinp,...
    Revelleinp,OmegaCainp,OmegaArinp,xCO2dryinp] = deal(BAlkinp);
F=(~isnan(PHic)); % if PHic = NaN, pH calculation was not performed or did not converge
[BAlkinp(F),OHinp(F), PAlkinp(F),SiAlkinp(F),AmmAlkinp(F),...
    HSAlkinp(F), Hfreeinp(F),HSO4inp(F),HFinp(F)] = CalculateAlkParts(PHic(F));
PAlkinp(F)                = PAlkinp(F)+PengCorrection(F);
Revelleinp(F)             = RevelleFactor(TAc(F)-PengCorrection(F), TCc(F));
[OmegaCainp(F),OmegaArinp(F),CAL] = CaSolubility(Sal(F), TempCi(F), Pdbari(F), TCc(F), PHic(F));
xCO2dryinp(~isnan(PCic),1) = PCic(~isnan(PCic),1)./VPFac(~isnan(PCic),1); % ' this assumes pTot = 1 atm
SIRinp = HCO3ic./(Hfreeinp.*1e6);

% % Just for reference, convert pH at input conditions to the other scales
pHicT = nan(ntps,1);
pHicS = nan(ntps,1);
pHicF = nan(ntps,1);
pHicN = nan(ntps,1);
[pHicT(F),pHicS(F),pHicF(F),pHicN(F)]=FindpHOnAllScales(PHic(F));

% Merge the Ks at input into an array. Ks at output will be glued to this later.
KIVEC=[K0 K1 K2 -log10(K1) -log10(K2) KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S];

% Calculate the constants for all samples at output conditions
Constants(TempCo,Pdbaro);

% Added by JM Epitalon
% For computing derivative with respect to Ks, one has to perturb the value of one K
% Requested perturbation is passed through global variables PertK and Perturb
if (~ isempty(PertK))
    switch PertK
        case {'K0'}
            K0 = K0 + Perturb;
        case {'K1'}
            K1 = K1 + Perturb;
        case {'K2'}
            K2 = K2 + Perturb;
        case {'KB'}
            KB = KB + Perturb;
        case {'KW'}
            KW = KW + Perturb;
        case {'BOR'}
            TB = TB + Perturb;
    end
end                  

% For output conditions, using conservative TA and TC, calculate pH, fCO2
% and pCO2, HCO3, CO3, and CO2
F=(~isnan(TAc) & ~isnan(TCc)); % i.e., do for all samples that have TA and TC values
PHoc=nan(ntps,1);
[CO3oc,HCO3oc,FCoc] = deal(PHoc);
PHoc(F) = CalculatepHfromTATC(TAc(F)-PengCorrection(F), TCc(F)); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
    FCoc(F) = CalculatefCO2fromTCpH(TCc(F), PHoc(F));
    [CO3oc(F),HCO3oc(F)] = CalculateCO3HCO3fromTCpH(TCc(F),PHoc(F));

% Generate the associated pCO2 value:
PCoc  = FCoc./FugFac;
% Generate the associated CO2 value:
CO2oc = FCoc.*K0;

% Calculate Other Params At Output Conditions:
BAlkout    = nan(ntps,1); % Generate empty vectors
[OHout,PAlkout,SiAlkout,AmmAlkout,HSAlkout,Hfreeout,HSO4out,HFout,...
    Revelleout,OmegaCaout,OmegaArout,xCO2dryout] = deal(BAlkout);
F=(~isnan(PHoc)); % if PHoc = NaN, pH calculation was not performed or did not converge
[BAlkout(F),OHout(F),PAlkout(F),SiAlkout(F),AmmAlkout(F),...
    HSAlkout(F), Hfreeout(F),HSO4out(F),HFout(F)] = CalculateAlkParts(PHoc(F));
PAlkout(F)                 = PAlkout(F)+PengCorrection(F);
Revelleout(F)              = RevelleFactor(TAc(F)-PengCorrection(F), TCc(F));
[OmegaCaout(F),OmegaArout(F),CAL] = CaSolubility(Sal(F), TempCo(F), Pdbaro(F), TCc(F), PHoc(F));
xCO2dryout(~isnan(PCoc),1)    = PCoc(~isnan(PCoc))./VPFac(~isnan(PCoc)); % ' this assumes pTot = 1 atm
SIRout = HCO3oc./(Hfreeout.*1e6);

% Just for reference, convert pH at output conditions to the other scales
pHocT = nan(ntps,1);
pHocS = nan(ntps,1);
pHocF = nan(ntps,1);
pHocN = nan(ntps,1);
[pHocT(F),pHocS(F),pHocF(F),pHocN(F)]=FindpHOnAllScales(PHoc(F));

KOVEC=[K0 K1 K2 -log10(K1) -log10(K2) KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S];
TVEC =[TB TF TS TP TSi TNH4 TH2S CAL];
% keyboard
% Saving data in array, 99 columns, as many rows as samples input
DATA=[TAc*1e6         TCc*1e6        PHic           PCic*1e6        FCic*1e6...
      HCO3ic*1e6      CO3ic*1e6      CO2ic*1e6      BAlkinp*1e6     OHinp*1e6...
      PAlkinp*1e6     SiAlkinp*1e6   AmmAlkinp*1e6  HSAlkinp*1e6    Hfreeinp*1e6... %%% Multiplied Hfreeinp *1e6, svh20100827
      Revelleinp      OmegaCainp     OmegaArinp     xCO2dryinp*1e6  SIRinp...
      PHoc            PCoc*1e6       FCoc*1e6       HCO3oc*1e6      CO3oc*1e6...
      CO2oc*1e6       BAlkout*1e6    OHout*1e6      PAlkout*1e6     SiAlkout*1e6...
      AmmAlkout*1e6   HSAlkout*1e6   Hfreeout*1e6   Revelleout      OmegaCaout... %%% Multiplied Hfreeout *1e6, svh20100827
      OmegaArout      xCO2dryout*1e6 SIRout         pHicT           pHicS...
      pHicF           pHicN          pHocT          pHocS           pHocF...
      pHocN           TEMPIN         TEMPOUT        PRESIN          PRESOUT...
      PAR1TYPE        PAR2TYPE       K1K2CONSTANTS  KSO4CONSTANT    KFCONSTANT...
      BORON           pHSCALEIN      SAL            PO4             SI...
      NH4             H2S            KIVEC          KOVEC           TVEC*1e6       fH];
DATA(isnan(DATA))=-999;

HEADERS={'TAlk';'TCO2';'pHin';'pCO2in';'fCO2in';'HCO3in';'CO3in';...
    'CO2in';'BAlkin';'OHin';'PAlkin';'SiAlkin';'AmmAlkin';'HSAlkin';...
    'Hfreein';'RFin';'OmegaCAin';'OmegaARin';'xCO2in';'SIRin';'pHout';...
    'pCO2out';'fCO2out';'HCO3out';'CO3out';'CO2out';'BAlkout';'OHout';...
    'PAlkout';'SiAlkout';'AmmAlkout';'HSAlkout';'Hfreeout';'RFout';'OmegaCAout';...
    'OmegaARout';'xCO2out';'SIRout';'pHinTOTAL';'pHinSWS';'pHinFREE';'pHinNBS';...
    'pHoutTOTAL';'pHoutSWS';'pHoutFREE';'pHoutNBS';'TEMPIN';'TEMPOUT';...
    'PRESIN';'PRESOUT';'PAR1TYPE';'PAR2TYPE';'K1K2CONSTANTS';'KSO4CONSTANT';... KSO4CONSTANTS => KSO4CONSTANT // MPH
    'KFCONSTANT';'BORON';'pHSCALEIN';'SAL';'PO4';'SI';'NH4';'H2S';'K0input';...
    'K1input';'K2input';'pK1input';'pK2input';'KWinput';'KBinput';'KFinput';...
    'KSinput';'KP1input';'KP2input';'KP3input';'KSiinput';'KNH4input';...
    'KH2Sinput';'K0output';'K1output';'K2output';'pK1output';'pK2output';...
    'KWoutput';'KBoutput';'KFoutput';'KSoutput';'KP1output';'KP2output';...
    'KP3output';'KSioutput';'KNH4output';'KH2Soutput';'TB';'TF';'TS';...
    'TP';'TSi';'TNH4';'TH2S';'TCal';'fH'};

NICEHEADERS={...
    '01 - TAlk             (umol/kgSW) ';
    '02 - TCO2             (umol/kgSW) ';
    '03 - pHin             ()          ';
    '04 - pCO2in           (uatm)      ';
    '05 - fCO2in           (uatm)      ';
    '06 - HCO3in           (umol/kgSW) ';
    '07 - CO3in            (umol/kgSW) ';
    '08 - CO2in            (umol/kgSW) ';
    '09 - BAlkin           (umol/kgSW) ';
    '10 - OHin             (umol/kgSW) ';
    '11 - PAlkin           (umol/kgSW) ';
    '12 - SiAlkin          (umol/kgSW) ';
    '13 - AmmAlkin         (umol/kgSW) ';
    '14 - HSAlkin          (umol/kgSW) ';
    '15 - Hfreein          (umol/kgSW) ';
    '16 - RevelleFactorin  ()          ';
    '17 - OmegaCain        ()          ';
    '18 - OmegaArin        ()          ';
    '19 - xCO2in           (ppm)       ';
    '20 - SIRin            ()          ';
    '21 - pHout            ()          ';
    '22 - pCO2out          (uatm)      ';
    '23 - fCO2out          (uatm)      ';
    '24 - HCO3out          (umol/kgSW) ';
    '25 - CO3out           (umol/kgSW) ';
    '26 - CO2out           (umol/kgSW) ';
    '27 - BAlkout          (umol/kgSW) ';
    '28 - OHout            (umol/kgSW) ';
    '29 - PAlkout          (umol/kgSW) ';
    '30 - SiAlkout         (umol/kgSW) ';
    '31 - AmmAlkout        (umol/kgSW) ';
    '32 - HSAlkout         (umol/kgSW) ';
    '33 - Hfreeout         (umol/kgSW) ';
    '34 - RevelleFactorout ()          ';
    '35 - OmegaCaout       ()          ';
    '36 - OmegaArout       ()          ';
    '37 - xCO2out          (ppm)       ';
    '38 - SIRout           ()          ';
    '39 - pHin (Total)     ()          ';
    '40 - pHin (SWS)       ()          ';
    '41 - pHin (Free)      ()          ';
    '42 - pHin (NBS )      ()          ';
    '43 - pHout(Total)     ()          ';
    '44 - pHout(SWS)       ()          ';
    '45 - pHout(Free)      ()          ';
    '46 - pHout(NBS )      ()          ';
    '47 - TEMPIN           (Deg C)     ';    
    '48 - TEMPOUT          (Deg C)     ';
    '49 - PRESIN           (dbar)      ';
    '50 - PRESOUT          (dbar)      ';
    '51 - PAR1TYPE         ()          ';
    '52 - PAR2TYPE         ()          ';
    '53 - K1K2CONSTANTS    ()          ';
    '54 - KSO4CONSTANT     ()          ';
    '55 - KFCONSTANT       ()          ';
    '56 - BORON            ()          ';
    '57 - pHSCALEIN        ()          ';
    '58 - SAL              (umol/kgSW) ';
    '59 - PO4              (umol/kgSW) ';
    '60 - SI               (umol/kgSW) ';
    '61	- NH4	           (umol/kgSW) ';
    '62	- H2S	           (umol/kgSW) ';
    '63 - K0input          ()          ';
    '64 - K1input          ()          ';
    '65 - K2input          ()          ';
    '66 - pK1input         ()          ';
    '67 - pK2input         ()          ';
    '68 - KWinput          ()          ';
    '69 - KBinput          ()          ';
    '70 - KFinput          ()          ';
    '71 - KSinput          ()          ';
    '72 - KP1input         ()          ';
    '73 - KP2input         ()          ';
    '74 - KP3input         ()          ';
    '75 - KSiinput         ()          ';
    '76 - KNH4input        ()          ';
    '77 - KH2Sinput        ()          ';  
    '78 - K0output         ()          ';
    '79 - K1output         ()          ';
    '80 - K2output         ()          ';
    '81 - pK1output        ()          ';
    '82 - pK2output        ()          ';
    '83 - KWoutput         ()          ';
    '84 - KBoutput         ()          ';
    '85 - KFoutput         ()          ';
    '86 - KSoutput         ()          ';
    '87 - KP1output        ()          ';
    '88 - KP2output        ()          ';
    '89 - KP3output        ()          ';
    '90 - KSioutput        ()          ';
    '91 - KNH4output       ()          ';
    '92 - KH2Soutput       ()          ';
    '93 - KCaoutput        ()          ';
    '94 - KAroutput        ()          ';
    '95 - TB               (umol/kgSW) ';
    '96 - TF               (umol/kgSW) ';
    '97 - TS               (umol/kgSW) ';
    '98 - TP               (umol/kgSW) ';
    '99 - TSi              (umol/kgSW) ';
    '100 - TNH4            (umol/kgSW) ';
    '101 - TH2S            (umol/kgSW) ';
    '102 - TCal            (umol/kgSW) ';
    '103 - fH              (umol/kgSW) '};

clear global F K2 KP3 Pdbari Sal TS VPFac ntps 
clear global FugFac KB KS Pdbaro T TSi BORON WhichKs pHScale 
clear global K KF KSi KNH4 KH2S PengCorrection TB TempCi WhoseKSO4 WhoseKF WhoseTB sqrSal 
clear global K0 KP1 KW RGasConstant TF TempCo fH 
clear global K1 KP2 Pbar RT TP TempK logTempK
	
end % end main function


%**************************************************************************
% Subroutines:
%**************************************************************************


function Constants(TempC,Pdbar)
global pHScale WhichKs WhoseKSO4 WhoseKF WhoseTB sqrSal Pbar RT;
global K0 fH FugFac VPFac ntps TempK logTempK;
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi CAL RGasConstant Sal p_opt;

% SUB Constants, version 04.01, 10-13-97, written by Ernie Lewis.
% Inputs: pHScale%, WhichKs%, WhoseKSO4%, Sali, TempCi, Pdbar
% Outputs: K0, K(), T(), fH, FugFac, VPFac
% This finds the Constants of the CO2 system in seawater or freshwater,
% corrects them for pressure, and reports them on the chosen pH scale.
% The process is as follows: the Constants (except KS, KF which stay on the
% free scale - these are only corrected for pressure) are
%       1) evaluated as they are given in the literature
%       2) converted to the SWS scale in mol/kg-SW or to the NBS scale
%       3) corrected for pressure
%       4) converted to the SWS pH scale in mol/kg-SW
%       5) converted to the chosen pH scale
%
%       PROGRAMMER'S NOTE: all logs are log base e
%       PROGRAMMER'S NOTE: all Constants are converted to the pH scale
%               pHScale% (the chosen one) in units of mol/kg-SW
%               except KS and KF are on the free scale
%               and KW is in units of (mol/kg-SW)^2
TempK    = TempC + 273.15;
RT       = RGasConstant.*TempK;
logTempK = log(TempK);
Pbar     = Pdbar ./ 10;

% Generate empty vectors for holding results
TB = nan(ntps,1);
TF = nan(ntps,1);
TS = nan(ntps,1);
CAL = nan(ntps,1);

% CalculateTB - Total Borate:
F=(WhichKs==8); % Pure water case.
if any(F)
    TB(F) = 0;
end
F=(WhichKs==6 | WhichKs==7);
if any(F)
    TB(F) = 0.0004106.*Sal(F)./35; % in mol/kg-SW
    % this is .00001173.*Sali
    % this is about 1% lower than Uppstrom's value
    % Culkin, F., in Chemical Oceanography,
    % ed. Riley and Skirrow, 1965:
    % GEOSECS references this, but this value is not explicitly
    % given here
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8); % All other cases
if any(F)
	FF=F&(WhoseTB==1); % If user opted for Uppstrom's values:
	if any(FF)
	    % Uppstrom, L., Deep-Sea Research 21:161-162, 1974:
	    % this is .000416.*Sali./35. = .0000119.*Sali
		% TB(FF) = (0.000232./10.811).*(Sal(FF)./1.80655); % in mol/kg-SW
	    TB(FF) =  0.0004157.*Sal(FF)./35; % in mol/kg-SW
	end
	FF=F&(WhoseTB==2); % If user opted for the Lee et al. values:
	if any(FF)
		% Lee, Kim, Byrne, Millero, Feely, Yong-Ming Liu. 2010.	
	 	% Geochimica Et Cosmochimica Acta 74 (6): 1801-1811.
		TB(FF) =  0.0004326.*Sal(FF)./35; % in mol/kg-SW
	end
end

% CalculateCAL - Total Calcium:
F=(WhichKs~=6 & WhichKs~=7);
    % Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
    % this is .010285.*Sali./35
    CAL(F) = 0.02128./40.087.*(Sal(F)./1.80655);
F=(WhichKs==6 | WhichKs==7);
    % *** CalculateCaforGEOSECS:
    % Culkin, F, in Chemical Oceanography, ed. Riley and Skirrow, 1965:
    % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982) 
    CAL(F) = 0.01026.*Sal(F)./35;

% CalculateTF;
% Riley, J. P., Deep-Sea Research 12:219-220, 1965:
% this is .000068.*Sali./35. = .00000195.*Sali
TF = (0.000067./18.998).*(Sal./1.80655); % in mol/kg-SW

% CalculateTS ;
% Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
% this is .02824.*Sali./35. = .0008067.*Sali
TS = (0.14./96.062).*(Sal./1.80655); % in mol/kg-SW

% CalculateK0:
% Weiss, R. F., Marine Chemistry 2:203-215, 1974.
TempK100  = TempK./100;
lnK0 = -60.2409 + 93.4517 ./ TempK100 + 23.3585 .* log(TempK100) + Sal .*...
    (0.023517 - 0.023656 .* TempK100 + 0.0047036 .* TempK100 .^2);
K0   = exp(lnK0);                 % this is in mol/kg-SW/atm
vCO2 = 32.3;                      % partial molal volume of CO2 (cm3 / mol)
                                  % from Weiss (1974, Appendix, paragraph 3)
if p_opt == 0
    prscorr = 1; % Set pressure correction to 1
elseif p_opt == 1
    prscorr = exp((-Pbar).*vCO2./RT); % Calculate pressure correction to K0
else
    disp('co2_press must be set to either 0 or 1'); % Display error message
end         
K0   = K0 .* prscorr; % this is in mol/kg-SW/atm

% CalculateIonS:
% This is from the DOE handbook, Chapter 5, p. 13/22, eq. 7.2.4:
IonS         = 19.924 .* Sal ./ (1000 - 1.005   .* Sal);

% CalculateKS:
lnKS   = nan(ntps,1); pKS  = nan(ntps,1); KS   = nan(ntps,1);
logKS0 = nan(ntps,1); logKSK0 = nan(ntps,1);
F=(WhoseKSO4==1);
if any(F)
    % Dickson, A. G., J. Chemical Thermodynamics, 22:113-127, 1990
    % The goodness of fit is .021.
    % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
    % TYPO on p. 121: the constant e9 should be e8.
    % This is from eqs 22 and 23 on p. 123, and Table 4 on p 121:
  lnKS(F) = -4276.1./TempK(F) + 141.328 - 23.093.*logTempK(F) +...             
      (-13856./TempK(F) + 324.57 - 47.986.*logTempK(F)).*sqrt(IonS(F)) +...     
      (35474./TempK(F) - 771.54 + 114.723.*logTempK(F)).*IonS(F) +...           
      (-2698./TempK(F)).*sqrt(IonS(F)).*IonS(F) + (1776./TempK(F)).*IonS(F).^2; 
	KS(F) = exp(lnKS(F))...            % this is on the free pH scale in mol/kg-H2O
        .* (1 - 0.001005 .* Sal(F));   % convert to mol/kg-SW
end
F=(WhoseKSO4==2);
if any(F)
    % Khoo et al, Analytical Chemistry, 49(1):29-34, 1977
    % KS was found by titrations with a hydrogen electrode
    % of artificial seawater containing sulfate (but without F)
    % at 3 salinities from 20 to 45 and artificial seawater NOT
    % containing sulfate (nor F) at 16 salinities from 15 to 45,
    % both at temperatures from 5 to 40 deg C.
    % KS is on the Free pH scale (inherently so).
    % It was given in mol/kg-H2O. I convert it to mol/kg-SW.
    % He finds log(beta) which = my pKS;
    % his beta is an association constant.
    % The rms error is .0021 in pKS, or about .5% in KS.
    % This is equation 20 on p. 33:
    pKS(F) = 647.59 ./ TempK(F) - 6.3451 + 0.019085.*TempK(F) - 0.5208.*sqrt(IonS(F));
    KS(F) = 10.^(-pKS(F))...          % this is on the free pH scale in mol/kg-H2O
        .* (1 - 0.001005.*Sal(F));    % convert to mol/kg-SW
end
F=(WhoseKSO4==3);
if any(F)
    % Waters and Millero, Marine Chemistry, 149: 8-22, 2013, with corrections from
    % Waters et al, Marine Chemistry, 165: 66-67, 2014
    logKS0(F) = 562.69486 - 102.5154.*logTempK(F) - 0.0001117033.*TempK(F).*TempK(F) + ...
        0.2477538.*TempK(F) - 13273.76./TempK(F);
    logKSK0(F) = (4.24666 - 0.152671.*TempK(F) + 0.0267059.*TempK(F).*logTempK(F) - 0.000042128.*TempK(F).*TempK(F)).*Sal(F).^0.5 + ...
        (0.2542181 - 0.00509534.*TempK(F) + 0.00071589.*TempK(F).*logTempK(F)).*Sal(F) + (-0.00291179 + 0.0000209968.*TempK(F)).*Sal(F).^1.5 + ...
        -0.0000403724.*Sal(F).^2;
    KS(F) = ((10.^(logKSK0(F))).*(10.^logKS0(F))) ... % this is on the free pH scale in mol/kg-H2O
        .* (1 - 0.001005.*Sal(F));                    % convert to mol/kg-SW
end

% CalculateKF:
KF = NaN(ntps, 1);  % added preallocation here and F-indexing below // MPH
F=(WhoseKF==1);
if any(F)
    % Dickson, A. G. and Riley, J. P., Marine Chemistry 7:89-99, 1979:
    lnKF = 1590.2./TempK - 12.641 + 1.525.*IonS.^0.5;
    KF(F)   = exp(lnKF(F))...                 % this is on the free pH scale in mol/kg-H2O
        .*(1 - 0.001005.*Sal(F));          % convert to mol/kg-SW
end
F=(WhoseKF==2);
if any(F)
    % Perez and Fraga 1987 (to be used for S: 10-40, T: 9-33)
    % P&F87 might actually be better than the fit of D&R79 above, which is based on only three salinities: [0 26.7 34.6]
    lnKF = 874./TempK - 9.68 + 0.111.*Sal.^0.5;
    KF(F)   = exp(lnKF(F));                   % this is on the free pH scale in mol/kg-SW
end

% CalculatepHScaleConversionFactors:
%       These are NOT pressure-corrected.
SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);
FREEtoTOT =  1 + TS./KS;

% CalculatefH
fH = nan(ntps,1);
% Use GEOSECS's value for cases 1,2,3,4,5 (and 6) to convert pH scales.
F=(WhichKs==8);
if any(F)
    fH(F) = 1; % this shouldn't occur in the program for this case
end
F=(WhichKs==7);
if any(F)
    fH(F) = 1.29 - 0.00204.*  TempK(F) + (0.00046 -...
        0.00000148.*TempK(F)).*Sal(F).*Sal(F);
    % Peng et al, Tellus 39B:439-458, 1987:
    % They reference the GEOSECS report, but round the value
    % given there off so that it is about .008 (1%) lower. It
    % doesn't agree with the check value they give on p. 456.
end
F=(WhichKs~=7 & WhichKs~=8);
if any(F)
    fH(F) = 1.2948 - 0.002036.*TempK(F) + (0.0004607 -...
        0.000001475.*TempK(F)).*Sal(F).^2;
    % Takahashi et al, Chapter 3 in GEOSECS Pacific Expedition,
    % v. 3, 1982 (p. 80);
end

% CalculateKB:
KB      = nan(ntps,1); logKB   = nan(ntps,1);
lnKBtop = nan(ntps,1); lnKB    = nan(ntps,1);
F=(WhichKs==8); % Pure water case
if any(F)
    KB(F) = 0;
end
F=(WhichKs==6 | WhichKs==7);
if any(F)
    % This is for GEOSECS and Peng et al.
    % Lyman, John, UCLA Thesis, 1957
    % fit by Li et al, JGR 74:5507-5525, 1969:
    logKB(F) = -9.26 + 0.00886.*Sal(F) + 0.01.*TempC(F);
    KB(F) = 10.^(logKB(F))...  % this is on the NBS scale
        ./fH(F);               % convert to the SWS scale
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8);
if any(F)
    % Dickson, A. G., Deep-Sea Research 37:755-766, 1990:
    lnKBtop(F) = -8966.9 - 2890.53.*sqrSal(F) - 77.942.*Sal(F) +...
        1.728.*sqrSal(F).*Sal(F) - 0.0996.*Sal(F).^2;
    lnKB(F) = lnKBtop(F)./TempK(F) + 148.0248 + 137.1942.*sqrSal(F) +...
        1.62142.*Sal(F) + (-24.4344 - 25.085.*sqrSal(F) - 0.2474.*...
        Sal(F)).*logTempK(F) + 0.053105.*sqrSal(F).*TempK(F);
    KB(F) = exp(lnKB(F))...    % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);         % convert to SWS pH scale
end

% CalculateKW:
lnKW = nan(ntps,1); KW = nan(ntps,1);
F=(WhichKs==7);
if any(F)
    % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
    lnKW(F) = 148.9802 - 13847.26./TempK(F) - 23.6521.*logTempK(F) +...
        (-79.2447 + 3298.72./TempK(F) + 12.0408.*logTempK(F)).*...
        sqrSal(F) - 0.019813.*Sal(F);
end
F=(WhichKs==8);
if any(F)
    % Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979
    % refit data of Harned and Owen, The Physical Chemistry of
    % Electrolyte Solutions, 1958
    lnKW(F) = 148.9802 - 13847.26./TempK(F) - 23.6521.*logTempK(F);
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8);
if any(F)
    % Millero, Geochemica et Cosmochemica Acta 59:661-677, 1995.
    % his check value of 1.6 umol/kg-SW should be 6.2
    lnKW(F) = 148.9802 - 13847.26./TempK(F) - 23.6521.*logTempK(F) +...
        (-5.977 + 118.67./TempK(F) + 1.0495.*logTempK(F)).*...
        sqrSal(F) - 0.01615.*Sal(F);
end
KW = exp(lnKW); % this is on the SWS pH scale in (mol/kg-SW)^2
F=(WhichKs==6);
if any(F)
    KW(F) = 0; % GEOSECS doesn't include OH effects
end

% CalculateKP1KP2KP3KSi:
KP1      = nan(ntps,1); KP2      = nan(ntps,1);
KP3      = nan(ntps,1); KSi      = nan(ntps,1);
lnKP1    = nan(ntps,1); lnKP2    = nan(ntps,1);
lnKP3    = nan(ntps,1); lnKSi    = nan(ntps,1);
F=(WhichKs==7);
if any(F)
    KP1(F) = 0.02;
    % Peng et al don't include the contribution from this term,
    % but it is so small it doesn't contribute. It needs to be
    % kept so that the routines work ok.
    % KP2, KP3 from Kester, D. R., and Pytkowicz, R. M.,
    % Limnology and Oceanography 12:243-252, 1967:
    % these are only for sals 33 to 36 and are on the NBS scale
    KP2(F) = exp(-9.039 - 1450./TempK(F))... % this is on the NBS scale
        ./fH(F);                          % convert to SWS scale
    KP3(F) = exp(4.466 - 7276./TempK(F))...  % this is on the NBS scale
        ./fH(F);                          % convert to SWS scale
    % Sillen, Martell, and Bjerrum,  Stability Constants of metal-ion complexes,
    % The Chemical Society (London), Special Publ. 17:751, 1964:
    KSi(F) = 0.0000000004...              % this is on the NBS scale
        ./fH(F);                          % convert to SWS scale
    
end
F=(WhichKs==6 | WhichKs==8);
if any(F)
    KP1(F) = 0; KP2(F) = 0; KP3(F) = 0; KSi(F) = 0;
    % Neither the GEOSECS choice nor the freshwater choice
    % include contributions from phosphate or silicate.
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8);
if any(F)
    % Yao and Millero, Aquatic Geochemistry 1:53-88, 1995
    % KP1, KP2, KP3 are on the SWS pH scale in mol/kg-SW.
    % KSi was given on the SWS pH scale in molal units.
    lnKP1(F) = -4576.752./TempK(F) + 115.54 - 18.453.*logTempK(F) + (-106.736./TempK(F) +...
        0.69171).*sqrSal(F) + (-0.65643./TempK(F) - 0.01844).*Sal(F);
    KP1(F) = exp(lnKP1(F));
    lnKP2(F) = -8814.715./TempK(F) + 172.1033 - 27.927.*logTempK(F) + (-160.34./TempK(F) +...
        1.3566).*sqrSal(F) + (0.37335./TempK(F) - 0.05778).*Sal(F);
    KP2(F) = exp(lnKP2(F));
    lnKP3(F) = -3070.75./TempK(F) - 18.126 + (17.27039./TempK(F) + 2.81197).*sqrSal(F) +...
        (-44.99486./TempK(F) - 0.09984).*Sal(F);
    KP3(F) = exp(lnKP3(F));
    lnKSi(F) = -8904.2./TempK(F) + 117.4 - 19.334.*logTempK(F) + (-458.79./TempK(F) +...
        3.5913).*sqrt(IonS(F)) + (188.74./TempK(F) - 1.5998).*IonS(F) +...
        (-12.1652./TempK(F) + 0.07871).*IonS(F).^2;
    KSi(F) = exp(lnKSi(F))...                % this is on the SWS pH scale in mol/kg-H2O
        .*(1 - 0.001005.*Sal(F));        % convert to mol/kg-SW
end

% Calculate KNH4 and KH2S: added by J. Sharp
KNH4           = nan(ntps,1); KH2S       = nan(ntps,1);
PKNH4expCW     = nan(ntps,1); lnKH2S     = nan(ntps,1);
F=(WhichKs==6 | WhichKs==7 | WhichKs==8); % GEOSECS or freshwater cases
if any(F)
    KNH4(F) = 0;
    KH2S(F) = 0;
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8); % All other cases
if any(F)
% Ammonia dissociation constant from Yao and Millero (1995)
%   KNH4(F) = (exp(-6285.33./TempK(F)+0.0001635.*TempK(F)-0.25444+...
%             (0.46532-123.7184./TempK(F)).*Sal(F).^0.5+(-0.01992+...
%             3.17556./TempK(F)).*Sal(F)))...
%             ./SWStoTOT(F);                    % convert to SWS pH scale
% Ammonia dissociation constant from Clegg and Whitfield (1995)
  PKNH4(F) = 9.244605-2729.33.*(1./298.15-1./TempK(F)) +...
          (0.04203362-11.24742./TempK(F)).*Sal(F).^0.25+...  % added missing (F) index on Sal // MPH
          (-13.6416+1.176949.*TempK(F).^0.5-...
          0.02860785.*TempK(F)+545.4834./TempK(F)).*Sal(F).^0.5+...
          (-0.1462507+0.0090226468.*TempK(F).^0.5-...
          0.0001471361.*TempK(F)+10.5425./TempK(F)).*Sal(F).^1.5+...
          (0.004669309-0.0001691742.*TempK(F).^0.5-...
          0.5677934./TempK(F)).*Sal(F).^2+...
          (-2.354039E-05+0.009698623./TempK(F)).*Sal(F).^2.5;
  KNH4(F)  = 10.^-PKNH4(F);                    % total scale, mol/kg-H2O
  KNH4(F)  = KNH4(F).*(1-0.001005.*Sal(F)); % mol/kg-SW
  KNH4(F)  = KNH4(F)./SWStoTOT(F);             % converts to SWS pH scale

% First hydrogen sulfide dissociation constant from Millero et al. (1988)
  KH2S(F)  = (exp(225.838-13275.3./TempK(F)-34.6435.*log(TempK(F))+...
              0.3449.*Sal(F).^0.5-0.0274.*Sal(F)))...
              ./SWStoTOT(F);                    % convert to SWS pH scale

end

% CalculateK1K2:
logK1    = nan(ntps,1); lnK1     = nan(ntps,1);
pK1      = nan(ntps,1); K1       = nan(ntps,1);
logK2    = nan(ntps,1); lnK2     = nan(ntps,1);
pK2      = nan(ntps,1); K2       = nan(ntps,1);
F=(WhichKs==1);
if any(F)
    % ROY et al, Marine Chemistry, 44:249-267, 1993
    % (see also: Erratum, Marine Chemistry 45:337, 1994
    % and Erratum, Marine Chemistry 52:183, 1996)
    % Typo: in the abstract on p. 249: in the eq. for lnK1* the
    % last term should have S raised to the power 1.5.
    % They claim standard deviations (p. 254) of the fits as
    % .0048 for lnK1 (.5% in K1) and .007 in lnK2 (.7% in K2).
    % They also claim (p. 258) 2s precisions of .004 in pK1 and
    % .006 in pK2. These are consistent, but Andrew Dickson
    % (personal communication) obtained an rms deviation of about
    % .004 in pK1 and .003 in pK2. This would be a 2s precision
    % of about 2% in K1 and 1.5% in K2.
    % T:  0-45  S:  5-45. Total Scale. Artificial sewater.
    % This is eq. 29 on p. 254 and what they use in their abstract:
    lnK1(F) = 2.83655 - 2307.1266./TempK(F) - 1.5529413.*logTempK(F) +...
        (-0.20760841 - 4.0484./TempK(F)).*sqrSal(F) + 0.08468345.*Sal(F) -...
        0.00654208.*sqrSal(F).*Sal(F);
    K1(F) = exp(lnK1(F))...            % this is on the total pH scale in mol/kg-H2O
        .*(1 - 0.001005.*Sal(F))...    % convert to mol/kg-SW
        ./SWStoTOT(F);                 % convert to SWS pH scale
    % This is eq. 30 on p. 254 and what they use in their abstract:
    lnK2(F) = -9.226508 - 3351.6106./TempK(F) - 0.2005743.*logTempK(F) +...
        (-0.106901773 - 23.9722./TempK(F)).*sqrSal(F) + 0.1130822.*Sal(F) -...
        0.00846934.*sqrSal(F).*Sal(F);
    K2(F) = exp(lnK2(F))...            % this is on the total pH scale in mol/kg-H2O
        .*(1 - 0.001005.*Sal(F))...    % convert to mol/kg-SW
        ./SWStoTOT(F);                 % convert to SWS pH scale
end
F=(WhichKs==2);
if any(F)
    % GOYET AND POISSON, Deep-Sea Research, 36(11):1635-1654, 1989
    % The 2s precision in pK1 is .011, or 2.5% in K1.
    % The 2s precision in pK2 is .02, or 4.5% in K2.
    % This is in Table 5 on p. 1652 and what they use in the abstract:
    pK1(F) = 812.27./TempK(F) + 3.356 - 0.00171.*Sal(F).*logTempK(F)...
        + 0.000091.*Sal(F).^2;
    K1(F) = 10.^(-pK1(F)); % this is on the SWS pH scale in mol/kg-SW
    % This is in Table 5 on p. 1652 and what they use in the abstract:
    pK2(F) = 1450.87./TempK(F) + 4.604 - 0.00385.*Sal(F).*logTempK(F)...
        + 0.000182.*Sal(F).^2;
    K2(F) = 10.^(-pK2(F)); % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==3);
if any(F)
    % HANSSON refit BY DICKSON AND MILLERO
    % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
    % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
    % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973
    % and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973.
    % on the SWS pH scale in mol/kg-SW.
    % Hansson gave his results on the Total scale (he called it
    % the seawater scale) and in mol/kg-SW.
    % Typo in DM on p. 1739 in Table 4: the equation for pK2*
    % for Hansson should have a .000132 *S^2
    % instead of a .000116 *S^2.
    % The 2s precision in pK1 is .013, or 3% in K1.
    % The 2s precision in pK2 is .017, or 4.1% in K2.
    % This is from Table 4 on p. 1739.
    pK1(F) = 851.4./TempK(F) + 3.237 - 0.0106.*Sal(F) + 0.000105.*Sal(F).^2;
    K1(F) = 10.^(-pK1(F)); % this is on the SWS pH scale in mol/kg-SW
    % This is from Table 4 on p. 1739.
    pK2(F) = -3885.4./TempK(F) + 125.844 - 18.141.*logTempK(F)...
        - 0.0192.*Sal(F) + 0.000132.*Sal(F).^2;
    K2(F) = 10.^(-pK2(F)); % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==4);
if any(F)
    % MEHRBACH refit BY DICKSON AND MILLERO
    % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
    % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
    % refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
    % on the SWS pH scale in mol/kg-SW.
    % Mehrbach et al gave results on the NBS scale.
    % The 2s precision in pK1 is .011, or 2.6% in K1.
    % The 2s precision in pK2 is .020, or 4.6% in K2.
	% Valid for salinity 20-40.
    % This is in Table 4 on p. 1739.
    pK1(F) = 3670.7./TempK(F) - 62.008 + 9.7944.*logTempK(F)...
             - 0.0118.*Sal(F) + 0.000116.*Sal(F).^2;
    K1(F) = 10.^(-pK1(F)); % this is on the SWS pH scale in mol/kg-SW
    % This is in Table 4 on p. 1739.
    pK2(F) = 1394.7./TempK(F) + 4.777 - 0.0184.*Sal(F) + 0.000118.*Sal(F).^2;
    K2(F) = 10.^(-pK2(F)); % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==5);
if any(F)
    % HANSSON and MEHRBACH refit BY DICKSON AND MILLERO
    % Dickson and Millero, Deep-Sea Research,34(10):1733-1743, 1987
    % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
    % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973,
    % Hansson, Acta Chemica Scandanavia, 27:931-944, 1973,
    % and Mehrbach et al, Limnol. Oceanogr.,18(6):897-907, 1973
    % on the SWS pH scale in mol/kg-SW.
    % Typo in DM on p. 1740 in Table 5: the second equation
    % should be pK2* =, not pK1* =.
    % The 2s precision in pK1 is .017, or 4% in K1.
    % The 2s precision in pK2 is .026, or 6% in K2.
	% Valid for salinity 20-40.
    % This is in Table 5 on p. 1740.
    pK1(F) = 845./TempK(F) + 3.248 - 0.0098.*Sal(F) + 0.000087.*Sal(F).^2;
    K1(F) = 10.^(-pK1(F)); % this is on the SWS pH scale in mol/kg-SW
    % This is in Table 5 on p. 1740.
    pK2(F) = 1377.3./TempK(F) + 4.824 - 0.0185.*Sal(F) + 0.000122.*Sal(F).^2;
    K2(F) = 10.^(-pK2(F)); % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==6 | WhichKs==7);
if any(F)
    % GEOSECS and Peng et al use K1, K2 from Mehrbach et al,
    % Limnology and Oceanography, 18(6):897-907, 1973.
	% I.e., these are the original Mehrbach dissociation constants.
    % The 2s precision in pK1 is .005, or 1.2% in K1.
    % The 2s precision in pK2 is .008, or 2% in K2.
    pK1(F) = - 13.7201 + 0.031334.*TempK(F) + 3235.76./TempK(F)...
        + 1.3e-5*Sal(F).*TempK(F) - 0.1032.*Sal(F).^0.5;
    K1(F) = 10.^(-pK1(F))...         % this is on the NBS scale
        ./fH(F);                     % convert to SWS scale
    pK2(F) = 5371.9645 + 1.671221.*TempK(F) + 0.22913.*Sal(F) + 18.3802.*log10(Sal(F))...
             - 128375.28./TempK(F) - 2194.3055.*log10(TempK(F)) - 8.0944e-4.*Sal(F).*TempK(F)...
             - 5617.11.*log10(Sal(F))./TempK(F) + 2.136.*Sal(F)./TempK(F); % pK2 is not defined for Sal=0, since log10(0)=-inf
    K2(F) = 10.^(-pK2(F))...         % this is on the NBS scale
        ./fH(F);                     % convert to SWS scale
end
F=(WhichKs==8);
if any(F)	
	% PURE WATER CASE
    % Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
    % K1 from refit data from Harned and Davis,
    % J American Chemical Society, 65:2030-2037, 1943.
    % K2 from refit data from Harned and Scholes,
    % J American Chemical Society, 43:1706-1709, 1941.
	% This is only to be used for Sal=0 water (note the absence of S in the below formulations)
    % These are the thermodynamic Constants:
    lnK1(F) = 290.9097 - 14554.21./TempK(F) - 45.0575.*logTempK(F);
    K1(F) = exp(lnK1(F));
    lnK2(F) = 207.6548 - 11843.79./TempK(F) - 33.6485.*logTempK(F);
    K2(F) = exp(lnK2(F));
end
F=(WhichKs==9);
if any(F)
    % From Cai and Wang 1998, for estuarine use.
	% Data used in this work is from:
	% K1: Mehrbach (1973) for S>15, for S<15: Mook and Keone (1975)
	% K2: Mehrbach (1973) for S>20, for S<20: Edmond and Gieskes (1970)
	% Sigma of residuals between fits and above data: ±0.015, +0.040 for K1 and K2, respectively.
	% Sal 0-40, Temp 0.2-30
  % Limnol. Oceanogr. 43(4) (1998) 657-668
	% On the NBS scale
	% Their check values for F1 don't work out, not sure if this was correctly published...
	F1 = 200.1./TempK(F) + 0.3220;
	pK1(F) = 3404.71./TempK(F) + 0.032786.*TempK(F) - 14.8435 - 0.071692.*F1.*Sal(F).^0.5 + 0.0021487.*Sal(F);
    K1(F)  = 10.^-pK1(F)...         % this is on the NBS scale
        ./fH(F);                    % convert to SWS scale (uncertain at low Sal due to junction potential);
	F2 = -129.24./TempK(F) + 1.4381;
	pK2(F) = 2902.39./TempK(F) + 0.02379.*TempK(F) - 6.4980 - 0.3191.*F2.*Sal(F).^0.5 + 0.0198.*Sal(F);
    K2(F)  = 10.^-pK2(F)...         % this is on the NBS scale
        ./fH(F);                    % convert to SWS scale (uncertain at low Sal due to junction potential); 
end
F=(WhichKs==10);
if any(F)
    % From Lueker, Dickson, Keeling, 2000
	% This is Mehrbach's data refit after conversion to the total scale, for comparison with their equilibrator work. 
    % Mar. Chem. 70 (2000) 105-119
    % Total scale and kg-sw
    pK1(F) = 3633.86./TempK(F)-61.2172+9.6777.*log(TempK(F))-0.011555.*Sal(F)+0.0001152.*Sal(F).^2;
	K1(F)  = 10.^-pK1(F)...           % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);                % convert to SWS pH scale
    pK2(F) = 471.78./TempK(F)+25.929 -3.16967.*log(TempK(F))-0.01781 .*Sal(F)+0.0001122.*Sal(F).^2;
	K2(F)  = 10.^-pK2(F)...           % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);                % convert to SWS pH scale
end
F=(WhichKs==11);
if any(F)
	% Mojica Prieto and Millero 2002. Geochim. et Cosmochim. Acta. 66(14) 2529-2540.
	% sigma for pK1 is reported to be 0.0056
	% sigma for pK2 is reported to be 0.010
	% This is from the abstract and pages 2536-2537
    pK1 =  -43.6977 - 0.0129037.*Sal(F) + 1.364e-4.*Sal(F).^2 + 2885.378./TempK(F) +  7.045159.*log(TempK(F));
    pK2 = -452.0940 + 13.142162.*Sal(F) - 8.101e-4.*Sal(F).^2 + 21263.61./TempK(F) + 68.483143.*log(TempK(F))...
				+ (-581.4428.*Sal(F) + 0.259601.*Sal(F).^2)./TempK(F) - 1.967035.*Sal(F).*log(TempK(F));
	K1(F) = 10.^-pK1; % this is on the SWS pH scale in mol/kg-SW
	K2(F) = 10.^-pK2; % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==12);
if any(F)
	% Millero et al., 2002. Deep-Sea Res. I (49) 1705-1723.
	% Calculated from overdetermined WOCE-era field measurements 
	% sigma for pK1 is reported to be 0.005
	% sigma for pK2 is reported to be 0.008
	% This is from page 1715
    pK1 =  6.359 - 0.00664.*Sal(F) - 0.01322.*TempC(F) + 4.989e-5.*TempC(F).^2;
    pK2 =  9.867 - 0.01314.*Sal(F) - 0.01904.*TempC(F) + 2.448e-5.*TempC(F).^2;
	K1(F) = 10.^-pK1; % this is on the SWS pH scale in mol/kg-SW
	K2(F) = 10.^-pK2; % this is on the SWS pH scale in mol/kg-SW
end
F=(WhichKs==13);
if any(F)
    % From Millero 2006 work on pK1 and pK2 from titrations
	% Millero, Graham, Huang, Bustos-Serrano, Pierrot. Mar.Chem. 100 (2006) 80-94.
    % S=1 to 50, T=0 to 50. On seawater scale (SWS). From titrations in Gulf Stream seawater.
	pK1_0 = -126.34048 + 6320.813./TempK(F) + 19.568224*log(TempK(F));
	A_1   = 13.4191*Sal(F).^0.5 + 0.0331.*Sal(F) - 5.33e-5.*Sal(F).^2;
	B_1   = -530.123*Sal(F).^0.5 - 6.103.*Sal(F);
	C_1   = -2.06950.*Sal(F).^0.5;
	pK1(F)= A_1 + B_1./TempK(F) + C_1.*log(TempK(F)) + pK1_0; % pK1 sigma = 0.0054
    K1(F) = 10.^-(pK1(F));
	pK2_0= -90.18333 + 5143.692./TempK(F) + 14.613358*log(TempK(F));	
	A_2   = 21.0894*Sal(F).^0.5 + 0.1248.*Sal(F) - 3.687e-4.*Sal(F).^2;
	B_2   = -772.483*Sal(F).^0.5 - 20.051.*Sal(F);
	C_2   = -3.3336.*Sal(F).^0.5;
	pK2(F)= A_2 + B_2./TempK(F) + C_2.*log(TempK(F)) + pK2_0; %pK2 sigma = 0.011
    K2(F) = 10.^-(pK2(F));
end
F=(WhichKs==14);
if any(F)
  % From Millero, 2010, also for estuarine use.
	% Marine and Freshwater Research, v. 61, p. 139-142.
	% Fits through compilation of real seawater titration results:
	% Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), Millero et al. (2006)
	% Constants for K's on the SWS;
	% This is from page 141
	pK10 = -126.34048 + 6320.813./TempK(F) + 19.568224.*log(TempK(F));
	% This is from their table 2, page 140.
	A1 = 13.4038.*Sal(F).^0.5 + 0.03206.*Sal(F) - 5.242e-5.*Sal(F).^2;
	B1 = -530.659.*Sal(F).^0.5 - 5.8210.*Sal(F);
	C1 = -2.0664*Sal(F).^0.5;
	pK1 = pK10 + A1 + B1./TempK(F) + C1.*log(TempK(F));
	K1(F) = 10.^-pK1;
	% This is from page 141
	pK20 =  -90.18333 + 5143.692./TempK(F) + 14.613358.*log(TempK(F));
	% This is from their table 3, page 140.
	A2 = 21.3728.*Sal(F).^0.5 + 0.1218.*Sal(F) - 3.688e-4.*Sal(F).^2;
	B2 = -788.289.*Sal(F).^0.5 - 19.189.*Sal(F);
	C2 = -3.374.*Sal(F).^0.5;
	pK2 = pK20 + A2 + B2./TempK(F) + C2.*log(TempK(F));
	K2(F) = 10.^-pK2;
end
F=(WhichKs==15);
% Added by J. C. Orr on 4 Dec 2016
if any(F)
    % From Waters, Millero, and Woosley, 2014
	% Mar. Chem., 165, 66-67, 2014
  % Corrigendum to "The free proton concentration scale for seawater pH".
	% Effectively, this is an update of Millero (2010) formulation (WhichKs==14)
	% Constants for K's on the SWS;
	pK10 = -126.34048 + 6320.813./TempK(F) + 19.568224.*log(TempK(F));
	A1 = 13.409160.*Sal(F).^0.5 + 0.031646.*Sal(F) - 5.1895e-5.*Sal(F).^2;
	B1 = -531.3642.*Sal(F).^0.5 - 5.713.*Sal(F);
	C1 = -2.0669166.*Sal(F).^0.5;
	pK1 = pK10 + A1 + B1./TempK(F) + C1.*log(TempK(F));
	K1(F) = 10.^-pK1;
	pK20 =  -90.18333 + 5143.692./TempK(F) + 14.613358.*log(TempK(F));
	A2 = 21.225890.*Sal(F).^0.5 + 0.12450870.*Sal(F) - 3.7243e-4.*Sal(F).^2;
	B2 = -779.3444.*Sal(F).^0.5 - 19.91739.*Sal(F);
	C2 = -3.3534679.*Sal(F).^0.5;
	pK2 = pK20 + A2 + B2./TempK(F) + C2.*log(TempK(F));
	K2(F) = 10.^-pK2;
end
F=(WhichKs==16);
% Added by J. D. Sharp on 9 Jul 2020
if any(F)
    % From Sulpis et al, 2020
	% Ocean Science Discussions, 16, 847-862
    % This study uses overdeterminations of the carbonate system to
    % iteratively fit K1 and K2
    pK1(F) = 8510.63./TempK(F)-172.4493+26.32996.*log(TempK(F))-0.011555.*Sal(F)+0.0001152.*Sal(F).^2;
	K1(F)  = 10.^-pK1(F)...           % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);                % convert to SWS pH scale
    pK2(F) = 4226.23./TempK(F)-59.4636+9.60817.*log(TempK(F))-0.01781 .*Sal(F)+0.0001122.*Sal(F).^2;
	K2(F)  = 10.^-pK2(F)...           % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);                % convert to SWS pH scale
end

F=(WhichKs==17);
% Added by J. D. Sharp on 15 Feb 2021
if any(F)
    % From Schockman & Byrne, 2021
	% Geochimica et Cosmochimica Acta, in press
    % This study uses spectrophotometric pH measurements to determine
    % K1*K2 with unprecedented precision, and presents a new
    % parameterization for K2 based on these determinations
    % K1 is taken from Waters, Millero, and Woosley, 2014, on the total pH scale:
    pK10 = -126.34048 + 6320.813./TempK(F) + 19.568224.*log(TempK(F));
	A1 = 13.568513.*Sal(F).^0.5 + 0.031645.*Sal(F) - 5.3834e-5.*Sal(F).^2;
	B1 = -539.2304.*Sal(F).^0.5 - 5.635.*Sal(F);
	C1 = -2.0901396.*Sal(F).^0.5;
	pK1 = pK10 + A1 + B1./TempK(F) + C1.*log(TempK(F));
	K1(F) = 10.^-pK1...               % this is on the total pH scale in mol/kg-sw
        ./SWStoTOT(F);                % convert to SWS pH scale
    % K2 is based on measurements of K1*K2:
    pK2 = 116.8067 - 3655.02./TempK(F) - 16.45817.*log(TempK(F)) + ...
        0.04523.*Sal(F) - 0.615.*Sal(F).^0.5 - 0.0002799.*Sal(F).^2 + ...
        4.969.*(Sal(F)./TempK(F));
    K2(F)  = 10.^-pK2...           % this is on the total pH scale in mol/kg-SW
        ./SWStoTOT(F);                % convert to SWS pH scale
end

%***************************************************************************
%CorrectKsForPressureNow:
% Currently: For WhichKs% = 1 to 7, all Ks (except KF and KS, which are on
%       the free scale) are on the SWS scale.
%       For WhichKs% = 6, KW set to 0, KP1, KP2, KP3, KSi don't matter.
%       For WhichKs% = 8, K1, K2, and KW are on the "pH" pH scale
%       (the pH scales are the same in this case); the other Ks don't matter.
%
%
% No salinity dependence is given for the pressure coefficients here.
% It is assumed that the salinity is at or very near Sali = 35.
% These are valid for the SWS pH scale, but the difference between this and
% the total only yields a difference of .004 pH units at 1000 bars, much
% less than the uncertainties in the values.
%****************************************************************************
% The sources used are:
% Millero, 1995:
%       Millero, F. J., Thermodynamics of the carbon dioxide system in the
%       oceans, Geochemica et Cosmochemica Acta 59:661-677, 1995.
%       See table 9 and eqs. 90-92, p. 675.
%       TYPO: a factor of 10^3 was left out of the definition of Kappa
%       TYPO: the value of R given is incorrect with the wrong units
%       TYPO: the values of the a's for H2S and H2O are from the 1983
%                values for fresh water
%       TYPO: the value of a1 for B(OH)3 should be +.1622
%        Table 9 on p. 675 has no values for Si.
%       There are a variety of other typos in Table 9 on p. 675.
%       There are other typos in the paper, and most of the check values
%       given don't check.
% Millero, 1992:
%       Millero, Frank J., and Sohn, Mary L., Chemical Oceanography,
%       CRC Press, 1992. See chapter 6.
%       TYPO: this chapter has numerous typos (eqs. 36, 52, 56, 65, 72,
%               79, and 96 have typos).
% Millero, 1983:
%       Millero, Frank J., Influence of pressure on chemical processes in
%       the sea. Chapter 43 in Chemical Oceanography, eds. Riley, J. P. and
%       Chester, R., Academic Press, 1983.
%       TYPO: p. 51, eq. 94: the value -26.69 should be -25.59
%       TYPO: p. 51, eq. 95: the term .1700t should be .0800t
%       these two are necessary to match the values given in Table 43.24
% Millero, 1979:
%       Millero, F. J., The thermodynamics of the carbon dioxide system
%       in seawater, Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
%       See table 5 and eqs. 7, 7a, 7b on pp. 1656-1657.
% Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
%       TYPO: the pressure dependence of K2 should have a 16.4, not 26.4
%       This matches the GEOSECS results and is in Edmond and Gieskes.
% Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic acid,
%       boric acid, and the pH of seawater, Limnology and Oceanography
%       13:403-417, 1968.
% Edmond, John M. and Gieskes, J. M. T. M., The calculation of the degree of
%       seawater with respect to calcium carbonate under in situ conditions,
%       Geochemica et Cosmochemica Acta, 34:1261-1291, 1970.
%****************************************************************************
% These references often disagree and give different fits for the same thing.
% They are not always just an update either; that is, Millero, 1995 may agree
%       with Millero, 1979, but differ from Millero, 1983.
% For WhichKs% = 7 (Peng choice) I used the same factors for KW, KP1, KP2,
%       KP3, and KSi as for the other cases. Peng et al didn't consider the
%       case of P different from 0. GEOSECS did consider pressure, but didn't
%       include Phos, Si, or OH, so including the factors here won't matter.
% For WhichKs% = 8 (freshwater) the values are from Millero, 1983 (for K1, K2,
%       and KW). The other aren't used (TB = TS = TF = TP = TSi = 0.), so
%       including the factors won't matter.
%****************************************************************************
%       deltaVs are in cm3/mole
%       Kappas are in cm3/mole/bar
%****************************************************************************

%CorrectK1K2KBForPressure:
deltaV    = nan(ntps,1); Kappa     = nan(ntps,1);
lnK1fac   = nan(ntps,1); lnK2fac   = nan(ntps,1);
lnKBfac   = nan(ntps,1);
F=(WhichKs==8);
if any(F)
    %***PressureEffectsOnK1inFreshWater:
    %               This is from Millero, 1983.
    deltaV(F)  = -30.54 + 0.1849 .*TempC(F) - 0.0023366.*TempC(F).^2;
    Kappa(F)   = (-6.22 + 0.1368 .*TempC(F) - 0.001233 .*TempC(F).^2)./1000;
    lnK1fac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    %***PressureEffectsOnK2inFreshWater:
    %               This is from Millero, 1983.
    deltaV(F)  = -29.81 + 0.115.*TempC(F) - 0.001816.*TempC(F).^2;
    Kappa(F)   = (-5.74 + 0.093.*TempC(F) - 0.001896.*TempC(F).^2)./1000;
    lnK2fac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    lnKBfac(F) = 0 ;%; this doesn't matter since TB = 0 for this case
end
F=(WhichKs==6 | WhichKs==7);
if any(F)
    %               GEOSECS Pressure Effects On K1, K2, KB (on the NBS scale)
    %               Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982 quotes
    %               Culberson and Pytkowicz, L and O 13:403-417, 1968:
    %               but the fits are the same as those in
    %               Edmond and Gieskes, GCA, 34:1261-1291, 1970
    %               who in turn quote Li, personal communication
    lnK1fac(F) = (24.2 - 0.085.*TempC(F)).*Pbar(F)./RT(F);
    lnK2fac(F) = (16.4 - 0.04 .*TempC(F)).*Pbar(F)./RT(F);
    %               Takahashi et al had 26.4, but 16.4 is from Edmond and Gieskes
    %               and matches the GEOSECS results
    lnKBfac(F) = (27.5 - 0.095.*TempC(F)).*Pbar(F)./RT(F);
end
F=(WhichKs~=6 & WhichKs~=7 & WhichKs~=8);
if any(F)
    %***PressureEffectsOnK1:
    %               These are from Millero, 1995.
    %               They are the same as Millero, 1979 and Millero, 1992.
    %               They are from data of Culberson and Pytkowicz, 1968.
    deltaV(F)  = -25.5 + 0.1271.*TempC(F);
    %                 'deltaV = deltaV - .151.*(Sali - 34.8); % Millero, 1979
    Kappa(F)   = (-3.08 + 0.0877.*TempC(F))./1000;
    %                 'Kappa = Kappa  - .578.*(Sali - 34.8)/1000.; % Millero, 1979
 	lnK1fac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    %               The fits given in Millero, 1983 are somewhat different.
    
    %***PressureEffectsOnK2:
    %               These are from Millero, 1995.
    %               They are the same as Millero, 1979 and Millero, 1992.
    %               They are from data of Culberson and Pytkowicz, 1968.
    deltaV(F)  = -15.82 - 0.0219.*TempC(F);
    %                  'deltaV = deltaV + .321.*(Sali - 34.8); % Millero, 1979
    Kappa(F)   = (1.13 - 0.1475.*TempC(F))./1000;
    %                 'Kappa = Kappa - .314.*(Sali - 34.8)./1000: % Millero, 1979
	lnK2fac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    %               The fit given in Millero, 1983 is different.
    %               Not by a lot for deltaV, but by much for Kappa. %
    
    %***PressureEffectsOnKB:
    %               This is from Millero, 1979.
    %               It is from data of Culberson and Pytkowicz, 1968.
    deltaV(F)  = -29.48 + 0.1622.*TempC(F) - 0.002608.*TempC(F).^2;
    %               Millero, 1983 has:
    %                 'deltaV = -28.56 + .1211.*TempCi - .000321.*TempCi.*TempCi
    %               Millero, 1992 has:
    %                 'deltaV = -29.48 + .1622.*TempCi + .295.*(Sali - 34.8)
    %               Millero, 1995 has:
    %                 'deltaV = -29.48 - .1622.*TempCi - .002608.*TempCi.*TempCi
    %                 'deltaV = deltaV + .295.*(Sali - 34.8); % Millero, 1979
    Kappa(F)   = -2.84./1000; % Millero, 1979
    %               Millero, 1992 and Millero, 1995 also have this.
    %                 'Kappa = Kappa + .354.*(Sali - 34.8)./1000: % Millero,1979
    %               Millero, 1983 has:
    %                 'Kappa = (-3 + .0427.*TempCi)./1000
    lnKBfac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
end

% CorrectKWForPressure:
lnKWfac   = nan(ntps,1);
F=(WhichKs==8);
if any(F)
    % PressureEffectsOnKWinFreshWater:
    %               This is from Millero, 1983.
    deltaV(F)  =  -25.6 + 0.2324.*TempC(F) - 0.0036246.*TempC(F).^2;
    Kappa(F)   = (-7.33 + 0.1368.*TempC(F) - 0.001233 .*TempC(F).^2)./1000;
 	lnKWfac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);

    %               NOTE the temperature dependence of KappaK1 and KappaKW
    %               for fresh water in Millero, 1983 are the same.
end
F=(WhichKs~=8);
if any(F)
    % GEOSECS doesn't include OH term, so this won't matter.
    % Peng et al didn't include pressure, but here I assume that the KW correction
    %       is the same as for the other seawater cases.
    % PressureEffectsOnKW:
    %               This is from Millero, 1983 and his programs CO2ROY(T).BAS.
    deltaV(F)  = -20.02 + 0.1119.*TempC(F) - 0.001409.*TempC(F).^2;
    %               Millero, 1992 and Millero, 1995 have:
    Kappa(F)   = (-5.13 + 0.0794.*TempC(F))./1000; % Millero, 1983
    %               Millero, 1995 has this too, but Millero, 1992 is different.
	lnKWfac(F) = (-deltaV(F) + 0.5.*Kappa(F).*Pbar(F)).*Pbar(F)./RT(F);
    %               Millero, 1979 does not list values for these.
end

% PressureEffectsOnKF:
%       This is from Millero, 1995, which is the same as Millero, 1983.
%       It is assumed that KF is on the free pH scale.
deltaV = -9.78 - 0.009.*TempC - 0.000942.*TempC.^2;
Kappa = (-3.91 + 0.054.*TempC)./1000;
lnKFfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKS:
%       This is from Millero, 1995, which is the same as Millero, 1983.
%       It is assumed that KS is on the free pH scale.
deltaV = -18.03 + 0.0466.*TempC + 0.000316.*TempC.^2;
Kappa = (-4.53 + 0.09.*TempC)./1000;
lnKSfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

% CorrectKP1KP2KP3KSiForPressure:
% These corrections don't matter for the GEOSECS choice (WhichKs% = 6) and
%       the freshwater choice (WhichKs% = 8). For the Peng choice I assume
%       that they are the same as for the other choices (WhichKs% = 1 to 5).
% The corrections for KP1, KP2, and KP3 are from Millero, 1995, which are the
%       same as Millero, 1983.
% PressureEffectsOnKP1:
deltaV = -14.51 + 0.1211.*TempC - 0.000321.*TempC.^2;
Kappa  = (-2.67 + 0.0427.*TempC)./1000;
lnKP1fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKP2:
deltaV = -23.12 + 0.1758.*TempC - 0.002647.*TempC.^2;
Kappa  = (-5.15 + 0.09  .*TempC)./1000;
lnKP2fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKP3:
deltaV = -26.57 + 0.202 .*TempC - 0.003042.*TempC.^2;
Kappa  = (-4.08 + 0.0714.*TempC)./1000;
lnKP3fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKSi:
%  The only mention of this is Millero, 1995 where it is stated that the
%    values have been estimated from the values of boric acid. HOWEVER,
%    there is no listing of the values in the table.
%    I used the values for boric acid from above.
deltaV = -29.48 + 0.1622.*TempC - 0.002608.*TempC.^2;
Kappa  = -2.84./1000;
lnKSifac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKNH4: added by J. Sharp
deltaV = -26.43 + 0.0889.*TempC - 0.000905.*TempC.^2;
Kappa  = (-5.03 + 0.0814.*TempC)./1000;
lnKNH4fac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;
% PressureEffectsOnKH2S: added by J. Sharp
deltaV = -11.07 - 0.009.*TempC - 0.000942.*TempC.^2;
Kappa  = (-2.89 + 0.054 .*TempC)./1000;
lnKH2Sfac = (-deltaV + 0.5.*Kappa.*Pbar).*Pbar./RT;

% CorrectKsForPressureHere:
K1fac  = exp(lnK1fac);  K1  = K1 .*K1fac;
K2fac  = exp(lnK2fac);  K2  = K2 .*K2fac;
KWfac  = exp(lnKWfac);  KW  = KW .*KWfac;
KBfac  = exp(lnKBfac);  KB  = KB .*KBfac;
KFfac  = exp(lnKFfac);  KF  = KF .*KFfac;
KSfac  = exp(lnKSfac);  KS  = KS .*KSfac;
KP1fac = exp(lnKP1fac); KP1 = KP1.*KP1fac;
KP2fac = exp(lnKP2fac); KP2 = KP2.*KP2fac;
KP3fac = exp(lnKP3fac); KP3 = KP3.*KP3fac;
KSifac = exp(lnKSifac); KSi = KSi.*KSifac;
KNH4fac= exp(lnKNH4fac);KNH4= KNH4.*KNH4fac; % added by J. Sharp
KH2Sfac= exp(lnKH2Sfac);KH2S= KH2S.*KH2Sfac; % added by J. Sharp

% CorrectpHScaleConversionsForPressure:
% fH has been assumed to be independent of pressure.
SWStoTOT  = (1 + TS./KS)./(1 + TS./KS + TF./KF);
FREEtoTOT =  1 + TS./KS;

%  The values KS and KF are already pressure-corrected, so the pH scale
%  conversions are now valid at pressure.

% FindpHScaleConversionFactor:
% this is the scale they will be put on
pHfactor  = nan(ntps,1);
F=(pHScale==1); %Total
pHfactor(F) = SWStoTOT(F);
F=(pHScale==2); %SWS, they are all on this now
pHfactor(F) = 1;
F=(pHScale==3); %pHfree
pHfactor(F) = SWStoTOT(F)./FREEtoTOT(F);
F=(pHScale==4); %pHNBS
pHfactor(F) = fH(F);

% ConvertFromSWSpHScaleToChosenScale:
K1   = K1.* pHfactor;  K2   = K2.* pHfactor;
KW   = KW.* pHfactor;  KB   = KB.* pHfactor;
KP1  = KP1.*pHfactor;  KP2  = KP2.*pHfactor;
KP3  = KP3.*pHfactor;  KSi  = KSi.*pHfactor;
KNH4 = KNH4.*pHfactor; KH2S = KH2S.*pHfactor;

% CalculateFugacityConstants:
% In previos versions of CO2SYS, the fugacity factor was calculated
% assuming pressure at one atmosphere, or close to it. Starting with
% v3.2.1, an option to use in situ pressure is provided.
%       Weiss, R. F., Marine Chemistry 2:203-215, 1974.
%       Delta and B in cm3/mol
FugFac=ones(ntps,1);
Delta = (57.7 - 0.118.*TempK);
b = -1636.75 + 12.0408.*TempK - 0.0327957.*TempK.^2 + 3.16528.*0.00001.*TempK.^3;
% For a mixture of CO2 and air at in situ pressure;
xc2 = 1; % assumed to be 1, though not strictly correct (xc2 = [1-xCO2]^2)
P1atm = 1.01325; % atmospheric pressure in bar
if p_opt == 0
    FugFac = exp((b + 2.*xc2.*Delta).*P1atm./RT); % FugFac at 1 atm
elseif p_opt == 1
    FugFac = exp((b + 2.*xc2.*Delta).*(P1atm+Pbar)./RT); % FugFac at in situ pressure
else
    disp('co2_press must be set to either 0 or 1'); % Display error message
end
F=(WhichKs==6 | WhichKs==7); % GEOSECS and Peng assume pCO2 = fCO2, or FugFac = 1
FugFac(F) = 1;

% CalculateVPFac:
% Weiss, R. F., and Price, B. A., Nitrous oxide solubility in water and
%       seawater, Marine Chemistry 8:347-359, 1980.
% They fit the data of Goff and Gratch (1946) with the vapor pressure
%       lowering by sea salt as given by Robinson (1954).
% This fits the more complicated Goff and Gratch, and Robinson equations
%       from 273 to 313 deg K and 0 to 40 Sali with a standard error
%       of .015%, about 5 uatm over this range.
% This may be on IPTS-29 since they didn't mention the temperature scale,
%       and the data of Goff and Gratch came before IPTS-48.
% The references are:
% Goff, J. A. and Gratch, S., Low pressure properties of water from -160 deg
%       to 212 deg F, Transactions of the American Society of Heating and
%       Ventilating Engineers 52:95-122, 1946.
% Robinson, Journal of the Marine Biological Association of the U. K.
%       33:449-455, 1954.
%       This is eq. 10 on p. 350.
%       This is in atmospheres.
VPWP = exp(24.4543 - 67.4509.*(100./TempK) - 4.8489.*log(TempK./100));
VPCorrWP = exp(-0.000544.*Sal);
VPSWWP = VPWP.*VPCorrWP;
VPFac = 1 - VPSWWP; % this assumes 1 atmosphere
end % end nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Core input functions, based on functions written by Ernie Lewis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=CalculatepHfromTATC(TAi, TCi)
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F;
% ' SUB CalculatepHfromTATC, version 04.01, 10-13-96, written by Ernie Lewis
% ' with modifications from Denis Pierrot.
% ' Inputs: TA, TC, K(), T()
% ' Output: pH
%
% ' This calculates pH from TA and TC using K1 and K2 by Newton's method.
% ' It tries to solve for the pH at which Residual = 0.
% ' The starting guess is determined by the method introduced by Munhoven
% ' (2013).
%
% ' Made this to accept vectors. It will continue iterating until all
% ' values in the vector are "abs(deltapH) < pHTol". SVH2007
% ' However, once a given abs(deltapH) is less than pHTol, that pH value
% ' will be locked in. This avoids erroneous contributions to results from
% ' other lines of input. JDS2020
K1F=K1(F);     K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
vl          = sum(F);  % VectorLength
% Find initital pH guess using method of Munhoven (2013)
pHGuess         = CalculatepHfromTATCMunhoven(TAi, TCi);
ln10            = log(10);
pH              = pHGuess;
pHTol           = 0.0001;  % tolerance for iterations end
deltapH(1:vl,1) = pHTol+1;
loopc=0;
nF=(abs(deltapH) > pHTol);
while any(nF)
    H         = 10.^(-pH);
    Denom     = (H.*H + K1F.*H + K1F.*K2F);
    CAlk      = TCi.*K1F.*(H + 2.*K2F)./Denom;
    BAlk      = TBF.*KBF./(KBF + H);
    OH        = KWF./H;
    PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
    PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
    PAlk      = TPF.*PhosTop./PhosBot;
    SiAlk     = TSiF.*KSiF./(KSiF + H);
    AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
    HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
    [~,~,pHfree,~] = FindpHOnAllScales(pH); % this converts pH to pHfree no matter the scale
    Hfree     = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree); % since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree); % since KF is on the free scale
    Residual  = TAi - CAlk - BAlk - OH - PAlk - SiAlk  - AmmAlk - HSAlk + Hfree + HSO4 + HF;
    % find Slope dTA/dpH;
    % (this is not exact, but keeps all important terms);
    Slope     = ln10.*(TCi.*K1F.*H.*(H.*H + K1F.*K2F + 4.*H.*K2F)./Denom./Denom + BAlk.*H./(KBF + H) + OH + H);
    deltapH   = Residual./Slope; %' this is Newton's method
    % ' to keep the jump from being too big:
    while any(abs(deltapH) > 1)
        FF=abs(deltapH)>1; deltapH(FF)=deltapH(FF)./2;
    end
    pH(nF) = pH(nF) + deltapH(nF);
    nF     = abs(deltapH) > pHTol;
    loopc=loopc+1;
 
    if loopc>10000
        Fr=find(abs(deltapH) > pHTol);
        pH(Fr)=NaN;  disp(['pH value did not converge for data on row(s): ' num2str((Fr)')]);
        deltapH=pHTol*0.9;
    end
end
varargout{1}=pH;
end % end nested function

function varargout=CalculatefCO2fromTCpH(TCx, pHx)
global K0 K1 K2 F
% ' SUB CalculatefCO2fromTCpH, version 02.02, 12-13-96, written by Ernie Lewis.
% ' Inputs: TC, pH, K0, K1, K2
% ' Output: fCO2
% ' This calculates fCO2 from TC and pH, using K0, K1, and K2.
H            = 10.^(-pHx);
fCO2x        = TCx.*H.*H./(H.*H + K1(F).*H + K1(F).*K2(F))./K0(F);
varargout{1} = fCO2x;
end % end nested function

function varargout=CalculateTCfromTApH(TAx, pHx)
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
K1F=K1(F);     K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
% ' SUB CalculateTCfromTApH, version 02.03, 10-10-97, written by Ernie Lewis.
% ' Inputs: TA, pH, K(), T()
% ' Output: TC
% ' This calculates TC from TA and pH.
H         = 10.^(-pHx);
BAlk      = TBF.*KBF./(KBF + H);
OH        = KWF./H;
PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
PAlk      = TPF.*PhosTop./PhosBot;
SiAlk     = TSiF.*KSiF./(KSiF + H);
AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
[~,~,pHfree,~] = FindpHOnAllScales(pHx); % this converts pH to pHfree no matter the scale
Hfree     = 10.^-pHfree; % this converts pHfree to Hfree
HSO4      = TSF./(1 + KSF./Hfree); %' since KS is on the free scale
HF        = TFF./(1 + KFF./Hfree); %' since KF is on the free scale
CAlk      = TAx - BAlk - OH - PAlk - SiAlk - AmmAlk - HSAlk + Hfree + HSO4 + HF;
TCxtemp   = CAlk.*(H.*H + K1F.*H + K1F.*K2F)./(K1F.*(H + 2.*K2F));
varargout{1} = TCxtemp;
end % end nested function

function varargout=CalculatepHfromTAfCO2(TAi, fCO2i)
global K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
% ' SUB CalculatepHfromTAfCO2, version 04.01, 10-13-97, written by Ernie
% ' Lewis with modifications from Denis Pierrot.
% ' Inputs: TA, fCO2, K0, K(), T()
% ' Output: pH
% ' This calculates pH from TA and fCO2 using K1 and K2 by Newton's method.
% ' It tries to solve for the pH at which Residual = 0.
% ' The starting guess is determined by the method introduced by Munhoven
% ' (2013) and extended by Humphreys et al. (2021).
%
% ' This will continue iterating until all values in the vector are
% ' "abs(deltapH) < pHTol"
% ' However, once a given abs(deltapH) is less than pHTol, that pH value
% ' will be locked in. This avoids erroneous contributions to results from
% ' other lines of input.
K0F=K0(F);     K1F=K1(F);     K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
vl         = sum(F); % vectorlength
% Find initital pH guess using method of Munhoven (2013)
CO2i       = fCO2i.*K0F; % Convert fCO2 to CO2
pHGuess    = CalculatepHfromTACO2Munhoven(TAi, CO2i);
ln10       = log(10);
pH         = pHGuess;
pHTol      = 0.0001; % tolerance
deltapH = pHTol+pH;
loopc=0;
nF=(abs(deltapH) > pHTol);
while any(nF)
    H         = 10.^(-pH);
    HCO3      = K0F.*K1F.*fCO2i./H;
    CO3       = K0F.*K1F.*K2F.*fCO2i./(H.*H);
    CAlk      = HCO3 + 2.*CO3;
    BAlk      = TBF.*KBF./(KBF + H);
    OH        = KWF./H;
    PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
    PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
    PAlk      = TPF.*PhosTop./PhosBot;
    SiAlk     = TSiF.*KSiF./(KSiF + H);
    AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
    HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
    [~,~,pHfree,~] = FindpHOnAllScales(pH); % this converts pH to pHfree no matter the scale
    Hfree     = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree); %' since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
    Residual  = TAi - CAlk - BAlk - OH - PAlk - SiAlk - AmmAlk - HSAlk + Hfree + HSO4 + HF;
    % '               find Slope dTA/dpH
    % '               (this is not exact, but keeps all important terms):
    Slope     = ln10.*(HCO3 + 4.*CO3 + BAlk.*H./(KBF + H) + OH + H);
    deltapH   = Residual./Slope; %' this is Newton's method
    % ' to keep the jump from being too big:
    while any(abs(deltapH) > 1)
        FF=abs(deltapH)>1; deltapH(FF)=deltapH(FF)./2;
    end
    pH(nF) = pH(nF) + deltapH(nF);
    nF     = abs(deltapH) > pHTol;
    loopc=loopc+1;
 
    if loopc>10000
        Fr=find(abs(deltapH) > pHTol);
        pH(Fr)=NaN;  disp(['pH value did not converge for data on row(s): ' num2str((Fr)')]);
        deltapH=pHTol*0.9;
    end
end
varargout{1}=pH;
end % end nested function

function varargout=CalculateTAfromTCpH(TCi, pHi)
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
% ' SUB CalculateTAfromTCpH, version 02.02, 10-10-97, written by Ernie Lewis.
% ' Inputs: TC, pH, K(), T()
% ' Output: TA
% ' This calculates TA from TC and pH.
K1F=K1(F);     K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
H         = 10.^(-pHi);
CAlk      = TCi.*K1F.*(H + 2.*K2F)./(H.*H + K1F.*H + K1F.*K2F);
BAlk      = TBF.*KBF./(KBF + H);
OH        = KWF./H;
PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
PAlk      = TPF.*PhosTop./PhosBot;
SiAlk     = TSiF.*KSiF./(KSiF + H);
AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
[~,~,pHfree,~] = FindpHOnAllScales(pHi); % this converts pH to pHfree no matter the scale
Hfree = 10.^-pHfree; % this converts pHfree to Hfree
HSO4      = TSF./(1 + KSF./Hfree);% ' since KS is on the free scale
HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
TActemp    = CAlk + BAlk + OH + PAlk + SiAlk + AmmAlk + HSAlk - Hfree - HSO4 - HF;
varargout{1}=TActemp;
end % end nested function

function varargout=CalculatepHfromTCfCO2(TCi, fCO2i)
global K0 K1 K2 F;
% ' SUB CalculatepHfromTCfCO2, version 02.02, 11-12-96, written by Ernie Lewis.
% ' Inputs: TC, fCO2, K0, K1, K2
% ' Output: pH
% ' This calculates pH from TC and fCO2 using K0, K1, and K2 by solving the
% '       quadratic in H: fCO2.*K0 = TC.*H.*H./(K1.*H + H.*H + K1.*K2).
% ' if there is not a real root, then pH is returned as missingn.
RR = K0(F).*fCO2i./TCi;
%       if RR >= 1
%          varargout{1}= missingn;
%          disp('nein!');return;
%       end
% check after sub to see if pH = missingn.
Discr = (K1(F).*RR).*(K1(F).*RR) + 4.*(1 - RR).*(K1(F).*K2(F).*RR);
H     = 0.5.*(K1(F).*RR + sqrt(Discr))./(1 - RR);
%       if (H <= 0)
%           pHctemp = missingn;
%       else
pHctemp = log(H)./log(0.1);
%       end
varargout{1}=pHctemp;
end % end nested function

function varargout=CalculateTCfrompHfCO2(pHi, fCO2i)
global K0 K1 K2 F;
% ' SUB CalculateTCfrompHfCO2, version 01.02, 12-13-96, written by Ernie Lewis.
% ' Inputs: pH, fCO2, K0, K1, K2
% ' Output: TC
% ' This calculates TC from pH and fCO2, using K0, K1, and K2.
H       = 10.^(-pHi);
TCctemp = K0(F).*fCO2i.*(H.*H + K1(F).*H + K1(F).*K2(F))./(H.*H);
varargout{1}=TCctemp;
end % end nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HCO3 input functions, modified from functions written by Ernie Lewis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=CalculateTAfrompHHCO3(pHi, HCO3i)
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
% ' SUB CalculateTAfrompHCO3, version 01.0, 3-19, added by J. Sharp
% ' Inputs: pH, HCO3, K(), T()
% ' Output: TA
% ' This calculates TA from pH and HCO3.
K1F=K1(F);     K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
H         = 10.^(-pHi);
CAlk      = HCO3i.*(2.*K2F./H + 1);
BAlk      = TBF.*KBF./(KBF + H);
OH        = KWF./H;
PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
PAlk      = TPF.*PhosTop./PhosBot;
SiAlk     = TSiF.*KSiF./(KSiF + H);
AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
[~,~,pHfree,~] = FindpHOnAllScales(pHi); % this converts pH to pHfree no matter the scale
Hfree = 10.^-pHfree; % this converts pHfree to Hfree
HSO4      = TSF./(1 + KSF./Hfree);% ' since KS is on the free scale
HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
TActemp     = CAlk + BAlk + OH + PAlk + SiAlk + AmmAlk + HSAlk - Hfree - HSO4 - HF;
varargout{1}=TActemp;
end % end nested function

function varargout=CalculatepHfromTAHCO3(TAi, HCO3i)
global K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
% ' SUB CalculatepHfromTAHCO3, version 01.0, 8-18, added by J. Sharp with
% ' modifications from Denis Pierrot.
% ' Inputs: TA, CO3, K0, K(), T()
% ' Output: pH
%
% ' This calculates pH from TA and CO3 using K1 and K2 by Newton's method.
% ' It tries to solve for the pH at which Residual = 0.
% ' The starting guess is determined by the method introduced by Munhoven
% ' (2013) and extended by Humphreys et al. (2021).
%
% ' This will continue iterating until all values in the vector are
% ' "abs(deltapH) < pHTol"
% ' However, once a given abs(deltapH) is less than pHTol, that pH value
% ' will be locked in. This avoids erroneous contributions to results from
% ' other lines of input.
K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
vl         = sum(F); % vectorlength
% Find initital pH guess using method of Munhoven (2013)
pHGuess    = CalculatepHfromTAHCO3Munhoven(TAi, HCO3i);
ln10       = log(10);
pH         = pHGuess;
pHTol      = 0.0001; % tolerance
deltapH    = pHTol+pH;
loopc=0;
nF=(abs(deltapH) > pHTol);
while any(nF)
    H         = 10.^(-pH);
    CAlk      = HCO3i.*(H+2.*K2F)./H;
    BAlk      = TBF.*KBF./(KBF + H);
    OH        = KWF./H;
    PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
    PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
    PAlk      = TPF.*PhosTop./PhosBot;
    SiAlk     = TSiF.*KSiF./(KSiF + H);
    AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
    HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
    [~,~,pHfree,~] = FindpHOnAllScales(pH); % this converts pH to pHfree no matter the scale
    Hfree = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree); %' since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
    Residual  = TAi - CAlk - BAlk - OH - PAlk - SiAlk - AmmAlk - HSAlk + Hfree + HSO4 + HF;
    % '               find Slope dTA/dpH
    % '               (this is not exact, but keeps all important terms):
    Slope = ln10 .* (2 .* HCO3i .* K2F ./ H + BAlk .* H ./ (KBF + H) + OH + H);
    deltapH   = Residual./Slope; %' this is Newton's method
    % ' to keep the jump from being too big:
    while any(abs(deltapH) > 1)
        FF=abs(deltapH)>1; deltapH(FF)=deltapH(FF)./2;
    end
    pH(nF) = pH(nF) + deltapH(nF);
    nF     = abs(deltapH) > pHTol;
    loopc=loopc+1;
 
    if loopc>10000
        Fr=find(abs(deltapH) > pHTol);
        pH(Fr)=NaN;  disp(['pH value did not converge for data on row(s): ' num2str((Fr)')]);
        deltapH=pHTol*0.9;
    end
end
varargout{1}=pH;
end % end nested function

function varargout=CalculatepHfromTCHCO3(TCi, HCO3i)
global K1 K2 F;
% ' SUB CalculatepHfromTCHCO3, version 01.0, 3-19, added by J. Sharp
% ' Inputs: TC, HCO3, K0, K1, K2
% ' Output: pH
% ' This calculates pH from TC and HCO3 using K1 and K2 by solving the
% '       quadratic in H: TC = HCO3i.*(H./K1 + 1 + K2./H).
% '       Therefore:      0  = H.*H./K1 + (1-TC/HCO3i).*H + K2.
% ' if there is not a real root, then pH is returned as missingn.
RR = TCi./HCO3i;
%       if RR >= 1
%          varargout{1}= missingn;
%          disp('nein!');return;
%       end
% check after sub to see if pH = missingn.
Discr = ((1-RR).*(1-RR) - 4.*(1./(K1(F))).*(K2(F)));
H     = 0.5.*((-(1-RR)) - sqrt(Discr))./(1./(K1(F))); % Subtraction
%       if (H <= 0)
%           pHctemp = missingn;
%       else
pHctemp = log(H)./log(0.1);
%       end
varargout{1}=pHctemp;
end % end nested function

function varargout=CalculatepHfromfCO2HCO3(fCO2i, HCO3i)
global K0 K1 F;
% ' SUB CalculatepHfromfCO2HCO3, version 01.0, 3-19, added by J. Sharp
% ' Inputs: fCO2, HCO3, K0, K1, K2
% ' Output: pH
% ' This calculates pH from fCO2 and HCO3, using K0, K1, and K2.
H            = (fCO2i.*K0(F).*K1(F))./HCO3i;  % removed incorrect (F) index from HCO3i // MPH
pHx          = -log10(H);
varargout{1} = pHx;
end % end nested function

function varargout=CalculatepHfCO2fromTCHCO3(TCx, HCO3x)
% Outputs pH fCO2, in that order
% SUB CalculatepHfCO2fromTCHCO3, version 01.0, 3-19, added by J. Sharp
% Inputs: pHScale%, WhichKs%, WhoseKSO4%, TC, HCO3, Sal, K(), T(), TempC, Pdbar
% Outputs: pH, fCO2
% This calculates pH and fCO2 from TC and HCO3 at output conditions.
pHx   = CalculatepHfromTCHCO3(TCx, HCO3x); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
fCO2x = CalculatefCO2fromTCpH(TCx, pHx);
varargout{1} = pHx;
varargout{2} = fCO2x;
end % end nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CO3 input functions, modified from functions written by Ernie Lewis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=CalculateTAfrompHCO3(pHi, CO3i)
global K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
% ' SUB CalculateTAfrompHCO3, version 01.0, 8-18, added by J. Sharp
% ' Inputs: pH, CO3, K(), T()
% ' Output: TA
% ' This calculates TA from pH and CO3.
K1F=K1(F);     K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
H         = 10.^(-pHi);
CAlk      = CO3i.*(H./K2F + 2);
BAlk      = TBF.*KBF./(KBF + H);
OH        = KWF./H;
PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
PAlk      = TPF.*PhosTop./PhosBot;
SiAlk     = TSiF.*KSiF./(KSiF + H);
AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
[~,~,pHfree,~] = FindpHOnAllScales(pHi); % this converts pH to pHfree no matter the scale
Hfree = 10.^-pHfree; % this converts pHfree to Hfree
HSO4      = TSF./(1 + KSF./Hfree);% ' since KS is on the free scale
HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
TActemp     = CAlk + BAlk + OH + PAlk + SiAlk + AmmAlk + HSAlk - Hfree - HSO4 - HF;
varargout{1}=TActemp;
end % end nested function

function varargout=CalculatepHfromTACO3(TAi, CO3i)
global K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F
% ' SUB CalculatepHfromTACO3, version 01.0, 8-18, added by J. Sharp with
% ' modifications from Denis Pierrot.
% ' Inputs: TA, CO3, K0, K(), T()
% ' Output: pH
%
% ' This calculates pH from TA and CO3 using K1 and K2 by Newton's method.
% ' It tries to solve for the pH at which Residual = 0.
% ' The starting guess is determined by the method introduced by Munhoven
% ' (2013) and extended by Humphreys et al. (2021).
%
% ' This will continue iterating until all values in the vector are
% ' "abs(deltapH) < pHTol"
% ' However, once a given abs(deltapH) is less than pHTol, that pH value
% ' will be locked in. This avoids erroneous contributions to results from
% ' other lines of input.
K2F=K2(F);     KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);
vl         = sum(F); % vectorlength
% Find initital pH guess using method of Munhoven (2013)
pHGuess    = CalculatepHfromTACO3Munhoven(TAi, CO3i);
ln10       = log(10);
pH         = pHGuess;
pHTol      = 0.0001; % tolerance
deltapH    = pHTol+pH;
loopc=0;
nF=(abs(deltapH) > pHTol);
while any(nF)
    H         = 10.^(-pH);
    CAlk      = CO3i.*(H+2.*K2F)./K2F;
    BAlk      = TBF.*KBF./(KBF + H);
    OH        = KWF./H;
    PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
    PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
    PAlk      = TPF.*PhosTop./PhosBot;
    SiAlk     = TSiF.*KSiF./(KSiF + H);
    AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
    HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
    [~,~,pHfree,~] = FindpHOnAllScales(pH); % this converts pH to pHfree no matter the scale
    Hfree = 10.^-pHfree; % this converts pHfree to Hfree
    HSO4      = TSF./(1 + KSF./Hfree); %' since KS is on the free scale
    HF        = TFF./(1 + KFF./Hfree);% ' since KF is on the free scale
    Residual  = TAi - CAlk - BAlk - OH - PAlk - SiAlk - AmmAlk - HSAlk + Hfree + HSO4 + HF;
    % '               find Slope dTA/dpH
    % '               (this is not exact, but keeps all important terms):
    Slope = ln10 .* (-CO3i .* H ./ K2F + BAlk .* H ./ (KBF + H) + OH + H);
    deltapH   = Residual./Slope; %' this is Newton's method
    % ' to keep the jump from being too big:
    while any(abs(deltapH) > 1)
        FF=abs(deltapH)>1; deltapH(FF)=deltapH(FF)./2;
    end
    pH(nF) = pH(nF) + deltapH(nF);
    nF     = abs(deltapH) > pHTol;
    loopc=loopc+1;
 
    if loopc>10000 
        Fr=find(abs(deltapH) > pHTol);
        pH(Fr)=NaN;  disp(['pH value did not converge for data on row(s): ' num2str((Fr)')]);
        deltapH=pHTol*0.9;
    end
end
varargout{1}=pH;
end % end nested function

function varargout=CalculatepHfromTCCO3(TCi, CO3i)
global K1 K2 F;
% ' SUB CalculatepHfromTCCO3, version 01.0, 8-18, added by J. Sharp
% ' Inputs: TC, CO3, K0, K1, K2
% ' Output: pH
% ' This calculates pH from TC and CO3 using K1 and K2 by solving the
% '       quadratic in H: TC = CO3i.*(H.*H/(K1.*K2) + H./K2 + 1).
% '       Therefore:      0  = H.*H/(K1.*K2) + H./K2 + (1-TC./CO3i).
% ' if there is not a real root, then pH is returned as missingn.
RR = TCi./CO3i;
%       if RR >= 1
%          varargout{1}= missingn;
%          disp('nein!');return;
%       end
% check after sub to see if pH = missingn.
Discr = ((1./K2(F)).*(1./K2(F)) - 4.*(1./(K1(F).*K2(F))).*(1-RR));
H     = 0.5.*((-1./K2(F)) + sqrt(Discr))./(1./(K1(F).*K2(F))); % Addition
%       if (H <= 0)
%           pHctemp = missingn;
%       else
pHctemp = log(H)./log(0.1);
%       end
varargout{1}=pHctemp;
end % end nested function

function varargout=CalculatepHfromfCO2CO3(fCO2i, CO3i)
global K0 K1 K2 F;
% ' SUB CalculatepHfromfCO2CO3, version 01.0, 8-18, added by J. Sharp
% ' Inputs: fCO2, CO3, K0, K1, K2
% ' Output: pH
% ' This calculates pH from fCO2 and CO3, using K0, K1, and K2.
H            = sqrt((fCO2i.*K0(F).*K1(F).*K2(F))./CO3i);    % removed incorrect (F) index from CO3i // MPH
pHx          = -log10(H);
varargout{1} = pHx;
end % end nested function

function varargout=CalculatepHfCO2fromTCCO3(TCx, CO3x)
% Outputs pH fCO2, in that order
% SUB CalculatepHfCO2fromTCCO3, version 01.0, 8-18, added by J. Sharp
% Inputs: pHScale%, WhichKs%, WhoseKSO4%, TC, CO3, Sal, K(), T(), TempC, Pdbar
% Outputs: pH, fCO2
% This calculates pH and fCO2 from TC and CO3 at output conditions.
pHx   = CalculatepHfromTCCO3(TCx, CO3x); % pH is returned on the scale requested in "pHscale" (see 'constants'...)
fCO2x = CalculatefCO2fromTCpH(TCx, pHx);
varargout{1} = pHx;
varargout{2} = fCO2x;
end % end nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CO3/HCO3 input function, modified from functions written by Ernie Lewis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=CalculatepHfromCO3HCO3(CO3x, HCO3x)
global K2 F
% ' SUB CalculatepHfromCO3HCO3, version 01.0, 3-19, added by J. Sharp
% ' Inputs: CO3, HCO3, K2
% ' Output: pH
% ' This calculates fCO2 from TC and pH, using K2.
H            = HCO3x.*K2(F)./CO3x;
pHx          = -log10(H);
varargout{1} = pHx;
end % end nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CO3/HCO3/CO2 output functions, modified from functions written by Ernie Lewis.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=CalculateCO3HCO3fromTCpH(TCx, pHx)
global K1 K2 F
% ' SUB CalculateCO3HCO3CO2fromTCpH, version 01.0, 3-19, added by J. Sharp
% ' Inputs: TC, pH, K1, K2
% ' Output: CO3, HCO3, CO2
% ' This calculates CO3, HCO3, and CO2 from TC and pH, using K1, and K2.
H            = 10.^(-pHx);
CO3x         = TCx.*K1(F).*K2(F)./(K1(F).*H + H.*H + K1(F).*K2(F));
HCO3x        = TCx.*K1(F).*H./(K1(F).*H + H.*H + K1(F).*K2(F));
varargout{1} = CO3x;
varargout{2} = HCO3x;
end % end nested function

function varargout=CalculateCO3fromTCpH(TCx, pHx)
global K1 K2 F
% ' SUB CalculateCO3CO2fromTCpH, version 01.0, 3-19, added by J. Sharp
% ' Inputs: TC, pH, K1, K2
% ' Output: CO3, CO2
% ' This calculates CO3 and CO2 from TC and pH, using K1, and K2.
H            = 10.^(-pHx);
CO3x         = TCx.*K1(F).*K2(F)./(K1(F).*H + H.*H + K1(F).*K2(F));
varargout{1} = CO3x;
end % end nested function

function varargout=CalculateHCO3fromTCpH(TCx, pHx)
global K1 K2 F
% ' SUB CalculateHCO3CO2fromTCpH, version 01.0, 3-19, added by J. Sharp
% ' Inputs: TC, pH, K1, K2
% ' Output: HCO3, CO2
% ' This calculates HCO3 and CO2 from TC and pH, using K1, and K2.
H            = 10.^(-pHx);
HCO3x        = TCx.*K1(F).*H./(K1(F).*H + H.*H + K1(F).*K2(F));
varargout{1} = HCO3x;
end % end nested function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial pH estimates via Munhoven (2013), Humphreys et al (2021).
% Added by J. Sharp (3-30-21)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=CalculatepHfromTATCMunhoven(TAi, TCi)
global K1 K2 KB TB F;
K1F=K1(F);     K2F=K2(F);     TBF =TB(F);    KBF=KB(F);
g0 = K1F.*K2F.*KBF.*(1-(2.*TCi+TBF)./TAi);
g1 = K1F.*(KBF.*(1-TBF./TAi-TCi./TAi)+K2F.*(1-2.*TCi./TAi));
g2 = KBF.*(1-TBF./TAi)+K1F.*(1-TCi./TAi);
% Determine g21min
g21min = g2.^2-3.*g1;
g21min_positive = g21min > 0;
sq21 = nan(size(TAi,1),1);
sq21(g21min_positive) = sqrt(g21min(g21min_positive));
sq21(~g21min_positive) = 0;
% Determine Hmin
Hmin = nan(size(TAi,1),1);
g2_positive = g2 >=0;
Hmin(~g2_positive) = (-g2(~g2_positive) + sq21(~g2_positive))./3;
Hmin(g2_positive) = -g1(g2_positive)./(g2(g2_positive) + sq21(g2_positive));
% Calculate initial pH
pHGuess = nan(size(TAi,1),1);
idx = TAi <= 0;
pHGuess(idx) = -log10(1e-3);
idx = TAi > 0 & TAi < 2.*TCi + TBF;
pHGuess(idx & g21min_positive) = ...
    -log10(Hmin(idx & g21min_positive) + ...
    sqrt(-(Hmin(idx & g21min_positive).^3 + g2(idx & g21min_positive).*Hmin(idx & g21min_positive).^2 + ...
    g1(idx & g21min_positive).*Hmin(idx & g21min_positive) + ...
    g0(idx & g21min_positive))./sq21(idx & g21min_positive)));
pHGuess(idx & ~g21min_positive) = -log10(1e-7);
idx = TAi >= 2.*TCi + TBF;
pHGuess(idx) = -log10(1e-10);
varargout{1}=pHGuess;
end

function varargout=CalculatepHfromTACO2Munhoven(TAi, CO2x)
global K1 K2 KB TB F;
K1F=K1(F);     K2F=K2(F);     TBF =TB(F);    KBF=KB(F);
g0 = -2.*K1F.*K2F.*KBF.*CO2x./TAi;
g1 = -K1F.*(2.*K2F.*CO2x+KBF.*CO2x)./TAi;
g2 = KBF-(TBF.*KBF+K1F.*CO2x)./TAi;
% Determine Hmin
g21min = g2.^2-3.*g1;
g21min_positive = g21min > 0;
sq21 = nan(size(TAi,1),1);
sq21(g21min_positive) = sqrt(g21min(g21min_positive));
sq21(~g21min_positive) = 0;
Hmin = nan(size(TAi,1),1);
g2_positive = g2 >=0;
Hmin(~g2_positive) = (-g2(~g2_positive) + sq21(~g2_positive))./3;
Hmin(g2_positive) = -g1(g2_positive)./(g2(g2_positive) + sq21(g2_positive));
% Calculate initial pH
pHGuess = nan(size(TAi,1),1);
idx = TAi <= 0;
pHGuess(idx) = -log10(1e-3);
idx = TAi > 0;
pHGuess(idx & g21min_positive) = ...
    -log10(Hmin(idx & g21min_positive) + ...
    sqrt(-(Hmin(idx & g21min_positive).^3 + g2(idx & g21min_positive).*Hmin(idx & g21min_positive).^2 + ...
    g1(idx & g21min_positive).*Hmin(idx & g21min_positive)+...
    g0(idx & g21min_positive))./sq21(idx & g21min_positive)));
pHGuess(idx & ~g21min_positive) = -log10(1e-7);
varargout{1}=pHGuess;
end

function varargout=CalculatepHfromTAHCO3Munhoven(TAi, HCO3x)
global K1 K2 KB TB F;
K1F=K1(F);     K2F=K2(F);     TBF =TB(F);    KBF=KB(F);
g0 = 2.*K2F.*KBF.*HCO3x;
g1 = KBF.*(HCO3x+TBF-TAi)+2.*K2F.*HCO3x;
g2 = HCO3x-TAi;
% Calculate initial pH
pHGuess = nan(size(TAi,1),1);
idx = TAi <= HCO3x;
pHGuess(idx) = -log10(1e-3);
idx = TAi > HCO3x;
pHGuess(idx) = ...
    -log10((-g1(idx)-sqrt(g1(idx).^2-4.*g0(idx).*g2(idx)))./(2.*g2(idx)));
varargout{1}=pHGuess;
end

function varargout=CalculatepHfromTACO3Munhoven(TAi, CO3x)
global K1 K2 KB TB F;
K1F=K1(F);     K2F=K2(F);     TBF =TB(F);    KBF=KB(F);
g0 = K2F.*KBF.*(2.*CO3x+TBF-TAi);
g1 = KBF.*CO3x+K2F.*(2.*CO3x-TAi);
g2 = CO3x;
% Calculate initial pH
pHGuess = nan(size(TAi,1),1);
idx = TAi <= 2.*CO3x+TBF;
pHGuess(idx) = -log10(1e-3);
idx = TAi > 2.*CO3x+TBF;
pHGuess(idx) = ...
    -log10((-g1(idx)+sqrt(g1(idx).^2-4.*g0(idx).*g2(idx)))./(2.*g2(idx)));
varargout{1}=pHGuess;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=RevelleFactor(TAi, TCi)
% global WhichKs;
% ' SUB RevelleFactor, version 01.03, 01-07-97, written by Ernie Lewis.
% ' Inputs: WhichKs%, TA, TC, K0, K(), T()
% ' Outputs: Revelle
% ' This calculates the Revelle factor (dfCO2/dTC)|TA/(fCO2/TC).
% ' It only makes sense to talk about it at pTot = 1 atm, but it is computed
% '       here at the given K(), which may be at pressure <> 1 atm. Care must
% '       thus be used to see if there is any validity to the number computed.
TC0 = TCi;
dTC = 0.00000001;% ' 0.01 umol/kg-SW (lower than prior versions of CO2SYS)
% ' Find fCO2 at TA, TC + dTC
TCi = TC0 + dTC;
pHc= CalculatepHfromTATC(TAi, TCi);
fCO2c= CalculatefCO2fromTCpH(TCi, pHc);
fCO2plus = fCO2c;
% ' Find fCO2 at TA, TC - dTC
TCi = TC0 - dTC;
pHc= CalculatepHfromTATC(TAi, TCi);
fCO2c= CalculatefCO2fromTCpH(TCi, pHc);
fCO2minus = fCO2c;
% CalculateRevelleFactor:
Revelle = (fCO2plus - fCO2minus)./dTC./((fCO2plus + fCO2minus)./TC0); % Corrected error pointed out by MP Humphreys (https://pyco2sys.readthedocs.io/en/latest/validate/)
varargout{1}=Revelle; % keyboard
end % end nested function


function varargout=CalculateAlkParts(pH)
global K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S;
global TB TF TS TP TSi TNH4 TH2S F;
% ' SUB CalculateAlkParts, version 01.03, 10-10-97, written by Ernie Lewis.
% ' Inputs: pH, TC, K(), T()
% ' Outputs: BAlk, OH, PAlk, SiAlk, Hfree, HSO4, HF
% ' This calculates the various contributions to the alkalinity.
% ' Though it is coded for H on the total pH scale, for the pH values occuring
% ' in seawater (pH > 6) it will be equally valid on any pH scale (H terms
% ' negligible) as long as the K Constants are on that scale.

KWF =KW(F);
KP1F=KP1(F);   KP2F=KP2(F);   KP3F=KP3(F);   TPF=TP(F);
TSiF=TSi(F);   KSiF=KSi(F);   TNH4F=TNH4(F); KNH4F=KNH4(F);
TH2SF=TH2S(F); KH2SF=KH2S(F); TBF =TB(F);    KBF=KB(F);
TSF =TS(F);    KSF =KS(F);    TFF =TF(F);    KFF=KF(F);

H         = 10.^(-pH);
BAlk      = TBF.*KBF./(KBF + H);
OH        = KWF./H;
PhosTop   = KP1F.*KP2F.*H + 2.*KP1F.*KP2F.*KP3F - H.*H.*H;
PhosBot   = H.*H.*H + KP1F.*H.*H + KP1F.*KP2F.*H + KP1F.*KP2F.*KP3F;
PAlk      = TPF.*PhosTop./PhosBot;
SiAlk     = TSiF.*KSiF./(KSiF + H);
AmmAlk    = TNH4F.*KNH4F./(KNH4F + H);
HSAlk     = TH2SF.*KH2SF./(KH2SF + H);
[~,~,pHfree,~] = FindpHOnAllScales(pH); % this converts pH to pHfree no matter the scale
Hfree = 10.^-pHfree; % this converts pHfree to Hfree
HSO4      = TSF./(1 + KSF./Hfree); %' since KS is on the free scale
HF        = TFF./(1 + KFF./Hfree); %' since KF is on the free scale

varargout{1} = BAlk;  varargout{2} = OH; varargout{3} = PAlk;
varargout{4} = SiAlk; varargout{5} = AmmAlk; varargout{6} = HSAlk;
varargout{7} = Hfree; varargout{8} = HSO4; varargout{9} = HF;
end % end nested function


function varargout=CaSolubility(Sal, TempC, Pdbar, TC, pH)
global K1 K2 TempK logTempK sqrSal Pbar RT WhichKs CAL ntps F
global PertK    % Id of perturbed K
global Perturb  % perturbation
% '***********************************************************************
% ' SUB CaSolubility, version 01.05, 05-23-97, written by Ernie Lewis.
% ' Inputs: WhichKs%, Sal, TempCi, Pdbari, TCi, pHi, K1, K2
% ' Outputs: OmegaCa, OmegaAr
% ' This calculates omega, the solubility ratio, for calcite and aragonite.
% ' This is defined by: Omega = [CO3--]*[Ca++]./Ksp,
% '       where Ksp is the solubility product (either KCa or KAr).
% '***********************************************************************
% ' These are from:
% ' Mucci, Alphonso, The solubility of calcite and aragonite in seawater
% '       at various salinities, temperatures, and one atmosphere total
% '       pressure, American Journal of Science 283:781-799, 1983.
% ' Ingle, S. E., Solubility of calcite in the ocean,
% '       Marine Chemistry 3:301-319, 1975,
% ' Millero, Frank, The thermodynamics of the carbonate system in seawater,
% '       Geochemica et Cosmochemica Acta 43:1651-1661, 1979.
% ' Ingle et al, The solubility of calcite in seawater at atmospheric pressure
% '       and 35%o salinity, Marine Chemistry 1:295-307, 1973.
% ' Berner, R. A., The solubility of calcite and aragonite in seawater in
% '       atmospheric pressure and 34.5%o salinity, American Journal of
% '       Science 276:713-730, 1976.
% ' Takahashi et al, in GEOSECS Pacific Expedition, v. 3, 1982.
% ' Culberson, C. H. and Pytkowicz, R. M., Effect of pressure on carbonic acid,
% '       boric acid, and the pHi of seawater, Limnology and Oceanography
% '       13:403-417, 1968.
% '***********************************************************************
Ca=CAL(F);
Ar=nan(sum(F),1);
KCa=nan(sum(F),1);
KAr=nan(sum(F),1);
TempKx=TempK(F);
logTempKx=logTempK(F);
sqrSalx=sqrSal(F);
Pbarx=Pbar(F);
RTx=RT(F);
FF=(WhichKs(F)~=6 & WhichKs(F)~=7);
if any(FF)
% (below here, F isn't used, since almost always all rows match the above criterium,
%  in all other cases the rows will be overwritten later on).
    % CalciteSolubility:
    % '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
    logKCa = -171.9065 - 0.077993.*TempKx(FF) + 2839.319./TempKx(FF);
    logKCa = logKCa + 71.595.*logTempKx(FF)./log(10);
    logKCa = logKCa + (-0.77712 + 0.0028426.*TempKx(FF) + 178.34./TempKx(FF)).*sqrSalx(FF);
    logKCa = logKCa - 0.07711.*Sal(FF) + 0.0041249.*sqrSalx(FF).*Sal(FF);
    % '       sd fit = .01 (for Sal part, not part independent of Sal)
    KCa(FF) = 10.^(logKCa);% ' this is in (mol/kg-SW)^2
    % AragoniteSolubility:
    % '       Mucci, Alphonso, Amer. J. of Science 283:781-799, 1983.
    logKAr = -171.945 - 0.077993.*TempKx(FF) + 2903.293./TempKx(FF);
    logKAr = logKAr + 71.595.*logTempKx(FF)./log(10);
    logKAr = logKAr + (-0.068393 + 0.0017276.*TempKx(FF) + 88.135./TempKx(FF)).*sqrSalx(FF);
    logKAr = logKAr - 0.10018.*Sal(FF) + 0.0059415.*sqrSalx(FF).*Sal(FF);
    % '       sd fit = .009 (for Sal part, not part independent of Sal)
    KAr(FF)    = 10.^(logKAr);% ' this is in (mol/kg-SW)^2
    % PressureCorrectionForCalcite:
    % '       Ingle, Marine Chemistry 3:301-319, 1975
    % '       same as in Millero, GCA 43:1651-1661, 1979, but Millero, GCA 1995
    % '       has typos (-.5304, -.3692, and 10^3 for Kappa factor)
    deltaVKCa = -48.76 + 0.5304.*TempC(FF);
    KappaKCa  = (-11.76 + 0.3692.*TempC(FF))./1000;
    lnKCafac  = (-deltaVKCa + 0.5.*KappaKCa.*Pbarx(FF)).*Pbarx(FF)./RTx(FF);
    KCa(FF)       = KCa(FF).*exp(lnKCafac);
    % PressureCorrectionForAragonite:
    % '       Millero, Geochemica et Cosmochemica Acta 43:1651-1661, 1979,
    % '       same as Millero, GCA 1995 except for typos (-.5304, -.3692,
    % '       and 10^3 for Kappa factor)
    deltaVKAr = deltaVKCa + 2.8;
    KappaKAr  = KappaKCa;
    lnKArfac  = (-deltaVKAr + 0.5.*KappaKAr.*Pbarx(FF)).*Pbarx(FF)./RTx(FF);
    KAr(FF)       = KAr(FF).*exp(lnKArfac);
end
FF=(WhichKs(F)==6 | WhichKs(F)==7);
if any(FF)
    % *** CalculateKCaforGEOSECS:
    % Ingle et al, Marine Chemistry 1:295-307, 1973 is referenced in
    % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
    % but the fit is actually from Ingle, Marine Chemistry 3:301-319, 1975)
    KCa(FF) = 0.0000001.*(-34.452 - 39.866.*Sal(FF).^(1./3) +...
        110.21.*log(Sal(FF))./log(10) - 0.0000075752.*TempKx(FF).^2);
    % this is in (mol/kg-SW)^2
    %
    % *** CalculateKArforGEOSECS:
    % Berner, R. A., American Journal of Science 276:713-730, 1976:
    % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
    KAr(FF) = 1.45.*KCa(FF);% ' this is in (mol/kg-SW)^2
    % Berner (p. 722) states that he uses 1.48.
    % It appears that 1.45 was used in the GEOSECS calculations
    %
    % *** CalculatePressureEffectsOnKCaKArGEOSECS:
    % Culberson and Pytkowicz, Limnology and Oceanography 13:403-417, 1968
    % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982
    % but their paper is not even on this topic).
    % The fits appears to be new in the GEOSECS report.
    % I can't find them anywhere else.
    KCa(FF) = KCa(FF).*exp((36   - 0.2 .*TempC(FF)).*Pbarx(FF)./RTx(FF));
    KAr(FF) = KAr(FF).*exp((33.3 - 0.22.*TempC(FF)).*Pbarx(FF)./RTx(FF));
end
% Added by JM Epitalon
% For computing derivative with respect to KCa or KAr, one has to perturb the value of one K
% Requested perturbation is passed through global variables PertK and Perturb
if (~ isempty(PertK))
    switch PertK
        case {'KSPA'}   % solubility Product for Aragonite
            KAr = KAr + Perturb;
        case {'KSPC'}   % for Calcite
            KCa = KCa + Perturb;
        case {'CAL'}   % for calcium concentration
            Ca  = Ca  + Perturb;
    end
end

% CalculateOmegasHere:
H = 10.^(-pH);
CO3 = TC.*K1(F).*K2(F)./(K1(F).*H + H.*H + K1(F).*K2(F));
varargout{1} = CO3.*Ca./KCa; % OmegaCa, dimensionless
varargout{2} = CO3.*Ca./KAr; % OmegaAr, dimensionless
varargout{3} = CAL ;
end % end nested function

function varargout=FindpHOnAllScales(pH)
global pHScale K T TS KS TF KF fH F ntps;
% ' SUB FindpHOnAllScales, version 01.02, 01-08-97, written by Ernie Lewis.
% ' Inputs: pHScale%, pH, K(), T(), fH
% ' Outputs: pHNBS, pHfree, pHTot, pHSWS
% ' This takes the pH on the given scale and finds the pH on all scales.
%  TS = T(3); TF = T(2);
%  KS = K(6); KF = K(5);% 'these are at the given T, S, P
TSx=TS(F); KSx=KS(F); TFx=TF(F); KFx=KF(F);fHx=fH(F);
FREEtoTOT = (1 + TSx./KSx); % ' pH scale conversion factor
SWStoTOT  = (1 + TSx./KSx)./(1 + TSx./KSx + TFx./KFx);% ' pH scale conversion factor
% keyboard
factor=nan(sum(F),1);
nF=pHScale(F)==1;  %'"pHtot"
factor(nF) = 0;
nF=pHScale(F)==2; % '"pHsws"
factor(nF) = -log(SWStoTOT(nF))./log(0.1);
nF=pHScale(F)==3; % '"pHfree"
factor(nF) = -log(FREEtoTOT(nF))./log(0.1);
nF=pHScale(F)==4;  %'"pHNBS"
factor(nF) = -log(SWStoTOT(nF))./log(0.1) + log(fHx(nF))./log(0.1);
pHtot  = pH    - factor;    % ' pH comes into this sub on the given scale
pHNBS  = pHtot - log(SWStoTOT) ./log(0.1) + log(fHx)./log(0.1);
pHfree = pHtot - log(FREEtoTOT)./log(0.1);
pHsws  = pHtot - log(SWStoTOT) ./log(0.1);
varargout{1}=pHtot;
varargout{2}=pHsws;
varargout{3}=pHfree;
varargout{4}=pHNBS;
end % end nested function





%
% errors for CO2SYSv3 -----------------------------------------------------
%

% errors()
% This subroutine propagates uncertainties for the marine carbonate chemistry calculations
% from errors (or uncertainties) on inputs:
%  - pair of carbonate system variables
%  - nutrient (silicate and phosphate) concentrations
%  - concentrations of ammonium and hydrogen sulfide
%  - temperature and salinity
%  - calcium concentration (optional)
% plus errors in dissociation constants pK0, pK1, pK2, pKb, pKw, pKspa, and pKspc as well as total boron
%
% It calls derivnum, which computes numerical derivatives, and then
% it applies error propagation using the method of moments.
% The latter is a general technique to estimate the 2nd moment of a variable z
% (variance or standard deviation) based on a 1st-order approximation to z.
%
% This subroutine has been modified from its original version (Orr et al.
% 2018) to allow for input variables of CO2, HCO3, and CO3, the inclusion
% of ammonium and hydrogen sulfide, compatibility with CO2SYS.m(v3), and am
% optional input of calcium concentration uncertainty.
%
%**************************************************************************
%
%  **** SYNTAX:
%  [err, headers, units] = errors(PAR1,PAR2,PAR1TYPE,PAR2TYPE,...
%                                 SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,...
%                                 NH4,H2S,ePAR1,ePAR2,eSAL,eTEMP,eSI,ePO4,...
%                                 eNH4,eH2S,epK,eBt,r,pHSCALEIN,K1K2CONSTANTS,...
%                                 KSO4CONSTANT,KFCONSTANT,BORON,eCAL(optional))
% 
%  **** SYNTAX EXAMPLES:
%  [Result]                = errors(2400,2200,1,2,35,10,10,0,0,15,1,0,0,2,2,0.01,0.01,0,0,0,0,0,0,0,1,4,1,1,1)
%  [Result,Headers]        = errors(2400,   8,1,3,35,25,5,0,3000,15,1,0,0,2,0.001,0,0,0,0,0,0,0,0,0,1,4,1,1,1)
%  [Result,Headers,Units]  = errors(500,    8,5,3,35,25,5,0,4000,15,1,0,0,2,0.001,0,0,0,0,0,0,'','',0,1,4,1,1,1)
%  [A]                     = errors(2400,2000:10:2400,1,2,35,10,10,0,0,15,2,0,0,2,0,0,0,0,0,0,0,'','',0,1,4,1,1,1)
%  [A]                     = errors(2400,2200,1,2,0:1:35,0,25,4200,0,15,1,0,0,2,2,0,0,0,0,0,0,'','',0,1,4,1,1,1)
%  epK = [0.002, 0.0075, 0.015, 0.01, 0.01, 0.02, 0.02];
%  eBt = 0.02;
%  [A, hdr, units]         = errors(2400,2200,1,2,35,0,25,0:100:4200,0,15,1,0,0,2,2,0,0,0,0,0,0,epK,eBt,0,1,4,1,1,1)
%  
%**************************************************************************
%
% INPUT:
%
%   - ePAR1, ePAR2   :  uncertainty of PAR1 and PAR2 of input pair of CO2 system variables
%                       * Same units as PAR1 & PAR2
%   - eS, eT         :  uncertainty of Salinity and Temperature (same units as S and T)
%   - ePO4, eSI      :  uncertainty of Phosphate and Silicate total concentrations (same units as PO4 and SI [umol/kg])
%   - eNH4, eH2S     :  uncertainty of Ammonia and Hydrogen Sulfide total concentrations (same units as NH4 and H2S [umol/kg])
%   - epK            :  uncertainty of all seven dissociation constants (a vector) [pK units]
%   - eBt            :  uncertainty of total boron, given as fractional relative error (eBt=0.02 is a 2% error)
%   - r              :  correlation coefficient between PAR1 AND PAR2 (typicaly 0)
%   - others         :  same as input for subroutine  CO2SYS() (version 3, scalar or vectors)
%   - eCAL           :  uncertainty of calcium [mmol/kg] determined by ratio with salinity
%
% All parameters may be scalars or vectors except epK and eBt.
%   * epK must be vector of 7 values : errors of [pK0, pK1, pK2, pKb, pKw, pKspa, pKspc]. 
%     These errors are assumed to be the same for all rows of data.
%     These 7 values are in pK units
%
%     if epK is empty (= ''), this routine specifies default values.
%     These default standard errors are :
%        pK0   :  0.002 
%        pK1   :  0.0075
%        pK2   :  0.015
%        pKb   :  0.01    boric acid
%        pKw   :  0.01    water dissociation
%        pKspa :  0.02    solubility product of Aragonite 
%        pKspc :  0.02    solubility product of Calcite
%
%   * eBt is a scalar real number, fractional relative error (between 0.00 and 1.00)
%     for TB, where the default is eBt=0.02. It is assumed to be the same
%     for all rows of data.
%
% In constrast, ePAR1, ePAR2, eS, eT, ePO4, eSI, eNH4, and eH2S
%   - if vectors, are errors associated with each data point
%   - if scalars, are one error value associated to all data points
% The same for parameter 'r'.
%
% If no value is input for eCAL, it will not be evaluated.
%
% If 'r' is nonzero with a value between -1.0 and 1.0, it indicates the correlation 
% between uncertainties of the input pair of carbonate system variables.
% By default, 'r' is zero. However, for some pairs the user may want to specify a
% different value. For example, measurements of pCO2 and pH are often anti-correlated.
% The same goes for two other pairs: 'CO2 and CO3' and 'pCO2 and
% CO3'. But even for these cases, care is needed when using non-zero values of 'r'.
% 
% When the user propagates errors for an individual
% measurement, 'r' should ALWAYS be zero if each member of the input pair is
% measured independently. In this case, we are interested in the
% correlation between the uncertainties in those measurements, not in
% the correlation between the measurments themselves. Uncertainties from
% those measurements are probably not correlated if they come from
% different instruments. Conversely, if users are interested in the
% error in the mean of a distribution of measurements (i.e., if they are
% propagating standard errors instead of standard deviations), one
% should then also account for the correlation between the measurements of
% the two variables of the input pair.
% 
% For input pairs where one member is pH, this 'errors' routine automatically
% inverses the sign of 'r'.
% That inversion is done because the associated derivatives are computed in terms of 
% the hydrogen ion concentration H+, not pH. Therefore for each of these 6
% flags, if the user wants to compute 'r' that should be done (1) using
% the H+ concentration instead of pH, and (2) the sign of that computed 'r'
% should be inversed when passing it as an argument to this routine.
% To express perfect anticorrelation with pH, the user should 
% use 'r=-1.0'. 
% 
%**************************************************************************
%
% OUTPUT: * an array containing uncertainty for the following variables
%           (one row per sample):
%         *  a cell-array containing crudely formatted headers
%
%    POS  PARAMETER         UNIT
%
%    01 - TAlk              (umol/kgSW)
%    02 - TCO2              (umol/kgSW)
%    03 - [H+] in           (nmol/kgSW)
%    04 - pCO2 in           (uatm)
%    05 - fCO2 in           (uatm)
%    06 - HCO3 in           (umol/kgSW)
%    07 - CO3 in            (umol/kgSW)
%    08 - CO2 in            (umol/kgSW)
%    09 - RF in             ()
%    10 - OmegaCa in        ()
%    11 - OmegaAr in        ()
%    12 - xCO2 in           (ppm)
%    13 - [H+] out          (nmol/kgSW)
%    14 - pCO2 out          (uatm)
%    15 - fCO2 out          (uatm)
%    16 - HCO3 out          (umol/kgSW)
%    17 - CO3 out           (umol/kgSW)
%    18 - CO2 out           (umol/kgSW)
%    19 - RF out            ()
%    20 - OmegaCa out       ()
%    21 - OmegaAr out       ()
%    22 - xCO2 out          (ppm)
%
% NOTE: User-specified uncertainties for the input arguments are provided
%       in addition to uncertainties in output arguments at both input and
%       output conditions.
%

function [total_error, headers, units] = ...
        errors (PAR1, PAR2, PAR1TYPE, PAR2TYPE, SAL, TEMPIN, TEMPOUT, PRESIN, PRESOUT, SI, PO4,...
                NH4, H2S, ePAR1, ePAR2, eSAL, eTEMP, eSI, ePO4, eNH4, eH2S, epK, eBt, r, ...
                pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON,eCAL)

    global K0 K1 K2 KW KB KF KS KP1 KP2 KP3 KSi KNH4 KH2S E;
    
    % Input conditioning
    % ------------------
    
    % Handle optional calcium error input
    if (nargin < 30), eCAL=0; end
    
    % Determine lengths of input vectors
    veclengths=[length(PAR1) length(PAR2) length(PAR1TYPE)...
                length(PAR2TYPE) length(SAL) length(TEMPIN)...
                length(TEMPOUT) length(PRESIN) length(PRESOUT)...
                length(SI) length(PO4) length(NH4) length(H2S)...
                length(ePAR1) length(ePAR2) length(eSAL) length(eTEMP)...
                length(eSI) length(ePO4) length(eNH4) length(eH2S)...
                length(r) length(pHSCALEIN) length(K1K2CONSTANTS)...
                length(KSO4CONSTANT) length(KFCONSTANT) length(BORON)...
                length(eCAL)];

    if length(unique(veclengths))>2
            disp(' '); disp('*** INPUT ERROR: Input vectors must all be of same length, or of length 1. ***'); disp(' '); return
    end

    % Make column vectors of all input vectors
    PAR1         =PAR1         (:);
    PAR2         =PAR2         (:);
    PAR1TYPE     =PAR1TYPE     (:);
    PAR2TYPE     =PAR2TYPE     (:);
    SAL          =SAL          (:);
    TEMPIN       =TEMPIN       (:);
    TEMPOUT      =TEMPOUT      (:);
    PRESIN       =PRESIN       (:);
    PRESOUT      =PRESOUT      (:);
    SI           =SI           (:);
    PO4          =PO4          (:);
    NH4          =NH4          (:);
    H2S          =H2S          (:);
    ePAR1        =ePAR1        (:);
    ePAR2        =ePAR2        (:);
    eSAL         =eSAL         (:);
    eTEMP        =eTEMP        (:);
    eSI          =eSI          (:);
    ePO4         =ePO4         (:);
    eNH4         =eNH4         (:);
    eH2S         =eH2S         (:);
    r            =r            (:);
    pHSCALEIN    =pHSCALEIN    (:);
    K1K2CONSTANTS=K1K2CONSTANTS(:);
    KSO4CONSTANT =KSO4CONSTANT (:);
    KFCONSTANT   =KFCONSTANT   (:);
    BORON        =BORON        (:);
    eCAL         =eCAL         (:);
    
    % Find the longest column vector:
    ntps = max(veclengths);

    % Populate column vectors
    PAR1(1:ntps,1)          = PAR1(:)          ;
    PAR2(1:ntps,1)          = PAR2(:)          ;
    PAR1TYPE(1:ntps,1)      = PAR1TYPE(:)      ;
    PAR2TYPE(1:ntps,1)      = PAR2TYPE(:)      ;
    SAL(1:ntps,1)           = SAL(:)           ;
    TEMPIN(1:ntps,1)        = TEMPIN(:)        ;
    TEMPOUT(1:ntps,1)       = TEMPOUT(:)       ;
    PRESIN(1:ntps,1)        = PRESIN(:)        ;
    PRESOUT(1:ntps,1)       = PRESOUT(:)       ;
    SI(1:ntps,1)            = SI(:)            ;
    PO4(1:ntps,1)           = PO4(:)           ;
    NH4(1:ntps,1)           = NH4(:)           ;
    H2S(1:ntps,1)           = H2S(:)           ;
    ePAR1(1:ntps,1)         = ePAR1(:)         ;
    ePAR2(1:ntps,1)         = ePAR2(:)         ;
    eSAL(1:ntps,1)          = eSAL(:)          ;
    eTEMP(1:ntps,1)         = eTEMP(:)         ;
    eSI(1:ntps,1)           = eSI(:)           ;
    ePO4(1:ntps,1)          = ePO4(:)          ;
    eNH4(1:ntps,1)          = eNH4(:)          ;
    eH2S(1:ntps,1)          = eH2S(:)          ;
    r(1:ntps,1)             = r(:)             ;
    pHSCALEIN(1:ntps,1)     = pHSCALEIN(:)     ;
    K1K2CONSTANTS(1:ntps,1) = K1K2CONSTANTS(:) ;
    KSO4CONSTANT(1:ntps,1)  = KSO4CONSTANT(:)  ;
    KFCONSTANT(1:ntps,1)    = KFCONSTANT(:)    ;
    BORON(1:ntps,1)         = BORON(:)         ;
    eCAL(1:ntps,1)          = eCAL(:)          ;

    % Exclude input variables equal to -999
    E = (PAR1~=-999 & PAR2~=-999);
    
    % Default values for epK
    if (isempty(epK))
        epK = [0.002, 0.0075, 0.015, 0.01, 0.01, 0.02, 0.02];
    else
        % Check validity of epK
        if (length(epK) == 1 && epK == 0)
            % this means that the caller does not want to account for errors on dissoc. constants
            epK = [0.0 0.0 0.0 0.0 0.0 0.0 0.0];
        elseif (length(epK) ~= 7)
            error ('invalid parameter epK: ', epK)
        end
    end
    
    % Default value for eBt (also check for incorrectly specified values
    if (isempty(eBt))
        eBt = 0.02;
    elseif ( ~(isscalar(eBt)) )
        error ('invalid parameter eBt (must be scalar): ')
    elseif ( isscalar(eBt))
        if (eBt < 0 || eBt > 1)
           error ('The "eBt" input argument is the fractional error. It must be between 0 and 1. Default is 0.02 (a 2% error).')
	end
    end

    % names of dissociation constants
    Knames = {'K0','K1','K2','Kb','Kw','Kspa', 'Kspc'};

    % Convert error on pH to error on [H+] concentration
    % in case where first input variable is pH
    isH = (E & PAR1TYPE == 3);
    
    pH = PAR1(isH);
    epH = ePAR1(isH);       % Error on pH
    H  = 10.^(-pH);         % H+ concentration
    r(isH) = -r(isH);       % Inverse sign of 'r' if PAR1 is pH

    % dpH = d(-log10[H])
    %     = d(- ln[H] / ln[10] )
    %     = -(1/ln[10]) * d (ln[H])
    %     = -(1/ln[10]) * (dH / H)
    % Thus dH = - ln[1O] * [H] dpH
    eH =  log(10) * (H .* epH);     % Removed the minus sign because all errors (sigmas) are positive by definition
    eH =  eH * 1e9            ;     % Convert from mol/kg to nmol/kg (to have same units as partial derivative)
    ePAR1(isH) = eH;

    % Same conversion for second variable
    isH = (E & PAR2TYPE == 3);
    pH = PAR2(isH);
    epH = ePAR2(isH);       % Error on pH
    H  = 10.^(-pH);         % H+ concentration
    r(isH) = -r(isH);       % Inverse sign of 'r' if PAR2 is pH

    eH =   log(10) * (H .* epH);
    eH =  eH * 1e9;
    ePAR2(isH) = eH;

    % Convert calcium error from mmol/kg to mol/kg
    eCAL = eCAL.*1e-3;

    % initialise total square error
    sq_err = zeros(ntps,1);
    
    % Contribution of PAR1 to squared standard error
    if (any(ePAR1 ~= 0.0))
        % Compute sensitivities (partial derivatives)
        [deriv1(E,:),~,~,headers_err,units_err] = derivnum ('PAR1',...
            PAR1(E),PAR2(E),PAR1TYPE(E),PAR2TYPE(E),SAL(E),TEMPIN(E),...
            TEMPOUT(E),PRESIN(E),PRESOUT(E),SI(E),PO4(E),NH4(E),H2S(E),...
            pHSCALEIN(E),K1K2CONSTANTS(E),KSO4CONSTANT(E),KFCONSTANT(E),BORON(E));
        err = bsxfun(@times,deriv1,ePAR1);
	 sq_err = bsxfun(@plus,err*0., sq_err);
     sq_err = sq_err + err .* err;
    end

    % Contribution of PAR2 to squared standard error
    if (any (ePAR2 ~= 0.0))
        % Compute sensitivities (partial derivatives)
        [deriv2(E,:),~,~,headers_err,units_err] = derivnum ('PAR2',...
            PAR1(E),PAR2(E),PAR1TYPE(E),PAR2TYPE(E),SAL(E),TEMPIN(E),...
            TEMPOUT(E),PRESIN(E),PRESOUT(E),SI(E),PO4(E),NH4(E),H2S(E),...
            pHSCALEIN(E),K1K2CONSTANTS(E),KSO4CONSTANT(E),KFCONSTANT(E),BORON(E));
        err = bsxfun(@times,deriv2,ePAR2);
	 sq_err = bsxfun(@plus,err*0., sq_err);
     sq_err = sq_err + err .* err;
    end

    % Contribution of covariance of PAR1 and PAR2 to squared standard error
    covariance  = nan(ntps,1);
    if (any (r ~= 0.0) && any (ePAR1 ~= 0.0) && any (ePAR2 ~= 0.0))
        % Compute covariance from correlation coeff. & std deviations
        covariance(E) = r(E) .* ePAR1(E) .* ePAR2(E);
        % Contribution to squared error
        err2 = bsxfun(@times,2 * deriv1 .* deriv2, covariance);
        sq_err = sq_err + err2;
    end
    
    % Contribution of Silion (total dissolved inorganic concentration) to squared standard error
    %
    % Remark : does not compute error where SI = 0 
    %          because computation of sensitivity to SI fails in that case
    %
    SI_valid = (E & SI ~= 0 & eSI ~= 0);
    if (any (SI_valid))
        % Compute sensitivities (partial derivatives)
        [deriv, ~, ~, headers_err, units_err] = derivnum ('sil',PAR1(SI_valid),PAR2(SI_valid),PAR1TYPE(SI_valid),PAR2TYPE(SI_valid),...
                   SAL(SI_valid),TEMPIN(SI_valid),TEMPOUT(SI_valid),PRESIN(SI_valid),PRESOUT(SI_valid),...
                   SI(SI_valid),PO4(SI_valid),NH4(SI_valid),H2S(SI_valid),pHSCALEIN(SI_valid),K1K2CONSTANTS(SI_valid),KSO4CONSTANT(SI_valid),...
                   KFCONSTANT(SI_valid),BORON(SI_valid));
        err = bsxfun(@times,deriv,eSI(SI_valid));
        new_size = [ntps size(err,2)];
	sq_err = zeros(new_size) + sq_err;
        sq_err(SI_valid,:) = sq_err(SI_valid,:) + err .* err;
    end

    % Contribution of Phosphorus (total dissoloved inorganic concentration) to squared standard error
    %
    % Remark : does not compute error where PO4 = 0 
    %          because computation of sensitivity to PO4 fails in that case
    %
    PO4_valid = (E & PO4 ~= 0 & ePO4 ~= 0);
    if (any (PO4_valid))
        % Compute sensitivities (partial derivatives)
        [deriv, ~, ~, headers_err, units_err] = derivnum ('phos',PAR1(PO4_valid),PAR2(PO4_valid),PAR1TYPE(PO4_valid),PAR2TYPE(PO4_valid),...
                   SAL(PO4_valid),TEMPIN(PO4_valid),TEMPOUT(PO4_valid),PRESIN(PO4_valid),PRESOUT(PO4_valid),...
                   SI(PO4_valid),PO4(PO4_valid),NH4(PO4_valid),H2S(PO4_valid),pHSCALEIN(PO4_valid),K1K2CONSTANTS(PO4_valid),KSO4CONSTANT(PO4_valid),...
                   KFCONSTANT(PO4_valid),BORON(PO4_valid));
        err = bsxfun(@times,deriv,ePO4(PO4_valid));
        new_size = [ntps size(err,2)];
	sq_err = zeros(new_size) + sq_err;
        sq_err(PO4_valid,:) = sq_err(PO4_valid,:) + err .* err;
    end
    
    % Contribution of Ammonia (total dissolved inorganic concentration) to squared standard error
    %
    % Remark : does not compute error where NH4 = 0 
    %          because computation of sensitivity to SI fails in that case
    %
    NH4_valid = (E & NH4 ~= 0 & eNH4 ~= 0);
    if (any (NH4_valid))
        % Compute sensitivities (partial derivatives)
        [deriv, ~, ~, headers_err, units_err] = derivnum ('amm',PAR1(NH4_valid),PAR2(NH4_valid),PAR1TYPE(NH4_valid),PAR2TYPE(NH4_valid),...
                   SAL(NH4_valid),TEMPIN(NH4_valid),TEMPOUT(NH4_valid),PRESIN(NH4_valid),PRESOUT(NH4_valid),...
                   SI(NH4_valid),PO4(NH4_valid),NH4(NH4_valid),H2S(NH4_valid),pHSCALEIN(NH4_valid),K1K2CONSTANTS(NH4_valid),KSO4CONSTANT(NH4_valid),...
                   KFCONSTANT(NH4_valid),BORON(NH4_valid));
        err = bsxfun(@times,deriv,eNH4(NH4_valid));
        new_size = [ntps size(err,2)];
	sq_err = zeros(new_size) + sq_err;
        sq_err(NH4_valid,:) = sq_err(NH4_valid,:) + err .* err;
    end

    % Contribution of Hydrogen Sulfide (total dissoloved inorganic concentration) to squared standard error
    %
    % Remark : does not compute error where H2S = 0 
    %          because computation of sensitivity to PO4 fails in that case
    %
    H2S_valid = (E & H2S ~= 0 & eH2S ~= 0);
    if (any (H2S_valid))
        % Compute sensitivities (partial derivatives)
        [deriv, ~, ~, headers_err, units_err] = derivnum ('hyd',PAR1(H2S_valid),PAR2(H2S_valid),PAR1TYPE(H2S_valid),PAR2TYPE(H2S_valid),...
                   SAL(H2S_valid),TEMPIN(H2S_valid),TEMPOUT(H2S_valid),PRESIN(H2S_valid),PRESOUT(H2S_valid),...
                   SI(H2S_valid),PO4(H2S_valid),NH4(H2S_valid),H2S(H2S_valid),pHSCALEIN(H2S_valid),K1K2CONSTANTS(H2S_valid),KSO4CONSTANT(H2S_valid),...
                   KFCONSTANT(H2S_valid),BORON(H2S_valid));
        err = bsxfun(@times,deriv,eH2S(H2S_valid));
        new_size = [ntps size(err,2)];
	sq_err = zeros(new_size) + sq_err;
        sq_err(H2S_valid,:) = sq_err(H2S_valid,:) + err .* err;
    end

    % Contribution of T (temperature) to squared standard error
    if (any (eTEMP ~= 0.0))
        % Compute sensitivities (partial derivatives)
        [deriv(E,:),~,~,headers_err,units_err] = derivnum ('T',...
            PAR1(E),PAR2(E),PAR1TYPE(E),PAR2TYPE(E),SAL(E),TEMPIN(E),...
            TEMPOUT(E),PRESIN(E),PRESOUT(E),SI(E),PO4(E),NH4(E),H2S(E),...
            pHSCALEIN(E),K1K2CONSTANTS(E),KSO4CONSTANT(E),KFCONSTANT(E),BORON(E));
        err = bsxfun(@times,deriv,eTEMP);
	 sq_err = err*0. + sq_err;
     sq_err = sq_err + err .* err;
    end

    % Contribution of S (salinity) to squared standard error
    if (any (eSAL ~= 0.0))
        % Compute sensitivities (partial derivatives)
        [deriv(E,:),~,~,headers_err,units_err] = derivnum ('S',...
            PAR1(E),PAR2(E),PAR1TYPE(E),PAR2TYPE(E),SAL(E),TEMPIN(E),...
            TEMPOUT(E),PRESIN(E),PRESOUT(E),SI(E),PO4(E),NH4(E),H2S(E),...
            pHSCALEIN(E),K1K2CONSTANTS(E),KSO4CONSTANT(E),KFCONSTANT(E),BORON(E));
        err = bsxfun(@times,deriv,eSAL);
	 sq_err = err*0. + sq_err;
     sq_err = sq_err + err .* err;
    end

    % Calculate dissociation constants
    data = CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON);
        
    % Calculate [Ca++]
    % '       Riley, J. P. and Tongudai, M., Chemical Geology 2:263-269, 1967:
    % '       this is .010285.*Sali./35
    Ca = 0.02128./40.087.*(SAL./1.80655);% ' in mol/kg-SW

    % Contribution of all pKi to squared standard error
    for i = 1:length(epK)

        % if error on Ki is given
        if (epK(i) ~= 0.0)
            % Select Ki
            switch i
                case 1
                  Ki_in = data(:,63);    % K0in
                  Ki_out = data(:,78);   % K0out
                case 2
                  Ki_in = data(:,64);    % K1in
                  Ki_out = data(:,79);   % K1out
                case 3
                  Ki_in = data(:,65);    % K2in
                  Ki_out = data(:,80);   % K2out
                case 4
                  Ki_in = data(:,69);    % KBin
                  Ki_out = data(:,84);   % KBout
                case 5
                  Ki_in = data(:,68);    % KWin
                  Ki_out = data(:,83);   % KWout
                case 6
                  % Recompute KAr from OmegaAr and ions [Ca++] and [CO3--] concentrations
                  OmegaAr = data(:,18);
                  CO3 = data(:,7) * 1e-6;
                  Ki_in = CO3.*Ca./OmegaAr;  % KspAin
		          Ki_in(isnan(Ki_in)) = 0;
                  OmegaAr = data(:,36);
                  CO3 = data(:,25) * 1e-6;
                  Ki_out = CO3.*Ca./OmegaAr; % KspAout
		          Ki_out(isnan(Ki_out)) = 0;
                case 7
                  % Recompute KCa from OmegaCa and ions [Ca++] and [CO3--] concentrations
                  OmegaCa = data(:,17);
                  CO3 = data(:,7) * 1e-6;
                  Ki_in = CO3.*Ca./OmegaCa;  % KspCin
                  Ki_in(isnan(Ki_in)) = 0;
                  OmegaCa = data(:,35);
                  CO3 = data(:,25) * 1e-6;
                  Ki_out = CO3.*Ca./OmegaCa; % KspCout
                  Ki_out(isnan(Ki_out)) = 0;
            end

            % compute error on Ki from that on pKi
            eKi_in  = - epK(i) * Ki_in * log(10);  % error on K at input conditions
            eKi_out = - epK(i) * Ki_out * log(10); % error on K at output conditions
            % Compute sensitivities (partial derivatives)
            [deriv(E,:),~,~,headers_err,units_err] = derivnum (cell2mat(Knames(1,i)),...
                       PAR1(E),PAR2(E),PAR1TYPE(E),PAR2TYPE(E),...
                       SAL(E),TEMPIN(E),TEMPOUT(E),PRESIN(E),PRESOUT(E),...
                       SI(E),PO4(E),NH4(E),H2S(E),pHSCALEIN(E),K1K2CONSTANTS(E),...
                       KSO4CONSTANT(E),KFCONSTANT(E),BORON(E));
            % Derivatives at input conditions are calculated with respect
            % to constants at input conditions, and derivatives at output
            % conditions are calculated with respect to constants at output
            % conditions. So, error terms are calculated using absolute
            % errors in K values at both input and output conditions.
            err(:, 1:12) = bsxfun(@times, deriv(:, 1:12), eKi_in);
            err(:, 13:22) = bsxfun(@times, deriv(:, 13:22), eKi_out);
	     sq_err = err*0. + sq_err;
         sq_err = sq_err + err .* err;
        end
    end

    % Contribution of Boron (total dissoloved boron concentration) to squared standard error
    if (eBt ~= 0)
        % Compute sensitivities (partial derivatives)
        [deriv(E,:),~,~,headers_err,units_err] = derivnum ('bor',...
            PAR1(E),PAR2(E),PAR1TYPE(E),PAR2TYPE(E),SAL(E),TEMPIN(E),...
            TEMPOUT(E),PRESIN(E),PRESOUT(E),SI(E),PO4(E),NH4(E),H2S(E),...
            pHSCALEIN(E),K1K2CONSTANTS(E),KSO4CONSTANT(E),KFCONSTANT(E),BORON(E));
     err = bsxfun(@times, deriv, (eBt*data(:,93)*1e-6)); % where TB = data(:,93) in umol B/kg
     new_size = [ntps size(err,2)];
	 sq_err = zeros(new_size) + sq_err;
     sq_err = sq_err + err .* err;
    end
    
    % Contribution of Calcium (total dissoloved calcium concentration) to squared standard error
    CAL_valid = (E & eCAL ~= 0);
    if (any (CAL_valid))
        % Compute sensitivities (partial derivatives)
        [deriv,~,~,headers_err,units_err] = derivnum ('cal',...
            PAR1(CAL_valid),PAR2(CAL_valid),PAR1TYPE(CAL_valid),PAR2TYPE(CAL_valid),SAL(CAL_valid),TEMPIN(CAL_valid),...
            TEMPOUT(CAL_valid),PRESIN(CAL_valid),PRESOUT(CAL_valid),SI(CAL_valid),PO4(CAL_valid),NH4(CAL_valid),H2S(CAL_valid),...
            pHSCALEIN(CAL_valid),K1K2CONSTANTS(CAL_valid),KSO4CONSTANT(CAL_valid),KFCONSTANT(CAL_valid),BORON(CAL_valid));
    err = bsxfun(@times,deriv,eCAL(CAL_valid));
    new_size = [ntps size(err,2)];
	sq_err = zeros(new_size) + sq_err;
    sq_err(CAL_valid,:) = sq_err(CAL_valid,:) + err .* err;
    end

    % Compute and return resulting total error (or uncertainty)
    total_error = sqrt (sq_err);
    total_error(isnan(total_error))=-999;
    headers = strcat('u(',headers_err,')');
    units = strcat('(',units_err,')');
end

% 
% derivnum ----------------------------------------------------------------
%

% derivnum()
% This subroutine computes partial derivatives of output carbonate variables 
% with respect to input variables (two), plus nutrients (two), ammonium and hydrogen sulfide,
% temperature and salinity, dissociation constants, and total boron.
%
% It uses central differences, introducing a small perturbation 
% (plus and minus of a delta) in one INPUT variable and computes the
% resulting induced change in each OUTPUT variable
%
% This subroutine has been modified from its original version (Orr et al.
% 2018) to allow for input variables of CO2, HCO3, and CO3, the inclusion
% of ammonium and hydrogen sulfide, and compatibility with CO2SYS.m(v3)
%
% After numerical tests, the small PERTURBATION (delta) is chosen
% as a fraction of a reference value as follows:
% * 1.e-3 (0.1%) for the equilibrium constants and total boron 
% * 1.e-6        for the pair of CO2 system input variables (PAR1, PAR2)
% * 1.e-4        for temperature and salinity
% * 1.e-3        for total dissolved inorganic P and Si (PO4, SI)
%
%**************************************************************************
%
%  **** SYNTAX:
%  [deriv, headers_der, units_der, headers_err, units_err]=...
%                  derivnum(VARID,PAR1,PAR2,PAR1TYPE,PAR2TYPE,...
%                           SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,...
%                           SI,PO4,NH4,H2S,pHSCALEIN,K1K2CONSTANTS,...
%                           KSO4CONSTANT,KFCONSTANT,BORON)
% 
%  **** SYNTAX EXAMPLES:
%  [der, headers, units] = derivnum('par1',2400,2200,1,2,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [der, headers, units] = derivnum('sit', 2400,   8,1,3,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [deriv, headers]      = derivnum(  'T',  500,   8,5,3,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [deriv]               = derivnum('S',2400,2000:10:2400,1,2,35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [deriv]               = derivnum('K0',2400,2200,1,2,0:1:35,0,25,4200,0,15,1,0,0,1,4,1,1,1)
%  [deriv]               = derivnum('K1',2400,2200,1,2,35,0,25,0:100:4200,0,15,1,0,0,1,4,1,1,1)
%  
%**************************************************************************
%
% INPUT:
%
%   - VARID     :   character string to select the input variable for which
%                   derivatives will be taken with respect to. 
%                   This variable appears in denominator of each resulting derivative.
%                   = variable length, case insensitive, character code
%                     case 'par1'  :  Parameter 1 of the input pair (This is TAlk if PAR1TYPE is 1)
%                     case 'par2'  :  Parameter 2 of the input pair (This is TAlk if PAR2TYPE is 1)
%                     case 'sil', 'silt', 'tsil', 'silicate', or 'sit'     : Silicate concentration
%                     case 'phos', 'phost', 'tphos', 'phosphate', or 'pt'  : Phosphate concentration
%                     case 'amm', 'nh4', 'NH4'                    : Ammonia concentration
%                     case 'hyd', 'h2s', 'H2S'                    : Hydrogen Sulfide concentration
%                     case 't', 'temp' or 'temperature'           : temperature
%                     case 's', 'sal', or 'salinity'              : salinity
%                     case 'K0','K1','K2','Kb','Kw','Kspa', 'Kspc': dissociation constants 
%                     case 'bor': total boron
%
%   - all others :  same list of input parameters as in CO2SYS() subroutine (version 3, scalar or vectors)
%
%**************************************************************************%
%
% OUTPUT: 3 arrays	
%         a) an array containing the derivative of the following parameter (one row per sample):
%         b) the corresponding cell-array containing crudely formatted headers
%         c) the corresponding cell-array containing the units
%
%    POS  PARAMETER         UNIT
%
%    01 - TAlk              (umol/kgSW)
%    02 - TCO2              (umol/kgSW)
%    03 - [H+] in           (umol/kgSW)
%    04 - pCO2 in           (uatm)
%    05 - fCO2 in           (uatm)
%    06 - HCO3 in           (umol/kgSW)
%    07 - CO3 in            (umol/kgSW)
%    08 - CO2 in            (umol/kgSW)
%    09 - RF in             ()
%    10 - OmegaCa in        ()
%    11 - OmegaAr in        ()
%    12 - xCO2 in           (ppm)
%    13 - [H+] out          ()
%    14 - pCO2 out          (uatm)
%    15 - fCO2 out          (uatm)
%    16 - HCO3 out          (umol/kgSW)
%    17 - CO3 out           (umol/kgSW)
%    18 - CO2 out           (umol/kgSW)
%    19 - RF out            ()
%    20 - OmegaCa out       ()
%    21 - OmegaAr out       ()
%    22 - xCO2 out          (ppm)
%
% * 'in'  refers to INPUT  conditions (TEMPIN, PRESIN) as for CO2SYS
%   'out' refers to OUTPUT conditions (TEMPOUT, PRESOUT)
%
%
% CAUTION: derivnum.m is NOT necessarily designed to take partial derivatives
%          of input variables, only computed variables relative to input
%          variables. However, those partial derivatives of input variables
%          are kept in the OUT conditions to maintain consistency
%          of the order of output between the different input pairs.
%          Generally, we advise not to use these derivatives of input
%          variables. However, in some cases their results appear accurate.
%          Use them at your own risk.
%

function [derivatives, headers, units, headers_err, units_err] = ...
        derivnum (VARID,PAR1,PAR2,PAR1TYPE,PAR2TYPE, SAL,TEMPIN, ...
                  TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S, ...
                  pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT, ...
                  KFCONSTANT,BORON)

    % For computing derivative with respect to Ks, one has to call CO2sys with a perturbed K
    % Requested perturbation is passed through the following global variables
    global PertK    % Id of perturbed K
    global Perturb  % perturbation

    % Input conditioning
    % ------------------

    VARID = upper(VARID);
    % Determine lengths of input vectors
    veclengths=[length(PAR1) length(PAR2) length(PAR1TYPE)...
                length(PAR2TYPE) length(SAL) length(TEMPIN)...
                length(TEMPOUT) length(PRESIN) length(PRESOUT)...
                length(SI) length(PO4) length(NH4) length(H2S)...
                length(pHSCALEIN) length(K1K2CONSTANTS)...
                length(KSO4CONSTANT) length(KFCONSTANT)...
                length(BORON)];

    if length(unique(veclengths))>2
            disp(' '); disp('*** INPUT ERROR: Input vectors must all be of same length, or of length 1. ***'); disp(' '); return
    end

    % Make column vectors of all input vectors
    PAR1         =PAR1         (:);
    PAR2         =PAR2         (:);
    PAR1TYPE     =PAR1TYPE     (:);
    PAR2TYPE     =PAR2TYPE     (:);
    SAL          =SAL          (:);
    TEMPIN       =TEMPIN       (:);
    TEMPOUT      =TEMPOUT      (:);
    PRESIN       =PRESIN       (:);
    PRESOUT      =PRESOUT      (:);
    SI           =SI           (:);
    PO4          =PO4          (:);
    NH4          =NH4          (:);
    H2S          =H2S          (:);
    pHSCALEIN    =pHSCALEIN    (:);
    K1K2CONSTANTS=K1K2CONSTANTS(:);
    KSO4CONSTANT =KSO4CONSTANT (:);
    KFCONSTANT   =KFCONSTANT   (:);
    BORON        =BORON        (:);

    % Find the longest column vector:
    ntps = max(veclengths);

    % Populate column vectors
    PAR1(1:ntps,1)          = PAR1(:)          ;
    PAR2(1:ntps,1)          = PAR2(:)          ;
    PAR1TYPE(1:ntps,1)      = PAR1TYPE(:)      ;
    PAR2TYPE(1:ntps,1)      = PAR2TYPE(:)      ;
    SAL(1:ntps,1)           = SAL(:)           ;
    TEMPIN(1:ntps,1)        = TEMPIN(:)        ;
    TEMPOUT(1:ntps,1)       = TEMPOUT(:)       ;
    PRESIN(1:ntps,1)        = PRESIN(:)        ;
    PRESOUT(1:ntps,1)       = PRESOUT(:)       ;
    SI(1:ntps,1)            = SI(:)            ;
    PO4(1:ntps,1)           = PO4(:)           ;
    SI(1:ntps,1)            = SI(:)            ;
    PO4(1:ntps,1)           = PO4(:)           ;
    NH4(1:ntps,1)           = NH4(:)           ;
    H2S(1:ntps,1)           = H2S(:)           ;
    pHSCALEIN(1:ntps,1)     = pHSCALEIN(:)     ;
    K1K2CONSTANTS(1:ntps,1) = K1K2CONSTANTS(:) ;
    KSO4CONSTANT(1:ntps,1)  = KSO4CONSTANT(:)  ;
    KFCONSTANT(1:ntps,1)    = KFCONSTANT(:)    ;
    BORON(1:ntps,1)         = BORON(:)         ;

    % BASELINE:
    % --------
    
    carb = CO2SYS(PAR1,PAR2,PAR1TYPE,PAR2TYPE,SAL,TEMPIN,TEMPOUT,PRESIN,PRESOUT,SI,PO4,NH4,H2S,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON);
    % Compute [H+] in µmol/KgSW
    if (ndims(carb) == 2)
        Hin = 10.^(-carb(:,3)) * 1.e6;
        Hout = 10.^(-carb(:,21)) * 1.e6;
        carb = horzcat(carb(:,1:2), Hin, carb(:,3:20), Hout, carb(:,21:end));
    else
        Hin = 10.^(-carb(3)) * 1.e6;
        Hout = 10.^(-carb(21)) * 1.e6;
        carb = horzcat(carb(1:2), Hin, carb(3:20), Hout, carb(21:end));
    end    

    % Compute two slightly different values for input
    % -----------------------------------------------

    % Default values : not perturbed
    % no change in PAR1 and PAR2
    PAR11 = PAR1; PAR12 = PAR1;
    PAR21 = PAR2; PAR22 = PAR2;
    % no change in Sil total and Phos total (except for d/dSi, d/dPO4)
    SI1 = SI;    SI2 = SI;
    PO41  = PO4; PO42  = PO4;
    % no change in Amm total and H2S total (except for d/dAmm, d/dH2S)
    NH41  = NH4; NH42 = NH4;
    H2S1  = H2S; H2S2  = H2S;
    % no change in T and S (except for d/dS and d/dT, respectively)
    TEMP1 = TEMPIN; TEMP2 = TEMPIN;
    SAL1 = SAL; SAL2 = SAL;
	% no change in TEMPOUT except for d/dT (see why further below)	
    TEMPOUT1 = TEMPOUT; TEMPOUT2 = TEMPOUT;

    % Create empty vector for abs_dx (absolute delta)
    abs_dx  = nan(ntps,1);

    % Flag for dissociation constant as perturbed variable 
    flag_dissoc_K = 0;   % False
    % names of 7 dissociation constants and variable for total boron 
    K_names = {'K0', 'K1', 'K2', 'KB', 'KW', 'KSPA', 'KSPC', 'BOR', 'CAL'};    

    % Units for derivatives and for errors
    units_at = {'umol';'umol';'nmol';'total scale';'uatm kg';'uatm kg';'umol';'umol';...
                 'umol';'kg';'kg';'kg';'ppm kg';...
                 'nmol';'uatm kg';'uatm kg';'umol';'umol';...
                 'umol';'kg';'kg';'kg';'ppm kg';
                };
 
    units_kg  = {'umol/kg';'umol/kg';'nmol/kg';'total scale';'uatm';'uatm';'umol/kg';'umol/kg';...
                 'umol/kg';' ';' ';' ';'ppm';...
                 'nmol/kg';'uatm';'uatm';'umol/kg';'umol/kg';...
                 'umol/kg';' ';' ';' ';'ppm';
                };

    units_k = units_kg;
 
    units_pco2 = units_kg;
    
    units_err = units_kg;
 
    switch VARID
        case K_names
            flag_dissoc_K = 1;
            % Approximate values for K0, K1, K2, Kb, Kspa, Kspc, and CAL
            % They will be used to compute an absolute perturbation
            % value on these constants
            K = [0.034, 1.2e-06, 8.3e-10, 2.1e-09, 6.1e-14, 6.7e-07, ...
                 4.3e-07, 0.0004157, 0.0103];
            % Choose value of absolute perturbation
            [is_in_K_names, index] = ismember(VARID, K_names);
            perturbation = K(index) * 1.e-3;   % 0.1 percent of Kx value
            abs_dx = 2 * perturbation;
            denom_headers = VARID;
            denom_units = ' ';
            units = units_k;
            
        case {'PAR1', 'VAR1'}      % PAR1 (first variable of input pair) is perturbed
            % Define a relative delta
            delta = 1.e-6;
          
            PAR1ref = nan(size(PAR1TYPE));
            units = units_kg;
            denom_headers = '<PAR1>';
            denom_units = '<PAR1 units>';
                t1=PAR1TYPE==1;
                  PAR1ref(t1) = 2300.; % umol/kg (global surface average, Orr et al., 2017)
                t2=PAR1TYPE==2;
                  PAR1ref(t2) = 2000.; % umol/kg (global surface average, Orr et al., 2017)
                t3=PAR1TYPE==3;
                  PAR1ref(t3) = 1.0e-8; % mol/kg (equivalent to pH=8.0)
                t4=PAR1TYPE==4;
                  PAR1ref(t4) = 400.;  % uatm
                t5=PAR1TYPE==5;
                  PAR1ref(t5) = 400.; % uatm
                t6=PAR1TYPE==6;
                  PAR1ref(t6) = 1790.; % umol/kg
                t7=PAR1TYPE==7;
                  PAR1ref(t7) = 200.; % umol/kg
                t8=PAR1TYPE==8;
                  PAR1ref(t8) = 10.; % umol/kg
      
           % cases where first input variable is pH
            F = (PAR1TYPE == 3);
            H = 10.^(-PAR1(F)) ; % [H+] in mol/kg
            % Change slightly [H+]
            H1 = H - PAR1ref(F)*delta;
            H2 = H + PAR1ref(F)*delta;
            PAR11(F) = -log10(H1) ;
            PAR12(F) = -log10(H2) ;
            abs_dx(F) = (H2 - H1) * 1e9; % now in nmol/kg
           
            G = ~F;
            % Change slightly PAR1
            PAR11(G) = PAR1(G) - PAR1ref(G)*delta;
            PAR12(G) = PAR1(G) + PAR1ref(G)*delta;
            abs_dx(G) = PAR12(G) - PAR11(G);

      case {'PAR2', 'VAR2'}    % PAR2 (second variable of input pair) is perturbed
            % Define a relative delta
            delta = 1.e-6;
            
            PAR2ref = nan(size(PAR2TYPE));
            units = units_kg;
            denom_headers = '<PAR2>';
            denom_units = '<PAR2 units>';
                t1=PAR2TYPE==1;
                  PAR2ref(t1) = 2300.; % umol/kg (global surface average, Orr et al., 2017)
                t2=PAR2TYPE==2;
                  PAR2ref(t2) = 2000.; % umol/kg (global surface average, Orr et al., 2017)
                t3=PAR2TYPE==3;
                  PAR2ref(t3) = 1.0e-8; % mol/kg (equivalent to pH=8.0)
                t4=PAR2TYPE==4;
                  PAR2ref(t4) = 400.;  % uatm
                t5=PAR2TYPE==5;
                  PAR2ref(t5) = 400.; % uatm
                t6=PAR2TYPE==6;
                  PAR2ref(t6) = 1790.; % umol/kg
                t7=PAR2TYPE==7;
                  PAR2ref(t7) = 200.; % umol/kg
                t8=PAR2TYPE==8;
                  PAR2ref(t8) = 10.; % umol/kg
            
            
            % cases where second input variable is pH
            F = (PAR2TYPE == 3);
            H = 10.^(-PAR2(F)) ; % H+ in mol/kg
            % Change slightly [H+]
            H1 = H - PAR2ref(F)*delta;
            H2 = H + PAR2ref(F)*delta;
            PAR21(F) = -log10(H1) ;
            PAR22(F) = -log10(H2) ;
            abs_dx(F) = (H2 - H1) * 1e9;

            G = ~F;
            % Change slightly PAR2
            PAR21(G) = PAR2(G) - PAR2ref(G)*delta;
            PAR22(G) = PAR2(G) + PAR2ref(G)*delta;
            abs_dx(G) = PAR22(G) - PAR21(G);

       case {'SIL', 'TSIL', 'SILT', 'SILICATE', 'SIT'}    % Sil total
            % Define a relative delta
            delta = 1.e-3;
            SIref = 7.5; % umol/kg (global surface average, Orr et al., 2018)
            % Change slightly SI
            SI1 = SI - SIref*delta;
            SI2 = SI + SIref*delta;
            abs_dx = SI2 - SI1;
            denom_headers = 'Sit';
            denom_units = 'umol';
            units = units_at;
       case {'PHOS', 'TPHOS', 'PHOST', 'PHOSPHATE', 'PT'}    % Phos total
            % Define a relative delta
            delta = 1.e-3;
            PO4ref = 0.5; % umol/kg (global surface average, Orr et al., 2018)
            % Change slightly PO4
            PO41 = PO4 - PO4ref*delta;
            PO42 = PO4 + PO4ref*delta;
            abs_dx = PO42 - PO41;
            denom_headers = 'Pt';
            denom_units = 'umol';
            units = units_at;
       case {'AMM','NH4'}    % Amm total
            % Define a relative delta
            delta = 1.e-3;
            NH4ref = 1;
            % Change slightly NH4
            NH41 = NH4 - NH4ref*delta;
            NH42 = NH4 + NH4ref*delta;
            abs_dx = NH42 - NH41;
            denom_headers = 'NH4t';
            denom_units = 'umol';
            units = units_at;
       case {'HYD','H2S'}    % H2S total
            % Define a relative delta
            delta = 1.e-3;
            H2Sref = 1;
            % Change slightly H2S
            H2S1 = H2S - H2Sref*delta;
            H2S2 = H2S + H2Sref*delta;
            abs_dx = H2S2 - H2S1;
            denom_headers = 'H2St';
            denom_units = 'umol';
            units = units_at;
       case {'T', 'TEMP', 'TEMPERATURE'}    % Temperature
            % Define a relative delta
            delta = 1.e-4;
            TEMPref = 18.; % global surface mean (C)
            % Change slightly temperature
            TEMP1 = TEMPIN - TEMPref*delta;
            TEMP2 = TEMPIN + TEMPref*delta;
            abs_dx = TEMP2 - TEMP1;
            % When computing d/dT	
            % TEMPOUT changes must be the same as TEMPIN changes, e.g, for	
            % derivnum 'OUT' & 'IN' results to be the same when TEMPOUT=TEMPIN	
            TEMPOUT1 = TEMPOUT - TEMPref*delta;	
            TEMPOUT2 = TEMPOUT + TEMPref*delta;
            denom_headers = 'T';
            denom_units = 'C';
            units = units_kg;
       case {'S', 'SAL', 'SALINITY'}    % Salinity
            % Define a relative delta
            delta = 1.e-4;
            SALref = 35.;
            % Change slightly salinity
            SAL1 = SAL - SALref*delta;
            SAL2 = SAL + SALref*delta;
            abs_dx = SAL2 - SAL1;
            denom_headers = 'S';
            denom_units = 'psu';
            units = units_kg;
    end
    
    % PERTURBATION:
    % -------------

    % Point 1: (one dissociation constant or PAR1, PAR2, T or S is somewhat smaller)
    % if perturbed variable is a dissociation constant
    if (flag_dissoc_K)
        PertK = upper(VARID);
        Perturb = -perturbation;
    end
    cdel1 = CO2SYS ( ...
        PAR11,PAR21,PAR1TYPE,PAR2TYPE,SAL1,TEMP1,TEMPOUT1,PRESIN,PRESOUT,SI1,PO41,NH41,H2S1,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON);
    % Compute [H+]
    if (ndims(cdel1) == 2)
        Hin = 10.^(-cdel1(:,3))   * 1.e9; % to show H+ results in nmol/kg
        Hout = 10.^(-cdel1(:,21)) * 1.e9;
        cdel1 = horzcat(cdel1(:,1:2), Hin, cdel1(:,3:20), Hout, cdel1(:,21:end));
    else
        Hin = 10.^(-cdel1(3))     * 1.e9;
        Hout = 10.^(-cdel1(21))   * 1.e9;
        cdel1 = horzcat(cdel1(1:2), Hin, cdel1(3:20), Hout, cdel1(21:end));
    end

    % Point 2: (one dissociation constant or PAR1, PAR2, T or S is somewhat bigger)
    % if perturbed variable is a dissociation constant
    if (flag_dissoc_K)
        PertK = upper(VARID);
        Perturb = perturbation;
    end
    cdel2 = CO2SYS ( ...
        PAR12,PAR22,PAR1TYPE,PAR2TYPE,SAL2,TEMP2,TEMPOUT2,PRESIN,PRESOUT,SI2,PO42,NH42,H2S2,pHSCALEIN,K1K2CONSTANTS,KSO4CONSTANT,KFCONSTANT,BORON);
    % Compute [H+]
    if (ndims(cdel2) == 2)
        % Computed variable H+ (does not affect other computed
        % variables, i.e., it is the numerator of the derivative)
        Hin = 10.^(-cdel2(:,3))   * 1.e9; % to show H+ results in nmol/kg 
        Hout = 10.^(-cdel2(:,21)) * 1.e9;
        cdel2 = horzcat(cdel2(:,1:2), Hin, cdel2(:,3:20), Hout, cdel2(:,21:end));
    else
        Hin = 10.^(-cdel2(3))     * 1.e9;
        Hout = 10.^(-cdel2(21))   * 1.e9;
        cdel2 = horzcat(cdel2(1:2), Hin, cdel2(3:20), Hout, cdel2(21:end));
    end

    % if perturbed variable is a dissociation constant
    if (flag_dissoc_K)
        PertK = ''; %% Return to normal
    end
    
    % Drop unnecessary columns
    % Note: columns (04 - pHin) and (20 - pHout) drop    
    % Keep only the following columns
    %    01 - TAlk              (umol/kgSW)
    %    02 - TCO2              (umol/kgSW)
    %    03 - [H+] in           (nmol/kgSW)  lastly added
    %    05 - pCO2 in           (uatm)
    %    06 - fCO2 in           (uatm)
    %    07 - HCO3 in           (umol/kgSW)
    %    08 - CO3 in            (umol/kgSW)
    %    09 - CO2 in            (umol/kgSW)
    %    17 - RF in             ()
    %    18 - OmegaCa in        ()
    %    19 - OmegaAr in        ()
    %    20 - xCO2 in           (ppm)
    %    22 - [H+] out          (nmol/kgSW)  lastly added
    %    24 - pCO2 out          (uatm)
    %    25 - fCO2 out          (uatm)
    %    26 - HCO3 out          (umol/kgSW)
    %    27 - CO3 out           (umol/kgSW)
    %    28 - CO2 out           (umol/kgSW)
    %    36 - RF out            ()
    %    37 - OmegaCa out       ()
    %    38 - OmegaAr out       ()
    %    39 - xCO2 out          (ppm)
    keep = [1 2 3 5 6 7 8 9 17 18 19 20 22 24 25 26 27 28 36 37 38 39];
    
    % We will drop also some column headers
    headers = {'TAlk';'TCO2';'Hin';'pHin';'pCO2in';'fCO2in';'HCO3in';'CO3in';...
        'CO2in';'RFin';'OmegaCAin';'OmegaARin';'xCO2in';...
        'Hout';'pCO2out';'fCO2out';'HCO3out';'CO3out';...
        'CO2out';'RFout';'OmegaCAout';'OmegaARout';'xCO2out'};
    headers_err = headers;
    %units = {'umol';'umol';'nmol';'total scale';'uatm';'uatm';'umol';'umol';...
    %    'umol';' ';' ';'ppm';...
    %    'nmol';'uatm';'uatm';'umol';'umol';...
    %    'umol';' ';' ';'ppm';
    %    };
    % Initially, keep all headers except 'pHin'
    keep_head =  [1:3 5:23];

%     **** This was previously uncommented to eliminate partial derivatives
%     **** of input variables
%
%     % if all parameter PAR1 are of the same type
%     if all(PAR1TYPE == PAR1TYPE(1))
%         % Determine column number of PAR1
%         if PAR1TYPE(1) <= 3 % TAlk, TCO2 or pH
%             % By design of CO2sys, PARTYPE is equal to column number
%             col_number = PAR1TYPE(1);
%         else
%             % Because there is an extra column: [H+]
%             col_number = PAR1TYPE(1) + 1;
%         end
%         % Exclude input parameters PAR1
%         A = (keep ~= col_number);
%         keep = keep (A);
%         A = (keep_head ~= col_number);
%         keep_head = keep_head (A);
%     end
%     % if all parameter PAR2 are of the same type
%     if all(PAR2TYPE == PAR2TYPE(1))
%         % Determine column number of PAR1
%         if PAR2TYPE(1) <= 3 % TAlk, TCO2 or pH
%             % By design of CO2sys, PARTYPE is equal to column number
%             col_number = PAR2TYPE(1);
%         else
%             % Because there is an extra column: [H+]
%             col_number = PAR2TYPE(1) + 1;
%         end
%         % Exclude input parameters PAR2
%         A = (keep ~= col_number);
%         keep = keep (A);
%         A = (keep_head ~= col_number);
%         keep_head = keep_head (A);
%     end
    
    cdel1 = cdel1(:,keep);
    cdel2 = cdel2(:,keep);
    
    headers = headers(keep_head);
    units = units(keep_head);
    
    headers_err = headers_err(keep_head);
    units_err = units_err(keep_head);
    
    % concatenate strings to give each header its proper form, a partial derivative
    headers = strcat('d', headers, '/', 'd', denom_headers);
    % concatenate in an analogous way for the units
    units   = strcat('(',units, '/', denom_units,')');
    
    % Centered difference
    dy = (cdel2 - cdel1);

    % Compute ratio dy/dx
    if (isscalar(abs_dx))
        derivatives = dy ./ abs_dx;
    else
        derivatives = bsxfun(@rdivide, dy, abs_dx);
    end
    
%     **** This was previously uncommented to mask derivatives of
%     **** parameters that are directly related to input variables
%
%     % Mask values that should not be used with NaN (e.g., dHout/dT when PAR1 or PAR2 is pH)
%     switch VARID
%         case {'T', 'TEMP', 'TEMPERATURE'}
%             % For PAR1TYPE or PAR2TYPE = 3 (pH is input) make dHout/dT value a NaN
%             F = (PAR1TYPE==3 | PAR2TYPE==3); % either CO2 system input variable is pH
%             [is_in_headers, idx] = ismember('dHout/dT', headers);
%             if any(is_in_headers)
%                 derivatives(F,idx) = NaN ;
%             end
%             
%             % For PAR1TYPE or PAR2TYPE = 4 or 5 (pCO2 or fCO2 is input) make relevant d/dT values NaNs
%             F = (PAR1TYPE==4 | PAR2TYPE==4); %when pCO2 or fCO2 is input var
%             % masknan = {'dfCO2in/dT' 'dxCO2in/dT' 'dpCO2out/dT' 'dfCO2out/dT' 'dxCO2out/dT'};
%             masknan = {'dpCO2out/dT'};
%             [is_in_headers, idx] = ismember(masknan, headers);
%             if any(is_in_headers)
%                 derivatives(F,idx) = NaN ;
%             end
%             
%             % For PAR1TYPE or PAR2TYPE = 5 (fCO2 is input) make relevant d/dT values NaNs
%             F = (PAR1TYPE==5 | PAR2TYPE==5); %when fCO2 is input var
%             % masknan = {'dfCO2in/dT' 'dxCO2in/dT' 'dpCO2out/dT' 'dfCO2out/dT' 'dxCO2out/dT'};
%             masknan = {'dfCO2out/dT'};
%             [is_in_headers, idx] = ismember(masknan, headers);
%             if any(is_in_headers)
%                 derivatives(F,idx) = NaN ;
%             end
%         
%         case {'PAR1', 'VAR1'}
%             % For pH-pCO2 or pH-fCO2 pair (PAR1 is pH) 
%             F = (PAR1TYPE==3 & (PAR2TYPE==4 | PAR2TYPE==5));
%             masknan = {'dHout/dH' 'dpCO2out/dH' 'dfCO2out/dH'};
%             [is_in_headers, idx] = ismember(masknan, headers);
%             if any(is_in_headers)
%                 derivatives(F,idx) = NaN ;
%             end
%         
%         case {'PAR2', 'VAR2'}
%             % For pCO2-pH or fCO2-pH pair (PAR2 is pH)
%             F = ((PAR1TYPE==4 | PAR1TYPE==5) & PAR2TYPE==3);
%             masknan = {'dHout/dH' 'dpCO2out/dH' 'dfCO2out/dH'};
%             [is_in_headers, idx] = ismember(masknan, headers);
%             if any(is_in_headers)
%                 derivatives(F,idx) = NaN ;
%             end
%     end
    
end




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
    % opt.printcsv and opt.fname
    if ~isfield(opt,'printcsv') || isbad(opt.printcsv)
        opt.printcsv = 0; % default off
    elseif opt.printcsv > 1 || opt.printcsv < 0
        if opt.printmes ~= 0
            fprintf('Invalid CSV opt chosen. Assuming opt.csv = 1\n');
        end
    else
        if ~isfield(opt,'fname') || isbad(opt.fname)
            opt.fname = 'QUODcarb_output.csv';
            if opt.printmes ~= 0
                fprintf(['Invalid CSV filename. Assuming opt.fname' ...
                    ' = ''QUODcarb_output.csv'' \n']);
            end
        end
    end
    % opt.co2press
    if ~isfield(opt,'co2press') || isbad(opt.co2press)
        opt.co2press = 1; % on
        if opt.printmes ~=0
            fprintf('No opt.co2press chosen. Assuming opt.co2press = 1 (on). \n');
        end
    end
    % opt.Revelle
    if ~isfield(opt,'Revelle') || isbad(opt.Revelle)
        opt.Revelle = 0;
        if opt.printmes ~= 0
            fprintf('No opt.Revelle chosen. Assuming opt.Revelle = 0 (off). \n');
        end
    end
    if ~isfield(opt,'mpk')
        mpk0 = @(T,S,P) 1;
        mpk1 = @(T,S,P) 1;
        mpk2 = @(T,S,P) 1;
        mpk = @(T,S,P) [mpk0(T,S,P);mpk1(T,S,P);mpk2(T,S,P);1;1;1;1;1;1;1;1;1;1;1;1;1;1];
        opt.mpk = mpk;
    end
    if ~isfield(opt,'gmpk')        
        gmpk0 = @(T,S,P) [0 0 0];
        gmpk1 = @(T,S,P) [0 0 0];
        gmpk2 = @(T,S,P) [0 0 0];
        gmpk = @(T,S,P) [gmpk0(T,S,P);...
                         gmpk1(T,S,P);...
                         gmpk2(T,S,P);...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0;...
                         0,0,0];
        opt.gmpk = gmpk;
    end
    if ~isfield(opt,'empk')
        empk0 = 1;
        empk1 = 1;
        empk2 = 1;
        empk = [empk0;empk1;empk2;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
        opt.empk = empk;
    end
    % opt.turnoff
end




function [pK,gpK,epK] = calc_pK(opt,T,S,P)
% base equations COPIED FROM co2sys.m Orr et al. (2018)  Github
% Originally from  van Heuven et al. (2011)
% Original co2sys is from Lewis and Wallace (1998)

% INPUT:
%   T  = Temp (deg C)
%   S  = Salinity
%   P  = pressure (dbar)

% OUTPUT:
%    pK  = [pK0;pK1;pK2;pKb;pKw;pKs;pKf;pKp1;pKp2;pKp3;pKsi;pKnh4;pKh2s;pp2f;pKar;pKca];
%           -log10 values of equilibrium constants
%   gpK  = [pK_T, pK_S, pK_P]; 
%           first derivatives (gradient of pK wrt T, S, P) 
%           (size: length(pK) x 3 )
%   epK  = [epK0;epK1;epK2;epKb;epKw;epKs;epKf;epKp1;epKp2;epKp3;epKsi;epKnh4;epKh2s;epp2f;epKar;epKca];
%           errors of pK (1 standard deviation)

    TK   = T + 273.15; % convert to Kelvin
    Rgas = 83.14462618; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
    RT   = Rgas * TK;
    RT_T = Rgas;
    
    Pbar    = P / 10; % convert from dbar to bar
    Pbar_P  = 1 / 10;
    A       = 19.924; B = 1000; C = 1.005;
    ions    = @(S) A * S / ( B - C * S); % from DOE handbook
    ions_S  = @(S) ( A * B) / ( B - C * S)^2;
    LOG10   = log(10);
    p       = @(x) -log10(x);
    q       = @(x) 10.^(-x);  % inverse p, i.e., a backward p
    dpdx    = @(x) -1 / (x * LOG10);        % p'
    dqdx    = @(x) -LOG10 * 10.^( -x );     % q'
    
    % corrections for pressure---------------------------------------------
    % Millero 1995, 1992, 1982, 1979; Takahashi et al. 1982;
    %   Culberson & Pytkowicz 1968; Edmond & Gieskes 1970.
    dV      = @(T,a) a(1) + a(2) * T +     a(3) * T^2; 
    dV_T    = @(T,a)        a(2)     + 2 * a(3) * T;
    Ka      = @(T,b) ( b(1) + b(2) * T ) / 1000;
    Ka_T    = @(T,b) b(2) / 1000;
    ppfac   = @(T,Pbar,a,b)... 
              -(( -dV(T,a)    * Pbar  + 0.5 * Ka(T,b)   * Pbar^2 ) / RT )   / ( LOG10 );
    ppfac_T = @(T,Pbar,a,b)...
              -(( -dV_T(T,a) * Pbar   + 0.5 * Ka_T(T,b) * Pbar^2 ) / RT )   / ( LOG10 ) ...
              -(-( -dV(T,a)   * Pbar  + 0.5 * Ka(T,b)   * Pbar^2 ) * RT_T / RT^2 ) / ( LOG10 );        
    ppfac_P = @(T,Pbar,Pbar_P,a,b)... 
              -(( -dV(T,a)         +    Ka(T,b)  * Pbar ) * Pbar_P / RT) / ( LOG10 ) ;
    
    % compute the pK's and their derivatives w.r.t. T,P,and S -------------
    [pp2f    , gpp2f , epp2f ] = calc_p2f(opt,T,RT,RT_T,Pbar,Pbar_P); 
    [pKs     , gpKs  , epKs  ] = calc_pKs(opt,T,S,Pbar,Pbar_P); 
    [pKf     , gpKf  , epKf  ] = calc_pKf(opt,T,S,Pbar,Pbar_P); 
    [pSWS2tot, gpSWS2tot     ] = calc_pSWS2tot(opt,S,pKs,gpKs,pKf,gpKf);
    [pfH     , gpfH  , epfH  ] = calc_pfH(opt,T,S);
    [pK0     , gpK0  , epK0  ] = calc_pK0(opt,T,RT,RT_T,S,Pbar,Pbar_P); 
    [pKb     , gpKb  , epKb  ] = calc_pKb(opt,T,S,Pbar,Pbar_P,pSWS2tot,gpSWS2tot,pfH,gpfH); 
    [pKw     , gpKw  , epKw  ] = calc_pKw(opt,T,S,Pbar,Pbar_P); 
    [pKp1    , gpKp1 , epKp1 ] = calc_pKp1(opt,T,S,Pbar,Pbar_P,pfH,gpfH); 
    [pKp2    , gpKp2 , epKp2 ] = calc_pKp2(opt,T,S,Pbar,Pbar_P,pfH,gpfH); 
    [pKp3    , gpKp3 , epKp3 ] = calc_pKp3(opt,T,S,Pbar,Pbar_P,pfH,gpfH); 
    [pKsi    , gpKsi , epKsi ] = calc_pKsi(opt,T,S,Pbar,Pbar_P,pfH,gpfH); 
    [pK1     , gpK1  , epK1  ] = calc_pK1(opt,T,S,Pbar,Pbar_P,pfH,gpfH,pSWS2tot,gpSWS2tot); 
    [pK2     , gpK2  , epK2  ] = calc_pK2(opt,T,S,Pbar,Pbar_P,pfH,gpfH,pSWS2tot,gpSWS2tot);
    [pKnh4   , gpKnh4, epKnh4] = calc_pKnh4(opt,T,S,Pbar,Pbar_P,pSWS2tot,gpSWS2tot);
    [pKh2s   , gpKh2s, epKh2s] = calc_pKh2s(opt,T,S,Pbar,Pbar_P,pSWS2tot,gpSWS2tot);
    [pKar    , gpKar , epKar ] = calc_pKar(opt,T,S,Pbar,Pbar_P,pfH,gpfH);
    [pKca    , gpKca , epKca ] = calc_pKca(opt,T,S,Pbar,Pbar_P,pfH,gpfH);


    % pressure correction for Ks (Millero, 1995) --------------------------
    a = [ -18.03; 0.0466; 0.000316 ];
    b = [-4.53; 0.09 ];
    pKs     = pKs + ppfac(T,Pbar,a,b);
    pKs_T   = gpKs(1); % isn't this needed?
    pKs_P   = gpKs(3);
    pKs_T   = pKs_T + ppfac_T(T,Pbar,a,b);
    pKs_P   = pKs_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKs(1) = pKs_T;
    gpKs(3) = pKs_P;

    % pressure correction for Kf (Millero, 1995) --------------------------
    a = [ -9.78; -0.009; -0.000942 ];
    b = [ -3.91; 0.054 ];       
    pKf     = pKf + ppfac(T,Pbar,a,b);
    pKf_T   = gpKf(1);
    pKf_P   = gpKf(3);
    pKf_T   = pKf_T + ppfac_T(T,Pbar,a,b);
    pKf_P   = pKf_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKf(1) = pKf_T;
    gpKf(3) = pKf_P;

    % pressure correction for Kb (Millero, 1979) --------------------------
    a = [ -29.48; 0.1622; -0.002608 ];
    b = [  -2.84;   0.0 ];     
    pKb     = pKb + ppfac(T,Pbar,a,b);
    pKb_T   = gpKb(1);
    pKb_P   = gpKb(3);
    pKb_T   = pKb_T + ppfac_T(T,Pbar,a,b);
    pKb_P   = pKb_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKb(1) = pKb_T;
    gpKb(3) = pKb_P;

    % pressure correction for Kw (Millero, 1983) --------------------------
    a = [ -20.02; 0.1119; -0.001409];
    b = [ -5.13; 0.0794 ];
    pKw     = pKw + ppfac(T,Pbar,a,b);
    pKw_T   = gpKw(1);
    pKw_P   = gpKw(3);
    pKw_T   = pKw_T + ppfac_T(T,Pbar,a,b);
    pKw_P   = pKw_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKw(1) = pKw_T;
    gpKw(3) = pKw_P;

    % pressure correction for Kp1 (Millero, 1995; same as Millero, 1983) --
    a = [ -14.51; 0.1211; -0.000321 ];
    b = [  -2.67; 0.0427 ];
    pKp1        = pKp1 + ppfac(T,Pbar,a,b);
    pKp1_T      = gpKp1(1);
    pKp1_P      = gpKp1(3);
    pKp1_T      = pKp1_T + ppfac_T(T,Pbar,a,b);
    pKp1_P      = pKp1_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKp1(1)    = pKp1_T;
    gpKp1(3)    = pKp1_P;
        
    % pressure correction for Kp2 (Millero, 1995; same as Millero, 1983) --
    a = [ -23.12; 0.1758; -0.002647 ];
    b = [ -5.15; 0.09 ]; 
    pKp2        = pKp2 + ppfac(T,Pbar,a,b);
    pKp2_T      = gpKp2(1);
    pKp2_P      = gpKp2(3);
    pKp2_T      = pKp2_T + ppfac_T(T,Pbar,a,b);
    pKp2_P      = pKp2_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKp2(1)    = pKp2_T;
    gpKp2(3)    = pKp2_P;
        
    % pressure correction for Kp3 (Millero, 1995; same as Millero, 1983) --
    a = [ -26.57; 0.202; -0.003042 ];
    b = [ -4.08; 0.0714 ];
    pKp3        = pKp3 + ppfac(T,Pbar,a,b);
    pKp3_T      = gpKp3(1);
    pKp3_P      = gpKp3(3);
    pKp3_T      = pKp3_T + ppfac_T(T,Pbar,a,b);
    pKp3_P      = pKp3_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKp3(1)    = pKp3_T;
    gpKp3(3)    = pKp3_P;

    % pressure correction for Ksi 
    % (Millero, 1995; used the values from boric acid)
    a = [ -29.48; 0.1622; -0.002608 ];
    b =[ -2.84; 0];     
    pKsi        = pKsi + ppfac(T,Pbar,a,b);
    pKsi_T      = gpKsi(1);
    pKsi_P      = gpKsi(3);
    pKsi_T      = pKsi_T + ppfac_T(T,Pbar,a,b);
    pKsi_P      = pKsi_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKsi(1)    = pKsi_T;
    gpKsi(3)    = pKsi_P;

    % pressure correction for K1 (Millero, 1995) --------------------------
    % only for opt.cK1K2 ~=6 & ~=7 & ~=8
    a = [-25.5; 0.1271; 0];
    b = [ -3.08; 0.0877 ];
    pK1     = pK1 + ppfac(T,Pbar,a,b);
    pK1_T   = gpK1(1);
    pK1_P   = gpK1(3);
    pK1_T   = pK1_T + ppfac_T(T,Pbar,a,b);
    pK1_P   = pK1_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpK1(1) = pK1_T;
    gpK1(3) = pK1_P;

    % pressure correction for K2 (Millero, 1995) --------------------------
    % only for opt.cK1K2 ~=6 & ~=7 & ~=8
    a = [ -15.82; -0.0219; 0 ];
    b = [ 1.13; -0.1475 ];
    pK2     = pK2 + ppfac(T,Pbar,a,b);
    pK2_T   = gpK2(1);
    pK2_P   = gpK2(3);
    pK2_T   = pK2_T + ppfac_T(T,Pbar,a,b);
    pK2_P   = pK2_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpK2(1) = pK2_T;
    gpK2(3) = pK2_P;

    % pressure correction for Knh4 (added to CO2SYSv3 by J. Sharp) --------
    a = [ -26.43; 0.0889; -0.000905 ];
    b = [ -5.03; 0.0814 ];
    pKnh4       = pKnh4 + ppfac(T,Pbar,a,b);
    pKnh4_T     = gpKnh4(1);
    pKnh4_P     = gpKnh4(3);
    pKnh4_T     = pKnh4_T + ppfac_T(T,Pbar,a,b);
    pKnh4_P     = pKnh4_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKnh4(1)   = pKnh4_T;
    gpKnh4(3)   = pKnh4_P;
        
    % pressure correction for Kh2s (added to CO2SYSv3 by J. Sharp) --------
    a = [ -11.07; -0.009; -0.000942 ];
    b = [ -2.89;  0.054 ]; 
    pKh2s       = pKh2s + ppfac(T,Pbar,a,b);
    pKh2s_T     = gpKh2s(1);
    pKh2s_P     = gpKh2s(3);
    pKh2s_T     = pKh2s_T + ppfac_T(T,Pbar,a,b);
    pKh2s_P     = pKh2s_P + ppfac_P(T,Pbar,Pbar_P,a,b);
    gpKh2s(1)   = pKh2s_T;
    gpKh2s(3)   = pKh2s_P;
        
    % correct pH scale conversion factors for pressure --------------------
    [pSWS2tot, gpSWS2tot, pFREE2tot, gpFREE2tot] = calc_pSWS2tot(opt,S,pKs,gpKs,pKf,gpKf);

    % find pH scale conversion factor -------------------------------------
    % (pressure corrected)
    if opt.phscale == 1 % pH_total
        phfac = pSWS2tot;
        gphfac = gpSWS2tot;
    elseif opt.phscale == 2 % pH_SWS, they are all on this now
        phfac = 0;
        gphfac = 0;
    elseif opt.phscale == 3 % pH_free
        phfac = pSWS2tot - pFREE2tot;
        gphfac = gpSWS2tot - gpFREE2tot;
    elseif opt.phscale == 4 % pH_NBS
        phfac = pfH;
        gphfac = gpfH;
    else 
        error('Need to input valid pH scale 1-4');
    end
    

    % convert from SWS to chosen pH scale ---------------------------------    
    pK1   = pK1   + phfac;  gpK1   = gpK1   + gphfac;
    pK2   = pK2   + phfac;  gpK2   = gpK2   + gphfac;
    pKw   = pKw   + phfac;  gpKw   = gpKw   + gphfac;
    pKb   = pKb   + phfac;  gpKb   = gpKb   + gphfac; 
    pKp1  = pKp1  + phfac;  gpKp1  = gpKp1  + gphfac;
    pKp2  = pKp2  + phfac;  gpKp2  = gpKp2  + gphfac;
    pKp3  = pKp3  + phfac;  gpKp3  = gpKp3  + gphfac;
    pKsi  = pKsi  + phfac;  gpKsi  = gpKsi  + gphfac;
    pKnh4 = pKnh4 + phfac;  gpKnh4 = gpKnh4 + gphfac;
    pKh2s = pKh2s + phfac;  gpKh2s = gpKh2s + gphfac;
    % pKar, pKca, pKs, and pKf do not need the conversion

   
    % ---------------------------------------------------------------------
    % output
    % ---------------------------------------------------------------------
    pK   = [  pK0;   pK1;   pK2;    pKb;    pKw;   pKs;   pKf;  pKp1; ...
             pKp2;  pKp3;  pKsi;  pKnh4;  pKh2s;  pp2f;  pKar;  pKca; pfH];

    pK = pK .* opt.mpk(T,S,P);
    
    gpK  = [ gpK0;  gpK1;  gpK2;   gpKb;   gpKw;  gpKs;  gpKf; gpKp1; ...
            gpKp2; gpKp3; gpKsi; gpKnh4; gpKh2s; gpp2f; gpKar; gpKca; gpfH];

    gpK = pK .* opt.gmpk(T,S,P) + gpK.*opt.mpk(T,S,P);
    
    epK  = [ epK0;  epK1;  epK2;   epKb;   epKw;  epKs;  epKf; epKp1; ...
            epKp2; epKp3; epKsi; epKnh4; epKh2s; epp2f; epKar; epKca; epfH];
    epK = epK.*opt.empk;
    
    % pK0   = 1,  pK1  = 2,  pK2  = 3,  pKb  = 4,  pKw  = 5,  pKs   = 6, 
    % pKf   = 7,  pKp1 = 8,  pKp2 = 9,  pKp3 = 10, pKsi = 11, pKnh4 = 12, 
    % pKh2s = 13, pp2f = 14, pKar = 15, pKca = 16, pfH  = 17 

    % ---------------------------------------------------------------------
    % subfunctions
    % ---------------------------------------------------------------------

    function [pp2f,gpp2f,epp2f] = calc_p2f(opt,T,RT,RT_T,Pbar,Pbar_P)
    % pCO2 to fCO2 conversion, Weiss 1974 Marine Chemistry, 2: 203-215
    % pg. 207 -> valid to within 0.1% (assuming 1 sigma -MF)
        TK = T + 273.15; % convert to Kelvin
        Pstd = 1.01325;
        delC = (57.7 - 0.118.*TK);
        delC_T = -0.118;

        a = [ -1636.75; 12.0408; -0.0327957; 3.16528e-5 ]; 
        
        f  = @(T,a) a(1) + a(2) * T +     a(3) * T^2 +     a(4) * T^3;
        df = @(T,a)        a(2)     + 2 * a(3) * T   + 3 * a(4) * T^2;
  
        if opt.co2press == 0
            % FugFac at 1 atm
            pp2f   = -( ( f(TK,a) + 2 * delC ) * Pstd / RT ) / LOG10;
            pp2f_T = -( ( df(TK,a) + 2 * delC_T ) * Pstd / RT - ...
                    ( f(TK,a) + 2 * delC ) * Pstd * RT_T / RT^2 ) / LOG10;
            pp2f_S = 0;
            pp2f_P = 0;

        elseif opt.co2press == 1
           % FugFac at in situ pressure
           pp2f   = -( ( f(TK,a) + 2 * delC ) * (Pstd+Pbar) / RT ) / LOG10;
           pp2f_T = -( ( df(TK,a) + 2 * delC_T ) * (Pstd+Pbar) / RT - ...
                   ( f(TK,a) + 2 * delC ) * (Pstd+Pbar) * RT_T / RT^2 ) / LOG10;
           pp2f_S = 0;
           pp2f_P = -( ( f(TK,a) + 2 * delC )*Pbar_P / RT ) / LOG10;
        end
        gpp2f = [pp2f_T, pp2f_S, pp2f_P]; % gradient of p(p2f);

        p2f = q(pp2f);
        ep2f = 0.001 * p2f; % 0.1% relative uncertainty
        my_abs = @(x) sqrt(x*x);
        epp2f = my_abs( p(p2f + ep2f) - pp2f );
    end
    
    function [pK0,gpK0,epK0] = calc_pK0(opt,T,RT,RT_T,S,Pbar,Pbar_P)
    % calculate K0, Weiss 1974, Marine Chemistry, 2: 203-215
    % "data show a root-mean-square deviation from the final fitted
    % equation of1.4*10^-4 mols/l*atm in K0, or about 0.3%"
    % (assuming 0.3% for 1 sigma -MF)
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
        
        pK0_P = 0; % no pressure correction at 1atm

        if opt.co2press == 1   % 1.01325 Bar = 1 atm
            vCO2 = 32.3; % partial molal volume of CO2 (cm3 / mol)
                         % from Weiss (1974, Appendix, paragraph 3)

            % pressure correction at in situ pressure
            ppfacK0   = - ( ((-Pbar)*vCO2) / RT ) / LOG10 ;
            ppfacK0_T = - ( ((-Pbar)*vCO2) * (-RT_T) / RT^2 ) / LOG10 ;
            ppfacK0_P = - ( (-vCO2 *Pbar_P / RT ) ) / LOG10;

            pK0   = pK0 + ppfacK0 ;
            pK0_T = pK0_T + ppfacK0_T;
            pK0_P = pK0_P + ppfacK0_P ;
        end
        % derivatives aka gradient
        gpK0 = [pK0_T, pK0_S, pK0_P];

        % Weiss (1974) reports 0.2 - 0.3% uncertainty on K0
        K0 = q(pK0);
        eK0 = K0 * 0.003; % 0.3% relative on K0
        my_abs = @(x) sqrt(x*x);
        epK0 = my_abs( p(K0 + eK0) - pK0 );
    end
    
    function [pKs,gpKs,epKs] = calc_pKs(opt,T,S,Pbar,Pbar_P)
        TK = T + 273.15; % convert to Kelvin
        IonS = ions(S);
        IonS_S = ions_S(S);
        sqrtIonS_S = 0.5*IonS_S/sqrt(IonS);

        % all calculated on free pH scale
        % stay on free pH scale (no conversion to SWS or total)
        if opt.KSO4 == 1
            % calculate Ks (Dickson 1990a)------------------------------
            % "goodness of fit: 0.021" (assuming 1 sigma -MF)
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
                   f(TK,a5) * IonS_S * 2.0 * IonS ) / ...
                    LOG10 - 0.001005 * dpdx(1 - 0.001005 * S );
        
            elnKs = 0.021; % from Dickson 1990a pg. 123
            my_abs = @(x) sqrt(x*x);
            epKs = my_abs( -elnKs/LOG10 ); % 0.021 on lnKs, -lnK/LOG10 converts to epK
        elseif opt.KSO4 == 2
            % calculate Ks (Khoo et al 1977)----------------------------
            % pg. 33 "the standard deviation from regression is 0.0021 in
            % log Beta_HSO4"
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

            epKs = 0.0021; % given in CO2SYS from Khoo et al 1977
        elseif opt.KSO4 == 3
            % calculate Ks (Waters and Millero, 2013) ---------------------
            % log(K^(.)_HSO4) = logKS0
            % log((K^(*)_HSO4)/(K^(.)_HSO4)) = logKSK0
            % KS = (10^logKSK0))*(10^logKS0)
            % (^ from Sharp's code)
            % eKs = 0.007, assuming 95% confidence (2 sigma) -MF
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
            logKSK0_T = ( df3(TK,a3) + df2(TK,a4) ) + ...
                          df3(TK,a5) + ...
                          df3(TK,a6) ;
            logKSK0_S = ( f3(TK,a3) + f2(TK,a4) ) * (-0.5/sqrt(S)) + ...
                        ( f3(TK,a6) * 0.5 * sqrt(S)) + ...
                        ( a7 * 2 * S );

            pKs_T = (-logKS0_T - logKSK0_T) ;
            pKs_S = (-logKSK0_S) - 0.001005 * dpdx(1-0.001005*S);

            Ks = q(pKs);
            eKs = 0.007/2; % QUODcarb uses 1sigma
            my_abs = @(x) sqrt(x*x);
            epKs = my_abs( p(Ks + eKs) - pKs );
        end
        pKs_P = 0.0;
        gpKs = [pKs_T, pKs_S, pKs_P];
    end
    
    function [pKf,gpKf,epKf] = calc_pKf(opt,T,S,Pbar,Pbar_P)
        TK = T + 273.15; % convert to Kelvin
        IonS = ions(S);
        IonS_S = ions_S(S);
        sqrtIonS_S = 0.5*IonS_S/sqrt(IonS);

        % all calculated on free pH scale
        % stay on free pH scale (no conversion to SWS or total)
        if opt.KF == 1
        % calculate Kf (Dickson and Riley 1979)----------------------------   
            a = 1590.2; b = -12.641; c = 1.525;

            pKf = -( a/TK + b + c * sqrt(IonS) ) ./LOG10 ...
                + p(1 - 0.001005 * S) ; % free pH scale
            pKf_T = -(-a/(TK.^2) ) / LOG10; 
            pKf_S =  -c * sqrtIonS_S / LOG10 + dpdx(1 - 0.001005 * S)...
                * (-0.001005) ;  
        elseif opt.KF == 2
        % calculate Perez and Fraga (1987)---------------------------------
        % (to be used for S: 10-40, T: 9-33) (from Sharp's code)
        % "observed experimental error was not greater than pm 0.05 units
        % of lnBeta_HF" pg. 163 (assume 1sigma -MF)
            a = 874;  b = -9.68;  c = 0.111;

            pKf = -( a/TK + b + c * sqrt(S) ) ./LOG10 ; % free pH scale
            pKf_T = -(-a/(TK.^2) ) / LOG10; 
            pKf_S = ( -c * (0.5/sqrt(S)) )/ LOG10 ;  
        end
        pKf_P = 0.0;
        gpKf = [pKf_T, pKf_S, pKf_P];

        elnKf = 0.05; % 0.05 on lnKf
        my_abs = @(x) sqrt(x*x);
        epKf = my_abs( -elnKf/ LOG10 ) ; %  -/LOG10 converts to epK 
        % none found in Dickson's, so take Perez and Fraga's
    end
    
        
    function [pKb,gpKb,epKb] = calc_pKb(opt,T,S,Pbar,Pbar_P,pSWS2tot,gpSWS2tot,pfH,gpfH)   
        TK = T + 273.15; % convert to Kelvin

        if opt.K1K2 == 6 || opt.K1K2 == 7
        % GEOSECS Peng et al.
        % Lyman, John, UCLA Thesis, 1957
        % fit by Li et al, JGR 74:5507-5525, 1969
        % (copied from Orr's code)
            pKb = -( -9.26 + 0.00886 * S + 0.01*T) ; % our 'p' is -log10
            pKb_T = -0.01;
            pKb_S = -0.00886;
            pKb_P = 0;

        % convert pKb from NBS scale to SWS scale
            pKb = pKb - pfH ;
            pKb_T = pKb_T - gpfH(1);
            pKb_S = pKb_S - gpfH(2);

        else
        % calculate Kb (Dickson 1990b)--------------------------------------
            a1 = [  -8966.9; -2890.53; -77.942; 1.728; -0.0996 ];
            a2 = [ 148.0248; 137.1942; 1.62142;  0.0 ;    0.0  ];
            a3 = [ -24.4344;  -25.085; -0.2474;  0.0 ;    0.0  ];
            a4 = [   0.0   ; 0.053105;   0.0  ;  0.0 ;    0.0  ];

            f  = @(S,a) a(1) + (a(2) .* sqrt(S)) + (a(3) .* S) + ...
                (a(4) .* S.^(1.5)) + (a(5) .* (S.^2));
            df = @(S,a)  0.5 .* a(2) .* S^(-0.5) + a(3) + ...
                1.5 .* a(4) .* sqrt(S) + 2 .* a(5) .* S;

            pKb   = -(  (f(S,a1) ./ TK) + ...
                f(S,a2) + (f(S,a3) .* log(TK)) + ...
                (f(S,a4) .* TK) ) ./ LOG10;
            pKb   = pKb - pSWS2tot; % convert from total to SWS pH scale
            pKb_T = -( (-f(S,a1) ./ (TK.^2) ) + ...
                ( f(S,a3) ./ TK ) +  f(S,a4) ) ./ LOG10;
            pKb_T = pKb_T - gpSWS2tot(1);
            pKb_S = -( ( df(S,a1) ./ TK ) + df(S,a2) + ...
                (df(S,a3) .* log(TK)) + (df(S,a4) .* TK) ) ./ LOG10;
            pKb_S = pKb_S - gpSWS2tot(2);
            pKb_P = 0;
        end
        gpKb = [pKb_T, pKb_S, pKb_P];

        elnKb = 0.004; % pg 764 Dickson (1990b) (assume 1 sigma -MF)
        my_abs = @(x) sqrt(x*x);
        epKb = my_abs( -elnKb/LOG10 ) ; % convert from lnKb to pKb with -/LOG10
        % none found in Li et al's paper
    end
    
    function [pKw,gpKw,epKw] = calc_pKw(opt,T,S,Pbar,Pbar_P)
        TK = T + 273.15; % convert to Kelvin

        if opt.K1K2 == 7
        % calculate Kw (Millero, 1979)-------------------------------------
        % temperature dependent parameters only
            a1 = [ 148.9802; -13847.26; -23.6521 ];
            a2 = [ -79.2447;   3298.72;  12.0408 ];
            a3 = [ -0.019813;      0.0;    0.0   ];

            f  = @(T,a)   a(1) + a(2) / T    + a(3) * log(T);
            df = @(T,a)        - a(2) / T^2  + a(3) / T;

            pKw   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
            pKw_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
            pKw_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;
            pKw_P = 0;

            elnKw = 0.0214; % combine Harned and Owen (1958) 0.0014 (table 1)
                            % plus Culberson and Pytkowicz (1973) 0.020 (table 3)
            my_abs = @(x) sqrt(x*x);
            epKw = my_abs( -elnKw/LOG10 ) ; % convert to epK with -/LOG10
        elseif opt.K1K2 == 8
        % calculate Kw (Millero 1979)--------------------------------------
        % refit data of Harned and Owen, 1958
            a1 = [ 148.9802; -13847.26; -23.6521 ];

            f  = @(T,a)   a(1) + a(2) / T    + a(3) * log(T);
            df = @(T,a)        - a(2) / T^2  + a(3) / T;

            pKw   = -(  f(TK,a1) ) / LOG10;
            pKw_T = -( df(TK,a1) ) / LOG10;
            pKw_S = -(   0   ) / LOG10;
            pKw_P = 0;
            elnKw = 0.0014; % table 1
            my_abs = @(x) sqrt(x*x);
            epKw = my_abs( -elnKw/LOG10 ) ; % convert to epK with -/LOG10            
        else
        % calculate Kw (Millero, 1995)-------------------------------------
            a1 = [ 148.9802; -13847.26; -23.6521 ];
            a2 = [   -5.977;    118.67;   1.0495 ];
            a3 = [ -0.01615;     0.0  ;     0.0  ];

            f  = @(T,a)   a(1) + a(2) / T    + a(3) * log(T);
            df = @(T,a)        - a(2) / T^2  + a(3) / T;

            pKw   = -( f(TK,a1) + ...
                       f(TK,a2) * sqrt(S) + ...
                       f(TK,a3) * S ) / LOG10;
            % pKw   =  f(TK,a2) * (S)  ;
            pKw_T = -( df(TK,a1) + ...
                       df(TK,a2) * sqrt(S) + ...
                       df(TK,a3) * S ) / LOG10;
            pKw_S = -( ( f(TK,a2) .* ( 0.5 ./ sqrt(S) ) ) + ...
                       f(TK,a3) ) ./ LOG10;
            % pKw_S = (  f(TK,a2)  );
            pKw_P = 0;

            elnKw = 0.01; % pg 670 Millero, 1995
            my_abs = @(x) sqrt(x*x);
            epKw = my_abs( -elnKw/LOG10 ) ; % convert to epK with -/LOG10
        end 

        gpKw = [pKw_T,pKw_S,pKw_P];
    end
    
    function [pKp1,gpKp1,epKp1] = calc_pKp1(opt,T,S,Pbar,Pbar_P,pfH,gpfH)
        TK = T + 273.15; % convert to Kelvin

        if opt.K1K2 == 7
        % calculate Kp1----------------------------------------------------
        % Peng et al don't include the contribution from this term,
        % but it is so small it doesn't contribute. It needs to be kept
        % so that the routines work ok. (From Orr's code)
            pKp1 = p(0.02) - pfH; % Kp1 = 0.02, convert NBS to SWS pH scale
            pKp1_T = 0 - gpfH(1); 
            pKp1_S = 0 - gpfH(2);        
            pKp1_P = 0;
        else
        % calculate Kp1 (Yao and Millero 1995)-----------------------------       
            a1 = [ -4576.752;    115.54;  -18.453 ];
            a2 = [  -106.736;   0.69171;     0.0  ];
            a3 = [  -0.65643;  -0.01844;     0.0  ];

            f  = @(T,a)   a(1) / T   + a(2) + a(3) * log(T);
            df = @(T,a) - a(1) / T^2        + a(3) / T;

            pKp1   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
            pKp1_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
            pKp1_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;
            pKp1_P = 0;
        end
        gpKp1 = [pKp1_T,pKp1_S,pKp1_P];

        epKp1 = 0.09; % pg 84 Yao and Millero, 1995 (assume 1 sigma -MF)
    end
    
    function [pKp2,gpKp2,epKp2] = calc_pKp2(opt,T,S,Pbar,Pbar_P,pfH,gpfH)
        TK = T + 273.15; % convert to Kelvin
        if opt.K1K2 == 7

        % calculate Kp2 (Kester and Pytkowicz, 1967)-----------------------
            pKp2 = -(-9.039 - 1450/TK) / LOG10  ; 
            pKp2 = pKp2 - pfH ; % convert from NBS to SWS pH scale

            pKp2_T = -( 1450 / TK^2 ) / LOG10 - gpfH(1) ;
            pKp2_S = -gpfH(2) ;
            pKp2_P = -gpfH(3) ;
        else
        % calculate Kp2 (Yao and Millero 1995)-----------------------------
            a1 = [ -8814.715;  172.1033; -27.927 ]; 
            a2 = [   -160.34;    1.3566;    0.0  ];
            a3 = [   0.37335;  -0.05778;    0.0  ];
        
            f  = @(T,a)   a(1) / T   + a(2) + a(3) * log(T);
            df = @(T,a) - a(1) / T^2        + a(3) / T;

            pKp2   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S)       +  f(TK,a3) * S ) / LOG10;
            pKp2_T = -( df(TK,a1) + df(TK,a2) * sqrt(S)       + df(TK,a3) * S ) / LOG10;
            pKp2_S = -(              f(TK,a2) * 0.5 / sqrt(S) +  f(TK,a3)     ) / LOG10;
            pKp2_P = 0 ;
        end
        gpKp2 = [pKp2_T, pKp2_S,pKp2_P];

        epKp2 = 0.03; % pg 84 Yao and Millero, 1995 (assume 1 sigma -MF)
    end
    
    function [pKp3,gpKp3,epKp3] = calc_pKp3(opt,T,S,Pbar,Pbar_P,pfH,gpfH)
        TK = T + 273.15; % convert to Kelvin

        if opt.K1K2 == 7
        % calculate Kp3 (Kester and Pytkowicz, 1967)-----------------------
            pKp3 = -(4.466 - 7276/TK) / LOG10;
            pKp3 = pKp3 - pfH ; % convert from NBS to SWS pH scale

            pKp3_T = -( 7276 / TK^2 ) / LOG10 - gpfH(1) ;
            pKp3_S = -gpfH(2) ;
            pKp3_P = -gpfH(3) ;    
        else
        % calculate Kp3 (Yao and Millero 1995)-----------------------------
            a1 = [    -3070.75; -18.126  ];
            a2 = [    17.27039;  2.81197 ];
            a3 = [   -44.99486; -0.09984 ];
                
            f  = @(T,a)   a(1) / T   + a(2);
            df = @(T,a) - a(1) / T^2;
        
            pKp3   = -(  f(TK,a1) +  f(TK,a2) * sqrt(S) +  f(TK,a3) * S ) / LOG10;
            pKp3_T = -( df(TK,a1) + df(TK,a2) * sqrt(S) + df(TK,a3) * S ) / LOG10;
            pKp3_S = -(              f(TK,a2) * 0.5 / sqrt(S) + f(TK,a3) ) / LOG10;
            pKp3_P = 0 ;   
        end
        gpKp3 = [pKp3_T, pKp3_S, pKp3_P]; 

        epKp3 = 0.2; % pg 84 Yao and Millero, 1995 (assume 1 sigma -MF)
    end

    
    function [pKsi,gpKsi,epKsi] = calc_pKsi(opt,T,S,Pbar,Pbar_P,pfH,gpfH)
        TK = T + 273.15; % convert to Kelvin
        IonS = ions(S);
        IonS_S = ions_S(S);
        sqrtIonS_S = 0.5 * IonS_S / sqrt(IonS);

        if opt.K1K2 == 7
        % calculate Ksi (Sillen, Martell, and Bjerrum, 1964)---------------
            % Ksi = 4e-10; from CO2SYS
            pKsi = p(4e-10) - pfH ; % also convert NBS scale to SWS pH scale
            pKsi_T = -gpfH(1);
            pKsi_S = -gpfH(2);
            pKsi_P = 0;
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
            pKsi_P = 0;
        end
        gpKsi = [pKsi_T, pKsi_S, pKsi_P];

        epKsi = 0.02; % pg 84 Yao and Millero, 1995 (assume 1 sigma -MF)
    end
    
    function [pK1,gpK1,epK1] = calc_pK1(opt,T,S,Pbar,Pbar_P,pfH,gpfH,pSWS2tot,gpSWS2tot)
        % calculate pK1 based on user choice input with opt.cK1K2
        TK = T + 273.15; % convert to Kelvin

        if opt.K1K2 == 1
        % Roy et al, Marine Chemistry, 44:249-267, 1993 -------------------
            % (see also: Erratum, Marine Chemistry 45:337, 1994
            % and Erratum, Marine Chemistry 52:183, 1996)
            % Typo: in the abstract on p. 249: in the eq. for lnK1* the
            % last term should have S raised to the power 1.5.
            % They claim standard deviations (p. 254) of the fits as
            % .0048 for lnK1 (.5% in K1) and .007 in lnK2 (.7% in K2).
            % They also claim (p. 258) 2s precisions of .004 in pK1 and
            % .006 in pK2. These are consistent, but Andrew Dickson
            % (personal communication) obtained an rms deviation of about
            % .004 in pK1 and .003 in pK2. This would be a 2s precision
            % of about 2% in K1 and 1.5% in K2.
            % T:  0-45  S:  5-45. Total Scale. Artificial seawater.
            % This is eq. 29 on p. 254 and what they use in their abstract:    
            a1 = [ 2.83655     ; -2307.1266; -1.5529413];
            a2 = [ -0.20760841 ; -4.0484   ;  0.0      ];
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
            pK1_P = 0;
            % pass all pK1 out of this function as SWS scale

            epK1 = 0.004/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 2 
        % Goyet and Poisson, Deep-Sea Research, 36(11):1635-1654, 1989 ----
            % The 2s precision in pK1 is .011, or 2.5% in K1.
            % The 2s precision in pK2 is .02, or 4.5% in K2.
            % This is in Table 5 on p. 1652 and what they use in the abstract:
            a = 812.27; b = 3.356; c = -0.00171; d = 0.000091;
            pK1 = a ./ TK + b + c .* S .* log(TK) + d .* (S.^2) ; % SWS scale
            pK1_T = - a ./ (TK.^2) + c .* S ./ TK ;
            pK1_S = c .* log(TK) + 2 .* d .* S ;
            pK1_P = 0;
            epK1 = 0.011/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 3
        % Hansson refit by Dickson and Millero, 1987 ----------------------
            % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973
            % and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973.
            % on the SWS pH scale in mol/kg-SW.
            % Hansson gave his results on the Total scale (he called it
            % the seawater scale) and in mol/kg-SW.
            % Typo in DM on p. 1739 in Table 4: the equation for pK2*
            % for Hansson should have a .000132 *S^2
            % instead of a .000116 *S^2.
            % The 2s precision in pK1 is .013, or 3% in K1.
            % The 2s precision in pK2 is .017, or 4.1% in K2.
            % This is from Table 4 on p. 1739.
            a = 851.4; b = 3.237; c = -0.0106; d = 0.000105;
            pK1 = a ./ TK + b + c .* S + d .* (S.^2) ; % SWS scale
            pK1_T = - a ./ (TK.^2) ;
            pK1_S = c + 2 .* d .* S ;
            pK1_P = 0;
            epK1 = 0.013/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 4
        % Mehrbach refit by Dickson and Millero, 1987 ---------------------
            % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
            % on the SWS pH scale in mol/kg-SW.
            % Mehrbach et al gave results on the NBS scale.
            % The 2s precision in pK1 is .011, or 2.6% in K1.
            % The 2s precision in pK2 is .020, or 4.6% in K2.
	        % Valid for salinity 20-40.
            % This is in Table 4 on p. 1739.
            a = 3670.7; b = -62.008; c = 9.7944; d = -0.0118; g = 0.000116;
            pK1   =  a ./ TK   + b + c .* log(TK) + d .* S + g .* S^2;
            pK1_T = -a ./ TK^2     + c ./ TK;
            pK1_S = d + 2.*g.*S;
            pK1_P = 0;
            epK1 = 0.011/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 5
        % Hansson and Mehrbach refit by Dickson and Millero, 1987 ---------
            % Dickson and Millero, Deep-Sea Research,34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973,
            % Hansson, Acta Chemica Scandanavia, 27:931-944, 1973,
            % and Mehrbach et al, Limnol. Oceanogr.,18(6):897-907, 1973
            % on the SWS pH scale in mol/kg-SW.
            % Typo in DM on p. 1740 in Table 5: the second equation
            % should be pK2* =, not pK1* =.
            % The 2s precision in pK1 is .017, or 4% in K1.
            % The 2s precision in pK2 is .026, or 6% in K2.
	        % Valid for salinity 20-40.
            % This is in Table 5 on p. 1740.    
            a = 845; b = 3.248; c = -0.0098; d = 0.000087;
            pK1   = a ./ TK + b + c .* S + d .* (S.^2) ;
            pK1_T = -a ./ (TK.^2) ;
            pK1_S = c + 2 .* d .* S ;
            pK1_P = 0;
            epK1 = 0.017/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 6 || opt.K1K2 == 7
        % GEOSECS and Peng et al use K1, K2 from Mehrbach et al, ----------
            % Limnology and Oceanography, 18(6):897-907, 1973.
	        % I.e., these are the original Mehrbach dissociation constants.
            % The 2s precision in pK1 is .005, or 1.2% in K1.
            % The 2s precision in pK2 is .008, or 2% in K2.

            a = -13.7201; b = 0.031334; c = 3235.76; 
            d = 1.3e-5; g = -0.1032;

            pK1 = a + b .* TK + c ./ TK + d .* S .* TK + g .* sqrt(S) ; % NBS scale
            pK1 = pK1 - pfH ; % convert from NBS to SWS pH scale
            pK1_T = ( b - c ./ (TK.^2) + d .* S ) - gpfH(1) ; % SWS scale
            pK1_S = ( d .* TK + 0.5 .* g ./ sqrt(S) ) - gpfH(2) ; % SWS scale
            pK1_P = 0;
            epK1 = 0.005/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 8
        % PURE WATER CASE -------------------------------------------------
            % Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
            % K1 from refit data from Harned and Davis,
            % J American Chemical Society, 65:2030-2037, 1943.
            % K2 from refit data from Harned and Scholes,
            % J American Chemical Society, 43:1706-1709, 1941.
	            % This is only to be used for Sal=0 water 
                % (note the absence of S in the below formulations)
            % (0-50 C)
            a = 290.9097; b = -14554.21; c = -45.0575;
            pK1 = - (a + b ./ TK + c .* log(TK) ) / LOG10;
            pK1_T = - ( -b ./ (TK.^2) + c ./ TK ) / LOG10;
            pK1_S = 0;
            pK1_P = 0;
            elnK1 = 0.0024;
            my_abs = @(x) sqrt(x*x);
            epK1 = my_abs( -elnK1/LOG10 ) ; % convert lnK to epK with -/LOG10

        elseif opt.K1K2 == 9
        % Cai and Wang 1998, for estuarine use ----------------------------
            % Data used in this work is from:
	        % K1: Mehrbach (1973) for S>15, for S<15: Mook and Keone (1975)
	        % K2: Mehrbach (1973) for S>20, for S<20: Edmond and Gieskes (1970)
	        % Sigma of residuals between fits and above data: ±0.015, 
            % +0.040 for K1 and K2, respectively.
	        % Sal 0-40, Temp 0.2-30
            % Limnol. Oceanogr. 43(4) (1998) 657-668    
            a = 3404.71; b = 0.032786; c = -14.8435; d = -0.071692; 
            g = 0.0021487; h = 200.1; l = 0.3220;
            f1 = h ./ TK + l;

            pK1 = a ./ TK + b .* TK + c + d .* f1 .* sqrt(S) + g .* S ;
            pK1 = pK1 - pfH; % convert from NBS scale to SWS scale
            pK1_T = ( -a ./ (TK.^2) + b + d .* sqrt(S) .* (-h ./ (TK.^2)) )...
                - gpfH(1) ; % Temp derivative convert to sws scale
            pK1_S = ( 0.5 .* d .* f1 ./ sqrt(S) + g ) - gpfH(2) ;
            pK1_P = 0;
            epK1 = 0.015;

        elseif opt.K1K2 == 10
        % Leuker, Dickson, Keeling, 2000 ----------------------------------
            % This is Mehrbach's data refit after conversion to the 
            % total scale, for comparison with their equilibrator work. 
            % Mar. Chem. 70 (2000) 105-119
            % rms deviation is 0.0055 in pK1 and 0.0100 in pK2
            a = 3633.86; b = -61.2172; c = 9.6777;
            d = -0.011555; g = 0.0001152;

            pK1 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2);
            pK1 = pK1 - pSWS2tot; % convert from total scale to SWS scale
            pK1_T = ( -a ./ (TK.^2) + c ./ TK ) - gpSWS2tot(1);
            pK1_S = ( d + 2 .* g .* S ) - gpSWS2tot(2);
            pK1_P = 0;
            epK1 = 0.0055;

        elseif opt.K1K2 == 11
        % Mojica Prieto and Millero, 2002 ---------------------------------
            % Geochim. et Cosmochim. Acta. 66(14) 2529-2540
            % sigma for pK1 is reported to be 0.0056
	        % sigma for pK2 is reported to be 0.010
	        % This is from the abstract and pages 2536-2537
            a = -43.6977; b = -0.0129037; c = 1.364e-4;
            d = 2885.378; g = 7.045159;

            pK1 = a + b .* S + c .* (S.^2) + d ./ TK + g .* log(TK) ;
            pK1_T = -d ./ (TK.^2) + g ./ TK ;
            pK1_S = b + 2 .* c .* S ;
            pK1_P = 0;
            epK1 = 0.0056;

        elseif opt.K1K2 == 12
        % Millero et al, 2002 ---------------------------------------------
            % Deep-Sea Res. I (49) 1705-1723.
	        % Calculated from overdetermined WOCE-era field measurements 
	        % sigma for pK1 is reported to be 0.005
	        % sigma for pK2 is reported to be 0.008
	        % This is from page 1716
            a = 6.359; b = -0.00664; c = -0.01322; d = 4.989e-5;
            pK1 = a + b .* S + c .* T + d .* (T.^2) ; % tempC
            pK1_T = c + 2 .* d .* T ;
            pK1_S = b ;
            pK1_P = 0;
            epK1 = 0.005;

        elseif opt.K1K2 == 13
        % Millero et al (2006) --------------------------------------------
            % Millero, Graham, Huang, Bustos-Serrano, Pierrot, 2006
            % Mar.Chem. 100 (2006) 80-94.
            % S=1 to 50, T=0 to 50. On seawater scale (SWS). 
            % From titrations in Gulf Stream seawater.
            % sigma pK1 = 0.0054, sigma pK2 = 0.011, from abstract (-MF)
            a1 = [-126.34048;  13.4191; 0.0331; -5.33e-5] ;
            a2 = [  6320.813; -530.123; -6.103;    0.0  ] ;
            a3 = [ 19.568224; -2.06950;   0.0 ;    0.0  ] ;
            
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK1 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK1_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK1_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);
            pK1_P = 0;
            epK1 = 0.0054; 

        elseif opt.K1K2 == 14
        % Millero 2010, for estuarine use ---------------------------------
            % Marine and Freshwater Research, v. 61, p. 139-142.
	        % Fits through compilation of real seawater titration results:
	        % Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), 
            % Millero et al. (2006)
	        % Constants for K's on the SWS; This is from page 141
            % sigma pK1 = 0.005 & sigma pK2 = 0.010 (-MF)
            a1 = [-126.34048;  13.4038; 0.03206; -5.242e-5] ; 
            a2 = [  6320.813; -530.659; -5.8210;    0.0   ] ;
            a3 = [ 19.568224;  -2.0664;   0.0  ;    0.0   ] ;
                        
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK1 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK1_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK1_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);
            pK1_P = 0;
            epK1 = 0.005; % pg 141

        elseif opt.K1K2 == 15
        % Waters, Millero, and Woosley, 2014 ------------------------------
            % Mar. Chem., 165, 66-67, 2014
            % Corrigendum to "The free proton concentration scale for seawater pH".
	        % Effectively, this is an update of Millero (2010) formulation (WhichKs==14)
            % sigma pK1 = 0.0055 & sigma pK2 = 0.0110 (-MF)
            a1 = [-126.34048;  13.409160; 0.031646; -5.1895e-5] ;
            a2 = [  6320.813;  -531.3642;   -5.713;    0.0    ] ;
            a3 = [ 19.568224; -2.0669166;   0.0   ;    0.0    ] ;

            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK1 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK1_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK1_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);
            pK1_P = 0;
            epK1 = 0.0055;

        elseif opt.K1K2 == 16
        % Sulpis et al, 2020 ----------------------------------------------
            % Ocean Science Discussions, 16, 847-862
            % This study uses overdeterminations of the carbonate system to
            % iteratively fit K1 and K2
            % Relative overall uncertainties ~2.5% (sigK/K) for both
            a = 8510.63; b = -172.4493; c = 26.32996; 
            d = -0.011555; g = 0.0001152; 

            pK1 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2) ;
            pK1 = pK1 - pSWS2tot; % convert from tot to SWS scale
            pK1_T = ( -a ./ (TK.^2) + c ./ TK ) - gpSWS2tot(1);
            pK1_S = ( d + 2 .* g .* S ) - gpSWS2tot(2);
            pK1_P = 0;
            K1 = q(pK1);
            eK1 = 0.025 * K1; % ~2.5 % uncertainty on K, pg 854
            my_abs = @(x) sqrt(x*x);
            epK1 = my_abs( p(K1 + eK1) - pK1 ); 

        elseif opt.K1K2 == 17
        % Schockman and Byrne, 2021 ---------------------------------------
            % Geochimica et Cosmochimica Acta, in press
            % This study uses spectrophotometric pH measurements to determine
            % K1*K2 with unprecedented precision, and presents a new
            % parameterization for K2 based on these determinations
            % K1 is taken from Waters, Millero, and Woosley, 2014, on the total pH scale:
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
            pK1_P = 0;
            epK1 = 0.0055; % same as Waters and Millero formulation

        end

        gpK1 = [pK1_T, pK1_S, pK1_P];     
    end
    
    function [pK2,gpK2,epK2] = calc_pK2(opt,T,S,Pbar,Pbar_P,pfH,gpfH,pSWS2tot,gpSWS2tot)
    % calculate pK2 based on user choice input with opt.cK1K2
        TK = T + 273.15; % convert to Kelvin

        if opt.K1K2 == 1
        % Roy et al, Marine Chemistry, 44:249-267, 1993 ------------------
            % (see also: Erratum, Marine Chemistry 45:337, 1994
            % and Erratum, Marine Chemistry 52:183, 1996)
            % Typo: in the abstract on p. 249: in the eq. for lnK1* the
            % last term should have S raised to the power 1.5.
            % They claim standard deviations (p. 254) of the fits as
            % .0048 for lnK1 (.5% in K1) and .007 in lnK2 (.7% in K2).
            % They also claim (p. 258) 2s precisions of .004 in pK1 and
            % .006 in pK2. These are consistent, but Andrew Dickson
            % (personal communication) obtained an rms deviation of about
            % .004 in pK1 and .003 in pK2. This would be a 2s precision
            % of about 2% in K1 and 1.5% in K2.
            % T:  0-45  S:  5-45. Total Scale. Artificial seawater.
            % This is eq. 29 on p. 254 and what they use in their abstract:
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
            pK2_P = 0;
            % pass all pK2 out of this function as SWS scale

            epK2 = 0.003/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 2 
        % Goyet and Poisson, Deep-Sea Research, 36(11):1635-1654, 1989 ----
            % The 2s precision in pK1 is .011, or 2.5% in K1.
            % The 2s precision in pK2 is .02, or 4.5% in K2.
            % This is in Table 5 on p. 1652 and what they use in the abstract:
            a = 1450.87; b = 4.604; c = -0.00385; d = 0.000182;
            pK2 = a ./ TK + b + c .* S .* log(TK) + d .* (S.^2) ; % SWS scale
            pK2_T = - a ./ (TK.^2) + c .* S ./ TK ;
            pK2_S = c .* log(TK) + 2 .* d .* S ;
            pK2_P = 0;
            epK2 = 0.02/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 3
        % Hansson refit by Dickson and Millero, 1987 ----------------------
            % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973
            % and Hansson, Acta Chemica Scandanavia, 27:931-944, 1973.
            % on the SWS pH scale in mol/kg-SW.
            % Hansson gave his results on the Total scale (he called it
            % the seawater scale) and in mol/kg-SW.
            % Typo in DM on p. 1739 in Table 4: the equation for pK2*
            % for Hansson should have a .000132 *S^2
            % instead of a .000116 *S^2.
            % The 2s precision in pK1 is .013, or 3% in K1.
            % The 2s precision in pK2 is .017, or 4.1% in K2.
            % This is from Table 4 on p. 1739.
            a = -3885.4; b = 125.844; c = -18.141; 
            d = -0.0192; g = 0.000132;

            pK2 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2) ; % SWS scale

            pK2_T = - a ./ (TK.^2) + c ./ TK ;
            pK2_S = d + 2 .* g .* S ;
            pK2_P = 0;
            epK2 = 0.017/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 4
        % Mehrbach refit by Dickson and Millero, 1987 ---------------------
            % Dickson and Millero, Deep-Sea Research, 34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Mehrbach et al, Limn Oc, 18(6):897-907, 1973
            % on the SWS pH scale in mol/kg-SW.
            % Mehrbach et al gave results on the NBS scale.
            % The 2s precision in pK1 is .011, or 2.6% in K1.
            % The 2s precision in pK2 is .020, or 4.6% in K2.
	        % Valid for salinity 20-40.
            % This is in Table 4 on p. 1739.
            a = 1394.7; b = 4.777; c = -0.0184; d = 0.000118;
            pK2 = a ./ TK + b + c .* S + d .* (S.^2);
            pK2_T = -a ./ (TK.^2);
            pK2_S = c + 2 .* d .* S;
            pK2_P = 0;
            epK2 = 0.020/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 5
        % Hansson and Mehrbach refit by Dickson and Millero, 1987 ---------
            % Dickson and Millero, Deep-Sea Research,34(10):1733-1743, 1987
            % (see also Corrigenda, Deep-Sea Research, 36:983, 1989)
            % refit data of Hansson, Deep-Sea Research, 20:461-478, 1973,
            % Hansson, Acta Chemica Scandanavia, 27:931-944, 1973,
            % and Mehrbach et al, Limnol. Oceanogr.,18(6):897-907, 1973
            % on the SWS pH scale in mol/kg-SW.
            % Typo in DM on p. 1740 in Table 5: the second equation
            % should be pK2* =, not pK1* =.
            % The 2s precision in pK1 is .017, or 4% in K1.
            % The 2s precision in pK2 is .026, or 6% in K2.
	        % Valid for salinity 20-40.
            % This is in Table 5 on p. 1740.    
            a = 1377.3; b = 4.824; c = -0.0185; d = 0.000122;
            pK2   = a ./ TK + b + c .* S + d .* (S.^2) ;
            pK2_T = -a ./ (TK.^2) ;
            pK2_S = c + 2 .* d .* S ;
            pK2_P = 0;
            epK2 = 0.026/2; % QUODcarb uses 1sigma
        
        elseif opt.K1K2 == 6 || opt.K1K2 == 7
        % GEOSECS and Peng et al use K1, K2 from Mehrbach et al, ----------
            % Limnology and Oceanography, 18(6):897-907, 1973.
	        % I.e., these are the original Mehrbach dissociation constants.
            % The 2s precision in pK1 is .005, or 1.2% in K1.
            % The 2s precision in pK2 is .008, or 2% in K2.
            a1 = [ 5371.9654; -128375.28;  1.671221 ] ;
            a2 = [  0.22913 ;    2.136  ; -8.0944e-4] ;
            a3 = [  18.3802 ;  -5617.11 ;     0.0   ] ;
            b = -2194.3055;

            f  = @(T,a) a(1) +  a(2) ./ T   + a(3) .* T ;
            df = @(T,a)    - a(2) ./ (T.^2) + a(3) ./ T ;
            
            pK2 = ( f(TK,a1) + f(TK,a2) .* S + f(TK,a3) .* log10(S) + ...
                b .* log10(TK) ) - pfH ; % convert from NBS to SWS scale

            pK2_T = ( df(TK,a1) + df(TK,a2) .* S + df(TK,a3) .* log10(S) + ...
                b ./ (TK .* LOG10) ) - gpfH(1) ;
            pK2_S = ( f(TK,a2) + f(TK,a3) ./ (S .* LOG10) ) - gpfH(2) ;
            pK2_P = 0;
            epK2 = 0.008/2; % QUODcarb uses 1sigma

        elseif opt.K1K2 == 8
        % PURE WATER CASE -------------------------------------------------
            % Millero, F. J., Geochemica et Cosmochemica Acta 43:1651-1661, 1979:
            % K1 from refit data from Harned and Davis,
            % J American Chemical Society, 65:2030-2037, 1943.
            % K2 from refit data from Harned and Scholes,
            % J American Chemical Society, 43:1706-1709, 1941.
	            % This is only to be used for Sal=0 water 
                % (note the absence of S in the below formulations)
            a = 207.6548; b = -11843.79; c = -33.6485;
            pK2 = - (a + b ./ TK + c .* log(TK) ) / LOG10;
            pK2_T = - ( -b ./ (TK.^2) + c ./ TK ) / LOG10;
            pK2_S = 0;
            pK2_P = 0;
            elnK2 = 0.0033;
            my_abs = @(x) sqrt(x*x);
            epK2 = my_abs( -elnK2/LOG10 ) ; % convert lnK to epK with -/LOG10

        elseif opt.K1K2 == 9
        % Cai and Wang 1998, for estuarine use ----------------------------
            % Data used in this work is from:
	        % K1: Mehrbach (1973) for S>15, for S<15: Mook and Keone (1975)
	        % K2: Mehrbach (1973) for S>20, for S<20: Edmond and Gieskes (1970)
	        % Sigma of residuals between fits and above data: ±0.015, 
            % +0.040 for K1 and K2, respectively.
	        % Sal 0-40, Temp 0.2-30
            % Limnol. Oceanogr. 43(4) (1998) 657-668    
            a = 2902.39; b = 0.02379; c = -6.4980; d = -0.3191; 
            g = 0.0198; h = -129.24; l = 1.4381;
            f1 = h ./ TK + l;

            pK2 = a ./ TK + b .* TK + c + d .* f1 .* sqrt(S) + g .* S ;
            pK2 = pK2 - pfH; % convert from NBS scale to SWS scale

            pK2_T = ( -a ./ (TK.^2) + b + d .* sqrt(S) .* (-h ./ (TK.^2)) )...
                - gpfH(1) ; % Temp derivative convert to sws scale
            pK2_S = ( 0.5 .* d .* f1 ./ sqrt(S) + g ) - gpfH(2) ;
            pK2_P = 0;
            epK2 = 0.040;

        elseif opt.K1K2 == 10
        % Leuker, Dickson, Keeling, 2000 ----------------------------------
            % This is Mehrbach's data refit after conversion to the 
            % total scale, for comparison with their equilibrator work. 
            % Mar. Chem. 70 (2000) 105-119
            % rms deviation is 0.0055 in pK1 and 0.0100 in pK2 (-MF)
            a = 471.78; b = 25.929; c = -3.16967;
            d = -0.01781; g = 0.0001122;

            pK2 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2);
            pK2 = pK2 - pSWS2tot; % convert from total scale to SWS scale

            pK2_T = ( -a ./ (TK.^2) + c ./ TK ) - gpSWS2tot(1);
            pK2_S = ( d + 2 .* g .* S ) - gpSWS2tot(2);
            pK2_P = 0;
            epK2 = 0.0100;

        elseif opt.K1K2 == 11
        % Mojica Prieto and Millero, 2002 ---------------------------------
            % Geochim. et Cosmochim. Acta. 66(14) 2529-2540
            % sigma for pK1 is reported to be 0.0056
	        % sigma for pK2 is reported to be 0.010
	        % This is from the abstract and pages 2536-2537
            a1 = [ -452.0940;  21263.61; 68.483143 ] ;
            a2 = [ 13.142162; -581.4428; -1.967035 ] ;
            a3 = [ -8.101e-4;  0.259601;    0.0    ] ;

            f = @(T,a) a(1) + a(2) ./ TK + a(3) .* log(TK) ;
            df = @(T,a) -a(2) ./ (TK.^2) + a(3) ./ TK ;

            pK2 = f(TK,a1) + f(TK,a2) .* S + f(TK,a3) .* (S.^2) ;

            pK2_T = df(TK,a1) + df(TK,a2) .* S + df(TK,a3) .* (S.^2) ;
            pK2_S = f(TK,a2) + 2 .* f(TK,a3) .* S ;
            pK2_P = 0;
            epK2 = 0.010;

        elseif opt.K1K2 == 12
        % Millero et al, 2002 ---------------------------------------------
            % Deep-Sea Res. I (49) 1705-1723.
	        % Calculated from overdetermined WOCE-era field measurements 
	        % sigma for pK1 is reported to be 0.005
	        % sigma for pK2 is reported to be 0.008
	        % This is from page 1715-1716
            a = 9.867; b = -0.01314; c = -0.01904; d = 2.448e-5;
            pK2 = a + b .* S + c .* T + d .* (T.^2) ; % tempC
            pK2_T = c + 2 .* d .* T ;
            pK2_S = b ;
            pK2_P = 0;
            epK2 = 0.008;

        elseif opt.K1K2 == 13
        % Millero et al (2006) --------------------------------------------
            % Millero, Graham, Huang, Bustos-Serrano, Pierrot, 2006
            % Mar.Chem. 100 (2006) 80-94.
            % S=1 to 50, T=0 to 50. On seawater scale (SWS). 
            % From titrations in Gulf Stream seawater.
            % sigma pK1 = 0.0054, sigma pK2 = 0.011, from abstract (-MF)
            a1 = [ -90.18333;  21.0894;  0.1248; -3.687e-4 ] ;
            a2 = [  5143.692; -772.483; -20.051;    0.0    ] ;
            a3 = [ 14.613358;  -3.3336;    0.0 ;    0.0    ] ;
            
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK2 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;

            pK2_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK2_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);
            pK2_P = 0;
            epK2 = 0.011;

        elseif opt.K1K2 == 14
        % Millero 2010, for estuarine use ---------------------------------
            % Marine and Freshwater Research, v. 61, p. 139-142.
	        % Fits through compilation of real seawater titration results:
	        % Mehrbach et al. (1973), Mojica-Prieto & Millero (2002), 
            % Millero et al. (2006)
	        % Constants for K's on the SWS; This is from page 141
            % sigma pK1 = 0.005 & sigma pK2 = 0.010 (-MF)
            a1 = [ -90.18333;  21.3728;  0.1218; -3.688e-4] ; 
            a2 = [  5143.692; -788.289; -19.189;    0.0   ] ;
            a3 = [ 14.613358;   -3.374;   0.0  ;    0.0   ] ;
                        
            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK2 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK2_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK2_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);
            pK2_P = 0;
            epK2 = 0.010; % pg 141

        elseif opt.K1K2 == 15
        % Waters, Millero, and Woosley, 2014 ------------------------------
            % Mar. Chem., 165, 66-67, 2014
            % Corrigendum to "The free proton concentration scale for seawater pH".
	        % Effectively, this is an update of Millero (2010) formulation (WhichKs==14)
            % sigma pK1 = 0.0055 & sigma pK2 = 0.0110 (-MF)
            a1 = [  -90.18333;  21.225890; 0.12450870; -3.7243e-4] ;
            a2 = [  5143.692;  -779.3444;  -19.91739;    0.0    ] ;
            a3 = [ 14.613358; -3.3534679;     0.0   ;    0.0    ] ;

            f  = @(S,a) a(1) + a(2) .* sqrt(S) + a(3) .* S + a(4) .* (S.^2);
            df = @(S,a)  0.5 .* a(2)./ sqrt(S) + a(3) + 2 .* a(4) .* S ;

            pK2 = f(S,a1)  + f(S,a2) ./ TK + f(S,a3) .* log(TK) ;
            pK2_T = -f(S,a2) ./ (TK.^2) + f(S,a3) ./ TK ;
            pK2_S = df(S,a1) + df(S,a2) ./ TK + df(S,a3) .* log(TK);
            pK2_P = 0;
            epK2 = 0.0110;

        elseif opt.K1K2 == 16
        % Sulpis et al, 2020 ----------------------------------------------
            % Ocean Science Discussions, 16, 847-862
            % This study uses overdeterminations of the carbonate system to
            % iteratively fit K1 and K2
            % Relative overall uncertainties ~2.5% (sigK/K) for both (-MF)
            a = 4226.23; b = -59.4636; c = 9.60817; 
            d = -0.01781; g = 0.0001122; 

            pK2 = a ./ TK + b + c .* log(TK) + d .* S + g .* (S.^2) ;
            pK2 = pK2 - pSWS2tot; % convert from tot to SWS scale

            pK2_T = ( -a ./ (TK.^2) + c ./ TK ) - gpSWS2tot(1);
            pK2_S = ( d + 2 .* g .* S ) - gpSWS2tot(2);
            pK2_P = 0;
            K2 = q(pK2);
            eK2 = 0.025 * K2; % ~2.5 % uncertainty in K, pg 854
            my_abs = @(x) sqrt(x*x);
            epK2 = my_abs( p(K2 + eK2) - pK2 ); 
            
        elseif opt.K1K2 == 17
        % Schockman and Byrne, 2021 ---------------------------------------
            % Geochimica et Cosmochimica Acta, in press
            % This study uses spectrophotometric pH measurements to determine
            % K1*K2 with unprecedented precision, and presents a new
            % parameterization for K2 based on these determinations
            % K1 is taken from Waters, Millero, and Woosley, 2014, on the total pH scale:
            a = 116.8067; b = -3655.02; c = -16.45817; 
            d = 0.04523; g = -0.615; h = -0.0002799; l = 4.969;

            pK2 = a + b ./ TK + c .* log(TK) + d .* S + g .* sqrt(S) + ...
                h .* (S.^2) + l .* (S ./ TK) ;
            pK2 = pK2 - pSWS2tot ; % convert from pHtot to SWS

            pK2_T = ( -b ./ (TK.^2) + c ./ TK - l .* S ./ (TK.^2) ) ...
                - gpSWS2tot(1) ;
            pK2_S = ( d + 0.5 .* g ./ sqrt(S) + 2 .* h .* S + l ./ TK ) ...
                - gpSWS2tot(2);
            pK2_P = 0;
            epK2 = 0.010; % in abstract

        end
        
        gpK2 = [pK2_T, pK2_S, pK2_P]; 
    end
    
    function [pKnh4, gpKnh4,epKnh4] = calc_pKnh4(opt,T,S,Pbar,Pbar_P,pSWS2tot,gpSWS2tot)
    % calcaulate pKnh4
    TK = T + 273.15; % convert to Kelvin
        if (opt.K1K2 == 6 || opt.K1K2 == 7 || opt.K1K2 == 8)
            % Knh4 = 0;
            % how to put in log space?
            % invalid K1K2 options for now

        else
        % Ammonia, added by Sharp et al 2021, from Clegg and Whitfield (1995)
        a = 9.244605; b = -2729.33; c = 1/298.15; d = 0.04203362; g = -11.24742; 

        f  = @(T,a) a(1)  + a(2) * sqrt(T) + a(3) * T + a(4) / T;
        df = @(T,a)    0.5 *a(2) / sqrt(T) + a(3)     - a(4) / T^2;

        a1 = [ -13.6416     ; 1.176949     ; -0.02860785  ; 545.4834    ];
        a2 = [ -0.1462507   ; 0.0090226468 ; -0.0001471361; 10.5425     ];
        a3 = [ 0.004669309  ; -0.0001691742;   0.0        ; -0.5677934  ];
        a4 = [ -2.354039e-05;   0.0        ;   0.0        ; 0.009698623 ];

        pKnh4 = a + b * (c - (1 / TK)) + (d + (g / TK)) * S^(0.25) + ...
                f(TK,a1) * S^0.5 + ...
                f(TK,a2) * S^1.5 + ...
                f(TK,a3) * S^2   + ...
                f(TK,a4) * S^2.5 + ...
                + p(1 - 0.001005 * S) ...
                - pSWS2tot; % convert from total scale to SWS pH scale 
        
        pKnh4_T = b / (TK^2) - g / (TK^2) * S^(0.25) + ...
                  df(TK,a1) * S^0.5 + ...
                  df(TK,a2) * S^1.5 + ...
                  df(TK,a3) * S^2   + ...
                  df(TK,a4) * S^2.5 ...
                  - gpSWS2tot(1); % convert from total scale to SWS pH scale

        pKnh4_S = 0.25 * (d + g / TK) * S^(-0.75) + ...
                  f(TK,a1) * 0.5 * S^-0.5 + ...
                  f(TK,a2) * 1.5 * S^0.5 + ...
                  f(TK,a3) * 2.0 * S   + ...
                  f(TK,a4) * 2.5 * S^1.5 - ...
                  0.001005 * dpdx(1 - 0.001005 * S) ...
                  - gpSWS2tot(2); % convert from total scale to SWS pH scale 
        pKnh4_P = 0;
        end
        gpKnh4 = [pKnh4_T, pKnh4_S, pKnh4_P];

        epKnh4 = 0.00017; % pg 2416 of Clegg and Whitefield (1995)
    end
    
    function [pKh2s,gpKh2s,epKh2s] = calc_pKh2s(opt,T,S,Pbar,Pbar_P,pSWS2tot,gpSWS2tot)
    % calculate pKh2s
    TK = T + 273.15; % convert to Kelvin
        if (opt.K1K2 == 6 || opt.K1K2 == 7 || opt.K1K2 == 8)
            % Kh2s = 0;
            % how to convert to p space?

        else
        % Millero et al (1988)
            a = 225.838; b = -13275.3; c = -34.6435;
            d = 0.3449; h = -0.0274;
        
            pKh2s = (- (a + b./TK + c.*log(TK) + d.*sqrt(S) + h.*S) ./ LOG10 ) ...
                   - pSWS2tot; % convert from total scale to SWS pH scale
            pKh2s_T = (- (-b./(TK.^2) + c./TK ) ./ LOG10 ) ...
                   - gpSWS2tot(1); % convert from total scale to SWS pH scale
            pKh2s_S = (- (0.5.*d./sqrt(S) + h) ./ LOG10 ) ...
                      - gpSWS2tot(2); % convert from total scale to SWS pH scale
            pKh2s_P = 0;
        end
        gpKh2s = [pKh2s_T, pKh2s_S, pKh2s_P];

        epKh2s = 0.033; % from Millero et al (1988), in abstract
    end
    
    function [pKar, gpKar,epKar] = calc_pKar(opt,T,S,Pbar,Pbar_P,pfH,gpfH)
    % Aragonite solubility
        TK = T + 273.15;
        Rgas = 83.14462618; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
        RT   = Rgas * TK;
        RT_T = Rgas;

        if opt.K1K2 == 6 || opt.K1K2 == 7
        % calculate Kar-Aragonite for GEOSECS -----------------------------
            % Berner, R. A., American Journal of Science 276:713-730, 1976:
            % (quoted in Takahashi et al, GEOSECS Pacific Expedition v. 3, 1982)
            a = 0.0000001; b = -34.452;      c = -39.866; 
            d = 110.21;    g = 0.0000075752;
            
            % Berner (p.722) states that he uses 1.48.
            % It appears that 1.45 was used in the GEOSECS calculations.
            Kar = 1.45 .* a .* ( b + c .* S^(1/3) + d .* log10(S) ...
                    + g .* (TK.^2) ) ;
            Kar_T = 1.45 .* a .* 2 .* g .* TK ; 
            Kar_S = 1.45 .* a .* ( (1/3) .* c .* S ^ (-2/3) + ...
                    d ./ (S .* LOG10) ) ;

            pKar = p(Kar); % - pfH ; % 'p' it and convert to SWS pH scale            
            pKar_T = dpdx(Kar) .* Kar_T ; %- gpfH(1) ; 
            pKar_S = dpdx(Kar) .* Kar_S ; %- gpfH(2) ;

            % pressure correction
            pKar = pKar - ((33.3 - 0.22 .* T) .* Pbar ./ RT ) ./ LOG10 ; % T = tempC
            pKar_T = pKar_T - (Pbar .* Rgas .* ...
                    (59.8730 .* T - 9155.988) ./ (RT)^2 ) ./ LOG10;
            pKar_P = - ((33.3 - 0.22 .* T) * Pbar_P ./ RT ) ./ LOG10;

        else
        % calculate Kar-Aragonite (Mucci, 1983) ---------------------------
            a1 = [ -171.945; -0.077993; 2903.293; 71.595] ;
            a2 = [-0.068393; 0.0017276;   88.135;   0.0 ] ;
            b = -0.10018; c = 0.0059415;

            f  = @(T,a) a(1) + a(2) .* T + a(3) ./ T + a(4) .* log10(T) ;
            df = @(T,a)      a(2) - a(3)./ (T.^2) + a(4) ./ (T .* LOG10); 

            log10Kar = f(TK,a1) + f(TK,a2) .* sqrt(S) + ...
                        b .* S + c .* S ^(3/2) ;
            pKar = -log10Kar ; % pK = -log10(K);

            pKar_T = - ( df(TK,a1) + df(TK,a2) .* sqrt(S) ) ;
            pKar_S = - ( 0.5 .* f(TK,a2) ./ sqrt(S) + ...
                        b + (3/2) .* c .* sqrt(S) ) ;

            % pressure correction
            d1 = [(-48.76 + 2.8); 0.5304; 0.0 ] ;
            d2 = [-11.76; 0.3692] ;

            pKar = pKar + ppfac(T,Pbar,d1,d2);
            pKar_T = pKar_T + ppfac_T(T,Pbar,d1,d2);
            pKar_P = ppfac_P(T,Pbar,Pbar_P,d1,d2);

            % pKa std = 0.03 (Mucci 1983)
            % std for Sal part = 0.009 (as given in CO2SYSv3)
        end
        gpKar = [pKar_T, pKar_S, pKar_P];

        epKar = 0.009; % from Mucci table 7
    end

    function [pKca, gpKca,epKca] = calc_pKca(opt,T,S,Pbar,Pbar_P,pfH,gpfH)
    % Calcite solubility
        TK = T + 273.15;
        Rgas = 83.14462618; % RgasConstant, ml bar-1 K-1 mol-1, DOEv2
        RT   = Rgas * TK;
        RT_T = Rgas;

        if opt.K1K2 == 6 || opt.K1K2 == 7
        % calculate Kca-Calcite (Berner, 1976) ----------------------------
            a = 0.0000001; b = -34.452;      c = -39.866; 
            d = 110.21;    g = 0.0000075752;

            Kca = a .* ( b + c .* S^(1/3) + d .* log10(S) ...
                    + g .* (TK.^2) ) ;
            Kca_T = a .* 2 .* g .* TK ; 
            Kca_S = a .* ( (1/3) .* c .* S ^ (-2/3) + ...
                    d ./ (S .* LOG10) ) ;

            pKca = p(Kca) ;            
            pKca_T = dpdx(Kca) .* Kca_T ; 
            pKca_S = dpdx(Kca) .* Kca_S ; 

            % pressure correction
            pKca = pKca - ((36 - 0.2 .* T) .* Pbar ./ RT ) ./ LOG10 ; % T = tempC
            pKca_T = pKca_T - (Pbar .* Rgas .* ...
                    (54.43 .* T - 9888.03) ./ (RT)^2 ) ./ LOG10;
            pKca_P = - ((36 - 0.2 .* T) *Pbar_P ./ RT ) ./ LOG10;

        else
        % calculate Kca-Calcite (Mucci, 1983) -----------------------------
            a1 = [-171.9065; -0.077993; 2839.319; 71.595] ;
            a2 = [ -0.77712; 0.0028426;   178.34;   0.0 ] ;
            b = -0.07711; c = 0.0041249;

            f  = @(T,a) a(1) + a(2) .* T + a(3) ./ T + a(4) .* log10(T) ;
            df = @(T,a)      a(2) - a(3)./ (T.^2) + a(4) ./ (T .* LOG10); 

            log10Kca = f(TK,a1) + f(TK,a2) .* sqrt(S) + ...
                        b .* S + c .* S ^(3/2) ;
            pKca = -log10Kca ; % pK = -log10(K);

            pKca_T = - ( df(TK,a1) + df(TK,a2) .* sqrt(S) );
            pKca_S = - ( 0.5 .* f(TK,a2) ./ sqrt(S) + ...
                        b + (3/2) .* c .* sqrt(S) ) ;

            % pressure correction
            d1 = [-48.76; 0.5304; 0.0 ] ;
            d2 = [-11.76; 0.3692] ;

            pKca = pKca + ppfac(T,Pbar,d1,d2);
            pKca_T = pKca_T + ppfac_T(T,Pbar,d1,d2);
            pKca_P = ppfac_P(T,Pbar,Pbar_P,d1,d2);
            
            % pKca std = 0.03 (Mucci 1983)
            % std for Sal part = 0.01 (as given in CO2SYSv3)
        end
        gpKca = [pKca_T, pKca_S, pKca_P];

        epKca = 0.010; % from Mucci 1983, table 7
    end
    
    function [pSWS2tot,gpSWS2tot,pFREE2tot,gpFREE2tot] = calc_pSWS2tot(opt,S,pKs,gpKs,pKf,gpKf)
    % pH scale conversion factors (not pressure corrected)----------- 

    % calculate TF (Riley 1965)--------------------------------------------
        TF = (0.000067./18.998).*(S./1.80655); % mol/kg-SW
        TF_S = (0.000067/18.998)*(1/1.80655); % derivative wrt S
        pTF = p(TF);
        pTF_S = dpdx(TF).*TF_S;
        
    % calculate TS (Morris & Riley 1966)-----------------------------------
        TS = (0.14./96.062).*(S./1.80655); % mol/kg-SW
        TS_S = (0.14/96.062)*(1/1.80655); % derivative wrt S
        pTS = p(TS);
        pTS_S = dpdx(TS).*TS_S;
        
        top = 1 + q(pTS - pKs) ;
        top_T = dqdx(pTS - pKs) * (-gpKs(1)) ;
        top_S = dqdx(pTS - pKs) * (pTS_S - gpKs(2)) ;
        top_P = dqdx(pTS - pKs) * (-gpKs(3)) ;
        bot = top + q(pTF - pKf) ;
        bot_T = top_T + dqdx(pTF - pKf) * (-gpKf(1)) ;
        bot_S = top_S + dqdx(pTF - pKf) * (pTF_S - gpKf(2));
        bot_P = top_P + dqdx(pTF - pKf) * ( -gpKf(3) );
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

    function [pfH,gpfH,epfH] = calc_pfH(opt,T,S)
        % fH = [H]/(1 + TS/Ks)
        TK = T + 273.15; % convert to Kelvin
        if opt.K1K2 == 8
            %fH = 1; % shouldn't occur
        elseif opt.K1K2 == 7
        % fH def'n: Peng et al, Tellus 39B: 439-458, 1987: ----------------
            % They reference the GEOSECS report, but round the value
            % given there off so that it is about .008 (1%) lower. It
            % doesn't agree with the check value they give on p. 456.
            fH   = ( 1.29 - 0.00204.*(TK) + ...
                (0.00046 - (0.00000148 * (TK))) * (S).^2 ) ;
            pfH = p(fH);
            gpfH_T = dpdx(pfH) - 0.00204 - 0.00000148 * (S.^2)  ;
            gpfH_S = dpdx(pfH) * 2 * (0.00046 - (0.00000148 * (TK))) * S ;
        else
        % fH def'n: Takahashi et al, Ch 3 in GEOSECS v.3, 1982 ------------
            fH = 1.2948 - 0.002036*TK + ...
                ( 0.0004607 - 0.000001475 * (TK) ) .* (S).^2 ;
            fH_T = -0.002036 + ( -0.000001475 ) .* (S).^2  ;
            fH_S = 2 * ( 0.0004607 - ( 0.000001475 .* (TK) ) ) .* (S) ;
            pfH = p(fH);
            gpfH_T = dpdx(fH) * fH_T;
            gpfH_S = dpdx(fH) * fH_S;
        end
        gpfH = [gpfH_T, gpfH_S, 0];
        % assumed independent of pressure
        
        efH = 0.005; % ± 0.005 on fH from Culberson, Pytkowicz, 
                     % and Hawley 1970 Journal of Marine Research
                     % assume 1 sigma -MF
        my_abs = @(x) sqrt(x*x);
        epfH = my_abs( p(fH + efH) - pfH ) ; 
    end 
    
end




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







