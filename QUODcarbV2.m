function [est,obs,iflag] = QUODcarbV2(obs,sys)
% [est,obs,iflag] = QUODcarb(obs,sys);
%
% OUTPUT:
%   est := posterior estimates of co2-system variables, equilibrium constants, and precisions
%   obs := same as input except that the pK's have been added
% iflag := 0 if solver converged to specified accuracy 
%          1 after reaching maximum number of iterations without convergging
% INPUT:
%   obs  := co2-system measured quantities with precisions 
%   sys  := struct with stuff for co2-system solver (initialized using mksys.m)

    isgood = @(thing) (~isempty(thing) & ~sum(isnan(thing)));
    p = sys.p;
    q = sys.q;
    w = @(x,e) p(1+e./x).^(-2); % convert x+/-e into precision for p(x)
    
    nv = size(sys.K,2);
    yobs = nan(nv,1);
    wobs = nan(nv,1);
    nTP = length(obs.m);

    % make sure all the required fields in the obs struct exist
    if (~isfield(obs, 'sal'))
        error('Need to provide salinity measurement.');
    else
        yobs(sys.isal) = obs.sal;
    end
    if (~isfield(obs, 'esal'));
        obs.esal = w(obs.sal,0.002); % std = 0.002 PSU
        fprintf('Warning: Assuming salinity uncertainty is 0.002 PSU');
    else
        wobs(sys.isal) = w(obs.sal,obs.esal); % std e -> w
    end
    if (~isfield(obs,'TC'))
        obs.TC = [];
        yobs(sys.iTC) = nan;
    else
        yobs(sys.iTC) = p((obs.TC)*1e-6); % convt to mol/kg
    end
    if (~isfield(obs,'eTC'))
        obs.eTC = [];
        wobs(sys.iTC) = nan;
    else
        wobs(sys.iTC) = w(obs.TC,obs.eTC); % std e -> w
    end 
    if(~isfield(obs,'TA'))
        obs.TA = [];
        yobs(sys.iTA) = nan;
    else
        yobs(sys.iTA) = p((obs.TA)*1e-6); % convt to mol/kg
    end
    if (~isfield(obs,'eTA'))
        obs.eTA = [];
        wobs(sys.iTA) = nan;
    else
        wobs(sys.iTA) = w(obs.TA,obs.eTA); % std e -> w
    end
    if (~isfield(obs,'m'))
        error('Need to provide temperature and pressure measurement.')
    end
    % create obs structure with fieldnames
    if (~isfield(obs.m(1),'epK0'))
        obs.m(1).epK0 = [];
    end
    if (~isfield(obs.m(1),'pK0'))
        obs.m(1).pK0 = [];
    end
    if (~isfield(obs.m(1),'epK1'))
        obs.m(1).epK1 = [];
    end
    if (~isfield(obs.m(1),'pK1'))
        obs.m(1).pK1 = [];
    end
    if (~isfield(obs.m(1),'epK2'))
        obs.m(1).epK2 = [];
    end
    if (~isfield(obs.m(1),'pK2'))
        obs.m(1).pK2 = [];
    end
    if (~isfield(obs.m(1),'epp2f'))
        obs.m(1).epp2f = [];
    end
    if (~isfield(obs.m(1),'pp2f'))
        obs.m(1).pp2f = [];
    end
    if (~isfield(obs.m(1),'efco2'))
        obs.m(1).efco2 = [];
    end
    if (~isfield(obs.m(1),'fco2'))
        obs.m(1).fco2 = [];
    end
    if (~isfield(obs.m(1),'eco2st'))
        obs.m(1).eco2st = [];
    end
    if (~isfield(obs.m(1),'co2st'))
        obs.m(1).co2st = [];
    end
    if (~isfield(obs.m(1),'ehco3'))
        obs.m(1).ehco3 = [];
    end
    if (~isfield(obs.m(1),'hco3'))
        obs.m(1).hco3 = [];
    end
    if (~isfield(obs.m(1),'eco3'))
        obs.m(1).eco3 = [];
    end
    if (~isfield(obs.m(1),'co3'))
        obs.m(1).co3 = [];
    end
    if (~isfield(obs.m(1),'eph'))
        obs.m(1).ph = [];
    end
    if (~isfield(obs.m(1),'ph'))
        obs.m(1).eph = [];
    end

    if (ismember('borate',sys.abr))
        if (~isfield(obs.m(1),'epKb'))
            obs.m(1).epKb = [];
        end
        if (~isfield(obs.m(1),'pKb'))
            obs.m(1).pKb = [];
        end
        if (~isfield(obs,'TB'))
            % Culkin, F., in Chemical Oceanography,
            % ed. Riley and Skirrow, 1965: ( copied from Orr's code )
            fprintf('Warning: Assuming TB = 0.0004106 * sal / 35  mol/kg-SW.\n');
            obs.TB = (0.0004106 * obs.sal / 35) * (1e6) ; % convt to µmol/kg
            yobs(sys.iTB) = p(obs.TB);
        else
            yobs(sys.iTB) = p(obs.TB*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs, 'eTB'))
            obs.eTB = w((obs.TB),(( 35 / 0.0004106 ) * obs.esal)*(1e6) ) ;
            wobs(sys.iTB) = obs.eTB;
        else
            wobs(sys.iTB) = w(obs.TB,obs.eTB);
        end
        if (~isfield(obs.m(1),'eboh4'))
            obs.m(1).eboh4 = [];
        end
        if (~isfield(obs.m(1),'boh4'))
            obs.m(1).boh4 = [];
        end
        if (~isfield(obs.m(1),'eboh3'))
            obs.m(1).eboh3 = [];
        end
        if (~isfield(obs.m(1),'boh3'))
            obs.m(1).boh3 = [];
        end
    end

    if (ismember('water',sys.abr))
        if (~isfield(obs.m(1), 'epKw'))
            obs.m(1).epKw = [];
        end
        if (~isfield(obs.m(1), 'pKw'))
            obs.m(1).pKw = [];
        end
        if (~isfield(obs.m(1),'eoh'))
            obs.m(1).eoh = [];
        end
        if (~isfield(obs.m(1),'oh'))
            obs.m(1).oh = [];
        end
    end

    if (ismember('sulfate',sys.abr))
        if (~isfield(obs.m(1), 'epKs'))
            obs.m(1).epKs = [];
        end
        if (~isfield(obs.m(1), 'pKs'))
            obs.m(1).pKs = [];
        end
        if (~isfield(obs, 'TS'))
            % Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
            % copied from Orr's code
            fprintf('Warning: Assuming TS = ( 0.14 / 96.062 ) * ( sal / 1.80655 ) mol/kg-SW.\n'); 
            obs.TS = ( 0.14 / 96.062 ) * ( obs.sal / 1.80655 )*(1e6); % convt to µmol/kg
            yobs(sys.iTS) = p(obs.TS);
        else
            yobs(sys.iTS) = p(obs.TS*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs, 'eTS'))
            obs.eTS = w(obs.TS,( ( 0.14 / 96.062 ) / 1.80655 ) * obs.esal*1e6);
            wobs(sys.iTS) = obs.eTS;
        else
            wobs(sys.iTS) = w(obs.TS,obs.eTS);
        end
        if (~isfield(obs.m(1), 'phf'))
            obs.m(1).phf = [];
        end
        if (~isfield(obs.m(1), 'ephf'))
            obs.m(1).ephf = [];
        end
        if (~isfield(obs.m(1), 'so4'))
            obs.m(1).so4 = [];
        end
        if (~isfield(obs.m(1), 'eso4'))
            obs.m(1).eso4 = [];
        end
        if (~isfield(obs.m(1), 'hso4'))
            obs.m(1).hso4 = [];
        end     
        if (~isfield(obs.m(1), 'ehso4'))
            obs.m(1).ehso4 = [];
        end     
    end

    if (ismember('fluoride',sys.abr))
        if (~isfield(obs.m(1), 'epKf'))
            obs.m(1).epKf = [];
        end
        if (~isfield(obs.m(1), 'pKf'))
            obs.m(1).pKf = [];
        end
        if (~isfield(obs, 'TF'))
            % Riley, J. P., Deep-Sea Research 12:219-220, 1965:
            % this is .000068.*Sali./35. = .00000195.*Sali
            fprintf('Warning: Assuming TF = ( 0.000067 / 18.998 ) * ( sal / 1.80655 ). mol/kg-SW. \n');
            obs.TF = ( 0.000067 / 18.998 ) * ( obs.sal / 1.80655 )*1e6; % convt to µmol/kg-SW
            yobs(sys.iTF) = p(obs.TF);
        else
            yobs(sys.iTF) = p(obs.TF*1e-6); % convt µmol/kg to mol/kg
        end
        if (~isfield(obs, 'eTF'))
            obs.eTF = w(obs.TF,( ( 0.000067 / 18.998 ) / 1.80655 ) * obs.esal*1e6);
            wobs(sys.iTF) = obs.eTF;
        else
            wobs(sys.iTF) = w(obs.TF,obs.eTF);
        end
        if (~isfield(obs.m(1), 'F'))
            obs.m(1).F = [];
        end
        if (~isfield(obs.m(1), 'eF'))
            obs.m(1).eF = [];
        end
        if (~isfield(obs.m(1), 'HF'))
            obs.m(1).HF = [];
        end
        if (~isfield(obs.m(1), 'eHF'))
            obs.m(1).eHF = [];
        end
    end

    if (ismember('phosphate',sys.abr))
        if (~isfield(obs.m(1), 'epK1p'))
            obs.m(1).epK1p = [];
        end
        if (~isfield(obs.m(1), 'pK1p'))
            obs.m(1).pK1p = [];
        end
        if (~isfield(obs.m(1), 'epK2p'))
            obs.m(1).epK2p = [];
        end
        if (~isfield(obs.m(1), 'pK2p'))
            obs.m(1).pK2p = [];
        end
        if (~isfield(obs.m(1), 'epK3p'))
            obs.m(1).epK3p = [];
        end
        if (~isfield(obs.m(1), 'pK3p'))
            obs.m(1).pK3p = [];
        end
        if (~isfield(obs, 'TP'))
            fprintf('Warning: Assuming TP = 1e-3 µmol/kg-SW.\n');
            obs.TP = 1e-3; % µmol/kg
            yobs(sys.iTP) = p(obs.TP*1e-6); % convt µmol/kg to mol/kg
        else
            yobs(sys.iTP) = p(obs.TP*1e-6); % convt µmol/kg to mol/kg 
        end
        if (~isfield(obs, 'eTP'))
            fprintf('Warning: Setting eTP = 1e-3 µmol/kg-SW.\n')
            obs.eTP = 1e-3; % µmol/kg
            wobs(sys.iTP) = w(obs.TP,obs.eTP);
        else
            wobs(sys.iTP) = w(obs.TP,obs.eTP);
        end
        if (~isfield(obs.m(1), 'h3po4'))
            obs.m(1).h3po4 = [];
        end
        if (~isfield(obs.m(1), 'eh3po4'))
            obs.m(1).eh3po4 = [];
        end
        if (~isfield(obs.m(1), 'h2po4'))
            obs.m(1).h2po4 = [];
        end
        if (~isfield(obs.m(1), 'eh2po4'))
            obs.m(1).eh2po4 = [];
        end
        if (~isfield(obs.m(1), 'hpo4'))
            obs.m(1).hpo4 = [];
        end
        if (~isfield(obs.m(1), 'ehpo4'))
            obs.m(1).ehpo4 = [];
        end
        if (~isfield(obs.m(1), 'po4'))
            obs.m(1).po4 = [];
        end
        if (~isfield(obs.m(1), 'epo4'))
            obs.m(1).epo4 = [];
        end
    end

    if (ismember('silicate',sys.abr))
        if (~isfield(obs.m(1), 'epKsi'))
            obs.m(1).epKsi = [];
        end
        if (~isfield(obs.m(1), 'pKsi'))
            obs.m(1).pKsi = [];
        end
        if (~isfield(obs, 'TSi'))
            fprintf('Warning: Assuming TSi = 1e-3  µmol/kg-SW.\n');
            obs.TSi = 1e-3; % µmol/kg
            yobs(sys.iTSi) = p(obs.TSi*1e-6); % convt µmol/kg to mol/kg
        else
            yobs(sys.iTSi) = p(obs.TSi*1e-6);
        end
        if (~isfield(obs, 'eTSi'))
            fprintf('Warning: Assuming eTSi = 1e-3  µmol/kg-SW.\n');
            obs.eTSi = 1e-3; % µmol/kg
            wobs(sys.iTSi) = w(obs.TSi,obs.eTSi);
        else
            wobs(sys.iTSi) = w(obs.TSi,obs.eTSi);
        end
        if (~isfield(obs.m(1), 'sioh4'))
            obs.m(1).sioh4 = [];
        end
        if (~isfield(obs.m(1), 'esioh4'))
            obs.m(1).esioh4 = [];
        end
        if (~isfield(obs.m(1), 'siooh3'))
            obs.m(1).siooh3 = [];
        end
        if (~isfield(obs.m(1), 'esiooh3'))
            obs.m(1).esiooh3 = [];
        end
    end

    if (ismember('ammonia',sys.abr))
        if (~isfield(obs.m(1), 'epKnh4'))
            obs.m(1).epKnh4 = [];
        end
        if (~isfield(obs.m(1), 'pKnh4'))
            obs.m(1).pKnh4 = [];
        end
        if (~isfield(obs, 'TNH3'))
            fprintf('Warning: Assuming TNH3 = 1e-3  µmol/kg-SW.\n');
            obs.TNH3 = 1e-3; % µmol/kg
            yobs(sys.iTNH3) = p(obs.TNH3*1e-6); % convt µmol/kg to mol/kg
        else
            yobs(sys.iTNH3) = p(obs.TNH3*1e-6);
        end
        if (~isfield(obs, 'epTNH3'))
            fprintf('Warning: Assuming eTNH3 = 1e-3  µmol/kg-SW.\n');
            obs.epTNH3 = 1e-3; % µmol/kg
            wobs(sys.iTNH3) = w(obs.TNH3,obs.epTNH3);
        else
            wobs(sys.iTNH3) = w(obs.TNH3,obs.epTNH3);
        end
        if (~isfield(obs.m(1), 'nh3'))
            obs.m(1).nh3 = [];
        end
        if (~isfield(obs.m(1), 'enh3'))
            obs.m(1).enh3 = [];
        end
        if (~isfield(obs.m(1), 'nh4'))
            obs.m(1).nh4 = [];
        end
        if (~isfield(obs.m(1), 'enh4'))
            obs.m(1).enh4 = [];
        end
    end

    if (ismember('sulfide',sys.abr))
        if (~isfield(obs.m(1), 'epKh2s'))
            obs.m(1).epKh2s = [];
        end
        if (~isfield(obs.m(1), 'pKh2s'))
            obs.m(1).pKh2s = [];
        end
        if (~isfield(obs, 'TH2S'))
            fprintf('Warning: Assuming TH2S = 1e-3  µmol/kg-SW.\n');
            obs.TH2S = 1e-3; % µmol/kg
            yobs(sys.iTH2S) = p(obs.TH2S*1e-6); % convt µmol/kg to mol/kg
        else
            yobs(sys.iTH2S) = p(obs.TH2S*1e-6);
        end
        if (~isfield(obs, 'eTH2S'))
            fprintf('Warning: Assuming eTH2S = 1e-3  µmol/kg-SW.\n');
            obs.eTH2S = 1e-3; % µmol/kg
            wobs(sys.iTH2S) = w(obs.TH2S,obs.eTH2S);
        else
            wobs(sys.iTH2S) = w(obs.TH2S,obs.eTH2S);
        end
        if (~isfield(obs.m(1), 'hs'))
            obs.m(1).hs = [];
        end
        if (~isfield(obs.m(1), 'ehs'))
            obs.m(1).ehs = [];
        end
        if (~isfield(obs.m(1), 'h2s'))
            obs.m(1).h2s = [];
        end
        if (~isfield(obs.m(1), 'eh2s'))
            obs.m(1).eh2s = [];
        end
    end
        
    nTP = length(obs.m);
    yobs(sys.iTC) = obs.TC;
    yobs(sys.iTA) = obs.TA;
    
    for i = 1:nTP % loop over all the pressure and temperature sensitive components
        yobs(sys.m(i).iT) = obs.m(i).T;
        yobs(sys.m(i).iP) = obs.m(i).P;
        
        wobs(sys.m(i).iT) = obs.m(i).eT;
        wobs(sys.m(i).iP) = obs.m(i).eP;
        
        [pK,gpK] = local_pK( obs.m(i).T, obs.sal, obs.m(i).P );
        pK0   = pK(1);  pK1  = pK(2);  pK2   = pK(3);  pKb   = pK(4);  
        pKw   = pK(5);  pKs  = pK(6);  pKf   = pK(7);  pK1p  = pK(8);  
        pK2p  = pK(9);  pK3p = pK(10); pKsi  = pK(11); pKnh4 = pK(12); 
        pKh2s = pK(13); pp2f = pK(14);        
        
        %
        % add "observations" for the equilibrium constants 
        % and transfer from obs struct to yobs and wobs vectors
        %
        if (ismember('carbonate',sys.abr))
            if (isgood(obs.m(i).epK0))   
                wobs(sys.m(i).iK0) = obs.m(i).epK0;
            else
                wobs(sys.m(i).iK0) = (0.002)^(-2);
                obs.m(i).epK0 = wobs(sys.m(i).iK0);
            end
            if (isgood(obs.m(i).pK0))
                yobs(sys.m(i).iK0) = obs.m(i).pK0;
            else
                yobs(sys.m(i).iK0) = pK0;
                obs.m(i).pK0 = pK0;
            end

            if (isgood(obs.m(i).epK1))
                wobs(sys.m(i).iK1) = obs.m(i).epK1;
            else
                wobs(sys.m(i).iK1) = (0.01).^(-2);  % wK1 = 1/(1 + (0.01/pKsys(2)))^2 ;
                obs.m(i).epK1 = wobs(sys.m(i).iK1);
            end
            if (isgood(obs.m(i).pK1))
                yobs(sys.m(i).iK1) = obs.m(i).pK1;
            else
                yobs(sys.m(i).iK1) = pK1;
                obs.m(i).pK1 = pK1;
            end            
            
            if (isgood(obs.m(i).epK2))
                wobs(sys.m(i).iK2) = obs.m(i).epK2;
            else
                wobs(sys.m(i).iK2) = (0.02).^(-2);  % wK2 = 1/(1 + (0.02/pKsys(3)))^2 ;
                obs.m(i).epK2 = wobs(sys.m(i).iK2);
            end
            
            if (isgood(obs.m(i).pK2))
                yobs(sys.m(i).iK2) = obs.m(i).pK2;
            else
                yobs(sys.m(i).iK2) = pK2;
                obs.m(i).pK2 = pK2;
            end
            
            if (isgood(obs.m(i).epp2f))
                wobs(sys.m(i).ip2f) = obs.m(i).epp2f;
            else
                wobs(sys.m(i).ip2f) = (0.001).^(-2);
                obs.m(i).epp2f = wobs(sys.m(i).ip2f);
            end
            
            if (isgood(obs.m(i).pp2f))
                yobs(sys.m(i).ip2f) = obs.m(i).pp2f;
            else
                yobs(sys.m(i).ip2f) = pp2f;
                obs.m(i).pp2f = pp2f;
            end
            
            if (isgood(obs.m(i).co2st))
                yobs(sys.m(i).ico2st) = p(obs.m(i).co2st*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).ico2st) = nan;
                obs.m(i).co2st = nan;
            end
            if (isgood(obs.m(i).eco2st))
                wobs(sys.m(i).ico2st) = w(obs.m(i).co2st,obs.m(i).eco2st);
            else
                wobs(sys.m(i).ico2st) = nan;
                obs.m(i).eco2st = nan;
            end
            if (isgood(obs.m(i).hco3))
                yobs(sys.m(i).ihco3) = p(obs.m(i).hco3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).ihco3) = nan;
                obs.m(i).hco3 = nan;
            end
            if (isgood(obs.m(i).ehco3))
                wobs(sys.m(i).ihco3) = w(obs.m(i).hco3,obs.m(i).ehco3);
            else
                wobs(sys.m(i).ihco3) = nan;
                obs.m(i).ehco3 = nan;
            end
            if (isgood(obs.m(i).co3))
                yobs(sys.m(i).ico3) = p(obs.m(i).co3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).ico3) = nan;
                obs.m(i).co3 = nan;
            end
            if (isgood(obs.m(i).eco3))
                wobs(sys.m(i).ico3) = w(obs.m(i).co3,obs.m(i).eco3);
            else
                wobs(sys.m(i).ico3) = nan;
                obs.m(i).eco3 = nan;
            end
            if (isgood(obs.m(i).fco2))
                yobs(sys.m(i).ifco3) = p(obs.m(i).fco2*1e-6); % convt µatm to atm
            else
                yobs(sys.m(i).ifco2) = nan;
                obs.m(i).fco2 = nan;
            end
            if (isgood(obs.m(i).efco2))
                wobs(sys.m(i).ifco2) = w(obs.m(i).fco2,obs.m(i).efco2);
            else
                wobs(sys.m(i).ifco2) = nan;
                obs.m(i).efco2 = nan;
            end
            if (isgood(obs.m(i).ph))
                yobs(sys.m(i).iph) = obs.m(i).ph;
            else
                yobs(sys.m(i).iph) = nan;
                obs.m(i).ph = nan;
            end
            if (isgood(obs.m(i).eph))
                wobs(sys.m(i).iph) = (obs.m(i).eph)^(-2);
            else
                wobs(sys.m(i).iph) = nan;
                obs.m(i).eph = nan;
            end        
        end
        
        if (ismember('borate',sys.abr))
            if (isgood(obs.m(i).epKb))
                wobs(sys.m(i).iKb) = obs.m(i).epKb;
            else
                wobs(sys.m(i).iKb) = (0.01).^(-2);  % wKb = 1/(1 + (0.01/pKsys(4)))^2 ;
                obs.m(i).epKb = wobs(sys.m(i).iKb);
            end
            if (isgood(obs.m(i).pKb))
                yobs(sys.m(i).iKb) = obs.m(i).pKb;
            else
                yobs(sys.m(i).iKb) = pKb;
                obs.m(i).pKb = pKb;
            end
            if (isgood(obs.m(i).boh3))
                yobs(sys.m(i).iboh3) = p(obs.m(i).boh3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).iboh3) = nan;
                obs.m(i).boh3 = nan;
            end
            if (isgood(obs.m(i).eboh3))
                wobs(sys.m(i).iboh3) = w(obs.m(i).boh3,obs.m(i).eboh3);
            else
                wobs(sys.m(i).iboh3) = nan;
                obs.m(i).eboh3 = nan;
            end
            if (isgood(obs.m(i).boh4))
                yobs(sys.m(i).iboh4) = p(obs.m(i).boh4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).iboh4) = nan;
                obs.m(i).boh4 = nan;
            end
            if (isgood(obs.m(i).eboh4))
                wobs(sys.m(i).iboh4) = w(obs.m(i).boh4,obs.m(i).eboh4);
            else
                wobs(sys.m(i).iboh4) = nan;
                obs.m(i).eboh4 = nan;
            end
        end
        
        if (ismember('water',sys.abr))
            if (isgood(obs.m(i).epKw))
                wobs(sys.m(i).iKw) = obs.m(i).epKw;
            else
                wobs(sys.m(i).iKw) = (0.01).^(-2);  % wKw = 1/(1 + (0.01/pKsys(5)))^2 ;
                obs.m(i).epKw = wobs(sys.m(i).iKw);
            end
            if (isgood(obs.m(i).pKw))
                wobs(sys.m(i).iKw) = obs.m(i).pKw;
            else
                yobs(sys.m(i).iKw) = pKw;
                obs.m(i).pKw = pKw;
            end
            if (isgood(obs.m(i).oh))
                yobs(sys.m(i).ioh) = p(obs.m(i).oh*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).ioh) = nan;
                obs.m(i).oh = nan;
            end
            if (isgood(obs.m(i).eoh))
                wobs(sys.m(i).ioh) = w(obs.m(i).oh,obs.m(i).eoh);
            else
                wobs(sys.m(i).ioh) = nan;
                obs.m(i).eoh = nan;
            end
            
        end
        
        if (ismember('sulfate',sys.abr))
            if (isgood(obs.m(i).epKs))
                wobs(sys.m(i).iKs) = obs.m(i).epKs;
            else
                wobs(sys.m(i).iKs) = (0.0021).^(-2); % wKs = 1/(1 + (0.0021/pKsys(6)))^2 ;
                obs.m(i).epKs = wobs(sys.m(i).iKs);
            end
            if (isgood(obs.m(i).pKs))
                yobs(sys.m(i).iKs) = obs.m(i).pKs;
            else
                yobs(sys.m(i).iKs) = pKs;
                obs.m(i).pKs = pKs;
            end
            if (isgood(obs.m(i).hso4))
                yobs(sys.m(i).ihso4) = p(obs.m(i).hso4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).ihso4) = nan;
                obs.m(i).hso4 = nan;
            end
            if (isgood(obs.m(i).ehso4))
                wobs(sys.m(i).ihso4) = w(obs.m(i).hso4,obs.m(i).ehso4);
            else
                wobs(sys.m(i).ihso4) = nan;
                obs.m(i).ehso4 = nan;
            end
            if (isgood(obs.m(i).so4))
                yobs(sys.m(i).iso4) = p(obs.m(i).so4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).iso4) = nan;
                obs.m(i).so4 = nan;
            end
            if (isgood(obs.m(i).eso4))
                wobs(sys.m(i).iso4) = w(obs.m(i).so4,obs.m(i).eso4);
            else
                wobs(sys.m(i).iso4) = nan;
                obs.m(i).eso4 = nan;
            end
            if (isgood(obs.m(i).phf))
                yobs(sys.m(i).iphf) = obs.m(i).phf; % hydrogen free
            else
                yobs(sys.m(i).iphf) = nan;
                obs.m(i).phf = nan;
            end
            if (isgood(obs.m(i).ephf))
                wobs(sys.m(i).iphf) = (obs.m(i).ephf)^(-2);
            else
                wobs(sys.m(i).iphf) = nan;
                obs.m(i).ephf = nan;
            end

        end
        
        if (ismember('fluoride',sys.abr))
            if (isgood(obs.m(i).epKf))
                wobs(sys.m(i).iKf) = obs.m(i).epKf;
            else
                wobs(sys.m(i).iKf) = (0.02).^(-2);   % wKF = 1/(p(1 + 0.02/KF))^2 ; % 2% relative error
                obs.m(i).epKf = wobs(sys.m(i).iKf);
            end
            if (isgood(obs.m(i).pKf))
                yobs(sys.m(i).iKf) = obs.m(i).pKf;
            else
                yobs(sys.m(i).iKf) = pKf;
                obs.m(i).pKf = pKf;
            end
            if (isgood(obs.m(i).F))
                yobs(sys.m(i).iF) = p(obs.m(i).F*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).iF) = nan;
                obs.m(i).F = nan;
            end
            if (isgood(obs.m(i).eF))
                wobs(sys.m(i).iF) = w(obs.m(i).F,obs.m(i).eF);
            else
                wobs(sys.m(i).iF) = nan;
                obs.m(i).eF = nan;
            end
            if (isgood(obs.m(i).HF))
                yobs(sys.m(i).iHF) = p(obs.m(i).HF*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).iHF) = nan;
                obs.m(i).HF = nan;
            end
            if (isgood(obs.m(i).eHF))
                wobs(sys.m(i).iHF) = w(ons.m(i).HF,obs.m(i).eHF);
            else
                wobs(sys.m(i).iHF) = nan;
                obs.m(i).eHF = nan;
            end
        end
        
        if (ismember('phosphate',sys.abr))
            if (isgood(obs.m(i).epK1p))
                wobs(sys.m(i).iK1p) = obs.m(i).epK1p;
            else
                wobs(sys.m(i).iK1p) = (0.09).^(-2);  % wK1p = 1/(1 + (0.09/pKsys(8)))^2 ;
                obs.m(i).epK1p = wobs(sys.m(i).iK1p);
            end
            if (isgood(obs.m(i).pK1p))
                yobs(sys.m(i).iK1p) = obs.m(i).pK1p;
            else
                yobs(sys.m(i).iK1p) = pK1p;
                obs.m(i).pK1p = pK1p;
            end
            if (isgood(obs.m(i).epK2p))
                wobs(sys.m(i).iK2p) = obs.m(i).epK2p;
            else
                wobs(sys.m(i).iK2p) = (0.03).^(-2);  % wK2p = 1/(1 + (0.03/pKsys(9)))^2 ;
                obs.m(i).epK2p = wobs(sys.m(i).iK2p);
            end
            if (isgood(obs.m(i).pK2p))
                yobs(sys.m(i).iK2p) = obs.m(i).pK2p;
            else
                yobs(sys.m(i).iK2p) = pK2p;
                obs.m(i).pK2p = pK2p;
            end
            if (isgood(obs.m(i).epK3p))
                wobs(sys.m(i).iK3p) = obs.m(i).epK3p;
            else
                wobs(sys.m(i).iK3p) = (0.02).^(-2);  % wK3p = 1/(1 + (0.02/pKsys(10)))^2 ;
                obs.m(i).epK3p = wobs(sys.m(i).iK3p);
            end
            if (isgood(obs.m(i).pK3p))
                yobs(sys.m(i).iK3p) = obs.m(i).pK3p;
            else
                yobs(sys.m(i).iK3p) = pK3p;
                obs.m(i).pK3p = pK3p;
            end
            if (isgood(obs.m(i).h3po4))
                yobs(sys.m(i).ih3po4) = p(obs.m(i).h3po4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).ih3po4) = nan;
                obs.m(i).h3po4 = nan;
            end
            if (isgood(obs.m(i).eh3po4))
                wobs(sys.m(i).ih3po4) = w(obs.m(i).h3po4,obs.m(i).eh3po4);
            else
                wobs(sys.m(i).ih3po4) = nan;
                obs.m(i).eh3po4 = nan;
            end
            if (isgood(obs.m(i).h2po4))
                yobs(sys.m(i).ih2po4) = p(obs.m(i).h2po4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).ih2po4) = nan;
                obs.m(i).h2po4 = nan;
            end
            if (isgood(obs.m(i).eh2po4))
                wobs(sys.m(i).ih2po4) = w(obs.m(i).h2po4,obs.m(i).eh2po4);
            else
                wobs(sys.m(i).ih2po4) = nan;
                obs.m(i).eh2po4 = nan;
            end
            if (isgood(obs.m(i).hpo4))
                yobs(sys.m(i).ihpo4) = p(obs.m(i).hpo4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).ihpo4) = nan;
                obs.m(i).hpo4 = nan;
            end
            if (isgood(obs.m(i).ehpo4))
                wobs(sys.m(i).ihpo4) = w(obs.m(i).hpo4,obs.m(i).ehpo4);
            else
                wobs(sys.m(i).ihpo4) = nan;
                obs.m(i).ehpo4 = nan;
            end
            if (isgood(obs.m(i).po4))
                yobs(sys.m(i).ipo4) = p(obs.m(i).po4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).ipo4) = nan;
                obs.m(i).po4 = nan;
            end
            if (isgood(obs.m(i).epo4))
                wobs(sys.m(i).ipo4) = w(obs.m(i).po4,obs.m(i).epo4);
            else
                wobs(sys.m(i).ipo4) = nan;
                obs.m(i).epo4 = nan;
            end

        end

        if (ismember('silicate',sys.abr))
            if (isgood(obs.m(i).epKsi))
                wobs(sys.m(i).iKsi) = obs.m(i).epKsi;
            else
                wobs(sys.m(i).iKsi) = (0.02).^(-2);  % wKSi = 1/(1 + (0.02/pKsys(11)))^2 ;
                obs.m(i).epKsi = wobs(sys.m(i).iKsi);
            end
            if (isgood(obs.m(i).pKsi))
                yobs(sys.m(i).iKsi) = obs.m(i).pKsi;
            else
                yobs(sys.m(i).iKsi) = pKsi;
                obs.m(i).pKsi = pKsi;
            end
            if (isgood(obs.m(i).sioh4))
                yobs(sys.m(i).isioh4) = p(obs.m(i).sioh4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).isioh4) = nan;
                obs.m(i).sioh4 = nan;
            end
            if (isgood(obs.m(i).esioh4))
                wobs(sys.m(i).isioh4) = w(obs.m(i).sioh4,obs.m(i).esioh4);
            else
                wobs(sys.m(i).isioh4) = nan;
                obs.m(i).esioh4 = nan;
            end
            if (isgood(obs.m(i).siooh3))
                yobs(sys.m(i).isiooh3) = p(obs.m(i).siooh3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).isiooh3) = nan;
                obs.m(i).siooh3 = nan;
            end
            if (isgood(obs.m(i).esiooh3))
                wobs(sys.m(i).isiooh3) = w(obs.m(i).siooh3,obs.m(i).esiooh3);
            else
                wobs(sys.m(i).isiooh3) = nan;
                obs.m(i).esiooh3 = nan;
            end
            
        end

        if (ismember('ammonia',sys.abr))
            if (isgood(obs.m(i).epKnh4))
                wobs(sys.m(i).iKnh4) = obs.m(i).epKnh4;
            else
                wobs(sys.m(i).iKnh4) = (0.00017).^(-2);  % wKnh4 = 1/(1 + (0.00017/pKsys(11)))^2 ;
                obs.m(i).epKnh4 = wobs(sys.m(i).iKnh4);
            end
            if (isgood(obs.m(i).pKnh4))
                yobs(sys.m(i).iKnh4) = obs.m(i).pKnh4;
            else
                yobs(sys.m(i).iKnh4) = pKnh4;
                obs.m(i).pKnh4 = pKnh4;
            end
            if (isgood(obs.m(i).nh4))
                yobs(sys.m(i).inh4) = p(obs.m(i).nh4*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).inh4) = nan;
                obs.m(i).nh4 = nan;
            end
            if (isgood(obs.m(i).enh4))
                wobs(sys.m(i).inh4) = w(obs.m(i).nh4,obs.m(i).enh4);
            else
                wobs(sys.m(i).inh4) = nan;
                obs.m(i).enh4 = nan;
            end
            if (isgood(obs.m(i).nh3))
                yobs(sys.m(i).inh3) = p(obs.m(i).nh3*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).inh3) = nan;
                obs.m(i).nh3 = nan;
            end
            if (isgood(obs.m(i).enh3))
                wobs(sys.m(i).inh3) = w(obs.m(i).nh3,obs.m(i).enh3);
            else
                wobs(sys.m(i).inh3) = nan;
                obs.m(i).enh3 = nan;
            end
        end

        if (ismember('sulfide',sys.abr))
            if (isgood(obs.m(i).epKh2s))
                wobs(sys.m(i).iKh2s) = obs.m(i).epKh2s;
            else
                wobs(sys.m(i).iKh2s) = (0.033).^(-2);  % wKh2s = 1/(1 + (0.033/pKsys(11)))^2 ;
                obs.m(i).epKh2s = wobs(sys.m(i).iKh2s);
            end
            if (isgood(obs.m(i).pKh2s))
                yobs(sys.m(i).iKh2s) = obs.m(i).pKh2s;
            else
                yobs(sys.m(i).iKh2s) = pKh2s;
                obs.m(i).pKh2s = pKh2s;
            end
            if (isgood(obs.m(i).h2s))
                yobs(sys.m(i).ih2s) = p(obs.m(i).h2s*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).ih2s) = nan;
                obs.m(i).h2s = nan;
            end
            if (isgood(obs.m(i).eh2s))
                wobs(sys.m(i).ih2s) = w(obs.m(i).h2s,obs.m(i).eh2s);
            else
                wobs(sys.m(i).ih2s) = nan;
                obs.m(i).eh2s = nan;
            end
            if (isgood(obs.m(i).hs))
                yobs(sys.m(i).ihs) = p(obs.m(i).hs*1e-6); % convt µmol/kg to mol/kg
            else
                yobs(sys.m(i).ihs) = nan;
                obs.m(i).hs = nan;
            end
            if (isgood(obs.m(i).ehs))
                wobs(sys.m(i).ihs) = w(obs.m(i).hs,obs.m(i).ehs);
            else
                wobs(sys.m(i).ihs) = nan;
                obs.m(i).ehs = nan;
            end
        end
    end
    gun = @(z) grad_limpco2(z,yobs,wobs,sys);
    z0 = init(yobs,sys);
    tol = 1e-7;
    %keyboard
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
    %
    % populate est
    %
    % ebar = @(j) [ q( z(j) - sigy(j) ), q( z(j) + sigy(j) ) ];
    ebar = @(j) (0.5 * ( q( z(j) - sigy(j) ) - q( z(j) + sigy(j) ) ) );
        
    est.sal = z(sys.isal);
    est.esal = sigy(sys.isal);
    est.TC = q(z(sys.iTC));
    est.eTC = ebar(sys.iTC);
    est.TA = q(z(sys.iTA));
    est.eTA = ebar(sys.iTA);
    if ismember('borate', sys.abr)
        est.TB = q(z(sys.iTB))*1e6; % convt mol/kg to µmol/kg
        est.eTB = ebar(sys.iTB)*1e6;
    end
    if ismember('sulfate', sys.abr)
        est.TS = q(z(sys.iTS))*1e6; % convt mol/kg to µmol/kg
        est.eTS = ebar(sys.iTS)*1e6;
    end
    if ismember('fluoride', sys.abr)
        est.TF = q(z(sys.iTF))*1e6; % convt mol/kg to µmol/kg
        est.eTF = ebar(sys.iTF)*1e6;
    end
    if ismember('phosphate', sys.abr)
        est.TP = q(z(sys.iTP))*1e6; % convt mol/kg to µmol/kg
        est.eTP = ebar(sys.iTP)*1e6;
    end
    if ismember('silicate', sys.abr)
        est.TSi = q(z(sys.iTSi))*1e6; % convt mol/kg to µmol/kg
        est.eTSi = ebar(sys.iTSi)*1e6;
    end
    if ismember('ammonia', sys.abr)
        est.TNH3 = q(z(sys.iTNH3))*1e6; % convt mol/kg to µmol/kg
        est.eTNH3 = ebar(sys.iTNH3)*1e6;
    end
    if ismember('sulfide', sys.abr)
        est.TH2S = q(z(sys.iTH2S))*1e6; % convt mol/kg to µmol/kg
        est.eTH2S = ebar(sys.iTH2S)*1e6;
    end
    
    for i = 1:nTP
        % thermodynamic state & errorbar
        est.m(i).T  = z(sys.m(i).iT);
        est.m(i).eT = sigy(sys.m(i).iT);
        est.m(i).P  = z(sys.m(i).iP);        
        est.m(i).eP = sigy(sys.m(i).iP);
        % chemical equilibrium & errorbar
        est.m(i).ph      = z(sys.m(i).iph);
        est.m(i).eph     = sigy(sys.m(i).iph);
        est.m(i).fco2    = q(z(sys.m(i).ifco2))*1e6; % convt atm to µatm
        est.m(i).efco2   = ebar(sys.m(i).ifco2)*1e6;
        est.m(i).pco2    = q(z(sys.m(i).ipco2))*1e6;% convt atm to µatm
        est.m(i).epco2   = ebar(sys.m(i).ipco2)*1e6;
        est.m(i).pco2st  = z(sys.m(i).ico2st);
        est.m(i).epco2st = sigy(sys.m(i).ico2st);
        est.m(i).hco3   = q(z(sys.m(i).ihco3))*1e6;
        est.m(i).ehco3  = ebar(sys.m(i).ihco3)*1e6;
        est.m(i).co3    = q(z(sys.m(i).ico3))*1e6;
        est.m(i).eco3   = ebar(sys.m(i).ico3)*1e6;          
        % equilibrium pK's & errorbar
        est.m(i).pp2f  = z(sys.m(i).ip2f);
        est.m(i).epp2f = sigy(sys.m(i).ip2f);
        est.m(i).pK0   = z(sys.m(i).iK0);
        est.m(i).epK0  = sigy(sys.m(i).iK0);
        est.m(i).pK1   = z(sys.m(i).iK1);
        est.m(i).epK1  = sigy(sys.m(i).iK1);
        est.m(i).pK2   = z(sys.m(i).iK2);
        est.m(i).epK2  = sigy(sys.m(i).iK2);
        if ismember('borate', sys.abr)
            %chemical equilibrium & errorbar
            est.m(i).boh4  = q(z(sys.m(i).iboh4))*1e6; % convt mol/kg to µmol/kg
            est.m(i).eboh4 = ebar(sys.m(i).iboh4)*1e6;
            est.m(i).boh3  = q(z(sys.m(i).iboh3))*1e6;
            est.m(i).eboh3 = ebar(sys.m(i).iboh3)*1e6;
            % equilibrium pK & errorbar
            est.m(i).pKb  = z(sys.m(i).iKb);
            est.m(i).epKb = sigy(sys.m(i).iKb);
        end
        if ismember('water', sys.abr)
            % chemical equilibrium & errorbar
            est.m(i).oh   = q(z(sys.m(i).ioh))*1e6; % convt 
            est.m(i).eoh  = ebar(sys.m(i).ioh);
            % equilibrium pK & errorbar
            est.m(i).pKw  = z(sys.m(i).iKw);
            est.m(i).epKw = sigy(sys.m(i).iKw); 
        end
        if ismember('sulfate', sys.abr)
            % chemical equilibrium & errorbar
            est.m(i).phf   = z(sys.m(i).iphf);
            est.m(i).ephf  = sigy(sys.m(i).iphf);
            est.m(i).so4   = q(z(sys.m(i).iso4))*1e6; % convt
            est.m(i).eso4  = ebar(sys.m(i).iso4)*1e6;
            est.m(i).hso4  = q(z(sys.m(i).ihso4))*1e6;
            est.m(i).ehso4 = ebar(sys.m(i).ihso4)*1e6;
            % equilibrium pK & errorbar
            est.m(i).pKs  = z(sys.m(i).iKs);
            est.m(i).epKs = sigy(sys.m(i).iKs);
        end
        if ismember('fluoride', sys.abr)
            % chemical equilibrium & errorbar
            est.m(i).F    = q(z(sys.m(i).iF))*1e6; % convt
            est.m(i).eF   = ebar(sys.m(i).iF)*1e6;
            est.m(i).HF   = q(z(sys.m(i).iHF))*1e6;
            est.m(i).eHF  = ebar(sys.m(i).iHF)*1e6;
            % equilibrium pK & errorbar
            est.m(i).pKf  = z(sys.m(i).iKf);
            est.m(i).epKf = sigy(sys.m(i).iKf);
        end
        if ismember('phosphate', sys.abr)
            % chemical equilibrium & errorbar
            est.m(i).po4    = q(z(sys.m(i).ipo4))*1e6; % convt
            est.m(i).epo4   = ebar(sys.m(i).ipo4)*1e6;
            est.m(i).hpo4   = q(z(sys.m(i).ihpo4))*1e6;
            est.m(i).ehpo4  = ebar(sys.m(i).ihpo4)*1e6;
            est.m(i).h2po4  = q(z(sys.m(i).ih2po4))*1e6;
            est.m(i).eh2po4 = ebar(sys.m(i).ih2po4)*1e6;
            est.m(i).h3po4  = q(z(sys.m(i).ih3po4))*1e6;
            est.m(i).eh3po4 = ebar(sys.m(i).ih3po4)*1e6;
            % equilibrium pK & errorbar
            est.m(i).pK1p  = z(sys.m(i).iK1p);
            est.m(i).epK1p = sigy(sys.m(i).iK1p);
            est.m(i).pK2p  = z(sys.m(i).iK2p);
            est.m(i).epK2p = sigy(sys.m(i).iK2p);
            est.m(i).pK3p  = z(sys.m(i).iK3p);
            est.m(i).epK3p = sigy(sys.m(i).iK3p);
        end
        if ismember('silicate', sys.abr)
            % chemical equilibrium & errorbar
            est.m(i).sioh4   = q(z(sys.m(i).isioh4))*1e6; % convt
            est.m(i).esioh4  = ebar(sys.m(i).isioh4)*1e6;
            est.m(i).siooh3  = q(z(sys.m(i).isiooh3))*1e6;
            est.m(i).esiooh3 = ebar(sys.m(i).isiooh3)*1e6;
            % equilibrium pK & errorbar
            est.m(i).pKsi  = z(sys.m(i).iKsi);
            est.m(i).epKsi = sigy(sys.m(i).iKsi);
        end
        if ismember('ammonia', sys.abr)
            % chemical equilibrium & errorbar
            est.m(i).nh3    = q(z(sys.m(i).inh3))*1e6; % convt
            est.m(i).enh3   = ebar(sys.m(i).inh3)*1e6;
            est.m(i).nh4    = q(z(sys.m(i).inh4))*1e6;
            est.m(i).enh4   = ebar(sys.m(i).inh4)*1e6;
            % equilibrium pK & errorbar
            est.m(i).pKnh4  = z(sys.m(i).iKnh4);
            est.m(i).epKnh4 = sigy(sys.m(i).iKnh4);
        end
        if ismember('sulfide', sys.abr)
            % chemical equilibrium & errorbar
            est.m(i).hs    = q(z(sys.m(i).ihs))*1e6; % convt
            est.m(i).ehs   = ebar(sys.m(i).ihs)*1e6;
            est.m(i).h2s   = q(z(sys.m(i).ih2s))*1e6;
            est.m(i).eh2s  = ebar(sys.m(i).ih2s)*1e6;
            % equilibrium pK & errorbar
            est.m(i).pKh2s  = z(sys.m(i).ipKh2s);
            est.m(i).epKh2s = sigy(sys.m(i).iKh2s);
        end
    end
end

function [g,H] = grad_limpco2(z,y,w,sys)
    I = eye(length(z));
    H = zeros(length(z),length(z));
    im = sqrt(-1);
    %test_g = zeros(length(z),1);
    for k = 1:length(z)
        [f,g] = limpco2( z + im * eps^3 * I(:,k), y, w, sys);
        %   test_g(k) = imag(f)/eps^3;
        H(k,:) = imag( g(:) ) / eps^3 ;
    end
    g = real( g(:) );
    % [ g, test_g ]
    %keyboard
end

function [f,g] = limpco2(z,y,w,sys)
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
    nTP = length(sys.m);
    if (ismember('sulfate',sys.abr)) % MF doesn't understand why this is here
        nlam = nlam+nTP; % one extra lagrange multiplier for each (T,P)-dependent free to total ph conversions
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
    
    nrk = size(K,1);
    zpK = zeros(nrk,1);
    zgpK = zeros(nrk,nv);
    f2t = [];
    for i = 1:nTP
        [pK, gpK] = local_pK( x(sys.m(i).iT), x(sys.isal), x(sys.m(i).iP) );
        iTSP = [ sys.m(i).iT, sys.isal, sys.m(i).iP];
        if (ismember('carbonate',sys.abr))
            zpK(sys.m(i).kK0)          = pK(1); 
            zgpK(sys.m(i).kK0, iTSP )  = gpK(1,:); % ∂T, ∂S, ∂P
        
            zpK(sys.m(i).kK1)          = pK(2);  
            zgpK(sys.m(i).kK1, iTSP )  = gpK(2,:); % ∂T, ∂S, ∂P       
        
            zpK(sys.m(i).kK2)          = pK(3); 
            zgpK(sys.m(i).kK2, iTSP )  = gpK(3,:); % ∂T, ∂S, ∂P  
        
            zpK(sys.m(i).kp2f)         = pK(14); 
            zgpK(sys.m(i).kp2f, iTSP ) = gpK(14,:); % ∂T, ∂S, ∂P 
        end
    
        if (ismember('borate',sys.abr))
            zpK(sys.m(i).kKb)          = pK(4); 
            zgpK(sys.m(i).kKb, iTSP ) = gpK(4,:); % ∂T, ∂S, ∂P
        end

        if (ismember('water',sys.abr))
            zpK(sys.m(i).kKw)          = pK(5); 
            zgpK(sys.m(i).kKw, iTSP )  = gpK(5,:); % ∂T, ∂S, ∂P  
        end

        if (ismember('sulfate',sys.abr))
            zpK(sys.m(i).kKs)          = pK(6); 
            zgpK(sys.m(i).kKs, iTSP )  = gpK(6,:); % ∂T, ∂S, ∂P  
            f2t = [f2t;sys.m(i).f2t(x)];
        end
        
        if (ismember('fluoride',sys.abr))
            zpK(sys.m(i).kKf)          = pK(7); 
            zgpK(sys.m(i).kKf, iTSP )  = gpK(7,:); % ∂T, ∂S, ∂P
        end
        
        if (ismember('phosphate',sys.abr))
            zpK(sys.m(i).kK1p)         = pK(8); 
            zgpK(sys.m(i).kK1p, iTSP ) = gpK(8,:); % ∂T, ∂S, ∂P 
            
            zpK(sys.m(i).kK2p)         = pK(9); 
            zgpK(sys.m(i).kK2p,iTSP )  = gpK(9,:); % ∂T, ∂S, ∂P 
            
            zpK(sys.m(i).kK3p)         = pK(10); 
            zgpK(sys.m(i).kK3p, iTSP ) = gpK(10,:); % ∂T, ∂S, ∂P 
        end

        if (ismember('silicate',sys.abr))
            zpK(sys.m(i).kKsi)         = pK(11); 
            zgpK(sys.m(i).kKsi, iTSP ) = gpK(11,:); % ∂T, ∂S, ∂P 
        end

        if (ismember('ammonia',sys.abr))
            zpK(sys.m(i).kKnh4)         = pK(12); 
            zgpK(sys.m(i).kKnh4, iTSP ) = gpK(12,:); % ∂T, ∂S, ∂P  
        end

        if (ismember('sulfide',sys.abr))
            zpK(sys.m(i).kKh2s)         = pK(13); 
            zgpK(sys.m(i).kKh2s, iTSP ) = gpK(13,:); % ∂T, ∂S, ∂P 
        end
    end
    
    % constraint equations
    if (ismember('sulfate',sys.abr))  
        c = [  M * q( x ); ... 
               (-K * x) + zpK;...
               f2t ] ; 
    else
        c = [  M * q( x ); ...
               (-K * x) + zpK  ] ; 
    end

    f = 0.5 *  e.' * W * e  + lam.' * c ;  % limp, method of lagrange multipliers
        
    % -(-1/2 sum of squares) + constraint eqns, minimize f => grad(f) = 0
    if ( nargout > 1 ) % compute the gradient
        if (ismember('sulfate',sys.abr))
            gf2t = zeros(nTP,nv);
            for i = 1:nTP
                gf2t(i,[ sys.m(i).iph, sys.m(i).iKs, sys.iTS, sys.m(i).iphf ] ) = sys.m(i).gf2t(x);
            end
            dcdx = [ M * diag( sys.dqdx( x ) ); ...
                     (-K + zgpK) ;...
                     gf2t ]; % constraint eqns wrt -log10(concentrations)
        else
            dcdx = [ M * diag( sys.dqdx( x ) ); ...
                     (-K + zgpK) ]; % c'
        end    
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
        if (ismember('hf',sys.variables))
            dhfdx2 = zeros(nc,nc);
            ii = [sys.iKs,sys.iTS];
            dhfdx2(ii,ii) = sys.ggf2t(x);
            gg = gg + lam(nr+1)*dhfdx2;
        end
        H = [  PP.'*W*PP + gg , dcdx.'  ; ...
               dcdx         , zeros(nlam)  ];
    end
end




function z0 = init(yobs,sys);
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
        h = 0.5*( ( gam - 1 ) * K1 + ( ( 1 - gam )^2 * K1^2 - 4 * K1 * K2 * ( 1 - 2 * gam ) ).^0.5 ) ;    
        hco3 =  h * alk / (h + 2 * K2 );
        co2st = h * hco3 / K1 ;
        co3 = 0.5 * ( alk - hco3 ) ;
        fco2 = co2st/K0;
        pco2 = fco2/p2f;
    
        y0(sys.m(i).iph)    = p(h);
        y0(sys.m(i).ihco3)  = p(hco3);
        y0(sys.m(i).ico2st) = p(co2st);
        y0(sys.m(i).ico3)   = p(co3);
        y0(sys.m(i).ifco2)  = p(fco2);
        y0(sys.m(i).ipco2)  = p(pco2);
    
        if (ismember('water',sys.abr))
            Kw = q(y0(sys.m(i).iKw));
            oh = Kw / h;
            y0(sys.m(i).ioh) = p(oh);
        end
        
        if (ismember('borate',sys.abr));
            Kb = q(y0(sys.m(i).iKb));
            TB = q(yobs(sys.iTB));
            boh4 = TB / ( 1 + h / Kb );
            boh3 = TB - boh4;        
            y0(sys.iTB)   = p(TB);
            y0(sys.m(i).iboh3) = p(boh3);
            y0(sys.m(i).iboh4) = p(boh4);
        end
    
        if (ismember('sulfate',sys.abr));
            Ks = q(y0(sys.m(i).iKs));
            TS = q(yobs(sys.iTS));
            hf = h / ( 1 + TS / Ks );
            hso4 = TS / ( 1 + Ks / hf);
            so4  = Ks * hso4 / hf;
            y0(sys.iTS)   = p(TS);
            y0(sys.m(i).iphf)   = p(hf);
            y0(sys.m(i).ihso4) = p(hso4);
            y0(sys.m(i).iso4)  = p(so4);
        end
    
        if (ismember('fluoride',sys.abr));
            Kf = q(y0(sys.m(i).iKf));
            TF = q(yobs(sys.iTF));
            HF = TF / ( 1 + Kf / h );
            F  = Kf * HF / h;
            y0(sys.iTF) = p(TF);
            y0(sys.m(i).iF)  = p(F);
            y0(sys.m(i).iHF) = p(HF);
        end
    
        if (ismember('phosphate',sys.abr));
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
    
        if (ismember('silicate',sys.abr));
            Ksi = q(y0(sys.m(i).iKsi));
            TSi = q(yobs(sys.iTSi));
            siooh3 = TSi / ( 1 + h / Ksi );
            sioh4  = TSi - siooh3;
            y0(sys.iTSi)    = p(TSi);
            y0(sys.m(i).isiooh3) = p(siooh3);
            y0(sys.m(i).isioh4)  = p(sioh4);
        end
        if (ismember('ammonia',sys.abr));
            error('Need to implement initialization for ammonia');
        end
        if (ismember('sulfide',sys.abr))
            error('Need to implement initialization for sulfide');
        end
    end
    % add the Lagrange multipliers
    nlam = size(sys.M,1)+size(sys.K,1);
    if(ismember('sulfate',sys.abr))
        nlam = nlam+nTP;
    end
    lam = zeros(nlam,1);
    z0 = [y0(:);lam(:)];
end

