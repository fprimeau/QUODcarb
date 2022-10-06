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

    isgood = @(thing) (~isempty(thing) & ~isnan(thing));
    p = sys.p;
    q = sys.q;
    w = @(x,e) p(1+e./x).^(-2); % convert x+/-e into precision for p(x)
    
    nv = size(sys.K,2);
    yobs = nan(nv,1);
    wobs = nan(nv,1);
    % make sure all the required fields in the obs struct exist
    if (~isfield(obs, 'sal'))
        error('Need to provide salinity measurement.');
    else
        yobs(sys.isal) = obs.sal;
    end
    
    if (~isfield(obs, 'wsal'));
        obs.wsal = (0.002).^(-2);
        fprintf('Warning: Assuming salinity uncertainty is 0.002 PSU');
    else
        wobs(sys.isal) = obs.wsal;
    end
    if (~isfield(obs,'TC'))
        obs.TC = [];
        yobs(sys.iTC) = nan;
    else
        yobs(sys.iTC) = p(obs.TC);
    end
    if (~isfield(obs,'wTC'))
        obs.wTC = [];
        wobs(sys.iTC) = nan;
    end 
    if(~isfield(obs,'TA'))
        obs.TA = [];
        yobs(sys.iTA) = nan;
    else
        yobs(sys.iTA) = p(obs.TA);
    end
    if (~isfield(obs,'wTA'))
        obs.wTA = [];
        wobs(sys.iTA) = nan;
    else
        wobs(sys.iTA) = obs.wTA;
    end
    if (~isfield(obs,'m'))
        error('Need to provide temperature and pressure measurement.')
    end
    
    if (~isfield(obs.m(1),'wK0'))
        obs.m(1).wK0 = [];
    end
    if (~isfield(obs.m(1),'pK0'))
        obs.m(1).pK0 = [];
    end
    if (~isfield(obs.m(1),'wK1'))
        obs.m(1).wK1 = [];
    end
    if (~isfield(obs.m(1),'pK1'))
        obs.m(1).pK1 = [];
    end
    if (~isfield(obs.m(1),'wK2'))
        obs.m(1).wK2 = [];
    end
    if (~isfield(obs.m(1),'pK2'))
        obs.m(1).pK2 = [];
    end
    if (~isfield(obs.m(1),'wp2f'))
        obs.m(1).wp2f = [];
    end
    if (~isfield(obs.m(1),'pp2f'))
        obs.m(1).pp2f = [];
    end
    if (~isfield(obs.m(1),'wfco2'))
        obs.m(1).wfco2 = [];
    end
    if (~isfield(obs.m(1),'fco2'))
        obs.m(1).fco2 = [];
    end
    if (~isfield(obs.m(1),'wco2st'))
        obs.m(1).wco2st = [];
    end
    if (~isfield(obs.m(1),'co2st'))
        obs.m(1).co2st = [];
    end
    if (~isfield(obs.m(1),'whco3'))
        obs.m(1).whco3 = [];
    end
    if (~isfield(obs.m(1),'hco3'))
        obs.m(1).hco3 = [];
    end
    if (~isfield(obs.m(1),'wco3'))
        obs.m(1).wco3 = [];
    end
    if (~isfield(obs.m(1),'co3'))
        obs.m(1).co3 = [];
    end
    if (~isfield(obs.m(1),'wph'))
        obs.m(1).ph = [];
    end
    if (~isfield(obs.m(1),'ph'))
        obs.m(1).wph = [];
    end


    if (ismember('borate',sys.abr))
        if (~isfield(obs.m(1),'wKb'))
            obs.m(1).wKb = [];
        end
        if (~isfield(obs.m(1),'pKb'))
            obs.m(1).pKb = [];
        end
        if (~isfield(obs,'TB'))
            % Culkin, F., in Chemical Oceanography,
            % ed. Riley and Skirrow, 1965: ( copied from Orr's code )
            fprintf('Warning: Assuming TB = 0.0004106 * sal / 35  mol/kg-SW.\n');
            obs.TB = 0.0004106 * obs.sal / 35 ;
            yobs(sys.iTB) = p(obs.TB);
        else
            yobs(sys.iTB) = p(obs.TB);
        end
        if (~isfield(obs, 'wTB'))
            obs.wTB = ( 35 / 0.0004106 )^2 * obs.wsal ;
            wobs(sys.iTB) = obs.wTB;
        else
            wobs(sys.iTB) = obs.wTB;
        end
        if (~isfield(obs.m(1),'wboh4'))
            obs.m(1).wboh4 = [];
        end
        if (~isfield(obs.m(1),'boh4'))
            obs.m(1).boh4 = [];
        end
        if (~isfield(obs.m(1),'wboh3'))
            obs.m(1).wboh3 = [];
        end
        if (~isfield(obs.m(1),'boh3'))
            obs.m(1).boh3 = [];
        end
    end

    if (ismember('water',sys.abr))
        if (~isfield(obs.m(1), 'wKw'))
            obs.m(1).wKw = [];
        end
        if (~isfield(obs.m(1), 'pKw'))
            obs.m(1).pKw = [];
        end
        if (~isfield(obs.m(1),'woh'));
            obs.m(1).woh = [];
        end
        if (~isfield(obs.m(1),'oh'));
            obs.m(1).oh = [];
        end
    end

    if (ismember('sulfate',sys.abr))
        if (~isfield(obs.m(1), 'wKs'))
            obs.m(1).wKs = [];
        end
        if (~isfield(obs.m(1), 'pKs'))
            obs.m(1).pKs = [];
        end
        if (~isfield(obs, 'TS'))
            % Morris, A. W., and Riley, J. P., Deep-Sea Research 13:699-705, 1966:
            % copied from Orr's code
            fprintf('Warning: Assuming TS = ( 0.14 / 96.062 ) * ( sal / 1.80655 ) mol/kg-SW.\n'); 
            obs.TS = ( 0.14 / 96.062 ) * ( obs.sal / 1.80655 );
            yobs(sys.iTS) = p(obs.TS);
        else
            yobs(sys.iTS) = p(obs.TS);
        end
        if (~isfield(obs, 'wTS'))
            obs.wTS = ( ( 0.14 / 96.062 ) / 1.80655 )^(-2) * obs.wsal;
            wobs(sys.iTS) = obs.wTS;
        else
            wobs(sys.iTS) = obs.wTS;
        end
        if (~isfield(obs.m(1), 'phf'))
            obs.m(1).phf = [];
        end
        if (~isfield(obs.m(1), 'wphf'))
            obs.m(1).wphf = [];
        end
        if (~isfield(obs.m(1), 'so4'))
            obs.m(1).so4 = [];
        end
        if (~isfield(obs.m(1), 'wso4'))
            obs.m(1).wso4 = [];
        end
        if (~isfield(obs.m(1), 'hso4'))
            obs.m(1).hso4 = [];
        end     
        if (~isfield(obs.m(1), 'whso4'))
            obs.m(1).whso4 = [];
        end     
    end

    if (ismember('fluoride',sys.abr))
        if (~isfield(obs.m(1), 'wKf'))
            obs.m(1).wKf = [];
        end
        if (~isfield(obs.m(1), 'pKf'))
            obs.m(1).pKf = [];
        end
        if (~isfield(obs, 'TF'))
            % Riley, J. P., Deep-Sea Research 12:219-220, 1965:
            % this is .000068.*Sali./35. = .00000195.*Sali
            fprintf('Warning: Assuming TF = ( 0.000067 / 18.998 ) * ( sal / 1.80655 ). mol/kg-SW. \n');
            obs.TF = ( 0.000067 / 18.998 ) * ( obs.sal / 1.80655 ); % in mol/kg-SW
            yobs(sys.iTF) = p(obs.TF);
        else
            yobs(sys.iTF) = p(obs.TF);
        end
        if (~isfield(obs, 'wTF'))
            obs.wTF = ( ( 0.000067 / 18.998 ) / 1.80655 )^(-2) * obs.wsal;
            wobs(sys.iTF) = obs.wTF;
        else
            wobs(sys.iTF) = obs.wTF;
        end
        if (~isfield(obs.m(1), 'F'))
            obs.m(1).F = [];
        end
        if (~isfield(obs.m(1), 'wF'))
            obs.m(1).wF = [];
        end
        if (~isfield(obs.m(1), 'HF'))
            obs.m(1).HF = [];
        end
        if (~isfield(obs.m(1), 'wHF'))
            obs.m(1).wHF = [];
        end
    end

    if (ismember('phosphate',sys.abr))
        if (~isfield(obs.m(1), 'wK1p'))
            obs.m(1).wK1p = [];
        end
        if (~isfield(obs.m(1), 'pK1p'))
            obs.m(1).pK1p = [];
        end
        if (~isfield(obs.m(1), 'wK2p'))
            obs.m(1).wK2p = [];
        end
        if (~isfield(obs.m(1), 'pK2p'))
            obs.m(1).pK2p = [];
        end
        if (~isfield(obs.m(1), 'wK3p'))
            obs.m(1).wK3p = [];
        end
        if (~isfield(obs.m(1), 'pK3p'))
            obs.m(1).pK3p = [];
        end
        if (~isfield(obs, 'TP'))
            fprintf('Warning: Assuming TP = 1e-9 mol/kg-SW.\n');
            obs.TP = 1e-9;
            yobs(sys.iTP) = p(obs.TP);
        else
            yobs(sys.iTP) = p(obs.TP);
        end
        if (~isfield(obs, 'wTP'))
            fprintf('Warning: Setting wTP = 1000 (mol/kg-SW)^-2.\n')
            obs.wTP = 1000;
            wobs(sys.iTP) = obs.wTP;
        else
            wobs(sys.iTP) = obs.wTP;
        end
        if (~isfield(obs.m(1), 'h3po4'))
            obs.m(1).h3po4 = [];
        end
        if (~isfield(obs.m(1), 'wh3po4'))
            obs.m(1).wh3po4 = [];
        end
        if (~isfield(obs.m(1), 'h2po4'))
            obs.m(1).h2po4 = [];
        end
        if (~isfield(obs.m(1), 'wh2po4'))
            obs.m(1).wh2po4 = [];
        end
        if (~isfield(obs.m(1), 'hpo4'))
            obs.m(1).hpo4 = [];
        end
        if (~isfield(obs.m(1), 'whpo4'))
            obs.m(1).whpo4 = [];
        end
        if (~isfield(obs.m(1), 'po4'))
            obs.m(1).po4 = [];
        end
        if (~isfield(obs.m(1), 'wpo4'))
            obs.m(1).wpo4 = [];
        end
    end

    if (ismember('silicate',sys.abr))
        if (~isfield(obs.m(1), 'wKsi'))
            obs.m(1).wKsi = [];
        end
        if (~isfield(obs.m(1), 'pKsi'))
            obs.m(1).pKsi = [];
        end
        if (~isfield(obs, 'TSi'));
            fprintf('Warning: Assuming TSi = 1e-9  mol/kg-SW.\n');
            obs.TSi = 1e-9;
            yobs(sys.iTSi) = p(obs.TSi);
        else
            yobs(sys.iTSi) = p(obs.TSi);
        end
        if (~isfield(obs, 'wTSi'))
            fprintf('Warning: Assuming wTSi = 1000  (mol/kg-SW)^-2.\n');
            obs.wTSi = 1000;
            wobs(sys.iTSi) = obs.wTSi;
        else
            wobs(sys.iTSi) = obs.wTSi;
        end
        if (~isfield(obs.m(1), 'sioh4'))
            obs.m(1).sioh4 = [];
        end
        if (~isfield(obs.m(1), 'wsioh4'))
            obs.m(1).wsioh4 = [];
        end
        if (~isfield(obs.m(1), 'siooh3'))
            obs.m(1).siooh3 = [];
        end
        if (~isfield(obs.m(1), 'wsiooh3'))
            obs.m(1).wsiooh3 = [];
        end
    end

    if (ismember('ammonia',sys.abr))
        if (~isfield(obs.m(1), 'wKnh4'))
            obs.m(1).wKnh4 = [];
        end
        if (~isfield(obs.m(1), 'pKnh4'))
            obs.m(1).pKnh4 = [];
        end
        if (~isfield(obs, 'TNH3'));
            fprintf('Warning: Assuming TNH3 = 1e-9  mol/kg-SW.\n');
            obs.TNH3 = 0;
            yobs(sys.iTNH3) = p(obs.TNH3);
        else
            yobs(sys.iTNH3) = p(obs.TNH3);
        end
        if (~isfield(obs, 'wTNH3'));
            fprintf('Warning: Assuming wTNH3 = 1000  (mol/kg-SW)^-2.\n');
            obs.wTNH3 = 1000;
            wobs(sys.iTNH3) = obs.wTNH3;
        else
            wobs(sys.iTNH3) = obs.wTNH3;
        end
        
        if (~isfield(obs.m(1), 'nh3'))
            obs.m(1).nh3 = [];
        end
        if (~isfield(obs.m(1), 'nh4'))
            obs.m(1).nh4 = [];
        end
    end

    if (ismember('sulfide',sys.abr))
        if (~isfield(obs.m(1), 'wKh2s'))
            obs.m(1).wKh2s = [];
        end
        if (~isfield(obs.m(1), 'pKh2s'))
            obs.m(1).pKh2s = [];
        end
        if (~isfield(obs, 'TH2S'))
            fprintf('Warning: Assuming TH2S = 1e-9  mol/kg-SW.\n');
            obs.TH2S = 0;
            yobs(sys.iTH2S) = p(obs.TH2S);
        else
            yobs(sys.iTH2S) = p(obs.TH2S);
        end
        if (~isfield(obs, 'wTH2S'))
            fprintf('Warning: Assuming wTH2S = 1000  (mol/kg-SW)^-2.\n');
            obs.wTH2S = 1000;
            wobs(sys.iTH2S) = obs.wTH2S;
        else
            wobs(sys.iTH2S) = obs.wTH2S;
        end
        if (~isfield(obs.m(1), 'hs'))
            obs.m(1).hs = [];
        end
        if (~isfield(obs.m(1), 'h2s'))
            obs.m(1).h2s = [];
        end
    end
        
    nTP = length(obs.m);
    yobs(sys.iTC) = obs.TC;
    yobs(sys.iTA) = obs.TA;
    
    for i = 1:nTP % loop over all the pressure and temperature sensitive components
        yobs(sys.m(i).iT) = obs.m(i).T;
        yobs(sys.m(i).iP) = obs.m(i).P;
        
        wobs(sys.m(i).iT) = obs.m(i).wT;
        wobs(sys.m(i).iP) = obs.m(i).wP;
        
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
            if (isgood(obs.m(i).wK0))
                wobs(sys.m(i).iK0) = obs.m(i).wK0;
            else
                wobs(sys.m(i).iK0) = (0.002).^(-2);
                obs.m(i).wK0 = wobs(sys.m(i).iK0);
            end
            if (isgood(obs.m(i).pK0))
                yobs(sys.m(i).iK0) = obs.m(i).pK0;
            else
                yobs(sys.m(i).iK0) = pK0;
                obs.m(i).pK0 = pK0;
            end

            if (isgood(obs.m(i).wK1))
                wobs(sys.m(i).iK1) = obs.m(i).wK1;
            else
                wobs(sys.m(i).iK1) = (0.01).^(-2);  % wK1 = 1/(1 + (0.01/pKsys(2)))^2 ;
                obs.m(i).wK1 = wobs(sys.m(i).iK1);
            end
            if (isgood(obs.m(i).pK1))
                yobs(sys.m(i).iK1) = obs.m(i).pK1;
            else
                yobs(sys.m(i).iK1) = pK1;
                obs.m(i).pK1 = pK1;
            end            
            
            if (isgood(obs.m(i).wK2))
                wobs(sys.m(i).iK2) = obs.m(i).wK2;
            else
                wobs(sys.m(i).iK2) = (0.02).^(-2);  % wK2 = 1/(1 + (0.02/pKsys(3)))^2 ;
                obs.m(i).wK2 = wobs(sys.m(i).iK2);
            end
            
            if (isgood(obs.m(i).pK2))
                yobs(sys.m(i).iK2) = obs.m(i).pK2;
            else
                yobs(sys.m(i).iK2) = pK2;
                obs.m(i).pK2 = pK2;
            end
            
            if (isgood(obs.m(i).wp2f))
                wobs(sys.m(i).ip2f) = obs.m(i).wp2f;
            else
                wobs(sys.m(i).ip2f) = (0.001).^(-2);
                obs.m(i).wp2f = wobs(sys.m(i).ip2f);
            end
            
            if (isgood(obs.m(i).pp2f))
                yobs(sys.m(i).ip2f) = obs.m(i).pp2f;
            else
                yobs(sys.m(i).ip2f) = pp2f;
                obs.m(i).pp2f = pp2f;
            end
            
            if (isgood(obs.m(i).co2st))
                yobs(sys.m(i).ico2st) = p(obs.m(i).co2st);
            else
                yobs(sys.m(i).ico2st) = nan;
            end
            if (isgood(obs.m(i).wco2st))
                wobs(sys.m(i).ico2st) = obs.m(i).wco2st;
            else
                wobs(sys.m(i).ico2st) = nan;
            end
            if (isgood(obs.m(i).hco3))
                yobs(sys.m(i).ihco3) = p(obs.m(i).hco3);
            else
                yobs(sys.m(i).ihco3) = nan;
            end
            if (isgood(obs.m(i).whco3))
                wobs(sys.m(i).ihco3) = obs.m(i).whco3;
            else
                wobs(sys.m(i).ihco3) = nan;
            end
            if (isgood(obs.m(i).co3))
                yobs(sys.m(i).ico3) = p(obs.m(i).co3);
            else
                yobs(sys.m(i).ico3) = nan;
            end
            if (isgood(obs.m(i).wco3))
                wobs(sys.m(i).ico3) = obs.m(i).wco3;
            else
                wobs(sys.m(i).ico3) = nan;
            end
            if (isgood(obs.m(i).fco2))
                yobs(sys.m(i).ifco3) = p(obs.m(i).fco2);
            else
                yobs(sys.m(i).ifco2) = nan;
            end
            if (isgood(obs.m(i).wfco2))
                wobs(sys.m(i).ifco2) = obs.m(i).wfco2;
            else
                wobs(sys.m(i).ifco2) = nan;
            end
            if (isgood(obs.m(i).ph))
                yobs(sys.m(i).iph) = obs.m(i).ph;
            else
                yobs(sys.m(i).iph) = nan;
            end
            if (isgood(obs.m(i).wph))
                wobs(sys.m(i).iph) = obs.m(i).wph;
            else
                wobs(sys.m(i).iph) = nan;
            end        
        end
        
        if (ismember('borate',sys.abr))
            if (isgood(obs.m(i).wKb))
                wobs(sys.m(i).iKb) = obs.m(i).wKb;
            else
                wobs(sys.m(i).iKb) = (0.01).^(-2);  % wKb = 1/(1 + (0.01/pKsys(4)))^2 ;
                obs.m(i).wKb = wobs(sys.m(i).iKb);
            end
            if (isgood(obs.m(i).pKb))
                yobs(sys.m(i).iKb) = obs.m(i).pKb;
            else
                yobs(sys.m(i).iKb) = pKb;
                obs.m(i).pKb = pKb;
            end
            if (isgood(obs.m(i).boh3))
                yobs(sys.m(i).iboh3) = p(obs.m(i).boh3);
            else
                yobs(sys.m(i).iboh3) = nan;
            end
            if (isgood(obs.m(i).wboh3))
                wobs(sys.m(i).iboh3) = obs.m(i).wboh3;
            else
                wobs(sys.m(i).iboh3) = nan;
            end
            if (isgood(obs.m(i).boh4))
                yobs(sys.m(i).iboh4) = p(obs.m(i).boh4);
            else
                yobs(sys.m(i).iboh4) = nan;
            end
            if (isgood(obs.m(i).wboh4))
                wobs(sys.m(i).iboh4) = obs.m(i).wboh4;
            else
                wobs(sys.m(i).iboh4) = nan;
            end
        end
        
        if (ismember('water',sys.abr))
            if (isgood(obs.m(i).wKw))
                wobs(sys.m(i).iKw) = obs.m(i).wKw;
            else
                wobs(sys.m(i).iKw) = (0.01).^(-2);  % wKw = 1/(1 + (0.01/pKsys(5)))^2 ;
                obs.m(i).wKw = wobs(sys.m(i).iKw);
            end
            if (isgood(obs.m(i).pKw))
                wobs(sys.m(i).iKw) = obs.m(i).pKw;
            else
                yobs(sys.m(i).iKw) = pKw;
                obs.m(i).pKw = pKw;
            end
            if (isgood(obs.m(i).oh))
                yobs(sys.m(i).ioh) = p(obs.m(i).oh);
            else
                yobs(sys.m(i).ioh) = nan;
            end
            if (isgood(obs.m(i).woh))
                wobs(sys.m(i).ioh) = obs.m(i).woh;
            else
                wobs(sys.m(i).ioh) = nan;
            end
            
        end
        
        if (ismember('sulfate',sys.abr))
            if (isgood(obs.m(i).wKs))
                wobs(sys.m(i).iKs) = obs.m(i).wKs;
            else
                wobs(sys.m(i).iKs) = (0.0021).^(-2); % wKs = 1/(1 + (0.0021/pKsys(6)))^2 ;
                obs.m(i).wKs = wobs(sys.m(i).iKs);
            end
            if (isgood(obs.m(i).pKs))
                yobs(sys.m(i).iKs) = obs.m(i).pKs;
            else
                yobs(sys.m(i).iKs) = pKs;
                obs.m(i).pKs = pKs;
            end
            if (isgood(obs.m(i).hso4))
                yobs(sys.m(i).ihso4) = p(obs.m(i).hso4);
            else
                yobs(sys.m(i).ihso4) = nan;
            end
            if (isgood(obs.m(i).whso4))
                wobs(sys.m(i).ihso4) = obs.m(i).whso4;
            else
                wobs(sys.m(i).ihso4) = nan;
            end
            if (isgood(obs.m(i).so4))
                yobs(sys.m(i).iso4) = p(obs.m(i).so4);
            else
                yobs(sys.m(i).iso4) = nan;
            end
            if (isgood(obs.m(i).wso4))
                wobs(sys.m(i).iso4) = obs.m(i).wso4;
            else
                wobs(sys.m(i).iso4) = nan;
            end
            if (isgood(obs.m(i).phf))
                yobs(sys.m(i).iphf) = obs.m(i).phf;
            else
                yobs(sys.m(i).iphf) = nan;
            end
            if (isgood(obs.m(i).wphf))
                wobs(sys.m(i).iphf) = obs.m(i).wphf;
            else
                wobs(sys.m(i).iphf) = nan;
            end

        end
        
        if (ismember('fluoride',sys.abr))
            if (isgood(obs.m(i).wKf))
                wobs(sys.m(i).iKf) = obs.m(i).wKf;
            else
                wobs(sys.m(i).iKf) = (0.02).^(-2);   % wKF = 1/(p(1 + 0.02/KF))^2 ; % 2% relative error
                obs.m(i).wKf = wobs(sys.m(i).iKf);
            end
            if (isgood(obs.m(i).pKf))
                yobs(sys.m(i).iKf) = obs.m(i).pKf;
            else
                yobs(sys.m(i).iKf) = pKf;
                obs.m(i).pKf = pKf;
            end
            if (isgood(obs.m(i).F))
                yobs(sys.m(i).iF) = p(obs.m(i).F);
            else
                yobs(sys.m(i).iF) = nan;
            end
            if (isgood(obs.m(i).wF))
                wobs(sys.m(i).iF) = obs.m(i).wF;
            else
                wobs(sys.m(i).iF) = nan;
            end
            if (isgood(obs.m(i).HF))
                yobs(sys.m(i).iHF) = p(obs.m(i).HF);
            else
                yobs(sys.m(i).iHF) = nan;
            end
            if (isgood(obs.m(i).wF))
                wobs(sys.m(i).iHF) = obs.m(i).wHF;
            else
                wobs(sys.m(i).iHF) = nan;
            end
        end
        
        if (ismember('phosphate',sys.abr))
            if (isgood(obs.m(i).wK1p))
                wobs(sys.m(i).iK1p) = obs.m(i).wK1p;
            else
                wobs(sys.m(i).iK1p) = (0.09).^(-2);  % wK1p = 1/(1 + (0.09/pKsys(8)))^2 ;
                obs.m(i).wK1p = wobs(sys.m(i).iK1p);
            end
            if (isgood(obs.m(i).pK1p))
                yobs(sys.m(i).iK1p) = obs.m(i).pK1p;
            else
                yobs(sys.m(i).iK1p) = pK1p;
                obs.m(i).pK1p = pK1p;
            end
            if (isgood(obs.m(i).wK2p))
                wobs(sys.m(i).iK2p) = obs.m(i).wK2p;
            else
                wobs(sys.m(i).iK2p) = (0.03).^(-2);  % wK2p = 1/(1 + (0.03/pKsys(9)))^2 ;
                obs.m(i).wK2p = wobs(sys.m(i).iK2p);
            end
            if (~isgood(obs.m(i).pK2p))
                yobs(sys.m(i).iK2p) = obs.m(i).pK2p;
            else
                yobs(sys.m(i).iK2p) = pK2p;
                obs.m(i).pK2p = pK2p;
            end
            if (isgood(obs.m(i).wK3p))
                wobs(sys.m(i).iK3p) = obs.m(i).wK3p;
            else
                wobs(sys.m(i).iK3p) = (0.02).^(-2);  % wK3p = 1/(1 + (0.02/pKsys(10)))^2 ;
                obs.m(i).wK3p = wobs(sys.m(i).iK3p);
            end
            if (isgood(obs.m(i).pK3p))
                yobs(sys.m(i).iK3p) = obs.m(i).pK3p;
            else
                yobs(sys.m(i).iK3p) = pK3p;
                obs.m(i).pK3p = pK3p;
            end
            if (isgood(obs.m(i).h3po4))
                yobs(sys.m(i).ih3po4) = p(obs.m(i).h3po4);
            else
                yobs(sys.m(i).ih3po4) = nan;
            end
            if (isgood(obs.m(i).wh3po4))
                wobs(sys.m(i).ih3po4) = obs.m(i).wh3po4;
            else
                wobs(sys.m(i).ih3po4) = nan;
            end
            if (isgood(obs.m(i).h2po4))
                yobs(sys.m(i).ih2po4) = p(obs.m(i).h2po4);
            else
                yobs(sys.m(i).ih2po4) = nan;
            end
            if (isgood(obs.m(i).wh2po4))
                wobs(sys.m(i).ih2po4) = obs.m(i).wh2po4;
            else
                wobs(sys.m(i).ih2po4) = nan;
            end
            if (isgood(obs.m(i).hpo4))
                yobs(sys.m(i).ihpo4) = p(obs.m(i).hpo4);
            else
                yobs(sys.m(i).ihpo4) = nan;
            end
            if (isgood(obs.m(i).whpo4))
                wobs(sys.m(i).ihpo4) = obs.m(i).whpo4;
            else
                wobs(sys.m(i).ihpo4) = nan;
            end
            if (isgood(obs.m(i).po4))
                yobs(sys.m(i).ipo4) = p(obs.m(i).po4);
            else
                yobs(sys.m(i).ipo4) = nan;
            end
            if (isgood(obs.m(i).wpo4))
                wobs(sys.m(i).ipo4) = obs.m(i).wpo4;
            else
                wobs(sys.m(i).ipo4) = nan;
            end

        end

        if (ismember('silicate',sys.abr))
            if (isgood(obs.m(i).wKsi))
                wobs(sys.m(i).iKsi) = obs.m(i).wKsi;
            else
                wobs(sys.m(i).iKsi) = (0.02).^(-2);  % wKSi = 1/(1 + (0.02/pKsys(11)))^2 ;
                obs.m(i).wKsi = wobs(sys.m(i).iKsi);
            end
            if (isgood(obs.m(i).pKsi))
                yobs(sys.m(i).iKsi) = obs.m(i).pKsi;
            else
                yobs(sys.m(i).iKsi) = pKsi;
                obs.m(i).pKsi = pKsi;
            end
            if (isgood(obs.m(i).sioh4))
                yobs(sys.m(i).isioh4) = p(obs.m(i).sioh4);
            else
                yobs(sys.m(i).isioh4) = nan;
            end
            if (isgood(obs.m(i).wsioh4))
                wobs(sys.m(i).isioh4) = obs.m(i).wsioh4;
            else
                wobs(sys.m(i).isioh4) = nan;
            end
            if (isgood(obs.m(i).siooh3))
                yobs(sys.m(i).isiooh3) = p(obs.m(i).siooh3);
            else
                yobs(sys.m(i).isiooh3) = nan;
            end
            if (isgood(obs.m(i).wsiooh3))
                wobs(sys.m(i).isiooh3) = obs.m(i).wsiooh3;
            else
                wobs(sys.m(i).isiooh3) = nan;
            end
            
        end

        if (ismember('ammonia',sys.abr))
            if (isgood(obs.m(i).wKnh4))
                wobs(sys.m(i).iKnh4) = obs.m(i).wKnh4;
            else
                wobs(sys.m(i).iKnh4) = (0.00017).^(-2);  % wKnh4 = 1/(1 + (0.00017/pKsys(11)))^2 ;
                obs.m(i).wKnh4 = wobs(sys.m(i).iKnh4);
            end
            if (isgood(obs.m(i).pKnh4))
                yobs(sys.m(i).iKnh4) = obs.m(i).pKnh4;
            else
                yobs(sys.m(i).iKnh4) = pKnh4;
                obs.m(i).pKnh4 = pKnh4;
            end
            if (isgood(obs.m(i).nh4))
                yobs(sys.m(i).inh4) = p(obs.m(i).nh4);
            else
                yobs(sys.m(i).inh4) = nan;
            end
            if (isgood(obs.m(i).wnh4))
                wobs(sys.m(i).inh4) = obs.m(i).wnh4;
            else
                wobs(sys.m(i).inh4) = nan;
            end
            if (isgood(obs.m(i).nh3))
                yobs(sys.m(i).inh3) = p(obs.m(i).nh3);
            else
                yobs(sys.m(i).inh3) = nan;
            end
            if (isgood(obs.m(i).wnh3))
                wobs(sys.m(i).inh3) = obs.m(i).wnh3;
            else
                wobs(sys.m(i).inh3) = nan;
            end
        end

        if (ismember('sulfide',sys.abr))
            if (isgood(obs.m(i).wKh2s))
                wobs(sys.m(i).iKh2s) = obs.m(i).wKh2s;
            else
                wobs(sys.m(i).iKh2s) = (0.033).^(-2);  % wKh2s = 1/(1 + (0.033/pKsys(11)))^2 ;
                obs.m(i).wKh2s = wobs(sys.m(i).iKh2s);
            end
            if (isgood(obs.m(i).pKh2s))
                yobs(sys.m(i).iKh2s) = obs.m(i).pKh2s;
            else
                yobs(sys.m(i).iKh2s) = pKh2s;
                obs.m(i).pKh2s = pKh2s;
            end
            if (isgood(obs.m(i).h2s))
                yobs(sys.m(i).ih2s) = p(obs.m(i).h2s);
            else
                yobs(sys.m(i).ih2s) = nan;
            end
            if (isgood(obs.m(i).wh2s))
                wobs(sys.m(i).ih2s) = obs.m(i).wh2s;
            else
                wobs(sys.m(i).ih2s) = nan;
            end
            if (isgood(obs.m(i).hs))
                yobs(sys.m(i).ihs) = p(obs.m(i).hs);
            else
                yobs(sys.m(i).ihs) = nan;
            end
            if (isgood(obs.m(i).whs))
                wobs(sys.m(i).ihs) = obs.m(i).whs;
            else
                wobs(sys.m(i).ihs) = nan;
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
    ebar = @(j) [ q( z(j) - sigy(j) ), q( z(j) + sigy(j) ) ];
        
    est.sal = z(sys.isal);
    est.esal = sigy(sys.isal);
    est.TC = q(z(sys.iTC));
    est.TA = q(z(sys.iTA));
    est.eTC = ebar(sys.iTC);
    est.eTA = ebar(sys.iTA);
    if ismember('borate', sys.abr)
        est.TB = q(z(sys.iTB));
        est.eTB = ebar(sys.iTB);
    end
    if ismember('sulfate', sys.abr)
        est.TS = q(z(sys.iTS));
        est.eTS = ebar(sys.iTS);
    end
    if ismember('fluoride', sys.abr)
        est.TF = q(z(sys.iTF));
        est.eTF = ebar(sys.iTF);
    end
    if ismember('phosphate', sys.abr)
        est.TP = q(z(sys.iTP));
        est.eTP = ebar(sys.iTP);
    end
    if ismember('silicate', sys.abr)
        est.TSi = q(z(sys.iTSi));
        est.eTSi = ebar(sys.iTSi);
    end
    if ismember('ammonia', sys.abr)
        est.TNH3 = q(z(sys.iTNH3));
        est.eTNH3 = ebar(sys.iTNH3);
    end
    if ismember('sulfide', sys.abr)
        est.TH2S = q(z(sys.iTH2S));
        est.eTH2S = ebar(sys.iTH2S);
    end
    
    for i = 1:nTP
        % thermodynamic state
        est.m(i).T = z(sys.m(i).iT);
        est.m(i).P = z(sys.m(i).iP);
        % thermodynamic state errorbar
        est.m(i).eT = sigy(sys.m(i).iT);
        est.m(i).eP = sigy(sys.m(i).iP);
        % chemical equilibrium
        est.m(i).ph    = z(sys.m(i).iph);
        est.m(i).fco2  = q(z(sys.m(i).ifco2));
        est.m(i).pco2  = q(z(sys.m(i).ipco2));
        est.m(i).pco2st = z(sys.m(i).ico2st);
        est.m(i).phco3  = z(sys.m(i).ihco3);
        est.m(i).pco3   = z(sys.m(i).ico3);
        % chemical equilibrium errorbar
        est.m(i).eph    = sigy(sys.m(i).iph);
        est.m(i).efco2  = ebar(sys.m(i).ifco2);
        est.m(i).epco2  = ebar(sys.m(i).ipco2);
        est.m(i).epco2st = sigy(sys.m(i).ico2st);
        est.m(i).ephco3  = sigy(sys.m(i).ihco3);
        est.m(i).epco3   = sigy(sys.m(i).ico3);        
        % equilibrium pK's
        est.m(i).pp2f = z(sys.m(i).ip2f);
        est.m(i).pK0  = z(sys.m(i).iK0);
        est.m(i).pK1  = z(sys.m(i).iK1);
        est.m(i).pK2  = z(sys.m(i).iK2);
        % pK errorbars
        est.m(i).epp2f = sigy(sys.m(i).ip2f);
        est.m(i).epK0  = sigy(sys.m(i).iK0);
        est.m(i).epK1  = sigy(sys.m(i).iK1);
        est.m(i).epK2  = sigy(sys.m(i).iK2);
        if ismember('borate', sys.abr)
            %chemical equilibrium
            est.m(i).pboh4 = z(sys.m(i).iboh4);
            est.m(i).pboh3 = z(sys.m(i).iboh3);
            % chemical equilibrium errorbar
            est.m(i).epboh4 = sigy(sys.m(i).iboh4);
            est.m(i).epboh3 = sigy(sys.m(i).iboh3);
            % equilibrium pK
            est.m(i).pKb = z(sys.m(i).iKb);
            % pK errorbar
            est.m(i).epKb = sigy(sys.m(i).iKb);
        end
        if ismember('water', sys.abr)
            % chemical equilibrium
            est.m(i).poh = z(sys.m(i).ioh); 
            % chemical equilibrium errorbar
            est.m(i).epoh = sigy(sys.m(i).ioh);
            % equilibrium pK
            est.m(i).pKw = z(sys.m(i).iKw);
            % pK errorbar
            est.m(i).epKw = sigy(sys.m(i).iKw); 
        end
        if ismember('sulfate', sys.abr)
            % chemical equilibrium
            est.m(i).phf = z(sys.m(i).iphf);
            est.m(i).pso4 = z(sys.m(i).iso4);
            est.m(i).phso4 = z(sys.m(i).ihso4);
            % chemical equilibrium errorbar
            est.m(i).ephf = sigy(sys.m(i).iphf);
            est.m(i).epso4 = sigy(sys.m(i).iso4);
            est.m(i).ephso4 = sigy(sys.m(i).ihso4);
            % equilibrium pK
            est.m(i).pKs = z(sys.m(i).iKs);
            % pK errorbar
            est.m(i).epKs = sigy(sys.m(i).iKs);
        end
        if ismember('fluoride', sys.abr)
            % chemical equilibrium
            est.m(i).pF = z(sys.m(i).iF);
            est.m(i).pHF = z(sys.m(i).iHF);
            % chemical equilibrium errorbar
            est.m(i).epF = sigy(sys.m(i).iF);
            est.m(i).epHF = sigy(sys.m(i).iHF);
            % equilibrium pK
            est.m(i).pKf = z(sys.m(i).iKf);
            % pK errorbar
            est.m(i).epKf = sigy(sys.m(i).iKf);
        end
        if ismember('phosphate', sys.abr)
            % chemical equilibrium
            est.m(i).ppo4   = z(sys.m(i).ipo4);
            est.m(i).phpo4  = z(sys.m(i).ihpo4);
            est.m(i).ph2po4 = z(sys.m(i).ih2po4);
            est.m(i).ph3po4 = z(sys.m(i).ih3po4);
            % chemical equilibrium errorbar
            est.m(i).eppo4   = sigy(sys.m(i).ipo4);
            est.m(i).ephpo4  = sigy(sys.m(i).ihpo4);
            est.m(i).eph2po4 = sigy(sys.m(i).ih2po4);
            est.m(i).eph3po4 = sigy(sys.m(i).ih3po4);
            % equilibrium pK
            est.m(i).pK1p = z(sys.m(i).iK1p);
            est.m(i).pK2p = z(sys.m(i).iK2p);
            est.m(i).pK3p = z(sys.m(i).iK3p);
            % pK errorbar
            est.m(i).epK1p = sigy(sys.m(i).iK1p);
            est.m(i).epK2p = sigy(sys.m(i).iK2p);
            est.m(i).epK3p = sigy(sys.m(i).iK3p);
        end
        if ismember('silicate', sys.abr)
            % chemical equilibrium
            est.m(i).psioh4 = z(sys.m(i).isioh4);
            est.m(i).ppsiooh3 = z(sys.m(i).isiooh3);
            % chemical equilibrium errorbar
            est.m(i).epsioh4 = sigy(sys.m(i).isioh4);
            est.m(i).epsiooh3 = sigy(sys.m(i).isiooh3);
            % equilibrium pK
            est.m(i).pKsi = z(sys.m(i).iKsi);
            % pK errorbar
            est.m(i).epKsi = sigy(sys.m(i).iKsi);
        end
        if ismember('ammonia', sys.abr)
            % chemical equilibrium
            est.m(i).pnh3 = z(sys.m(i).inh3);
            est.m(i).pnh4 = z(sys.m(i).inh4);
            % chemical equilibrium errorbar
            est.m(i).epnh3 = sigy(sys.m(i).inh3);
            est.m(i).epnh4= sigy(sys.m(i).inh4);
            % equilibrium pK
            est.m(i).pKnh4 = z(sys.m(i).iKnh4);
            % pK errorbar
            est.m(i).epKnh4 = sigy(sys.m(i).iKnh4);
        end
        if ismember('sulfide', sys.abr)
            % chemical equilibrium
            est.m(i).phs = z(sys.m(i).ihs);
            est.m(i).ph2s = z(sys.m(i).ih2s);
            % chemical equilibrium errorbar
            est.m(i).ephs = sigy(sys.m(i).ihs);
            est.m(i).eph2s = sigy(sys.m(i).ih2s);
            % equilibrium pK
            est.m(i).pKh2s = z(sys.m(i).ipKh2s);
            % pK errorbar
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

