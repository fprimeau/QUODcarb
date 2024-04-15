function PrintCSV(est,obs,iflag,opt)
% subfunction of QUODcarb to print results to a CSV
% uses filename opt.fname 
%               opt.fname = 'QUODcarb_output.csv'(default)

    if opt.printcsv == 1
        nD = length(obs);
        for i = 1:nD
            if i == 1
                fname = fopen(opt.fname,'w');
                % make column headers
                make_headers(est(i),opt,fname); 
            end
            % fill one row with data
            parse_CSV(est(i),obs(i),iflag(i),opt,fname);
        end 
    end
end

function make_headers(est,opt,fid)
% Print the column headers  
    nTP = length(est.tp);
    % row 1
    fprintf(fid, '%s, ', '  ');
    fprintf(fid, 'est = output, ');
    fprintf(fid, '%s, ', '  ');

    fn = fieldnames(est);
    fnl = length(fn)-1;
    for i = 5:6:fnl % temperature independent totals
        fprintf(fid, '%s, ', '  '); % est
        fprintf(fid, '%s, ', '  '); % est.e
    end
    for j = 1:nTP
        fnj = fieldnames(est.tp);
        for i = 1:4:8 % T, P
            fprintf(fid, 'tp(%i), ', j); % est
            fprintf(fid, '%s, ', '  ');  % est.e
        end
        for i = 9:6:39 % fco2, pco2, hco3, co3, co2*, ph
            fprintf(fid, 'tp(%i), ', j); % est
            fprintf(fid, '%s, ', '  ');  % est.e
        end
        for i = 45:6:75 % ph_free, ph_tot, ph_sws, ph_nbs, fH, p2f
            fprintf(fid, 'tp(%i), ', j); % est
            fprintf(fid, '%s, ', '  ');  % est.e
        end
        for i = 79:6:length(fnj) % all the rest
            fprintf(fid, 'tp(%i), ', j); % est
            fprintf(fid, '%s, ', '  ');  % est.e
        end
    end

    fprintf(fid,'\n'); % finish first row

    % row 2
    fprintf(fid, '%s, ','iflag');

    fprintf(fid,'est.sal, '); % sal is first
    fprintf(fid,'est.esal, '); % esal

    for i = 5:6:fnl % temperature independent totals
        fprintf(fid,'est.%s, ',fn{i} );
        fprintf(fid,'est.e%s, ',fn{i} ); % e
    end
    for j = 1:nTP
        fnj = fieldnames(est.tp(j));
        for i = 1:4:8 % T, P
            fprintf(fid,'est.%s, ',  fnj{i});
            fprintf(fid,'est.e%s, ', fnj{i}); % e = error
        end
        for i = 9:6:39 % fco2, pco2, hco3, co3, co2*, ph
            fprintf(fid,'est.%s, ',  fnj{i});
            fprintf(fid,'est.e%s, ', fnj{i}); % e = error
        end
        for i = 45:6:75 % ph_free, ph_tot, ph_sws, ph_nbs, fH, p2f
            fprintf(fid,'est.%s, ',  fnj{i});
            fprintf(fid,'est.e%s, ', fnj{i}); % e = error
        end
        for i = 79:6:288 % all the rest (except Revelle if on)
            fprintf(fid,'est.%s, ',  fnj{i});
            fprintf(fid,'est.e%s, ', fnj{i}); % e = error
        end
        if opt.Revelle == 1
            fprintf(fid,'est.%s, ', fnj{end-1}); % Revelle
            fprintf(fid,'est.%s, ', fnj{end}); % dpfco2dpTA
        end    
    end

    fprintf(fid,'\n'); % finish second row

    % row 3
    fprintf(fid,'%s, ','(0=good)');
       
    fprintf(fid,'%s, ', '(PSU)'); % est.sal
    fprintf(fid,'%s, ', '(PSU)'); % est.esal
    
    for i = 3:6:fnl % temperature independent totals
        fprintf(fid,'%s, ', '(umol/kg)');
        fprintf(fid,'%s, ', '(umol/kg)');
    end
    for j = 1:nTP
        fprintf(fid, '%s, ', 'deg C'); % est.T
        fprintf(fid, '%s, ', 'deg C');
        fprintf(fid, '%s, ', 'dbar'); % est.P
        fprintf(fid, '%s, ', 'dbar');
        for i = 9:6:15 % fco2, pco2
            fprintf(fid, '%s, ', '(uatm)');
            fprintf(fid, '%s, ', '(uatm)');
        end
        for i = 21:6:33 % hco3, co3, co2*;
            fprintf(fid, '%s, ', '(umol/kg)');
            fprintf(fid, '%s, ', '(umol/kg)');
        end
        fprintf(fid, '%s, ', '(p units)'); % ph, log10 unitless
        fprintf(fid, '%s, ', '  ');
        for i = 45:6:63 % ph_free, ph_tot, ph_sws, ph_nbs
            fprintf(fid, '%s, ', '(p units)');      % log10 unitless
            fprintf(fid, '%s, ', '  ');
        end
        for i = 69:6:75 % fH, p2f aka FugFac
            fprintf(fid, '%s, ', '(unitless)'); 
            fprintf(fid, '%s, ', '  ');
        end
        for i = 79:6:91 % pK0, pK1, pK2
            fprintf(fid, '%s, ', '(p units)');      % log10 unitless
            fprintf(fid, '%s, ', '  ');
        end
        fprintf(fid, '%s, ', '(umol/kg)'); % oh (97)
        fprintf(fid, '%s, ', '(umol/kg)');
        fprintf(fid, '%s, ', '(p units)'); % pKw (103)
        fprintf(fid, '%s, ', '  ');
        for i = 1:3 % (boh3, boh4, pKb)(so4, hso4, pKs)(F, HF, pKf)
            fprintf(fid, '%s, ', '(umol/kg)'); % est 1st
            fprintf(fid, '%s, ', '(umol/kg)');
            fprintf(fid, '%s, ', '(umol/kg)'); % est 2nd
            fprintf(fid, '%s, ', '(umol/kg)');
            fprintf(fid, '%s, ', '(p units)'); % est.pK
            fprintf(fid, '%s, ', '  ');
        end
        for i = 1:4 % po4, hpo4, h2po4, h3po4
            fprintf(fid, '%s, ', '(umol/kg)');
            fprintf(fid, '%s, ', '(umol/kg)');
        end
        for i = 1:3 % pKp1, pKp2, pKp3
            fprintf(fid, '%s, ', '(p units)');
            fprintf(fid, '%s, ', '  ');
        end
        for i = 1:3 % (sioh4, siooh3, pKsi)(nh3, nh4, pKnh4)(HS, H2S, pKh2s)
            fprintf(fid, '%s, ', '(umol/kg)'); % est 1st
            fprintf(fid, '%s, ', '(umol/kg)');
            fprintf(fid, '%s, ', '(umol/kg)'); % est 2nd
            fprintf(fid, '%s, ', '(umol/kg)');
            fprintf(fid, '%s, ', '(p units)'); % est.pK
            fprintf(fid, '%s, ', '  ');
        end
        fprintf(fid, '%s, ', '(umol/kg)'); % ca
        fprintf(fid, '%s, ', '(umol/kg)');
        for i = 1:2 % (Omega_Ar, pKar)(Omega_Ca, pKca)
            fprintf(fid, '%s, ', '(umol/kg)'); % est 1st
            fprintf(fid, '%s, ', '(umol/kg)');
            fprintf(fid, '%s, ', '(p units)'); % est.pK
            fprintf(fid, '%s, ', '  ');
        end
        if opt.Revelle == 1
            fprintf(fid, '%s, ', '(unitless)'); % Revelle
            fprintf(fid, '%s, ', '  ');         % dpfco2dpTA
        end
    end

    fprintf(fid,'\n'); % finish third row

end

function parse_CSV(varargin)
% Print the output
    est     = varargin{1}; % posterior marginal precision
    obs     = varargin{2}; % measurement
    iflag   = varargin{3}; % Newton solver convergence flag: 0 converged, 1 not converged
    opt     = varargin{4}; % options
    fid   = varargin{5}; % filename with fopen command

    nTP     = length(est.tp);

    fprintf(fid,'%i, ',iflag);

    % salinity
    fprintf(fid,'%f, ', est.sal);
    fprintf(fid,'%f, ', est.esal);
    % TC (DIC)
    fprintf(fid,'%f, ', est.TC);
    fprintf(fid,'%f, ', est.eTC);
    % TA
    fprintf(fid,'%f, ', est.TA);
    fprintf(fid,'%f, ', est.eTA);
    % TB borate
    fprintf(fid,'%f, ', est.TB);
    fprintf(fid,'%f, ', est.eTB);
    % TS sulfate
    fprintf(fid,'%f, ', est.TS);
    fprintf(fid,'%f, ', est.eTS);
    % TF fluoride
    fprintf(fid,'%f, ', est.TF);
    fprintf(fid,'%f, ', est.eTF);
    % TP phosphate
    fprintf(fid,'%f, ', est.TP);
    fprintf(fid,'%f, ', est.eTP);
    % TSi silicate
    fprintf(fid,'%f, ', est.TSi);
    fprintf(fid,'%f, ', est.eTSi);
    % TNH4 nitrate
    fprintf(fid,'%f, ', est.TNH4);
    fprintf(fid,'%f, ', est.eTNH4);
    % TH2S sulfide
    fprintf(fid,'%f, ', est.TH2S);
    fprintf(fid,'%f, ', est.eTH2S);
    % TCa calcium solubility
    fprintf(fid,'%f, ', est.TCa);
    fprintf(fid,'%f, ', est.eTCa);

    for j = 1:nTP
        % Temp
        fprintf(fid,'%f, ', est.tp(j).T);
        fprintf(fid,'%f, ', est.tp(j).eT);
        % pres
        fprintf(fid,'%f, ', est.tp(j).P);
        fprintf(fid,'%f, ', est.tp(j).eP);

        % fco2
        fprintf(fid,'%f, ', est.tp(j).fco2);
        fprintf(fid,'%f, ', est.tp(j).efco2);
        % pco2
        fprintf(fid,'%f, ', est.tp(j).pco2);
        fprintf(fid,'%f, ', est.tp(j).epco2);
        % hco3
        fprintf(fid,'%f, ', est.tp(j).hco3);
        fprintf(fid,'%f, ', est.tp(j).ehco3);
        % co3
        fprintf(fid,'%f, ', est.tp(j).co3);
        fprintf(fid,'%f, ', est.tp(j).eco3);
        % co2* (co2st)
        fprintf(fid,'%f, ', est.tp(j).co2st);
        fprintf(fid,'%f, ', est.tp(j).eco2st);

        % ph
        fprintf(fid,'%f, ', est.tp(j).ph);
        fprintf(fid,'%f, ', est.tp(j).eph);
        % ph_free
        fprintf(fid,'%f, ', est.tp(j).ph_free);
        fprintf(fid,'%f, ', est.tp(j).eph_free);
        % ph_tot
        fprintf(fid,'%f, ', est.tp(j).ph_tot);
        fprintf(fid,'%f, ', est.tp(j).eph);
        % ph_sws
        fprintf(fid,'%f, ', est.tp(j).ph_sws);
        fprintf(fid,'%f, ', est.tp(j).eph);
        % ph_nbs
        fprintf(fid,'%f, ', est.tp(j).ph_nbs);
        fprintf(fid,'%f, ', est.tp(j).eph);

         % fH = activity coefficient
        fprintf(fid,'%f, ', est.tp(j).fH);
        fprintf(fid,'%f, ', est.tp(j).efH);
        % pp2f
        fprintf(fid,'%f, ', est.tp(j).pp2f);
        fprintf(fid,'%f, ', est.tp(j).epp2f);

        % pK0
        fprintf(fid,'%f, ', est.tp(j).pK0);
        fprintf(fid,'%f, ', est.tp(j).epK0);
        % pK1
        fprintf(fid,'%f, ', est.tp(j).pK1);
        fprintf(fid,'%f, ', est.tp(j).epK1);
        % pK2
        fprintf(fid,'%f, ', est.tp(j).pK2);
        fprintf(fid,'%f, ', est.tp(j).epK2);

        % oh
        fprintf(fid,'%f, ', est.tp(j).oh);
        fprintf(fid,'%f, ', est.tp(j).eoh);
        % pKw = [h][oh]
        fprintf(fid,'%f, ', est.tp(j).pKw);
        fprintf(fid,'%f, ', est.tp(j).epKw);

        % boh4
        fprintf(fid,'%f, ', est.tp(j).boh4);
        fprintf(fid,'%f, ', est.tp(j).eboh4);
        % boh3
        fprintf(fid,'%f, ', est.tp(j).boh3);
        fprintf(fid,'%f, ', est.tp(j).eboh3);
        % pKb = [h][boh4]/[boh3]
        fprintf(fid,'%f, ', est.tp(j).pKb);
        fprintf(fid,'%f, ', est.tp(j).epKb);

        % so4
        fprintf(fid,'%f, ', est.tp(j).so4);
        fprintf(fid,'%f, ', est.tp(j).eso4);
        % hso4
        fprintf(fid,'%f, ', est.tp(j).hso4);
        fprintf(fid,'%f, ', est.tp(j).ehso4);
        % pKs  = [hf][so4]/[hso4]
        fprintf(fid,'%f, ', est.tp(j).pKs);
        fprintf(fid,'%f, ', est.tp(j).epKs);

        % [F]
        fprintf(fid,'%f, ', est.tp(j).F);
        fprintf(fid,'%f, ', est.tp(j).eF);
        % [HF] hydrogen fluoride
        fprintf(fid,'%f, ', est.tp(j).HF);
        fprintf(fid,'%f, ', est.tp(j).eHF);
        % pKf = [h][F]/[HF]
        fprintf(fid,'%f, ', est.tp(j).pKf);
        fprintf(fid,'%f, ', est.tp(j).epKf);

        % po4
        fprintf(fid,'%f, ', est.tp(j).po4);
        fprintf(fid,'%f, ', est.tp(j).epo4);
        % hpo4
        fprintf(fid,'%f, ', est.tp(j).hpo4);
        fprintf(fid,'%f, ', est.tp(j).ehpo4);
        % h2po4
        fprintf(fid,'%f, ', est.tp(j).h2po4);
        fprintf(fid,'%f, ', est.tp(j).eh2po4);
        % h3po4
        fprintf(fid,'%f, ', est.tp(j).h3po4);
        fprintf(fid,'%f, ', est.tp(j).eh3po4);
        % pKp1 = [h][h2po4]/[h3po4]
        fprintf(fid,'%f, ', est.tp(j).pKp1);
        fprintf(fid,'%f, ', est.tp(j).epKp1);
        % pKp2 = [h][hpo4]/[h2po4]
        fprintf(fid,'%f, ', est.tp(j).pKp2);
        fprintf(fid,'%f, ', est.tp(j).epKp2);
        % pKp3 = [h][po4]/[hpo4]
        fprintf(fid,'%f, ', est.tp(j).pKp3);
        fprintf(fid,'%f, ', est.tp(j).epKp3);

        % sioh4
        fprintf(fid,'%f, ', est.tp(j).sioh4);
        fprintf(fid,'%f, ', est.tp(j).esioh4);
        % siooh3
        fprintf(fid,'%f, ', est.tp(j).siooh3);
        fprintf(fid,'%f, ', est.tp(j).esiooh3);
        % pKSi = [h][siooh3]/[sioh4]
        fprintf(fid,'%f, ', est.tp(j).pKsi);
        fprintf(fid,'%f, ', est.tp(j).epKsi);

        % nh3
        fprintf(fid,'%f, ', est.tp(j).nh3);
        fprintf(fid,'%f, ', est.tp(j).enh3);
        % nh4
        fprintf(fid,'%f, ', est.tp(j).nh4);
        fprintf(fid,'%f, ', est.tp(j).enh4);
        % pKnh4 = [h][nh3]/[nh4]
        fprintf(fid,'%f, ', est.tp(j).pKnh4);
        fprintf(fid,'%f, ', est.tp(j).epKnh4);

        % hs
        fprintf(fid,'%f, ', est.tp(j).HS);
        fprintf(fid,'%f, ', est.tp(j).eHS);
        % h2s
        fprintf(fid,'%f, ', est.tp(j).H2S);
        fprintf(fid,'%f, ', est.tp(j).eH2S);
        % pKh2s = [h][hs]/[h2s]
        fprintf(fid,'%f, ', est.tp(j).pKh2s);
        fprintf(fid,'%f, ', est.tp(j).epKh2s);

        % ca
        fprintf(fid,'%f, ', est.tp(j).ca);
        fprintf(fid,'%f, ', est.tp(j).eca);
        % OmegaAr
        fprintf(fid,'%f, ', est.tp(j).OmegaAr);
        fprintf(fid,'%f, ', est.tp(j).eOmegaAr);
        % OmegaCa
        fprintf(fid,'%f, ', est.tp(j).OmegaCa);
        fprintf(fid,'%f, ', est.tp(j).eOmegaCa);
        % pKar = [ca][co3]/[omegaAr]
        fprintf(fid,'%f, ', est.tp(j).pKar);
        fprintf(fid,'%f, ', est.tp(j).epKar);
        % pKca = [ca][co3]/[omegaCa]
        fprintf(fid,'%f, ', est.tp(j).pKca); 
        fprintf(fid,'%f, ', est.tp(j).epKca);

        if opt.Revelle == 1
            fprintf(fid,'%f, ', est.tp(j).Revelle);
            fprintf(fid,'%f, ', est.tp(j).dpfco2dpTA);
        end
    end

fprintf(fid,'\n'); % end row of data
out = [];

end % function

