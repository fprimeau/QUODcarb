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

function make_headers(est,opt,fname)
% Print the column headers
    
    nTP = length(est.tp);
    % row 1
    fprintf(fname, '%s, ', '  ');
    fprintf(fname,'obs = input, ');
    fprintf(fname,'NaN = not input, ');
    fprintf(fname,'est = output, ');
    fprintf(fname, '%s, ', '  ');

    fn = fieldnames(est);
    fnl = length(fn)-1;
    for i = 5:6:fnl % temperature independent totals
        fprintf(fname, '%s, ', '  '); % obs
        fprintf(fname, '%s, ', '  '); % obs.e
        fprintf(fname, '%s, ', '  '); % est
        fprintf(fname, '%s, ', '  '); % est.e
    end
    for j = 1:nTP
        fnj = fieldnames(est.tp);
        for i = 1:4:8 % T, P
            fprintf(fname, 'tp(%i), ', j); % obs
            fprintf(fname, '%s, ', '  ');  % obs.e
            fprintf(fname, 'tp(%i), ', j); % est
            fprintf(fname, '%s, ', '  ');  % est.e
        end
        for i = 9:6:39 % fco2, pco2, hco3, co3, co2*, ph
            fprintf(fname, 'tp(%i), ', j); % obs
            fprintf(fname, '%s, ', '  ');  % obs.e
            fprintf(fname, 'tp(%i), ', j); % est
            fprintf(fname, '%s, ', '  ');  % est.e
        end
        for i = 45:6:75 % ph_free, ph_tot, ph_sws, ph_nbs, fH, p2f
            fprintf(fname, 'tp(%i), ', j); % est
            fprintf(fname, '%s, ', '  ');  % est.e
        end
        for i = 79:6:length(fnj) % all the rest
            fprintf(fname, 'tp(%i), ', j); % est
            fprintf(fname, '%s, ', '  ');  % est.e
        end
    end

    fprintf(fname,'\n'); % finish first row

    % row 2
    fprintf(fname, '%s, ','iflag');

    fprintf(fname,'obs.sal, '); % sal is first
    fprintf(fname,'obs.esal, '); % esal
    fprintf(fname,'est.sal, ');
    fprintf(fname,'est.esal, ');

    for i = 5:6:fnl % temperature independent totals
        fprintf(fname,'obs.%s, ',fn{i} );
        fprintf(fname,'obs.e%s, ',fn{i} ); % e
        fprintf(fname,'est.%s, ',fn{i} );
        fprintf(fname,'est.e%s, ',fn{i} ); % e
    end
    for j = 1:nTP
        fnj = fieldnames(est.tp(j));
        for i = 1:4:8 % T, P
            fprintf(fname,'obs.%s, ',  fnj{i});
            fprintf(fname,'obs.e%s, ', fnj{i});
            fprintf(fname,'est.%s, ',  fnj{i});
            fprintf(fname,'est.e%s, ', fnj{i}); % e = error
        end
        for i = 9:6:39 % fco2, pco2, hco3, co3, co2*, ph
            fprintf(fname,'obs.%s, ',  fnj{i});
            fprintf(fname,'obs.e%s, ', fnj{i});
            fprintf(fname,'est.%s, ',  fnj{i});
            fprintf(fname,'est.e%s, ', fnj{i}); % e = error
        end
        for i = 45:6:75 % ph_free, ph_tot, ph_sws, ph_nbs, fH, p2f
            fprintf(fname,'est.%s, ',  fnj{i});
            fprintf(fname,'est.e%s, ', fnj{i}); % e = error
        end
        for i = 79:6:288 % all the rest (except Revelle if on)
            fprintf(fname,'est.%s, ',  fnj{i});
            fprintf(fname,'est.e%s, ', fnj{i}); % e = error
        end
        if opt.Revelle == 1
            fprintf(fname,'est.%s, ', fnj{end-1}); % Revelle
            fprintf(fname,'est.%s, ', fnj{end}); % dpfco2dpTA
        end    
    end

    fprintf(fname,'\n'); % finish second row

    % row 3
    fprintf(fname,'%s, ','(0=good)');
       
    fprintf(fname,'%s, ', '(PSU)'); % obs.sal
    fprintf(fname,'%s, ', '(PSU)'); % obs.esal
    fprintf(fname,'%s, ', '(PSU)'); % est.sal
    fprintf(fname,'%s, ', '(PSU)'); % est.esal
    
    for i = 3:6:fnl % temperature independent totals
        fprintf(fname,'%s, ', '(umol/kg)');
        fprintf(fname,'%s, ', '(umol/kg)');
        fprintf(fname,'%s, ', '(umol/kg)');
        fprintf(fname,'%s, ', '(umol/kg)');
    end
    for j = 1:nTP
        fprintf(fname, '%s, ', 'deg C'); % obs.T
        fprintf(fname, '%s, ', 'deg C');
        fprintf(fname, '%s, ', 'deg C'); % est.T
        fprintf(fname, '%s, ', 'deg C');
        fprintf(fname, '%s, ', 'dbar'); % obs.P
        fprintf(fname, '%s, ', 'dbar');
        fprintf(fname, '%s, ', 'dbar'); % est.P
        fprintf(fname, '%s, ', 'dbar');
        for i = 9:6:15 % fco2, pco2
            fprintf(fname, '%s, ', '(uatm)');
            fprintf(fname, '%s, ', '(uatm)');
            fprintf(fname, '%s, ', '(uatm)');
            fprintf(fname, '%s, ', '(uatm)');
        end
        for i = 21:6:33 % hco3, co3, co2*
            fprintf(fname, '%s, ', '(umol/kg)');
            fprintf(fname, '%s, ', '(umol/kg)');
            fprintf(fname, '%s, ', '(umol/kg)');
            fprintf(fname, '%s, ', '(umol/kg)');
        end
        fprintf(fname, '%s, ', '(p units)'); % ph, log10 unitless
        fprintf(fname, '%s, ', '  ');
        fprintf(fname, '%s, ', '(p units)');
        fprintf(fname, '%s, ', '  ');
        for i = 45:6:63 % ph_free, ph_tot, ph_sws, ph_nbs
            fprintf(fname, '%s, ', '(p units)');      % log10 unitless
            fprintf(fname, '%s, ', '  ');
        end
        for i = 69:6:75 % fH, p2f aka FugFac
            fprintf(fname, '%s, ', '(unitless)'); 
            fprintf(fname, '%s, ', '  ');
        end
        for i = 79:6:91 % pK0, pK1, pK2
            fprintf(fname, '%s, ', '(p units)');      % log10 unitless
            fprintf(fname, '%s, ', '  ');
        end
        fprintf(fname, '%s, ', '(umol/kg)'); % oh (97)
        fprintf(fname, '%s, ', '(umol/kg)');
        fprintf(fname, '%s, ', '(p units)'); % pKw (103)
        fprintf(fname, '%s, ', '  ');
        for i = 1:3 % (boh3, boh4, pKb)(so4, hso4, pKs)(F, HF, pKf)
            fprintf(fname, '%s, ', '(umol/kg)'); % est 1st
            fprintf(fname, '%s, ', '(umol/kg)');
            fprintf(fname, '%s, ', '(umol/kg)'); % est 2nd
            fprintf(fname, '%s, ', '(umol/kg)');
            fprintf(fname, '%s, ', '(p units)'); % est.pK
            fprintf(fname, '%s, ', '  ');
        end
        for i = 1:4 % po4, hpo4, h2po4, h3po4
            fprintf(fname, '%s, ', '(umol/kg)');
            fprintf(fname, '%s, ', '(umol/kg)');
        end
        for i = 1:3 % pKp1, pKp2, pKp3
            fprintf(fname, '%s, ', '(p units)');
            fprintf(fname, '%s, ', '  ');
        end
        for i = 1:3 % (sioh4, siooh3, pKsi)(nh3, nh4, pKnh4)(HS, H2S, pKh2s)
            fprintf(fname, '%s, ', '(umol/kg)'); % est 1st
            fprintf(fname, '%s, ', '(umol/kg)');
            fprintf(fname, '%s, ', '(umol/kg)'); % est 2nd
            fprintf(fname, '%s, ', '(umol/kg)');
            fprintf(fname, '%s, ', '(p units)'); % est.pK
            fprintf(fname, '%s, ', '  ');
        end
        fprintf(fname, '%s, ', '(umol/kg)'); % ca
        fprintf(fname, '%s, ', '(umol/kg)');
        for i = 1:2 % (Omega_Ar, pKar)(Omega_Ca, pKca)
            fprintf(fname, '%s, ', '(umol/kg)'); % est 1st
            fprintf(fname, '%s, ', '(umol/kg)');
            fprintf(fname, '%s, ', '(p units)'); % est.pK
            fprintf(fname, '%s, ', '  ');
        end
        if opt.Revelle == 1
            fprintf(fname, '%s, ', '(unitless)'); % Revelle
            fprintf(fname, '%s, ', '  ');
            % fprintf(fname, '%s, ', '(unitless)'); % dpfco2dpTA
            % fprintf(fname, '%s, ', '  ');
        end
    end

    fprintf(fname,'\n'); % finish third row

end


function parse_CSV(varargin)
% Print the output

    est     = varargin{1}; % posterior marginal precision
    obs     = varargin{2}; % measurement
    iflag   = varargin{3}; % Newton solver convergence flag: 0 converged, 1 not converged
    opt     = varargin{4}; % options
    fname   = varargin{5}; % filename with fopen command

    nTP     = length(est.tp);

    fprintf(fname,'%i, ',iflag);

    % salinity
    fprintf(fname,'%f, ', obs.sal);
    fprintf(fname,'%f, ', obs.esal);
    fprintf(fname,'%f, ', est.sal);
    fprintf(fname,'%f, ', est.esal);

    % TC (DIC)
    fprintf(fname,'%f, ', obs.TC);
    fprintf(fname,'%f, ', obs.eTC);
    fprintf(fname,'%f, ', est.TC);
    fprintf(fname,'%f, ', est.eTC);

    % TA
    fprintf(fname,'%f, ', obs.TA);
    fprintf(fname,'%f, ', obs.eTA);
    fprintf(fname,'%f, ', est.TA);
    fprintf(fname,'%f, ', est.eTA);

    % TB borate
    fprintf(fname,'%f, ', obs.TB);
    fprintf(fname,'%f, ', obs.eTB);
    fprintf(fname,'%f, ', est.TB);
    fprintf(fname,'%f, ', est.eTB);

    % TS sulfate
    fprintf(fname,'%f, ', obs.TS);
    fprintf(fname,'%f, ', obs.eTS);
    fprintf(fname,'%f, ', est.TS);
    fprintf(fname,'%f, ', est.eTS);

    % TF fluoride
    fprintf(fname,'%f, ', obs.TF);
    fprintf(fname,'%f, ', obs.eTF);
    fprintf(fname,'%f, ', est.TF);
    fprintf(fname,'%f, ', est.eTF);

    % TP phosphate
    fprintf(fname,'%f, ', obs.TP);
    fprintf(fname,'%f, ', obs.eTP);
    fprintf(fname,'%f, ', est.TP);
    fprintf(fname,'%f, ', est.eTP);

    % TSi silicate
    fprintf(fname,'%f, ', obs.TSi);
    fprintf(fname,'%f, ', obs.eTSi);
    fprintf(fname,'%f, ', est.TSi);
    fprintf(fname,'%f, ', est.eTSi);

    % TNH4 nitrate
    fprintf(fname,'%f, ', obs.TNH4);
    fprintf(fname,'%f, ', obs.eTNH4);
    fprintf(fname,'%f, ', est.TNH4);
    fprintf(fname,'%f, ', est.eTNH4);

    % TH2S sulfide
    fprintf(fname,'%f, ', obs.TH2S);
    fprintf(fname,'%f, ', obs.eTH2S);
    fprintf(fname,'%f, ', est.TH2S);
    fprintf(fname,'%f, ', est.eTH2S);

    % TCa calcium solubility
    fprintf(fname,'%f, ', obs.TCa);
    fprintf(fname,'%f, ', obs.eTCa);
    fprintf(fname,'%f, ', est.TCa);
    fprintf(fname,'%f, ', est.eTCa);

    for j = 1:nTP
        % Temp
        fprintf(fname,'%f, ', obs.tp(j).T);
        fprintf(fname,'%f, ', obs.tp(j).eT);
        fprintf(fname,'%f, ', est.tp(j).T);
        fprintf(fname,'%f, ', est.tp(j).eT);
        % pres
        fprintf(fname,'%f, ', obs.tp(j).P);
        fprintf(fname,'%f, ', obs.tp(j).eP);
        fprintf(fname,'%f, ', est.tp(j).P);
        fprintf(fname,'%f, ', est.tp(j).eP);

        % fco2
        fprintf(fname,'%f, ', obs.tp(j).fco2);
        fprintf(fname,'%f, ', obs.tp(j).efco2);
        fprintf(fname,'%f, ', est.tp(j).fco2);
        fprintf(fname,'%f, ', est.tp(j).efco2);
        % pco2
        fprintf(fname,'%f, ', obs.tp(j).pco2);
        fprintf(fname,'%f, ', obs.tp(j).epco2);
        fprintf(fname,'%f, ', est.tp(j).pco2);
        fprintf(fname,'%f, ', est.tp(j).epco2);
        % hco3
        fprintf(fname,'%f, ', obs.tp(j).hco3);
        fprintf(fname,'%f, ', obs.tp(j).ehco3);
        fprintf(fname,'%f, ', est.tp(j).hco3);
        fprintf(fname,'%f, ', est.tp(j).ehco3);
        % co3
        fprintf(fname,'%f, ', obs.tp(j).co3);
        fprintf(fname,'%f, ', obs.tp(j).eco3);
        fprintf(fname,'%f, ', est.tp(j).co3);
        fprintf(fname,'%f, ', est.tp(j).eco3);
        % co2* (co2st)
        fprintf(fname,'%f, ', obs.tp(j).co2st);
        fprintf(fname,'%f, ', obs.tp(j).eco2st);
        fprintf(fname,'%f, ', est.tp(j).co2st);
        fprintf(fname,'%f, ', est.tp(j).eco2st);

        % ph
        fprintf(fname,'%f, ', obs.tp(j).ph);
        fprintf(fname,'%f, ', obs.tp(j).eph);
        fprintf(fname,'%f, ', est.tp(j).ph);
        fprintf(fname,'%f, ', est.tp(j).eph);
        % ph_free
        fprintf(fname,'%f, ', est.tp(j).ph_free);
        fprintf(fname,'%f, ', est.tp(j).eph_free);
        % ph_tot
        fprintf(fname,'%f, ', est.tp(j).ph_tot);
        fprintf(fname,'%f, ', est.tp(j).eph);
        % ph_sws
        fprintf(fname,'%f, ', est.tp(j).ph_sws);
        fprintf(fname,'%f, ', est.tp(j).eph);
        % ph_nbs
        fprintf(fname,'%f, ', est.tp(j).ph_nbs);
        fprintf(fname,'%f, ', est.tp(j).eph);

         % fH = activity coefficient
        fprintf(fname,'%f, ', est.tp(j).fH);
        fprintf(fname,'%f, ', est.tp(j).efH);
        % pp2f
        fprintf(fname,'%f, ', est.tp(j).pp2f);
        fprintf(fname,'%f, ', est.tp(j).epp2f);

        % pK0
        fprintf(fname,'%f, ', est.tp(j).pK0);
        fprintf(fname,'%f, ', est.tp(j).epK0);
        % pK1
        fprintf(fname,'%f, ', est.tp(j).pK1);
        fprintf(fname,'%f, ', est.tp(j).epK1);
        % pK2
        fprintf(fname,'%f, ', est.tp(j).pK2);
        fprintf(fname,'%f, ', est.tp(j).epK2);

        % oh
        fprintf(fname,'%f, ', est.tp(j).oh);
        fprintf(fname,'%f, ', est.tp(j).eoh);
        % pKw = [h][oh]
        fprintf(fname,'%f, ', est.tp(j).pKw);
        fprintf(fname,'%f, ', est.tp(j).epKw);

        % boh4
        fprintf(fname,'%f, ', est.tp(j).boh4);
        fprintf(fname,'%f, ', est.tp(j).eboh4);
        % boh3
        fprintf(fname,'%f, ', est.tp(j).boh3);
        fprintf(fname,'%f, ', est.tp(j).eboh3);
        % pKb = [h][boh4]/[boh3]
        fprintf(fname,'%f, ', est.tp(j).pKb);
        fprintf(fname,'%f, ', est.tp(j).epKb);

        % so4
        fprintf(fname,'%f, ', est.tp(j).so4);
        fprintf(fname,'%f, ', est.tp(j).eso4);
        % hso4
        fprintf(fname,'%f, ', est.tp(j).hso4);
        fprintf(fname,'%f, ', est.tp(j).ehso4);
        % pKs  = [hf][so4]/[hso4]
        fprintf(fname,'%f, ', est.tp(j).pKs);
        fprintf(fname,'%f, ', est.tp(j).epKs);

        % [F]
        fprintf(fname,'%f, ', est.tp(j).F);
        fprintf(fname,'%f, ', est.tp(j).eF);
        % [HF] hydrogen fluoride
        fprintf(fname,'%f, ', est.tp(j).HF);
        fprintf(fname,'%f, ', est.tp(j).eHF);
        % pKf = [h][F]/[HF]
        fprintf(fname,'%f, ', est.tp(j).pKf);
        fprintf(fname,'%f, ', est.tp(j).epKf);

        % po4
        fprintf(fname,'%f, ', est.tp(j).po4);
        fprintf(fname,'%f, ', est.tp(j).epo4);
        % hpo4
        fprintf(fname,'%f, ', est.tp(j).hpo4);
        fprintf(fname,'%f, ', est.tp(j).ehpo4);
        % h2po4
        fprintf(fname,'%f, ', est.tp(j).h2po4);
        fprintf(fname,'%f, ', est.tp(j).eh2po4);
        % h3po4
        fprintf(fname,'%f, ', est.tp(j).h3po4);
        fprintf(fname,'%f, ', est.tp(j).eh3po4);
        % pKp1 = [h][h2po4]/[h3po4]
        fprintf(fname,'%f, ', est.tp(j).pKp1);
        fprintf(fname,'%f, ', est.tp(j).epKp1);
        % pKp2 = [h][hpo4]/[h2po4]
        fprintf(fname,'%f, ', est.tp(j).pKp2);
        fprintf(fname,'%f, ', est.tp(j).epKp2);
        % pKp3 = [h][po4]/[hpo4]
        fprintf(fname,'%f, ', est.tp(j).pKp3);
        fprintf(fname,'%f, ', est.tp(j).epKp3);

        % sioh4
        fprintf(fname,'%f, ', est.tp(j).sioh4);
        fprintf(fname,'%f, ', est.tp(j).esioh4);
        % siooh3
        fprintf(fname,'%f, ', est.tp(j).siooh3);
        fprintf(fname,'%f, ', est.tp(j).esiooh3);
        % pKSi = [h][siooh3]/[sioh4]
        fprintf(fname,'%f, ', est.tp(j).pKsi);
        fprintf(fname,'%f, ', est.tp(j).epKsi);

        % nh3
        fprintf(fname,'%f, ', est.tp(j).nh3);
        fprintf(fname,'%f, ', est.tp(j).enh3);
        % nh4
        fprintf(fname,'%f, ', est.tp(j).nh4);
        fprintf(fname,'%f, ', est.tp(j).enh4);
        % pKnh4 = [h][nh3]/[nh4]
        fprintf(fname,'%f, ', est.tp(j).pKnh4);
        fprintf(fname,'%f, ', est.tp(j).epKnh4);

        % hs
        fprintf(fname,'%f, ', est.tp(j).HS);
        fprintf(fname,'%f, ', est.tp(j).eHS);
        % h2s
        fprintf(fname,'%f, ', est.tp(j).H2S);
        fprintf(fname,'%f, ', est.tp(j).eH2S);
        % pKh2s = [h][hs]/[h2s]
        fprintf(fname,'%f, ', est.tp(j).pKh2s);
        fprintf(fname,'%f, ', est.tp(j).epKh2s);

        % ca
        fprintf(fname,'%f, ', est.tp(j).ca);
        fprintf(fname,'%f, ', est.tp(j).eca);
        % OmegaAr
        fprintf(fname,'%f, ', est.tp(j).OmegaAr);
        fprintf(fname,'%f, ', est.tp(j).eOmegaAr);
        % OmegaCa
        fprintf(fname,'%f, ', est.tp(j).OmegaCa);
        fprintf(fname,'%f, ', est.tp(j).eOmegaCa);
        % pKar = [ca][co3]/[omegaAr]
        fprintf(fname,'%f, ', est.tp(j).pKar);
        fprintf(fname,'%f, ', est.tp(j).epKar);
        % pKca = [ca][co3]/[omegaCa]
        fprintf(fname,'%f, ', est.tp(j).pKca);
        fprintf(fname,'%f, ', est.tp(j).epKca);

        if opt.Revelle == 1
            fprintf(fname,'%f, ', est.tp(j).Revelle);
            % fprintf(fname,'%f, ', nan); % est.tp(j).eRevelle);
            fprintf(fname,'%f, ', est.tp(j).dpfco2dpTA);
            % fprintf(fname,'%f, ', nan); % est.tp(j).edpfco2dpTA);
        end
    end

fprintf(fname,'\n'); % end row of data
out = [];


end % function


% input = find_inputs(obs); % find which are input

% function input = find_inputs(obs)
% % find which variables in obs struct are input
% % assume same types of inputs from first datapoint thru last
%     isgood = @(thing) ( ~isempty(thing) & ~sum(isnan(thing)) );
% 
%     input = {'sal'}; % 'sal' must be input
%     %fno = fieldnames(obs); % 23x1
%     %fntpo = fieldnames(obs.tp); % 92x1
% 
%     if isgood(obs.TC)
%         input = [input, 'TC'];
%     elseif isgood(obs.TA)
%         input = [input, 'TA'];
%     elseif isgood(obs.TB)
%         input = [input,'TB'];
%     elseif isgood(obs.TS)
%         input = [input,'TS'];
%     elseif isgood(obs.TF)
%         input = [input,'TF'];
%     elseif isgood(obs.TP)
%         input = [input,'TP'];
%     elseif isgood(obs.TSi)
%         input = [input,'TSi'];
%     elseif isgood(obs.TNH4)
%         input = [input,'TNH4'];
%     elseif isgood(obs.TH2S)
%         input = [input,'TH2S'];
%     elseif isgood(obs.TCa)
%         input = [input,'TCa'];
%     end
% 
%     nTP = length(obs.tp);
%     for i = 1:nTP
%         input.tp(i) = ['T','P']; % always has a T & P
%         if isgood(obs.tp(i).fco2)
%             input.tp(i) = [input.tp(i),'fco2'];
%         elseif isgood(obs.tp(i).pco2)
%             input.tp(i) = [input.tp(i),'pco2'];
%         elseif isgood(obs.tp(i).co2st)
%             input.tp(i) = [input.tp(i),'co2st'];
%         elseif isgood(obs.tp(i).co3)
%             input.tp(i) = [input.tp(i),'co3'];
%         elseif isgood(obs.tp(i).hco3)
%             input.tp(i) = [input.tp(i),'hco3'];
%         elseif isgood(obs.tp(i).ph)
%             input.tp(i) = [input.tp(i),'ph'];
%         elseif isgood(obs.tp(i).ph_free)
%             input.tp(i) = [input.tp(i),'ph_free'];
%         elseif isgood(obs.tp(i).ph_tot)
%             input.tp(i) = [input.tp(i),'ph_tot'];
%         elseif isgood(obs.tp(i).ph_sws)
%             input.tp(i) = [input.tp(i),'ph_sws'];
%         elseif isgood(obs.tp(i).ph_nbs)
%             input.tp(i) = [input.tp(i),'ph_nbs'];
%         end
%     end
% end











