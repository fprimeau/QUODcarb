
% v7/printCSVv8

function PrintCSVv8(opt,est,obs,iflag,fid,nD)
    if opt.printcsv == 1
          for i = 1:nD
              if i == 1
                  fid = fopen(opt.fid,'w');
                  parse_CSV(est(i),fid,opt); % make column headers
              end
              % fill one row with data
              parse_CSV(opt,est(i),obs(i),iflag(i),fid);
          end % for i = 1:nD
    end % if opt.printcsv == 1
end % function

function parse_CSV(varargin)
    
    if (nargin == 3)
        %
        % Print the column headers
        %
        est = varargin{1};
        fid = varargin{2};
        opt = varargin{3};

        nTP = length(est.tp);

        fprintf(fid, '%s, ', '  ');

        fn = fieldnames(est);
        fnl = length(fn)-1;

        fprintf(fid,'obs.sal, '); % sal is first
        fprintf(fid,'obs.esal, '); % esal
        fprintf(fid,'est.sal, '); 
        fprintf(fid,'est.esal, ');

        for i = 5:6:fnl
            fprintf(fid,'obs.%s, ',fn{i} );
            fprintf(fid,'obs.e%s, ',fn{i} ); % e
            fprintf(fid,'est.%s, ',fn{i} );
            fprintf(fid,'est.e%s, ',fn{i} ); % e
        end

        fnm = fieldnames(est.tp); % into T/P dependent
        for j = 1:nTP
            fnj = fieldnames(est.tp(j));
            for i = 1:4:8 % T, P  
                fprintf(fid,'obs.%s, ',  fnj{i});
                fprintf(fid,'obs.e%s, ', fnj{i});
                fprintf(fid,'est.%s, ',  fnj{i});
                fprintf(fid,'est.e%s, ', fnj{i}); % e = error
            end
            for i = 9:6:51 % ph, ph_free, ph_tot, ph_sws, ph_nbs, pfH, fco2, pco2
                fprintf(fid,'obs.%s, ',  fnj{i});
                fprintf(fid,'obs.e%s, ', fnj{i});
                fprintf(fid,'est.%s, ',  fnj{i});
                fprintf(fid,'est.e%s, ', fnj{i}); % e = error
            end
            for i = 57 % hco3 
                fprintf(fid,'est.%s, ',  fnj{i});
                fprintf(fid,'est.e%s, ', fnj{i}); % e = error
            end
            for i = 63 % co3
                fprintf(fid,'obs.%s, ',  fnj{i});
                fprintf(fid,'obs.e%s, ', fnj{i});
                fprintf(fid,'est.%s, ',  fnj{i});
                fprintf(fid,'est.e%s, ', fnj{i}); % e = error
            end
            for i = 69:6:75 % pco2st thru p2f
                fprintf(fid,'est.%s, ',  fnj{i});
                fprintf(fid,'est.e%s, ', fnj{i}); % e = error
            end
            for i = 81:6:93 % pK0, pK1, pK2
                fprintf(fid,'obs.%s, ',  fnj{i});
                fprintf(fid,'obs.e%s, ', fnj{i});
                fprintf(fid,'est.%s, ',  fnj{i});
                fprintf(fid,'est.e%s, ', fnj{i}); % e = error
            end
            for i = 99:6:length(fnj)
                if (    (strcmp(fnj(i),'oh'))       || ...  % oh
                        (strcmp(fnj(i),'boh4'))     || ...  % boh4
                        (strcmp(fnj(i),'boh3'))     || ...  % boh3
                        (strcmp(fnj(i),'so4'))      || ...  % so4
                        (strcmp(fnj(i),'hso4'))     || ...  % hso4
                        (strcmp(fnj(i),'F'))        || ...  % F
                        (strcmp(fnj(i),'HF'))       || ...  % HF
                        (strcmp(fnj(i),'po4'))      || ...  % po4
                        (strcmp(fnj(i),'hpo4'))     || ...  % hpo4
                        (strcmp(fnj(i),'h2po4'))    || ...  % h2po4
                        (strcmp(fnj(i),'h3po4'))    || ...  % h3po4
                        (strcmp(fnj(i),'sioh4'))    || ...  % sioh4
                        (strcmp(fnj(i),'siooh3'))   || ...  % siooh3
                        (strcmp(fnj(i),'nh3'))      || ...  % nh3
                        (strcmp(fnj(i),'nh4'))      || ...  % nh4
                        (strcmp(fnj(i),'hs'))       || ...  % hs
                        (strcmp(fnj(i),'h2s'))      || ...  % h2s
                        (strcmp(fnj(i),'ca'))       || ...  % ca
                        (strcmp(fnj(i),'OmegaAr'))  || ...  % OmegaAr
                        (strcmp(fnj(i),'OmegaCa'))      )   % OmegaCa
                    fprintf(fid,'est.%s, ',  fnj{i});
                    fprintf(fid,'est.e%s, ', fnj{i}); % e = error
                elseif ( (strcmp(fnj(i),'pKw'))     || ...  % pKw
                        (strcmp(fnj(i),'pKb'))      || ...  % pKb
                        (strcmp(fnj(i),'pKs'))      || ...  % pKs
                        (strcmp(fnj(i),'pKf'))      || ...  % pKf
                        (strcmp(fnj(i),'pK1p'))     || ...  % pK1p
                        (strcmp(fnj(i),'pK2p'))     || ...  % pK2p
                        (strcmp(fnj(i),'pK3p'))     || ...  % pK3p
                        (strcmp(fnj(i),'pKsi'))     || ...  % pKsi
                        (strcmp(fnj(i),'pKnh4'))    || ...  % pKnh4
                        (strcmp(fnj(i),'pKh2s'))    || ...  % pKh2s
                        (strcmp(fnj(i),'pKar'))     || ...  % pKar
                        (strcmp(fnj(i),'pKca'))         )   % pKca
                    fprintf(fid,'obs.%s, ',  fnj{i});
                    fprintf(fid,'obs.e%s, ', fnj{i});
                    fprintf(fid,'est.%s, ',  fnj{i});
                    fprintf(fid,'est.e%s, ', fnj{i}); % e = error
                end
            end
        end

        fprintf(fid,'\n'); % finish first line

        % second line, for tp(1), tp(2) etc. so mostly blank
        fprintf(fid, '%s, ','iflag');
        fprintf(fid, '%s, ', '  ');
        fprintf(fid, '%s, ', '  ');
        fprintf(fid, '%s, ', '  ');
        fprintf(fid, '%s, ', '  ');

        for i = 5:6:fnl % temperature independent
            fprintf(fid, '%s, ', '  ');
            fprintf(fid, '%s, ', '  ');
            fprintf(fid, '%s, ', '  ');
            fprintf(fid, '%s, ', '  ');
        end
        for j = 1:nTP
            fnj = fieldnames(est.tp(j));
            for i = 1:4:8
                fprintf(fid, 'tp(%i), ', j);
                fprintf(fid, '%s, ', '  ');
                fprintf(fid, 'tp(%i), ', j);
                fprintf(fid, '%s, ', '  ');
            end
            for i = 9:6:51
                fprintf(fid, 'tp(%i), ', j);
                fprintf(fid, '%s, ', '  ');
                fprintf(fid, 'tp(%i), ', j);
                fprintf(fid, '%s, ', '  ');
            end
            for i = 57
                fprintf(fid, 'tp(%i), ', j);
                fprintf(fid, '%s, ', '  ');
            end
            for i = 63
                fprintf(fid, 'tp(%i), ', j);
                fprintf(fid, '%s, ', '  ');
                fprintf(fid, 'tp(%i), ', j);
                fprintf(fid, '%s, ', '  ');
            end
            for i = 69:6:75
                fprintf(fid, 'tp(%i), ', j);
                fprintf(fid, '%s, ', '  ');
            end
            for i = 81:6:93
                fprintf(fid, 'tp(%i), ', j);
                fprintf(fid, '%s, ', '  ');
                fprintf(fid, 'tp(%i), ', j);
                fprintf(fid, '%s, ', '  ');
            end
            for i = 99:6:length(fnj)
                if (    (strcmp(fnj(i),'oh'))       || ...  % oh
                        (strcmp(fnj(i),'boh4'))     || ...  % boh4
                        (strcmp(fnj(i),'boh3'))     || ...  % boh3
                        (strcmp(fnj(i),'so4'))      || ...  % so4
                        (strcmp(fnj(i),'hso4'))     || ...  % hso4
                        (strcmp(fnj(i),'F'))        || ...  % F
                        (strcmp(fnj(i),'HF'))       || ...  % HF
                        (strcmp(fnj(i),'po4'))      || ...  % po4
                        (strcmp(fnj(i),'hpo4'))     || ...  % hpo4
                        (strcmp(fnj(i),'h2po4'))    || ...  % h2po4
                        (strcmp(fnj(i),'h3po4'))    || ...  % h3po4
                        (strcmp(fnj(i),'sioh4'))    || ...  % sioh4
                        (strcmp(fnj(i),'siooh3'))   || ...  % siooh3
                        (strcmp(fnj(i),'nh3'))      || ...  % nh3
                        (strcmp(fnj(i),'nh4'))      || ...  % nh4
                        (strcmp(fnj(i),'hs'))       || ...  % hs
                        (strcmp(fnj(i),'h2s'))      || ...  % h2s
                        (strcmp(fnj(i),'ca'))       || ...  % ca
                        (strcmp(fnj(i),'OmegaAr'))  || ...  % OmegaAr
                        (strcmp(fnj(i),'OmegaCa'))      )   % OmegaCa
                    fprintf(fid, 'tp(%i), ', j);
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnj(i),'pKw'))     || ...  % pKw
                        (strcmp(fnj(i),'pKb'))      || ...  % pKb
                        (strcmp(fnj(i),'pKs'))      || ...  % pKs
                        (strcmp(fnj(i),'pKf'))      || ...  % pKf
                        (strcmp(fnj(i),'pK1p'))     || ...  % pK1p
                        (strcmp(fnj(i),'pK2p'))     || ...  % pK2p
                        (strcmp(fnj(i),'pK3p'))     || ...  % pK3p
                        (strcmp(fnj(i),'pKsi'))     || ...  % pKsi
                        (strcmp(fnj(i),'pKnh4'))    || ...  % pKnh4
                        (strcmp(fnj(i),'pKh2s'))    || ...  % pKh2s
                        (strcmp(fnj(i),'pKar'))     || ...  % pKar
                        (strcmp(fnj(i),'pKca'))         )   % pKca
                    fprintf(fid, 'tp(%i), ', j);
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, 'tp(%i), ', j);
                    fprintf(fid, '%s, ', '  ');
                end
            end
        end

        fprintf(fid,'\n'); % finish second line

        % third line, units
        fprintf(fid,'%s, ','0=good');
        
        fprintf(fid,'%s, ', '(PSU)');
        fprintf(fid,'%s, ', '(PSU)');
        fprintf(fid,'%s, ', '(PSU)');
        fprintf(fid,'%s, ', '(PSU)');
        for i = 3:6:fnl
            if ( (i == 57 )  ) % TCal in mol/kg
                fprintf(fid,'%s, ', '(mol/kg)');
                fprintf(fid,'%s, ', '(mol/kg)');
                fprintf(fid,'%s, ', '(mol/kg)');
                fprintf(fid,'%s, ', '(mol/kg)');
            else
                fprintf(fid,'%s, ', '(umol/kg)');
                fprintf(fid,'%s, ', '(umol/kg)');
                fprintf(fid,'%s, ', '(umol/kg)');
                fprintf(fid,'%s, ', '(umol/kg)');
            end
        end
        for j = 1:nTP
            fnj = fieldnames(est.tp(j));
            fprintf(fid, '%s, ', 'deg C'); % T
            fprintf(fid, '%s, ', 'deg C');
            fprintf(fid, '%s, ', 'deg C');
            fprintf(fid, '%s, ', 'deg C');
            fprintf(fid, '%s, ', 'dbar'); % P
            fprintf(fid, '%s, ', 'dbar');
            fprintf(fid, '%s, ', 'dbar');
            fprintf(fid, '%s, ', 'dbar');
            for i = 9:6:length(fnj)
                if  ( (strcmp(fnj(i),'ph'))         || ...  % ph (chosen scale)
                        (strcmp(fnj(i),'ph_free'))  || ...  % ph_free
                        (strcmp(fnj(i),'ph_tot'))   || ...  % ph_tot
                        (strcmp(fnj(i),'ph_sws'))   || ...  % ph_sws
                        (strcmp(fnj(i),'ph_nbs'))   || ...  % ph_nbs
                        (strcmp(fnj(i),'pfH')) ) % )        % pfH
                    fprintf(fid, '%s, ', '(p units)');      % log10 unitless
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnj(i),'fco2')) || ...     % fco2
                        (strcmp(fnj(i),'pco2'))  )          % pco2
                    fprintf(fid, '%s, ', '(uatm)');
                    fprintf(fid, '%s, ', '(uatm)');
                    fprintf(fid, '%s, ', '(uatm)');
                    fprintf(fid, '%s, ', '(uatm)');
                elseif ( strcmp(fnj(i),'hco3'))             % hco3
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif (strcmp(fnj(i),'co3'))               % co3
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnj(i),'pco2st')) || ...   % pco2st
                        (strcmp(fnj(i),'pp2f'))  )          % pp2f
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnj(i),'pK0')) || ...      % pK0
                        (strcmp(fnj(i),'pK1')) || ...       % pK1
                        (strcmp(fnj(i),'pK2')) )            % pK2
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnj(i),'oh')) )            % oh
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnj(i),'pKw')) )           % pKw
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnj(i),'boh4')) || ...     % boh4
                        (strcmp(fnj(i),'boh3')) )           % boh3
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnj(i),'pKb')) )           % pKb
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( strcmp(fnj(i),'so4'))              % so4
                    fprintf(fid, '%s, ', '(mol/kg)');
                    fprintf(fid, '%s, ', '(mol/kg)');
                elseif (strcmp(fnj(i),'hso4'))              % hso4
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif (strcmp(fnj(i),'pKs'))               % pKs
                    fprintf(fid, '%s, ', '(p units)'); 
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnj(i),'F')) || ...        % F
                        (strcmp(fnj(i),'HF'))  )            % HF
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnj(i),'pKf')) )           % pKf
                    fprintf(fid, '%s, ', '(p units)'); 
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnj(i),'po4')) || ...      % po4
                        (strcmp(fnj(i),'hpo4')) || ...      % hpo4
                        (strcmp(fnj(i),'h2po4')) || ...     % h2po4
                        (strcmp(fnj(i),'h3po4')) )          % h3po4
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnj(i),'pK1p')) || ...     % pK1p
                        (strcmp(fnj(i),'pK2p')) || ...      % pK2p
                        (strcmp(fnj(i),'pK3p')) )           % pK3p
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnj(i),'sioh4')) || ...    % sioh4
                        (strcmp(fnj(i),'siooh3')) )         % siooh3
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnj(i),'pKsi'))  )         % pKsi
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnj(i),'nh3')) || ...      % nh3
                        (strcmp(fnj(i),'nh4')) )            % nh4
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnj(i),'pKnh4')) )         % pKnh4
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnj(i),'hs')) || ...       % hs
                        (strcmp(fnj(i),'h2s')) )            % h2s
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnj(i),'pKh2s')) )         % pKh2s
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnj(i),'ca')) )            % ca ( = TCal)
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnj(i),'OmegaAr')) || ...  % OmegaAr
                        (strcmp(fnj(i),'OmegaCa')) )        % OmegaCa
                    fprintf(fid, '%s, ', '   ');            % dimensionless
                    fprintf(fid, '%s, ', '   ');
                elseif ( (strcmp(fnj(i),'pKar')) || ...     % pKar
                        (strcmp(fnj(i),'pKca')) )           % pKca
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                end
            end % for i = 9:6:length(fnj)
        end % for j = 1:nTP
    
        fprintf(fid,'\n'); % 'return' i.e. finish third line

    else
        % Print the output

        opt     = varargin{1}; % posterior estimate
        est     = varargin{2}; % posterior marginal precision
        obs     = varargin{3}; % measurement
        iflag   = varargin{4}; % Newton solver convergence flag: 0 converged, 1 not converged
        fid     = varargin{5}; % file id
          
        nTP = length(est.tp);

        fprintf(fid,'%i, ',iflag);
        
        % salinity
        fprintf(fid,'%.4g, ', obs.sal); 
        fprintf(fid,'%.4g, ', obs.esal);        
        fprintf(fid,'%.4g, ', est.sal);
        fprintf(fid,'%.4g, ', est.esal);

        % TC (DIC)
        fprintf(fid,'%.4g, ', obs.TC); 
        fprintf(fid,'%.4g, ', obs.eTC);
        fprintf(fid,'%.4g, ', est.TC);
        fprintf(fid,'%.4g, ', est.eTC); 
        
        % TA
        fprintf(fid,'%.4g, ', obs.TA);
        fprintf(fid,'%.4g, ', obs.eTA);
        fprintf(fid,'%.4g, ', est.TA);
        fprintf(fid,'%.4g, ', est.eTA);
        
        % Kb = [h][boh4]/[boh3]
        fprintf(fid,'%f, ', obs.TB);
        fprintf(fid,'%f, ', obs.eTB);
        fprintf(fid,'%f, ', est.TB);
        fprintf(fid,'%f, ', est.eTB);

        % Ks  = [hf][so4]/[hso4]
        fprintf(fid,'%f, ', obs.TS);
        fprintf(fid,'%f, ', obs.eTS);
        fprintf(fid,'%f, ', est.TS);
        fprintf(fid,'%f, ', est.eTS);

        % Kf = [h][F]/[HF]
        fprintf(fid,'%f, ', obs.TF);
        fprintf(fid,'%f, ', obs.eTF);
        fprintf(fid,'%f, ', est.TF);
        fprintf(fid,'%f, ', est.eTF);

        if ismember('phosphate',opt.abr)
            % K1p = [h][h2po4]/[h3po4]
            % K2p = [h][hpo4]/[h2po4]
            % K3p = [h][po4]/[hpo4]            
            fprintf(fid,'%f, ', obs.TP);
            fprintf(fid,'%f, ', obs.eTP);
            fprintf(fid,'%f, ', est.TP);
            fprintf(fid,'%f, ', est.eTP);
        end
        if ismember('silicate',opt.abr)
            % KSi = [h][siooh3]/[sioh4]    
            fprintf(fid,'%f, ', obs.TSi);
            fprintf(fid,'%f, ', obs.eTSi);
            fprintf(fid,'%f, ', est.TSi);
            fprintf(fid,'%f, ', est.eTSi);
        end
        if ismember('ammonia',opt.abr)
            % Knh4 = [h][nh3]/[nh4]
            fprintf(fid,'%f, ', obs.TNH3);
            fprintf(fid,'%f, ', obs.eTNH3);
            fprintf(fid,'%f, ', est.TNH3);
            fprintf(fid,'%f, ', est.eTNH3);
        end
        if ismember('sulfide',opt.abr)
            % Kh2s = [h][hs]/[h2s]
            fprintf(fid,'%f, ', obs.TH2S);
            fprintf(fid,'%f, ', obs.eTH2S);
            fprintf(fid,'%f, ', est.TH2S);
            fprintf(fid,'%f, ', est.eTH2S);
        end
        if ismember('solubility',opt.abr)
            % Ksp = [co3][ca]/Omega
            fprintf(fid,'%f, ', obs.TCal);
            fprintf(fid,'%f, ', obs.eTCal);
            fprintf(fid,'%f, ', est.TCal);
            fprintf(fid,'%f, ', est.eTCal);
        end
        for j = 1:nTP
            % Temp
            fprintf(fid,'%f, ', obs.tp(j).T); 
            fprintf(fid,'%f, ', obs.tp(j).eT);
            fprintf(fid,'%f, ', est.tp(j).T); 
            fprintf(fid,'%f, ', est.tp(j).eT);
            % pres
            fprintf(fid,'%f, ', obs.tp(j).P); 
            fprintf(fid,'%f, ', obs.tp(j).eP);
            fprintf(fid,'%f, ', est.tp(j).P); 
            fprintf(fid,'%f, ', est.tp(j).eP);
            % ph
            fprintf(fid,'%f, ', obs.tp(j).ph); 
            fprintf(fid,'%f, ', obs.tp(j).eph);
            fprintf(fid,'%f, ', est.tp(j).ph); 
            fprintf(fid,'%f, ', est.tp(j).eph);
            % ph_free
            fprintf(fid,'%f, ', obs.tp(j).ph_free); 
            fprintf(fid,'%f, ', obs.tp(j).eph_free);
            fprintf(fid,'%f, ', est.tp(j).ph_free); 
            fprintf(fid,'%f, ', est.tp(j).eph_free);
            % ph_tot
            fprintf(fid,'%f, ', obs.tp(j).ph_tot); 
            fprintf(fid,'%f, ', obs.tp(j).eph);
            fprintf(fid,'%f, ', est.tp(j).ph_tot); 
            fprintf(fid,'%f, ', est.tp(j).eph);
            % keyboard
             % ph_sws
            fprintf(fid,'%f, ', obs.tp(j).ph_sws); 
            fprintf(fid,'%f, ', obs.tp(j).eph);
            fprintf(fid,'%f, ', est.tp(j).ph_sws); 
            fprintf(fid,'%f, ', est.tp(j).eph);
            % ph_nbs
            fprintf(fid,'%f, ', obs.tp(j).ph_nbs); 
            fprintf(fid,'%f, ', obs.tp(j).eph);
            fprintf(fid,'%f, ', est.tp(j).ph_nbs); 
            fprintf(fid,'%f, ', est.tp(j).eph);
            % pfH
            fprintf(fid,'%f, ', obs.tp(j).pfH); 
            fprintf(fid,'%f, ', obs.tp(j).epfH);
            fprintf(fid,'%f, ', est.tp(j).pfH); 
            fprintf(fid,'%f, ', est.tp(j).epfH);
            % fco2
            fprintf(fid,'%f, ', obs.tp(j).fco2); 
            fprintf(fid,'%f, ', obs.tp(j).efco2);
            fprintf(fid,'%f, ', est.tp(j).fco2); 
            fprintf(fid,'%f, ', est.tp(j).efco2);
            % pco2
            fprintf(fid,'%f, ', obs.tp(j).pco2); 
            fprintf(fid,'%f, ', obs.tp(j).epco2);
            fprintf(fid,'%f, ', est.tp(j).pco2); 
            fprintf(fid,'%f, ', est.tp(j).epco2);
            % hco3
            fprintf(fid,'%f, ', est.tp(j).hco3); 
            fprintf(fid,'%f, ', est.tp(j).ehco3);
            % co3
            fprintf(fid,'%f, ', obs.tp(j).co3); 
            fprintf(fid,'%f, ', obs.tp(j).eco3);
            fprintf(fid,'%f, ', est.tp(j).co3); 
            fprintf(fid,'%f, ', est.tp(j).eco3);            
            % pco2* (co2st)
            fprintf(fid,'%f, ', est.tp(j).pco2st); 
            fprintf(fid,'%f, ', est.tp(j).epco2st);
            % pp2f
            fprintf(fid,'%f, ', est.tp(j).pp2f); 
            fprintf(fid,'%f, ', est.tp(j).epp2f);
            % pK0
            fprintf(fid,'%f, ', obs.tp(j).pK0); 
            fprintf(fid,'%f, ', obs.tp(j).epK0);
            fprintf(fid,'%f, ', est.tp(j).pK0); 
            fprintf(fid,'%f, ', est.tp(j).epK0);
            % pK1
            fprintf(fid,'%f, ', obs.tp(j).pK1); 
            fprintf(fid,'%f, ', obs.tp(j).epK1);
            fprintf(fid,'%f, ', est.tp(j).pK1); 
            fprintf(fid,'%f, ', est.tp(j).epK1);
            % pK2
            fprintf(fid,'%f, ', obs.tp(j).pK2); 
            fprintf(fid,'%f, ', obs.tp(j).epK2);
            fprintf(fid,'%f, ', est.tp(j).pK2); 
            fprintf(fid,'%f, ', est.tp(j).epK2); 
            % oh
            fprintf(fid,'%f, ', est.tp(j).oh);
            fprintf(fid,'%f, ', est.tp(j).eoh);
            % pKw = [h][oh]
            fprintf(fid,'%f, ', obs.tp(j).pKw);
            fprintf(fid,'%f, ', obs.tp(j).epKw);
            fprintf(fid,'%f, ', est.tp(j).pKw);
            fprintf(fid,'%f, ', est.tp(j).epKw);

            % boh4
            fprintf(fid,'%f, ', est.tp(j).boh4);
            fprintf(fid,'%f, ', est.tp(j).eboh4);
            % boh3
            fprintf(fid,'%f, ', est.tp(j).boh3);
            fprintf(fid,'%f, ', est.tp(j).eboh3);
            % pKb = [h][boh4]/[boh3]
            fprintf(fid,'%f, ', obs.tp(j).pKb);
            fprintf(fid,'%f, ', obs.tp(j).epKb);
            fprintf(fid,'%f, ', est.tp(j).pKb);
            fprintf(fid,'%f, ', est.tp(j).epKb);

            % so4
            fprintf(fid,'%f, ', est.tp(j).so4);
            fprintf(fid,'%f, ', est.tp(j).eso4);
            % hso4
            fprintf(fid,'%f, ', est.tp(j).hso4);
            fprintf(fid,'%f, ', est.tp(j).ehso4);
            % pKs  = [hf][so4]/[hso4]
            fprintf(fid,'%f, ', obs.tp(j).pKs);
            fprintf(fid,'%f, ', obs.tp(j).epKs);
            fprintf(fid,'%f, ', est.tp(j).pKs);
            fprintf(fid,'%f, ', est.tp(j).epKs);

            % [F]
            fprintf(fid,'%f, ', est.tp(j).F);
            fprintf(fid,'%f, ', est.tp(j).eF);
            % [HF] hydrogen fluoride
            fprintf(fid,'%f, ', est.tp(j).HF);
            fprintf(fid,'%f, ', est.tp(j).eHF);
            % pKf = [h][F]/[HF]
            fprintf(fid,'%f, ', obs.tp(j).pKf);
            fprintf(fid,'%f, ', obs.tp(j).epKf);
            fprintf(fid,'%f, ', est.tp(j).pKf);
            fprintf(fid,'%f, ', est.tp(j).epKf);

            if ismember('phosphate',opt.abr)
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
                % pK1p = [h][h2po4]/[h3po4]
                fprintf(fid,'%f, ', obs.tp(j).pK1p); 
                fprintf(fid,'%f, ', obs.tp(j).epK1p);
                fprintf(fid,'%f, ', est.tp(j).pK1p); 
                fprintf(fid,'%f, ', est.tp(j).epK1p);                
                % pK2p = [h][hpo4]/[h2po4]
                fprintf(fid,'%f, ', obs.tp(j).pK2p); 
                fprintf(fid,'%f, ', obs.tp(j).epK2p);
                fprintf(fid,'%f, ', est.tp(j).pK2p); 
                fprintf(fid,'%f, ', est.tp(j).epK2p);                
                % pK3p = [h][po4]/[hpo4]
                fprintf(fid,'%f, ', obs.tp(j).pK3p); 
                fprintf(fid,'%f, ', obs.tp(j).epK3p);
                fprintf(fid,'%f, ', est.tp(j).pK3p); 
                fprintf(fid,'%f, ', est.tp(j).epK3p);                                              
            end
            if ismember('silicate',opt.abr)
                % sioh4
                fprintf(fid,'%f, ', est.tp(j).sioh4); 
                fprintf(fid,'%f, ', est.tp(j).esioh4);                
                % siooh3
                fprintf(fid,'%f, ', est.tp(j).siooh3); 
                fprintf(fid,'%f, ', est.tp(j).esiooh3);                
                % pKSi = [h][siooh3]/[sioh4]
                fprintf(fid,'%f, ', obs.tp(j).pKsi); 
                fprintf(fid,'%f, ', obs.tp(j).epKsi); 
                fprintf(fid,'%f, ', est.tp(j).pKsi); 
                fprintf(fid,'%f, ', est.tp(j).epKsi);                              
            end
            if ismember('ammonia',opt.abr)
                % nh3
                fprintf(fid,'%f, ', est.tp(j).nh3); 
                fprintf(fid,'%f, ', est.tp(j).enh3);                
                % nh4
                fprintf(fid,'%f, ', est.tp(j).nh4); 
                fprintf(fid,'%f, ', est.tp(j).enh4);                
                % pKnh4 = [h][nh3]/[nh4]
                fprintf(fid,'%f, ', obs.tp(j).pKnh4); 
                fprintf(fid,'%f, ', obs.tp(j).epKnh4);
                fprintf(fid,'%f, ', est.tp(j).pKnh4); 
                fprintf(fid,'%f, ', est.tp(j).epKnh4);                
            end
            if ismember('sulfide',opt.abr)                
                % hs
                fprintf(fid,'%f, ', est.tp(j).hs); 
                fprintf(fid,'%f, ', est.tp(j).ehs);                
                % h2s
                fprintf(fid,'%f, ', est.tp(j).h2s); 
                fprintf(fid,'%f, ', est.tp(j).eh2s);                
                % pKh2s = [h][hs]/[h2s]
                fprintf(fid,'%f, ', obs.tp(j).pKh2s); 
                fprintf(fid,'%f, ', obs.tp(j).epKh2s);
                fprintf(fid,'%f, ', est.tp(j).pKh2s); 
                fprintf(fid,'%f, ', est.tp(j).epKh2s);                
            end
            if ismember('solubility',opt.abr)
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
                fprintf(fid,'%f, ', obs.tp(j).pKar); 
                fprintf(fid,'%f, ', obs.tp(j).epKar);
                fprintf(fid,'%f, ', est.tp(j).pKar); 
                fprintf(fid,'%f, ', est.tp(j).epKar);                              
                % pKca = [ca][co3]/[omegaCa]
                fprintf(fid,'%f, ', obs.tp(j).pKca); 
                fprintf(fid,'%f, ', obs.tp(j).epKca);
                fprintf(fid,'%f, ', est.tp(j).pKca); 
                fprintf(fid,'%f, ', est.tp(j).epKca);
            end
        end
        
        fprintf(fid,'\n');
        out = [];

    end % if (nargin == 3) else

end % function




