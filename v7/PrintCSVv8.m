

% v7/printCSVv8

function out = PrintCSVv8(varargin)
    
    if (nargin == 3)
        %
        % Print the column headers
        %
        est = varargin{1};
        fid = varargin{2};
        opt = varargin{3};

        nTP = length(est.m);

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

        fnm = fieldnames(est.m); % into T/P dependent
        for j = 1:nTP
            fnj = fieldnames(est.m(j));
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
                        (strcmp(fnm(i),'hpo4'))     || ...  % hpo4
                        (strcmp(fnm(i),'h2po4'))    || ...  % h2po4
                        (strcmp(fnm(i),'h3po4'))    || ...  % h3po4
                        (strcmp(fnm(i),'sioh4'))    || ...  % sioh4
                        (strcmp(fnm(i),'siooh3'))   || ...  % siooh3
                        (strcmp(fnm(i),'nh3'))      || ...  % nh3
                        (strcmp(fnm(i),'nh4'))      || ...  % nh4
                        (strcmp(fnm(i),'hs'))       || ...  % hs
                        (strcmp(fnm(i),'h2s'))      || ...  % h2s
                        (strcmp(fnm(i),'ca'))       || ...  % ca
                        (strcmp(fnm(i),'OmegaAr'))  || ...  % OmegaAr
                        (strcmp(fnm(i),'OmegaCa'))      )   % OmegaCa
                    fprintf(fid,'est.%s, ',  fnj{i});
                    fprintf(fid,'est.e%s, ', fnj{i}); % e = error
                elseif ( (strcmp(fnm(i),'pKw'))     || ...  % pKw
                        (strcmp(fnm(i),'pKb'))      || ...  % pKb
                        (strcmp(fnm(i),'pKs'))      || ...  % pKs
                        (strcmp(fnm(i),'pKf'))      || ...  % pKf
                        (strcmp(fnm(i),'pK1p'))     || ...  % pK1p
                        (strcmp(fnm(i),'pK2p'))     || ...  % pK2p
                        (strcmp(fnm(i),'pK3p'))     || ...  % pK3p
                        (strcmp(fnm(i),'pKsi'))     || ...  % pKsi
                        (strcmp(fnm(i),'pKnh4'))    || ...  % pKnh4
                        (strcmp(fnm(i),'pKh2s'))    || ...  % pKh2s
                        (strcmp(fnm(i),'pKar'))     || ...  % pKar
                        (strcmp(fnm(i),'pKca'))         )   % pKca
                    fprintf(fid,'obs.%s, ',  fnj{i});
                    fprintf(fid,'obs.e%s, ', fnj{i});
                    fprintf(fid,'est.%s, ',  fnj{i});
                    fprintf(fid,'est.e%s, ', fnj{i}); % e = error
                end
            end
        end

        fprintf(fid,'\n'); % finish first line

        % second line, for m(1), m(2) etc. so mostly blank
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
            fnj = fieldnames(est.m(j));
            for i = 1:4:8
                fprintf(fid, 'm(%i), ', j);
                fprintf(fid, '%s, ', '  ');
                fprintf(fid, 'm(%i), ', j);
                fprintf(fid, '%s, ', '  ');
            end
            for i = 9:6:51
                fprintf(fid, 'm(%i), ', j);
                fprintf(fid, '%s, ', '  ');
                fprintf(fid, 'm(%i), ', j);
                fprintf(fid, '%s, ', '  ');
            end
            for i = 57
                fprintf(fid, 'm(%i), ', j);
                fprintf(fid, '%s, ', '  ');
            end
            for i = 63
                fprintf(fid, 'm(%i), ', j);
                fprintf(fid, '%s, ', '  ');
                fprintf(fid, 'm(%i), ', j);
                fprintf(fid, '%s, ', '  ');
            end
            for i = 69:6:75
                fprintf(fid, 'm(%i), ', j);
                fprintf(fid, '%s, ', '  ');
            end
            for i = 81:6:93
                fprintf(fid, 'm(%i), ', j);
                fprintf(fid, '%s, ', '  ');
                fprintf(fid, 'm(%i), ', j);
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
                        (strcmp(fnm(i),'hpo4'))     || ...  % hpo4
                        (strcmp(fnm(i),'h2po4'))    || ...  % h2po4
                        (strcmp(fnm(i),'h3po4'))    || ...  % h3po4
                        (strcmp(fnm(i),'sioh4'))    || ...  % sioh4
                        (strcmp(fnm(i),'siooh3'))   || ...  % siooh3
                        (strcmp(fnm(i),'nh3'))      || ...  % nh3
                        (strcmp(fnm(i),'nh4'))      || ...  % nh4
                        (strcmp(fnm(i),'hs'))       || ...  % hs
                        (strcmp(fnm(i),'h2s'))      || ...  % h2s
                        (strcmp(fnm(i),'ca'))       || ...  % ca
                        (strcmp(fnm(i),'OmegaAr'))  || ...  % OmegaAr
                        (strcmp(fnm(i),'OmegaCa'))      )   % OmegaCa
                    fprintf(fid, 'm(%i), ', j);
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnm(i),'pKw'))     || ...  % pKw
                        (strcmp(fnm(i),'pKb'))      || ...  % pKb
                        (strcmp(fnm(i),'pKs'))      || ...  % pKs
                        (strcmp(fnm(i),'pKf'))      || ...  % pKf
                        (strcmp(fnm(i),'pK1p'))     || ...  % pK1p
                        (strcmp(fnm(i),'pK2p'))     || ...  % pK2p
                        (strcmp(fnm(i),'pK3p'))     || ...  % pK3p
                        (strcmp(fnm(i),'pKsi'))     || ...  % pKsi
                        (strcmp(fnm(i),'pKnh4'))    || ...  % pKnh4
                        (strcmp(fnm(i),'pKh2s'))    || ...  % pKh2s
                        (strcmp(fnm(i),'pKar'))     || ...  % pKar
                        (strcmp(fnm(i),'pKca'))         )   % pKca
                    fprintf(fid, 'm(%i), ', j);
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, 'm(%i), ', j);
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
            fnm = fieldnames(est.m(j));
            fprintf(fid, '%s, ', 'deg C'); % T
            fprintf(fid, '%s, ', 'deg C');
            fprintf(fid, '%s, ', 'deg C');
            fprintf(fid, '%s, ', 'deg C');
            fprintf(fid, '%s, ', 'dbar'); % P
            fprintf(fid, '%s, ', 'dbar');
            fprintf(fid, '%s, ', 'dbar');
            fprintf(fid, '%s, ', 'dbar');
            for i = 9:6:length(fnm)
                if  ( (strcmp(fnm(i),'ph'))         || ...  % ph (chosen scale)
                        (strcmp(fnm(i),'ph_free'))  || ...  % ph_free
                        (strcmp(fnm(i),'ph_tot'))   || ...  % ph_tot
                        (strcmp(fnm(i),'ph_sws'))   || ...  % ph_sws
                        (strcmp(fnm(i),'ph_nbs'))   || ...  % ph_nbs
                        (strcmp(fnm(i),'pfH')) ) % )        % pfH
                    fprintf(fid, '%s, ', '(p units)');      % log10 unitless
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnm(i),'fco2')) || ...     % fco2
                        (strcmp(fnm(i),'pco2'))  )          % pco2
                    fprintf(fid, '%s, ', '(uatm)');
                    fprintf(fid, '%s, ', '(uatm)');
                    fprintf(fid, '%s, ', '(uatm)');
                    fprintf(fid, '%s, ', '(uatm)');
                elseif ( strcmp(fnm(i),'hco3'))             % hco3
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif (strcmp(fnm(i),'co3'))               % co3
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnm(i),'pco2st')) || ...   % pco2st
                        (strcmp(fnm(i),'pp2f'))  )          % pp2f
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnm(i),'pK0')) || ...      % pK0
                        (strcmp(fnm(i),'pK1')) || ...       % pK1
                        (strcmp(fnm(i),'pK2')) )            % pK2
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnm(i),'oh')) )            % oh
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnm(i),'pKw')) )           % pKw
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnm(i),'boh4')) || ...     % boh4
                        (strcmp(fnm(i),'boh3')) )           % boh3
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnm(i),'pKb')) )           % pKb
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( strcmp(fnm(i),'so4'))              % so4
                    fprintf(fid, '%s, ', '(mol/kg)');
                    fprintf(fid, '%s, ', '(mol/kg)');
                elseif (strcmp(fnm(i),'hso4'))              % hso4
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif (strcmp(fnm(i),'pKs'))               % pKs
                    fprintf(fid, '%s, ', '(p units)'); 
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnm(i),'F')) || ...        % F
                        (strcmp(fnm(i),'HF'))  )            % HF
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnm(i),'pKf')) )           % pKf
                    fprintf(fid, '%s, ', '(p units)'); 
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnm(i),'po4')) || ...      % po4
                        (strcmp(fnm(i),'hpo4')) || ...      % hpo4
                        (strcmp(fnm(i),'h2po4')) || ...     % h2po4
                        (strcmp(fnm(i),'h3po4')) )          % h3po4
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnm(i),'pK1p')) || ...     % pK1p
                        (strcmp(fnm(i),'pK2p')) || ...      % pK2p
                        (strcmp(fnm(i),'pK3p')) )           % pK3p
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnm(i),'sioh4')) || ...    % sioh4
                        (strcmp(fnm(i),'siooh3')) )         % siooh3
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnm(i),'pKsi'))  )         % pKsi
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnm(i),'nh3')) || ...      % nh3
                        (strcmp(fnm(i),'nh4')) )            % nh4
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnm(i),'pKnh4')) )         % pKnh4
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnm(i),'hs')) || ...       % hs
                        (strcmp(fnm(i),'h2s')) )            % h2s
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnm(i),'pKh2s')) )         % pKh2s
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                elseif ( (strcmp(fnm(i),'ca')) )            % ca ( = TCal)
                    fprintf(fid, '%s, ', '(umol/kg)');
                    fprintf(fid, '%s, ', '(umol/kg)');
                elseif ( (strcmp(fnm(i),'OmegaAr')) || ...  % OmegaAr
                        (strcmp(fnm(i),'OmegaCa')) )        % OmegaCa
                    fprintf(fid, '%s, ', '   ');            % dimensionless
                    fprintf(fid, '%s, ', '   ');
                elseif ( (strcmp(fnm(i),'pKar')) || ...     % pKar
                        (strcmp(fnm(i),'pKca')) )           % pKca
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '(p units)');
                    fprintf(fid, '%s, ', '  ');
                end
            end % for i = 9:6:length(fnm)
        end % for j = 1:nTP
    
        fprintf(fid,'\n'); % 'return' i.e. finish third line

    else
        % Print the output

        opt     = varargin{1}; % posterior estimate
        est     = varargin{2}; % posterior marginal precision
        obs     = varargin{3}; % measurement
        iflag   = varargin{4}; % Newton solver convergence flag: 0 converged, 1 not converged
        fid     = varargin{5}; % file id
          
        nTP = length(est.m);

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
            fprintf(fid,'%f, ', obs.m(j).T); 
            fprintf(fid,'%f, ', obs.m(j).eT);
            fprintf(fid,'%f, ', est.m(j).T); 
            fprintf(fid,'%f, ', est.m(j).eT);
            % pres
            fprintf(fid,'%f, ', obs.m(j).P); 
            fprintf(fid,'%f, ', obs.m(j).eP);
            fprintf(fid,'%f, ', est.m(j).P); 
            fprintf(fid,'%f, ', est.m(j).eP);
            % ph
            fprintf(fid,'%f, ', obs.m(j).ph); 
            fprintf(fid,'%f, ', obs.m(j).eph);
            fprintf(fid,'%f, ', est.m(j).ph); 
            fprintf(fid,'%f, ', est.m(j).eph);
            % ph_free
            fprintf(fid,'%f, ', obs.m(j).ph_free); 
            fprintf(fid,'%f, ', obs.m(j).eph_free);
            fprintf(fid,'%f, ', est.m(j).ph_free); 
            fprintf(fid,'%f, ', est.m(j).eph_free);
            % ph_tot
            fprintf(fid,'%f, ', obs.m(j).ph_tot); 
            fprintf(fid,'%f, ', obs.m(j).eph);
            fprintf(fid,'%f, ', est.m(j).ph_tot); 
            fprintf(fid,'%f, ', est.m(j).eph);
            % keyboard
             % ph_sws
            fprintf(fid,'%f, ', obs.m(j).ph_sws); 
            fprintf(fid,'%f, ', obs.m(j).eph);
            fprintf(fid,'%f, ', est.m(j).ph_sws); 
            fprintf(fid,'%f, ', est.m(j).eph);
            % ph_nbs
            fprintf(fid,'%f, ', obs.m(j).ph_nbs); 
            fprintf(fid,'%f, ', obs.m(j).eph);
            fprintf(fid,'%f, ', est.m(j).ph_nbs); 
            fprintf(fid,'%f, ', est.m(j).eph);
            % pfH
            fprintf(fid,'%f, ', obs.m(j).pfH); 
            fprintf(fid,'%f, ', obs.m(j).epfH);
            fprintf(fid,'%f, ', est.m(j).pfH); 
            fprintf(fid,'%f, ', est.m(j).epfH);
            % fco2
            fprintf(fid,'%f, ', obs.m(j).fco2); 
            fprintf(fid,'%f, ', obs.m(j).efco2);
            fprintf(fid,'%f, ', est.m(j).fco2); 
            fprintf(fid,'%f, ', est.m(j).efco2);
            % pco2
            fprintf(fid,'%f, ', obs.m(j).pco2); 
            fprintf(fid,'%f, ', obs.m(j).epco2);
            fprintf(fid,'%f, ', est.m(j).pco2); 
            fprintf(fid,'%f, ', est.m(j).epco2);
            % hco3
            fprintf(fid,'%f, ', est.m(j).hco3); 
            fprintf(fid,'%f, ', est.m(j).ehco3);
            % co3
            fprintf(fid,'%f, ', obs.m(j).co3); 
            fprintf(fid,'%f, ', obs.m(j).eco3);
            fprintf(fid,'%f, ', est.m(j).co3); 
            fprintf(fid,'%f, ', est.m(j).eco3);            
            % pco2* (co2st)
            fprintf(fid,'%f, ', est.m(j).pco2st); 
            fprintf(fid,'%f, ', est.m(j).epco2st);
            % pp2f
            fprintf(fid,'%f, ', est.m(j).pp2f); 
            fprintf(fid,'%f, ', est.m(j).epp2f);
            % pK0
            fprintf(fid,'%f, ', obs.m(j).pK0); 
            fprintf(fid,'%f, ', obs.m(j).epK0);
            fprintf(fid,'%f, ', est.m(j).pK0); 
            fprintf(fid,'%f, ', est.m(j).epK0);
            % pK1
            fprintf(fid,'%f, ', obs.m(j).pK1); 
            fprintf(fid,'%f, ', obs.m(j).epK1);
            fprintf(fid,'%f, ', est.m(j).pK1); 
            fprintf(fid,'%f, ', est.m(j).epK1);
            % pK2
            fprintf(fid,'%f, ', obs.m(j).pK2); 
            fprintf(fid,'%f, ', obs.m(j).epK2);
            fprintf(fid,'%f, ', est.m(j).pK2); 
            fprintf(fid,'%f, ', est.m(j).epK2); 
            % oh
            fprintf(fid,'%f, ', est.m(j).oh);
            fprintf(fid,'%f, ', est.m(j).eoh);
            % pKw = [h][oh]
            fprintf(fid,'%f, ', obs.m(j).pKw);
            fprintf(fid,'%f, ', obs.m(j).epKw);
            fprintf(fid,'%f, ', est.m(j).pKw);
            fprintf(fid,'%f, ', est.m(j).epKw);

            % boh4
            fprintf(fid,'%f, ', est.m(j).boh4);
            fprintf(fid,'%f, ', est.m(j).eboh4);
            % boh3
            fprintf(fid,'%f, ', est.m(j).boh3);
            fprintf(fid,'%f, ', est.m(j).eboh3);
            % pKb = [h][boh4]/[boh3]
            fprintf(fid,'%f, ', obs.m(j).pKb);
            fprintf(fid,'%f, ', obs.m(j).epKb);
            fprintf(fid,'%f, ', est.m(j).pKb);
            fprintf(fid,'%f, ', est.m(j).epKb);

            % so4
            fprintf(fid,'%f, ', est.m(j).so4);
            fprintf(fid,'%f, ', est.m(j).eso4);
            % hso4
            fprintf(fid,'%f, ', est.m(j).hso4);
            fprintf(fid,'%f, ', est.m(j).ehso4);
            % pKs  = [hf][so4]/[hso4]
            fprintf(fid,'%f, ', obs.m(j).pKs);
            fprintf(fid,'%f, ', obs.m(j).epKs);
            fprintf(fid,'%f, ', est.m(j).pKs);
            fprintf(fid,'%f, ', est.m(j).epKs);

            % [F]
            fprintf(fid,'%f, ', est.m(j).F);
            fprintf(fid,'%f, ', est.m(j).eF);
            % [HF] hydrogen fluoride
            fprintf(fid,'%f, ', est.m(j).HF);
            fprintf(fid,'%f, ', est.m(j).eHF);
            % pKf = [h][F]/[HF]
            fprintf(fid,'%f, ', obs.m(j).pKf);
            fprintf(fid,'%f, ', obs.m(j).epKf);
            fprintf(fid,'%f, ', est.m(j).pKf);
            fprintf(fid,'%f, ', est.m(j).epKf);

            if ismember('phosphate',opt.abr)
                % po4
                fprintf(fid,'%f, ', est.m(j).po4); 
                fprintf(fid,'%f, ', est.m(j).epo4);                
                % hpo4
                fprintf(fid,'%f, ', est.m(j).hpo4); 
                fprintf(fid,'%f, ', est.m(j).ehpo4);                
                % h2po4
                fprintf(fid,'%f, ', est.m(j).h2po4); 
                fprintf(fid,'%f, ', est.m(j).eh2po4);                
                % h3po4
                fprintf(fid,'%f, ', est.m(j).h3po4); 
                fprintf(fid,'%f, ', est.m(j).eh3po4);                
                % pK1p = [h][h2po4]/[h3po4]
                fprintf(fid,'%f, ', obs.m(j).pK1p); 
                fprintf(fid,'%f, ', obs.m(j).epK1p);
                fprintf(fid,'%f, ', est.m(j).pK1p); 
                fprintf(fid,'%f, ', est.m(j).epK1p);                
                % pK2p = [h][hpo4]/[h2po4]
                fprintf(fid,'%f, ', obs.m(j).pK2p); 
                fprintf(fid,'%f, ', obs.m(j).epK2p);
                fprintf(fid,'%f, ', est.m(j).pK2p); 
                fprintf(fid,'%f, ', est.m(j).epK2p);                
                % pK3p = [h][po4]/[hpo4]
                fprintf(fid,'%f, ', obs.m(j).pK3p); 
                fprintf(fid,'%f, ', obs.m(j).epK3p);
                fprintf(fid,'%f, ', est.m(j).pK3p); 
                fprintf(fid,'%f, ', est.m(j).epK3p);                                              
            end
            if ismember('silicate',opt.abr)
                % sioh4
                fprintf(fid,'%f, ', est.m(j).sioh4); 
                fprintf(fid,'%f, ', est.m(j).esioh4);                
                % siooh3
                fprintf(fid,'%f, ', est.m(j).siooh3); 
                fprintf(fid,'%f, ', est.m(j).esiooh3);                
                % pKSi = [h][siooh3]/[sioh4]
                fprintf(fid,'%f, ', obs.m(j).pKsi); 
                fprintf(fid,'%f, ', obs.m(j).epKsi); 
                fprintf(fid,'%f, ', est.m(j).pKsi); 
                fprintf(fid,'%f, ', est.m(j).epKsi);                              
            end
            if ismember('ammonia',opt.abr)
                % nh3
                fprintf(fid,'%f, ', est.m(j).nh3); 
                fprintf(fid,'%f, ', est.m(j).enh3);                
                % nh4
                fprintf(fid,'%f, ', est.m(j).nh4); 
                fprintf(fid,'%f, ', est.m(j).enh4);                
                % pKnh4 = [h][nh3]/[nh4]
                fprintf(fid,'%f, ', obs.m(j).pKnh4); 
                fprintf(fid,'%f, ', obs.m(j).epKnh4);
                fprintf(fid,'%f, ', est.m(j).pKnh4); 
                fprintf(fid,'%f, ', est.m(j).epKnh4);                
            end
            if ismember('sulfide',opt.abr)                
                % hs
                fprintf(fid,'%f, ', est.m(j).hs); 
                fprintf(fid,'%f, ', est.m(j).ehs);                
                % h2s
                fprintf(fid,'%f, ', est.m(j).h2s); 
                fprintf(fid,'%f, ', est.m(j).eh2s);                
                % pKh2s = [h][hs]/[h2s]
                fprintf(fid,'%f, ', obs.m(j).pKh2s); 
                fprintf(fid,'%f, ', obs.m(j).epKh2s);
                fprintf(fid,'%f, ', est.m(j).pKh2s); 
                fprintf(fid,'%f, ', est.m(j).epKh2s);                
            end
            if ismember('solubility',opt.abr)
                % ca
                fprintf(fid,'%f, ', est.m(j).ca); 
                fprintf(fid,'%f, ', est.m(j).eca);   
                % OmegaAr
                fprintf(fid,'%f, ', est.m(j).OmegaAr); 
                fprintf(fid,'%f, ', est.m(j).eOmegaAr);    
                % OmegaCa
                fprintf(fid,'%f, ', est.m(j).OmegaCa); 
                fprintf(fid,'%f, ', est.m(j).eOmegaCa);   
                % pKar = [ca][co3]/[omegaAr]
                fprintf(fid,'%f, ', obs.m(j).pKar); 
                fprintf(fid,'%f, ', obs.m(j).epKar);
                fprintf(fid,'%f, ', est.m(j).pKar); 
                fprintf(fid,'%f, ', est.m(j).epKar);                              
                % pKca = [ca][co3]/[omegaCa]
                fprintf(fid,'%f, ', obs.m(j).pKca); 
                fprintf(fid,'%f, ', obs.m(j).epKca);
                fprintf(fid,'%f, ', est.m(j).pKca); 
                fprintf(fid,'%f, ', est.m(j).epKca);
            end
        end
        
        fprintf(fid,'\n');
        out = [];

    end % if (nargin == 3) else

end % function




