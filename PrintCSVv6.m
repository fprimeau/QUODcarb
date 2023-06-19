
% printCSVv6

function out = PrintCSVv6(varargin)
    
    if (nargin == 3)
        %
        % Print the column headers
        %
        est = varargin{1};
        fid = varargin{2};
        opt = varargin{3};

        nTP = length(est.m);

        fprintf(fid,'%s, ','iflag');

        fn = fieldnames(est);
        fnl = length(fn)-1;

        fprintf(fid,'obs_sal, '); % sal is first
        fprintf(fid,'obs_esal, '); % esal
        fprintf(fid,'est_sal, '); 
        fprintf(fid,'est_esal, ');

        for i = 5:6:fnl
            fprintf(fid,'obs_%s, ',fn{i} );
            fprintf(fid,'obs_e%s, ',fn{i} ); % e
            fprintf(fid,'est_%s, ',fn{i} );
            fprintf(fid,'est_e%s, ',fn{i} ); % e
        end

        fnm = fieldnames(est.m); % into T/P dependent
        for j = 1:nTP
            for i = 1:4:8 % until pH
                fnj = fieldnames(est.m(j));
                fprintf(fid,'obs_%s_m(%g), ', fnj{i},j);
                fprintf(fid,'obs_e%s, ', fnj{i}); % error
                fprintf(fid,'est_%s_m(%g), ', fnj{i},j);
                fprintf(fid,'est_e%s, ', fnj{i}); % error
            end
            fprintf(fid,'obs_ph, '); % obs_ph
            fprintf(fid,'obs_eph, '); % eobs_ph
            fprintf(fid,'est_ph, ');
            fprintf(fid,'est_eph, ');
            if ismember('sulfate',opt.abr) && ismember('fluoride',opt.abr)
                fprintf(fid,'est_ph_tot, ');
                fprintf(fid,'est_ph_sws, ');
                fprintf(fid,'est_ph_free, ');
                fprintf(fid,'est_ph_nbs, ');
                for i = 15:6:length(fnm) % fco2 to p2f
                    fnj = fieldnames(est.m(j));
                    fprintf(fid,'obs_%s, ', fnj{i});
                    fprintf(fid,'obs_e%s, ', fnj{i}); % error
                    fprintf(fid,'est_%s, ', fnj{i});
                    fprintf(fid,'est_e%s, ', fnj{i}); % error
                end
            else
                for i = 11:6:length(fnm) % fco2 to p2f
                    fnj = fieldnames(est.m(j));
                    fprintf(fid,'obs_%s, ', fnj{i});
                    fprintf(fid,'obs_e%s, ', fnj{i}); % error
                    fprintf(fid,'est_%s, ', fnj{i});
                    fprintf(fid,'est_e%s, ', fnj{i}); % error
                end
            end
        end

        fprintf(fid,'\n');
        fprintf(fid,'%s, ','0=good');
        
        fprintf(fid,'%s, ', '(PSU)');
        fprintf(fid,'%s, ', '(PSU)');
        fprintf(fid,'%s, ', '(PSU)');
        fprintf(fid,'%s, ', '(PSU)');
        for i = 3:6:fnl
            if ( (i == 21 )  ) % TS in mol/kg
                fprintf(fid,'%s, ', '(mol/kg)');
                fprintf(fid,'%s, ', '(mol/kg)');
                fprintf(fid,'%s, ', '(mol/kg)');
                fprintf(fid,'%s, ', '(mol/kg)');
            elseif ( (i == 57 )  ) % TCal in mol/kg
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

            fprintf(fid, '%s, ', 'deg C');
            fprintf(fid, '%s, ', 'deg C');
            fprintf(fid, '%s, ', 'deg C');
            fprintf(fid, '%s, ', 'deg C');
            fprintf(fid, '%s, ', 'dbar');
            fprintf(fid, '%s, ', 'dbar');
            fprintf(fid, '%s, ', 'dbar');
            fprintf(fid, '%s, ', 'dbar');
            % next 2 are ph
            fprintf(fid, '%s, ', '(p units)'); % log10 unitless
            fprintf(fid, '%s, ', '  ');
            fprintf(fid, '%s, ', '(p units)');
            fprintf(fid, '%s, ', '  ');
            if ismember('sulfate',opt.abr) && ismember('fluoride',opt.abr)
                % 4 more ph outputs
                fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                fprintf(fid, '%s, ', '  ');
                fprintf(fid, '%s, ', '(p units)');
                fprintf(fid, '%s, ', '  ');
                
                for i = 15:6:length(fnm)
                    if ( (strcmp(fnm(i),'fco2')) || ...
                            (strcmp(fnm(i),'pco2'))  )
                        %fco2 & efco2 , pco2 & epco2
                        fprintf(fid, '%s, ', '(uatm)');
                        fprintf(fid, '%s, ', '(uatm)');
                        fprintf(fid, '%s, ', '(uatm)');
                        fprintf(fid, '%s, ', '(uatm)');
                    elseif ( (strcmp(fnm(i),'hco3')) || ...
                            (strcmp(fnm(i),'co3')) )
                        % hco3 % ehco3 , co3 & eco3
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pco2st')) || ...
                            (strcmp(fnm(i),'pp2f')) || ...
                            (strcmp(fnm(i),'pK0')) || ...
                            (strcmp(fnm(i),'pK1')) || ...
                            (strcmp(fnm(i),'pK2')) )
                        % pH , pco2st , pp2f , pK0 , pK1 , pK2
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'oh')) )
                        % oh & eoh
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pKw')) )
                        % pKw & epKw
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'boh4')) || ...
                            (strcmp(fnm(i),'boh3')) )
                        % boh4 & eboh4 , boh3 & eboh3
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pKb')) )
                        % pKb & epKb
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif (strcmp(fnm(i),'hf'))
                        % hf & ehf 
                        fprintf(fid, '%s, ', '(umol/kg)'); 
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'so4')) || ...
                            (strcmp(fnm(i),'hso4')) )
                        % so4 & eso4 , hso4 % ehso4
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif (strcmp(fnm(i),'pKs')) 
                        % phf & ephf , pKs & epKs
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'F')) || ...
                            (strcmp(fnm(i),'HF'))  )
                        % F & eF , HF & eHF
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pKf')) )
                        % pKf & epKf
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'po4')) || ...
                            (strcmp(fnm(i),'hpo4')) || ...
                            (strcmp(fnm(i),'h2po4')) || ...
                            (strcmp(fnm(i),'h3po4')) )
                        % po4 & epo4 , hpo4 % ehpo4 , h2po4 & eh2po4 , h3po4 & eh3po4
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pK1p')) || ...
                            (strcmp(fnm(i),'pK2p')) || ...
                            (strcmp(fnm(i),'pK3p')) )
                        % pK1p , pK2p , pK3p
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'sioh4')) || ...
                            (strcmp(fnm(i),'siooh3')) )
                        % sioh4 & esioh4 , siooh3 & esiooh3
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pKsi'))  )
                        % pKsi & epKsi
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'nh3')) || ...
                            (strcmp(fnm(i),'nh4')) )
                        % nh3 & enh3 , nh4 & enh4
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pKnh4')) )
                        % pKnh4 & epKnh4
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'hs')) || ...
                            (strcmp(fnm(i),'h2s')) )
                        % hs & ehs , h2s & eh2s
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pKh2s')) )
                        % pKh2s & epKh2s
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'ca')) )
                        % ca and eca ( = TCal)
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'OmegaAr')) || ...
                            (strcmp(fnm(i),'OmegaCa')) )
                        % OmegaAr & eOmegaAr, OmegaCa & eOmegaCa
                        fprintf(fid, '%s, ', '   '); % dimensionless
                        fprintf(fid, '%s, ', '   ');
                        fprintf(fid, '%s, ', '   ');
                        fprintf(fid, '%s, ', '   ');
                    elseif ( (strcmp(fnm(i),'pKar')) || ...
                            (strcmp(fnm(i),'pKca')) )
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    end
                end
            else
               for i = 11:6:length(fnm)
                    if ( (strcmp(fnm(i),'fco2')) || ...
                            (strcmp(fnm(i),'pco2'))  )
                        %fco2 & efco2 , pco2 & epco2
                        fprintf(fid, '%s, ', '(uatm)');
                        fprintf(fid, '%s, ', '(uatm)');
                        fprintf(fid, '%s, ', '(uatm)');
                        fprintf(fid, '%s, ', '(uatm)');
                    elseif ( (strcmp(fnm(i),'hco3')) || ...
                            (strcmp(fnm(i),'co3')) )
                        % hco3 % ehco3 , co3 & eco3
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pco2st')) || ...
                            (strcmp(fnm(i),'pp2f')) || ...
                            (strcmp(fnm(i),'pK0')) || ...
                            (strcmp(fnm(i),'pK1')) || ...
                            (strcmp(fnm(i),'pK2')) )
                        % pH , pco2st , pp2f , pK0 , pK1 , pK2
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'oh')) )
                        % oh & eoh
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pKw')) )
                        % pKw & epKw
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'boh4')) || ...
                            (strcmp(fnm(i),'boh3')) )
                        % boh4 & eboh4 , boh3 & eboh3
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pKb')) )
                        % pKb & epKb
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif (strcmp(fnm(i),'hf'))
                        % hf & ehf 
                        fprintf(fid, '%s, ', '(umol/kg)'); 
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'so4')) || ...
                            (strcmp(fnm(i),'hso4')) )
                        % so4 & eso4 , hso4 % ehso4
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif (strcmp(fnm(i),'pKs')) 
                        % phf & ephf , pKs & epKs
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'F')) || ...
                            (strcmp(fnm(i),'HF'))  )
                        % F & eF , HF & eHF
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pKf')) )
                        % pKf & epKf
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'po4')) || ...
                            (strcmp(fnm(i),'hpo4')) || ...
                            (strcmp(fnm(i),'h2po4')) || ...
                            (strcmp(fnm(i),'h3po4')) )
                        % po4 & epo4 , hpo4 % ehpo4 , h2po4 & eh2po4 , h3po4 & eh3po4
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pK1p')) || ...
                            (strcmp(fnm(i),'pK2p')) || ...
                            (strcmp(fnm(i),'pK3p')) )
                        % pK1p , pK2p , pK3p
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'sioh4')) || ...
                            (strcmp(fnm(i),'siooh3')) )
                        % sioh4 & esioh4 , siooh3 & esiooh3
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pKsi'))  )
                        % pKsi & epKsi
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'nh3')) || ...
                            (strcmp(fnm(i),'nh4')) )
                        % nh3 & enh3 , nh4 & enh4
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pKnh4')) )
                        % pKnh4 & epKnh4
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'hs')) || ...
                            (strcmp(fnm(i),'h2s')) )
                        % hs & ehs , h2s & eh2s
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'pKh2s')) )
                        % pKh2s & epKh2s
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    elseif ( (strcmp(fnm(i),'ca')) )
                        % ca and eca ( = TCal)
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                        fprintf(fid, '%s, ', '(umol/kg)');
                    elseif ( (strcmp(fnm(i),'OmegaAr')) || ...
                            (strcmp(fnm(i),'OmegaCa')) )
                        % OmegaAr & eOmegaAr, OmegaCa & eOmegaCa
                        fprintf(fid, '%s, ', '   '); % dimensionless
                        fprintf(fid, '%s, ', '   ');
                        fprintf(fid, '%s, ', '   ');
                        fprintf(fid, '%s, ', '   ');
                    elseif ( (strcmp(fnm(i),'pKar')) || ...
                            (strcmp(fnm(i),'pKca')) )
                        fprintf(fid, '%s, ', '(p units)'); % log10 unitless
                        fprintf(fid, '%s, ', '  ');
                        fprintf(fid, '%s, ', '(p units)');
                        fprintf(fid, '%s, ', '  ');
                    end
                end
            end

        end

        fprintf(fid,'\n');

    else
        %
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
        
        if ismember('borate',opt.abr)
            % Kb = [h][boh4]/[boh3]    
            fprintf(fid,'%f, ', obs.TB);
            fprintf(fid,'%f, ', obs.eTB);
            fprintf(fid,'%f, ', est.TB);
            fprintf(fid,'%f, ', est.eTB);            
        end
        if ismember('sulfate',opt.abr)
            % Ks  = [hf][so4]/[hso4]   
            fprintf(fid,'%f, ', obs.TS);
            fprintf(fid,'%f, ', obs.eTS);
            fprintf(fid,'%f, ', est.TS);
            fprintf(fid,'%f, ', est.eTS);
        end
        if ismember('fluoride',opt.abr)
            % Kf = [h][F]/[HF]    
            fprintf(fid,'%f, ', obs.TF);
            fprintf(fid,'%f, ', obs.eTF);
            fprintf(fid,'%f, ', est.TF);
            fprintf(fid,'%f, ', est.eTF);
        end
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
            % pH
            fprintf(fid,'%f, ', obs.m(j).ph); 
            fprintf(fid,'%f, ', obs.m(j).eph);
            fprintf(fid,'%f, ', est.m(j).ph); 
            fprintf(fid,'%f, ', est.m(j).eph);
            if ismember('sulfate',opt.abr) && ismember('fluoride',opt.abr)
                fprintf(fid,'%f, ', est.m(j).ph_tot); 
                fprintf(fid,'%f, ', est.m(j).ph_sws);
                fprintf(fid,'%f, ', est.m(j).ph_free); 
                fprintf(fid,'%f, ', est.m(j).ph_nbs);
            end
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
            fprintf(fid,'%f, ', obs.m(j).hco3); 
            fprintf(fid,'%f, ', obs.m(j).ehco3);
            fprintf(fid,'%f, ', est.m(j).hco3); 
            fprintf(fid,'%f, ', est.m(j).ehco3);
            % co3
            fprintf(fid,'%f, ', obs.m(j).co3); 
            fprintf(fid,'%f, ', obs.m(j).eco3);
            fprintf(fid,'%f, ', est.m(j).co3); 
            fprintf(fid,'%f, ', est.m(j).eco3);            
            % pco2* (co2st)
            fprintf(fid,'%f, ', obs.m(j).co2st); 
            fprintf(fid,'%f, ', obs.m(j).eco2st);
            fprintf(fid,'%f, ', est.m(j).pco2st); 
            fprintf(fid,'%f, ', est.m(j).epco2st);
            % pp2f
            fprintf(fid,'%f, ', obs.m(j).pp2f); 
            fprintf(fid,'%f, ', obs.m(j).epp2f);
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
            fprintf(fid,'%f, ', obs.m(j).oh);
            fprintf(fid,'%f, ', obs.m(j).eoh);
            fprintf(fid,'%f, ', est.m(j).oh);
            fprintf(fid,'%f, ', est.m(j).eoh);
            % pKw = [h][oh]
            fprintf(fid,'%f, ', obs.m(j).pKw);
            fprintf(fid,'%f, ', obs.m(j).epKw);
            fprintf(fid,'%f, ', est.m(j).pKw);
            fprintf(fid,'%f, ', est.m(j).epKw);
            if ismember('borate',opt.abr)
                % boh4
                fprintf(fid,'%f, ', obs.m(j).boh4); 
                fprintf(fid,'%f, ', obs.m(j).eboh4);
                fprintf(fid,'%f, ', est.m(j).boh4); 
                fprintf(fid,'%f, ', est.m(j).eboh4);
                % boh3
                fprintf(fid,'%f, ', obs.m(j).boh3); 
                fprintf(fid,'%f, ', obs.m(j).eboh3);
                fprintf(fid,'%f, ', est.m(j).boh3); 
                fprintf(fid,'%f, ', est.m(j).eboh3);
                % pKb = [h][boh4]/[boh3]
                fprintf(fid,'%f, ', obs.m(j).pKb); 
                fprintf(fid,'%f, ', obs.m(j).epKb);
                fprintf(fid,'%f, ', est.m(j).pKb); 
                fprintf(fid,'%f, ', est.m(j).epKb);              
            end
            if ismember('sulfate',opt.abr)
                % hf (Hfree)
                fprintf(fid,'%f, ', obs.m(j).phf); 
                fprintf(fid,'%f, ', obs.m(j).ephf);
                fprintf(fid,'%f, ', est.m(j).phf); 
                fprintf(fid,'%f, ', est.m(j).ephf);
                % so4
                fprintf(fid,'%f, ', obs.m(j).so4); 
                fprintf(fid,'%f, ', obs.m(j).eso4);
                fprintf(fid,'%f, ', est.m(j).so4); 
                fprintf(fid,'%f, ', est.m(j).eso4);
                % hso4
                fprintf(fid,'%f, ', obs.m(j).hso4); 
                fprintf(fid,'%f, ', obs.m(j).ehso4);
                fprintf(fid,'%f, ', est.m(j).hso4); 
                fprintf(fid,'%f, ', est.m(j).ehso4);
                % pKs  = [hf][so4]/[hso4]
                fprintf(fid,'%f, ', obs.m(j).pKs); 
                fprintf(fid,'%f, ', obs.m(j).epKs);
                fprintf(fid,'%f, ', est.m(j).pKs); 
                fprintf(fid,'%f, ', est.m(j).epKs);                                               
            end
            if ismember('fluoride',opt.abr)
                % [F]
                fprintf(fid,'%f, ', obs.m(j).F); 
                fprintf(fid,'%f, ', obs.m(j).eF);
                fprintf(fid,'%f, ', est.m(j).F); 
                fprintf(fid,'%f, ', est.m(j).eF);                
                % [HF] hydrogen fluoride
                fprintf(fid,'%f, ', obs.m(j).HF); 
                fprintf(fid,'%f, ', obs.m(j).eHF);
                fprintf(fid,'%f, ', est.m(j).HF); 
                fprintf(fid,'%f, ', est.m(j).eHF);                
                % pKf = [h][F]/[HF]
                fprintf(fid,'%f, ', obs.m(j).pKf); 
                fprintf(fid,'%f, ', obs.m(j).epKf);
                fprintf(fid,'%f, ', est.m(j).pKf); 
                fprintf(fid,'%f, ', est.m(j).epKf);                
            end
            if ismember('phosphate',opt.abr)
                % po4
                fprintf(fid,'%f, ', obs.m(j).po4); 
                fprintf(fid,'%f, ', obs.m(j).epo4);
                fprintf(fid,'%f, ', est.m(j).po4); 
                fprintf(fid,'%f, ', est.m(j).epo4);                
                % hpo4
                fprintf(fid,'%f, ', obs.m(j).hpo4); 
                fprintf(fid,'%f, ', obs.m(j).ehpo4);
                fprintf(fid,'%f, ', est.m(j).hpo4); 
                fprintf(fid,'%f, ', est.m(j).ehpo4);                
                % h2po4
                fprintf(fid,'%f, ', obs.m(j).h2po4); 
                fprintf(fid,'%f, ', obs.m(j).eh2po4);
                fprintf(fid,'%f, ', est.m(j).h2po4); 
                fprintf(fid,'%f, ', est.m(j).eh2po4);                
                % h3po4
                fprintf(fid,'%f, ', obs.m(j).h3po4); 
                fprintf(fid,'%f, ', obs.m(j).eh3po4);
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
                fprintf(fid,'%f, ', obs.m(j).sioh4); 
                fprintf(fid,'%f, ', obs.m(j).esioh4);
                fprintf(fid,'%f, ', est.m(j).sioh4); 
                fprintf(fid,'%f, ', est.m(j).esioh4);                
                % siooh3
                fprintf(fid,'%f, ', obs.m(j).siooh3); 
                fprintf(fid,'%f, ', obs.m(j).esiooh3);
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
                fprintf(fid,'%f, ', obs.m(j).nh3); 
                fprintf(fid,'%f, ', obs.m(j).enh3);
                fprintf(fid,'%f, ', est.m(j).nh3); 
                fprintf(fid,'%f, ', est.m(j).enh3);                
                % nh4
                fprintf(fid,'%f, ', obs.m(j).nh4); 
                fprintf(fid,'%f, ', obs.m(j).enh4);
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
                fprintf(fid,'%f, ', obs.m(j).hs); 
                fprintf(fid,'%f, ', obs.m(j).ehs);
                fprintf(fid,'%f, ', est.m(j).hs); 
                fprintf(fid,'%f, ', est.m(j).ehs);                
                % h2s
                fprintf(fid,'%f, ', obs.m(j).h2s); 
                fprintf(fid,'%f, ', obs.m(j).eh2s);
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
                fprintf(fid,'%f, ', obs.m(j).ca); 
                fprintf(fid,'%f, ', obs.m(j).eca);
                fprintf(fid,'%f, ', est.m(j).ca); 
                fprintf(fid,'%f, ', est.m(j).eca);   
                % OmegaAr
                fprintf(fid,'%f, ', obs.m(j).OmegaAr); 
                fprintf(fid,'%f, ', obs.m(j).eOmegaAr);
                fprintf(fid,'%f, ', est.m(j).OmegaAr); 
                fprintf(fid,'%f, ', est.m(j).eOmegaAr);    
                % OmegaCa
                fprintf(fid,'%f, ', obs.m(j).OmegaCa); 
                fprintf(fid,'%f, ', obs.m(j).eOmegaCa);
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
    end
end

