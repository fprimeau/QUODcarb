function out = PrintCSVv2(varargin)
    
    if (nargin == 3)
        %
        % Print the column headers
        %
        sys = varargin{1};
        est = varargin{2};
        fid = varargin{3};

        nTP = length(est.m);

        fprintf(fid,'%s, ','iflag');

        fn = fieldnames(est);
        for i = 1:2:(length(fn)-1) 
            fprintf(fid,'obs_%s, ',fn{i} );
            fprintf(fid,'obs_%s, ',fn{i+1} ); % e
            fprintf(fid,'est_%s, ',fn{i} );
            fprintf(fid,'est_%s, ',fn{i+1} ); % e
        end

        fnm = fieldnames(est.m);
        for j = 1:nTP
            for i = 1:2:(length(fnm))
                fnj = fieldnames(est.m(j));
                fprintf(fid,'obs_%s, ', fnj{i});
                fprintf(fid,'obs_%s, ', fnj{i+1});
                fprintf(fid,'est_%s, ', fnj{i});
                fprintf(fid,'est_%s, ', fnj{i+1});
            end
        end

        fprintf(fid,'\n');
        fprintf(fid,'%s, ','0=good');
        
        for i = 1:(length(fn)-1)
            if ( ( i == 1 ) || (i == 2 ) )
                fprintf(fid,'%s, ', '(PSU)');
                fprintf(fid,'%s, ', '(PSU)');
            else
                fprintf(fid,'%s, ', '(mol/kg)');
                fprintf(fid,'%s, ', '(mol/kg)');
            end
        end
        for j = 1:nTP
            for i = 1:(length(fnm))
                %fnj = fieldnames(est.m(j));
                if ( (i == 1) || (i == 3) ) % T & eT
                    fprintf(fid, '%s, ', 'deg C');
                    fprintf(fid, '%s, ', 'deg C');
                    fprintf(fid, '%s, ', 'deg C');
                    fprintf(fid, '%s, ', 'deg C');
                elseif ( (i == 2) || (i == 4) ) % P & eP
                    fprintf(fid, '%s, ', 'dbar');
                    fprintf(fid, '%s, ', 'dbar');
                    fprintf(fid, '%s, ', 'dbar');
                    fprintf(fid, '%s, ', 'dbar');
                elseif ( (i == 6) || (i == 7) || ... % fco2 & pco2
                         (i == 12) || (i == 13) )  % efco2 & epco2
                    fprintf(fid, '%s, ', '(uatm)');
                    fprintf(fid, '%s, ', '(uatm)');
                    fprintf(fid, '%s, ', '(uatm)');
                    fprintf(fid, '%s, ', '(uatm)');
                else
                    fprintf(fid, '%s, ', '(unitless)'); % log10 unitless
                    fprintf(fid, '%s, ', '  ');
                    fprintf(fid, '%s, ', '  '); 
                    fprintf(fid, '%s, ', '  ');
                end
            end
        end

        fprintf(fid,'\n');

    else
        %
        % Print the output

        sys     = varargin{1}; % posterior estimate
        est     = varargin{2}; % posterior marginal precision
        obs     = varargin{3}; % measurement
        iflag   = varargin{4}; % Newton solver convergence flag: 0 converged, 1 not converged
        fid     = varargin{5}; % file id
          
        nTP = length(est.m);

        fprintf(fid,'%i, ',iflag);
        
        % salinity
        fprintf(fid,'%.4g, ', obs.sal); 
        fprintf(fid,'%.4g, ', obs.wsal);        
        fprintf(fid,'%.4g, ', est.sal);
        fprintf(fid,'%.4g, ', est.esal);

        % TC (DIC)
        fprintf(fid,'%.4g, ', obs.TC); 
        fprintf(fid,'%.4g, ', obs.wTC);
        fprintf(fid,'%.4g, ', est.TC);
        fprintf(fid,'%.4g, ', est.eTC); 
        
        % TA
        fprintf(fid,'%.4g, ', obs.TA);
        fprintf(fid,'%.4g, ', obs.wTA);
        fprintf(fid,'%.4g, ', est.TA);
        fprintf(fid,'%.4g, ', est.eTA);
        
        if ismember('borate',sys.abr)
            % Kb = [h][boh4]/[boh3]    
            fprintf(fid,'%f, ', obs.TB);
            fprintf(fid,'%f, ', obs.wTB);
            fprintf(fid,'%f, ', est.TB);
            fprintf(fid,'%f, ', est.eTB);            
        end
        if ismember('sulfate',sys.abr)
            % Ks  = [hf][so4]/[hso4]   
            fprintf(fid,'%f, ', obs.TS);
            fprintf(fid,'%f, ', obs.wTS);
            fprintf(fid,'%f, ', est.TS);
            fprintf(fid,'%f, ', est.eTS);
        end
        if ismember('fluoride',sys.abr)
            % Kf = [h][F]/[HF]    
            fprintf(fid,'%f, ', obs.TF);
            fprintf(fid,'%f, ', obs.wTF);
            fprintf(fid,'%f, ', est.TF);
            fprintf(fid,'%f, ', est.eTF);
        end
        if ismember('phosphate',sys.abr)
            % K1p = [h][h2po4]/[h3po4]
            % K2p = [h][hpo4]/[h2po4]
            % K3p = [h][po4]/[hpo4]            
            fprintf(fid,'%f, ', obs.TP);
            fprintf(fid,'%f, ', obs.wTP);
            fprintf(fid,'%f, ', est.TP);
            fprintf(fid,'%f, ', est.eTP);
        end
        if ismember('silicate',sys.abr)
            % KSi = [h][siooh3]/[sioh4]    
            fprintf(fid,'%f, ', obs.TSi);
            fprintf(fid,'%f, ', obs.wTSi);
            fprintf(fid,'%f, ', est.TSi);
            fprintf(fid,'%f, ', est.eTSi);
        end
        if ismember('ammonia',sys.abr)
            % Knh4 = [h][nh3]/[nh4]
            fprintf(fid,'%f, ', obs.TNH3);
            fprintf(fid,'%f, ', obs.wTNH3);
            fprintf(fid,'%f, ', est.TNH3);
            fprintf(fid,'%f, ', est.eTNH3);
        end
        if ismember('sulfide',sys.abr)
            % Kh2s = [h][hs]/[h2s]
            fprintf(fid,'%f, ', obs.TH2S);
            fprintf(fid,'%f, ', obs.wTH2S);
            fprintf(fid,'%f, ', est.TH2S);
            fprintf(fid,'%f, ', est.eTH2S);
        end
        for j = 1:nTP
            % Temp
            fprintf(fid,'%f, ', obs.m(j).T); 
            fprintf(fid,'%f, ', obs.m(j).wT);
            fprintf(fid,'%f, ', est.m(j).T); 
            fprintf(fid,'%f, ', est.m(j).eT);
            % pres
            fprintf(fid,'%f, ', obs.m(j).P); 
            fprintf(fid,'%f, ', obs.m(j).wP);
            fprintf(fid,'%f, ', est.m(j).P); 
            fprintf(fid,'%f, ', est.m(j).eP);
            % pH
            fprintf(fid,'%f, ', obs.m(j).ph); 
            fprintf(fid,'%f, ', obs.m(j).wph);
            fprintf(fid,'%f, ', est.m(j).ph); 
            fprintf(fid,'%f, ', est.m(j).eph);
            % fco2
            fprintf(fid,'%f, ', obs.m(j).fco2); 
            fprintf(fid,'%f, ', obs.m(j).wpco2);
            fprintf(fid,'%f, ', est.m(j).fco2); 
            fprintf(fid,'%f, ', est.m(j).efco2);
            % pco2
            fprintf(fid,'%f, ', obs.m(j).pco2); 
            fprintf(fid,'%f, ', obs.m(j).wp2f);
            fprintf(fid,'%f, ', est.m(j).pco2); 
            fprintf(fid,'%f, ', est.m(j).epco2);
            % co2* (co2st)
            fprintf(fid,'%f, ', obs.m(j).co2st); 
            fprintf(fid,'%f, ', obs.m(j).wco2st);
            fprintf(fid,'%f, ', est.m(j).pco2st); 
            fprintf(fid,'%f, ', est.m(j).epco2st);
            % hco3
            fprintf(fid,'%f, ', obs.m(j).hco3); 
            fprintf(fid,'%f, ', obs.m(j).whco3);
            fprintf(fid,'%f, ', est.m(j).phco3); 
            fprintf(fid,'%f, ', est.m(j).ephco3);
            % co3
            fprintf(fid,'%f, ', obs.m(j).co3); 
            fprintf(fid,'%f, ', obs.m(j).wco3);
            fprintf(fid,'%f, ', est.m(j).pco3); 
            fprintf(fid,'%f, ', est.m(j).epco3);
            % p2f
            fprintf(fid,'%f, ', obs.m(j).pp2f); 
            fprintf(fid,'%f, ', obs.m(j).wp2f);
            fprintf(fid,'%f, ', est.m(j).pp2f); 
            fprintf(fid,'%f, ', est.m(j).epp2f);
            % K0
            fprintf(fid,'%f, ', obs.m(j).pK0); 
            fprintf(fid,'%f, ', obs.m(j).wK0);
            fprintf(fid,'%f, ', est.m(j).pK0); 
            fprintf(fid,'%f, ', est.m(j).epK0);
            % K1
            fprintf(fid,'%f, ', obs.m(j).pK1); 
            fprintf(fid,'%f, ', obs.m(j).wK1);
            fprintf(fid,'%f, ', est.m(j).pK1); 
            fprintf(fid,'%f, ', est.m(j).epK1);
            % K2
            fprintf(fid,'%f, ', obs.m(j).pK2); 
            fprintf(fid,'%f, ', obs.m(j).wK2);
            fprintf(fid,'%f, ', est.m(j).pK2); 
            fprintf(fid,'%f, ', est.m(j).epK2);                        
            if ismember('borate',sys.abr)
                % boh4
                fprintf(fid,'%f, ', obs.m(j).boh4); 
                fprintf(fid,'%f, ', obs.m(j).wboh4);
                fprintf(fid,'%f, ', est.m(j).pboh4); 
                fprintf(fid,'%f, ', est.m(j).epboh4);
                % boh3
                fprintf(fid,'%f, ', obs.m(j).boh3); 
                fprintf(fid,'%f, ', obs.m(j).wboh3);
                fprintf(fid,'%f, ', est.m(j).pboh3); 
                fprintf(fid,'%f, ', est.m(j).epboh3);
                % Kb = [h][boh4]/[boh3]
                fprintf(fid,'%f, ', obs.m(j).pKb); 
                fprintf(fid,'%f, ', obs.m(j).wKb);
                fprintf(fid,'%f, ', est.m(j).pKb); 
                fprintf(fid,'%f, ', est.m(j).epKb);              
            end
            if ismember('water',sys.abr)
                % oh
                fprintf(fid,'%f, ', obs.m(j).oh); 
                fprintf(fid,'%f, ', obs.m(j).woh);
                fprintf(fid,'%f, ', est.m(j).poh); 
                fprintf(fid,'%f, ', est.m(j).epoh);
                % Kw = [h][oh]
                fprintf(fid,'%f, ', obs.m(j).pKw); 
                fprintf(fid,'%f, ', obs.m(j).wKw);
                fprintf(fid,'%f, ', est.m(j).pKw); 
                fprintf(fid,'%f, ', est.m(j).epKw);               
            end
            if ismember('sulfate',sys.abr)
                % hf (Hfree)
                fprintf(fid,'%f, ', obs.m(j).phf); 
                fprintf(fid,'%f, ', obs.m(j).wphf);
                fprintf(fid,'%f, ', est.m(j).phf); 
                fprintf(fid,'%f, ', est.m(j).ephf);
                % so4
                fprintf(fid,'%f, ', obs.m(j).so4); 
                fprintf(fid,'%f, ', obs.m(j).wso4);
                fprintf(fid,'%f, ', est.m(j).pso4); 
                fprintf(fid,'%f, ', est.m(j).epso4);
                % hso4
                fprintf(fid,'%f, ', obs.m(j).hso4); 
                fprintf(fid,'%f, ', obs.m(j).whso4);
                fprintf(fid,'%f, ', est.m(j).phso4); 
                fprintf(fid,'%f, ', est.m(j).ephso4);
                % Ks  = [hf][so4]/[hso4]
                fprintf(fid,'%f, ', obs.m(j).pKs); 
                fprintf(fid,'%f, ', obs.m(j).wKs);
                fprintf(fid,'%f, ', est.m(j).pKs); 
                fprintf(fid,'%f, ', est.m(j).epKs);                                               
            end
            if ismember('fluoride',sys.abr)
                % [F]
                fprintf(fid,'%f, ', obs.m(j).F); 
                fprintf(fid,'%f, ', obs.m(j).wF);
                fprintf(fid,'%f, ', est.m(j).pF); 
                fprintf(fid,'%f, ', est.m(j).epF);                
                % [HF] hydrogen fluoride
                fprintf(fid,'%f, ', obs.m(j).HF); 
                fprintf(fid,'%f, ', obs.m(j).wHF);
                fprintf(fid,'%f, ', est.m(j).pHF); 
                fprintf(fid,'%f, ', est.m(j).epHF);                
                % Kf = [h][F]/[HF]
                fprintf(fid,'%f, ', obs.m(j).pKf); 
                fprintf(fid,'%f, ', obs.m(j).wKf);
                fprintf(fid,'%f, ', est.m(j).pKf); 
                fprintf(fid,'%f, ', est.m(j).epKf);                
            end
            if ismember('phosphate',sys.abr)
                % po4
                fprintf(fid,'%f, ', obs.m(j).po4); 
                fprintf(fid,'%f, ', obs.m(j).wpo4);
                fprintf(fid,'%f, ', est.m(j).ppo4); 
                fprintf(fid,'%f, ', est.m(j).eppo4);                
                % hpo4
                fprintf(fid,'%f, ', obs.m(j).hpo4); 
                fprintf(fid,'%f, ', obs.m(j).whpo4);
                fprintf(fid,'%f, ', est.m(j).phpo4); 
                fprintf(fid,'%f, ', est.m(j).ephpo4);                
                % h2po4
                fprintf(fid,'%f, ', obs.m(j).h2po4); 
                fprintf(fid,'%f, ', obs.m(j).wh2po4);
                fprintf(fid,'%f, ', est.m(j).ph2po4); 
                fprintf(fid,'%f, ', est.m(j).eph2po4);                
                % h3po4
                fprintf(fid,'%f, ', obs.m(j).h3po4); 
                fprintf(fid,'%f, ', obs.m(j).wh3po4);
                fprintf(fid,'%f, ', est.m(j).ph3po4); 
                fprintf(fid,'%f, ', est.m(j).eph3po4);                
                % K1p = [h][h2po4]/[h3po4]
                fprintf(fid,'%f, ', obs.m(j).pK1p); 
                fprintf(fid,'%f, ', obs.m(j).wK1p);
                fprintf(fid,'%f, ', est.m(j).pK1p); 
                fprintf(fid,'%f, ', est.m(j).epK1p);                
                % K2p = [h][hpo4]/[h2po4]
                fprintf(fid,'%f, ', obs.m(j).pK2p); 
                fprintf(fid,'%f, ', obs.m(j).wK2p);
                fprintf(fid,'%f, ', est.m(j).pK2p); 
                fprintf(fid,'%f, ', est.m(j).epK2p);                
                % K3p = [h][po4]/[hpo4]
                fprintf(fid,'%f, ', obs.m(j).pK3p); 
                fprintf(fid,'%f, ', obs.m(j).wK3p);
                fprintf(fid,'%f, ', est.m(j).pK3p); 
                fprintf(fid,'%f, ', est.m(j).epK3p);                                              
            end
            if ismember('silicate',sys.abr)
                % sioh4
                fprintf(fid,'%f, ', obs.m(j).sioh4); 
                fprintf(fid,'%f, ', obs.m(j).wsioh4);
                fprintf(fid,'%f, ', est.m(j).psioh4); 
                fprintf(fid,'%f, ', est.m(j).epsioh4);                
                % siooh3
                fprintf(fid,'%f, ', obs.m(j).siooh3); 
                fprintf(fid,'%f, ', obs.m(j).wsiooh3);
                fprintf(fid,'%f, ', est.m(j).psiooh3); 
                fprintf(fid,'%f, ', est.m(j).epsiooh3);                
                % KSi = [h][siooh3]/[sioh4]
                fprintf(fid,'%f, ', obs.m(j).pKsi); 
                fprintf(fid,'%f, ', obs.m(j).wKsi); 
                fprintf(fid,'%f, ', est.m(j).pKsi); 
                fprintf(fid,'%f, ', est.m(j).epKsi);                              
            end
            if ismember('ammonia',sys.abr)
                % nh3
                fprintf(fid,'%f, ', obs.m(j).nh3); 
                fprintf(fid,'%f, ', obs.m(j).wnh3);
                fprintf(fid,'%f, ', est.m(j).pnh3); 
                fprintf(fid,'%f, ', est.m(j).epnh3);                
                % nh4
                fprintf(fid,'%f, ', obs.m(j).nh4); 
                fprintf(fid,'%f, ', obs.m(j).wnh4);
                fprintf(fid,'%f, ', est.m(j).pnh4); 
                fprintf(fid,'%f, ', est.m(j).epnh4);                
                % Knh4 = [h][nh3]/[nh4]
                fprintf(fid,'%f, ', obs.m(j).Knh4); 
                fprintf(fid,'%f, ', obs.m(j).wKnh4);
                fprintf(fid,'%f, ', est.m(j).pKnh4); 
                fprintf(fid,'%f, ', est.m(j).epKnh4);                
            end
            if ismember('sulfide',sys.abr)                
                % hs
                fprintf(fid,'%f, ', obs.m(j).hs); 
                fprintf(fid,'%f, ', obs.m(j).whs);
                fprintf(fid,'%f, ', est.m(j).phs); 
                fprintf(fid,'%f, ', est.m(j).ephs);                
                % h2s
                fprintf(fid,'%f, ', obs.m(j).h2s); 
                fprintf(fid,'%f, ', obs.m(j).wh2s);
                fprintf(fid,'%f, ', est.m(j).ph2s); 
                fprintf(fid,'%f, ', est.m(j).eph2s);                
                % Kh2s = [h][hs]/[h2s]
                fprintf(fid,'%f, ', obs.m(j).Kh2s); 
                fprintf(fid,'%f, ', obs.m(j).wKh2s);
                fprintf(fid,'%f, ', est.m(j).pKh2s); 
                fprintf(fid,'%f, ', est.m(j).epKh2s);                
            end
        end
        
        %p = sys.p; % -log10(x)
        %q = sys.q; % 10^(-x)
        

        %for i = 1:(length(fn)-1)
        %    fprintf(fid,'%f, ', )

%         for i = 1:nv
%             if ( ( i == ih ) | ( i == ihf ) | ( i == iK0 )  | ( i == iK1 )  | ( i == iK2  ) | ( i == iKw )  | ( i == iKb ) | ...
%                  ( i == iKs) | ( i == iKf ) | ( i == iK1p ) | ( i == iK2p ) | ( i == iK3p ) | ( i == iKsi )|(i == iKnh3)|(i==iKh2s) )
%                 conv = 1;
%                 fprintf(fid, '%f, %f, ', conv * yobs(i), conv / sqrt( wobs(i) ) );
%                 fprintf(fid, '%f, %f, ', conv * y(i), conv / sqrt( w(i) ) );
%             elseif ( (i == ifco2) | (i == ipco2) ) 
%                 % observed
%                 conv = 1e6;
%                 yo   = q( yobs(i) );                
%                 yu = q( yobs(i) - 1 / sqrt( wobs(i) ) );
%                 yl = q( yobs(i) + 1 / sqrt( wobs(i) ) );
%                 sigo  = 0.5 * ( yu + yl );
%                 % posterior
%                 yp = q( y(i) );
%                 yu = q( y(i) - 1 / sqrt( w(i) ) );
%                 yl = q( y(i) + 1 / sqrt( w(i) ) );
%                 sigp  = 0.5 * ( yu - yl );
% 
%                 fprintf(fid,'%f, %f, ', conv * yo, conv * sigo);
%                 fprintf(fid,'%f, %f, ', conv * yp, conv * sigp);
%             
%             elseif ( ( i == iT ) | (i == iS ) | ( i == iP) )
%                 % observed
%                 yo   = yobs(i) ;
%                 yu = ( yobs(i) - 1 / sqrt( wobs(i) ) );
%                 yl = ( yobs(i) + 1 / sqrt( wobs(i) ) );
%                 sigo  = 0.5 * ( yu - yl );
%                 % posterior
%                 yp = ( y(i) );
%                 yu = ( y(i) - 1 / sqrt( w(i) ) );
%                 yl = ( y(i) + 1 / sqrt( w(i) ) );
%                 sigp  = 0.5 * ( yu - yl );
% 
%                 fprintf(fid,'%e, %e, ', yo, sigo);
%                 fprintf(fid,'%e, %e, ', yp, sigp);                
%             
%             else
%                 % observed
%                 conv = 1e6;
%                 yo   = q( yobs(i) );
%                 yu = q( yobs(i) - 1 / sqrt( wobs(i) ) );
%                 yl = q( yobs(i) + 1 / sqrt( wobs(i) ) );
%                 sigo  = 0.5 * ( yu - yl );
%                 % posterior
%                 yp = q( y(i) );
%                 yu = q( y(i) - 1 / sqrt( w(i) ) );
%                 yl = q( y(i) + 1 / sqrt( w(i) ) );
%                 sigp  = 0.5 * ( yu - yl );
% 
%                 fprintf(fid,'%e, %e, ', conv * yo, conv * sigo);
%                 fprintf(fid,'%e, %e, ', conv * yp, conv * sigp);
%                 
%             end
%         end
        
        fprintf(fid,'\n');
%         out = [];
    end
end

