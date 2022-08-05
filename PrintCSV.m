function out = PrintCSV(varargin)
    
    if (nargin == 2);
        %
        % Print the column headers
        %
        sys = varargin{1};
        fid = varargin{2};
        
        if (ismember('fco2',sys.variables))
            ifco2 = sys.ifco2;
        else
            ifco2 = 0;
        end

        if (ismember('pco2',sys.variables))
            ipco2 = sys.ipco2;
        else
            ipco2 = 0;
        end

        if (ismember('p2f',sys.variables))
            ip2f = sys.ip2f;
        else
            ip2f = 0;
        end   

        if (ismember('h',sys.variables))
            iph = sys.ih;
        else
            iph = 0;
        end
        if (ismember('hf',sys.variables))
            iphf = sys.ihf;
        else
            iphf = 0;
        end
        if (ismember('K0',sys.variables))
            ipK0 = sys.iK0;
        else
            ipK0 = 0;
        end
        if (ismember('K1',sys.variables))
            ipK1 = sys.iK1;
        else
            ipK1 = 0;
        end            
        if (ismember('K2',sys.variables))
            ipK2 = sys.iK2;
        else
            ipK2 = 0;
        end            
        if (ismember('Kw',sys.variables))
            ipKw = sys.iKw;
        else
            ipKw = 0;
        end            
        if (ismember('Kb',sys.variables))
            ipKb = sys.iKb;
        else
            ipKb = 0;
        end            
        if (ismember('KS',sys.variables))
            ipKS = sys.iKS;
        else
            ipKS = 0;
        end            
        if (ismember('KF',sys.variables))
            ipKF = sys.iKF;
        else
            ipKF = 0;
        end            
        if (ismember('K1p',sys.variables))
            ipK1p = sys.iK1p;
        else
            ipK1p = 0;
        end            
        if (ismember('K2p',sys.variables))
            ipK2p = sys.iK2p;
        else
            ipK2p = 0;
        end            
        if (ismember('K3p',sys.variables))
            ipK3p = sys.iK3p;
        else
            ipK3p = 0;
        end            
        if (ismember('KSi',sys.variables))
            ipKSi = sys.iKSi;
        else
            ipKSi = 0;
        end      
        
        nv = length(sys.variables);
        fprintf(fid,'%s, ','iflag');
        for i = 1:nv
            if ( ( i == iph  ) | ( i == iphf ) | ( i == ipK0 )  | ( i == ipK1 )  | ( i == ipK2 )  | ( i == ipKw ) | ( i == ipKb) | ...
                 ( i == ipKS ) | ( i == ipKF ) | ( i == ipK1p ) | ( i == ipK2p ) | ( i == ipK3p ) | ( i == ipKSi ) )
                fprintf(fid, 'obs_p%s, sig_obs_p%s, ',sys.variables{i},sys.variables{i});
                fprintf(fid, 'post_p%s, sig_post_p%s, ',sys.variables{i},sys.variables{i});
            else
                fprintf(fid,'obs_%s, sig_obs_%s, ',sys.variables{i},sys.variables{i});
                fprintf(fid,'post_%s, sig_post_%s, ',sys.variables{i},sys.variables{i});
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,'%s, ','0=good');
        for i = 1:nv
            if ( ( i == iph )  | ( i == iphf ) | ( i == ipK0 )  | ( i == ipK1 )  | ( i == ipK2 )  | ( i == ipKw ) | ( i == ipKb) | ...
                 ( i == ipKS ) | ( i == ipKF ) | ( i == ipK1p ) | ( i == ipK2p ) | ( i == ipK3p ) | ( i == ipKSi ) | ( i == ip2f) )
                fprintf(fid, '%s, %s, ', ' ', ' ');
                fprintf(fid, '%s, %s, ', ' ', ' ');
            elseif ( ( i == ipco2 ) | ( i == ifco2) )
                fprintf(fid,'%s, %s, ','(uatm)','(uatm)');
                fprintf(fid,'%s, %s, ','(uatm)','(uatm)');            
            else
                fprintf(fid,'%s, %s, ','(umol/kgSW)','(umol/kgSW)');
                fprintf(fid,'%s, %s, ','(umol/kgSW)','(umol/kgSW)');
            end
        end
        fprintf(fid,'\n');
    else
        %
        % Print the output
        %
        y     = varargin{1}; % posterior estimate
        w     = varargin{2}; % posterior marginal precision
        yobs  = varargin{3}; % measurement
        wobs  = varargin{4}; % measurement precision
        sys   = varargin{5}; % utility structure
        iflag = varargin{6}; % Newton solver convergence flag: 0 converged, 1 not converged
        fid   = varargin{7}; % file id
        if (ismember('fco2',sys.variables))
            ifco2 = sys.ifco2;
        else
            ifco2 = 0;
        end
        if (ismember('pco2',sys.variables))
            ipco2 = sys.ipco2;
        else
            ipco2 = 0;
        end        
        if (ismember('h',sys.variables))
            ih = sys.ih;
        else
            ih = 0;
        end
        if (ismember('hf',sys.variables))
            ihf = sys.ihf;
        else
            ihf = 0;
        end
        if (ismember('K0',sys.variables))
            iK0 = sys.iK0;
        else
            iK0 = 0;
        end
        if (ismember('K1',sys.variables))
            iK1 = sys.iK1;
        else
            iK1 = 0;
        end            
        if (ismember('K2',sys.variables))
            iK2 = sys.iK2;
        else
            iK2 = 0;
        end            
        if (ismember('Kw',sys.variables))
            iKw = sys.iKw;
        else
            iKw = 0;
        end            
        if (ismember('Kb',sys.variables))
            iKb = sys.iKb;
        else
            iKb = 0;
        end            
        if (ismember('KS',sys.variables))
            iKS = sys.iKS;
        else
            iKS = 0;
        end            
        if (ismember('KF',sys.variables))
            iKF = sys.iKF;
        else
            iKF = 0;
        end            
        if (ismember('K1p',sys.variables))
            iK1p = sys.iK1p;
        else
            iK1p = 0;
        end            
        if (ismember('K2p',sys.variables))
            iK2p = sys.iK2p;
        else
            iK2p = 0;
        end            
        if (ismember('K3p',sys.variables))
            iK3p = sys.iK3p;
        else
            iK3p = 0;
        end            
        if (ismember('KSi',sys.variables))
            iKSi = sys.iKSi;
        else
            iKSi = 0;
        end      
        
        
        p = sys.p; % -log10(x)
        q = sys.q; % 10^(-x)
        
        nv = length(sys.variables);
        fprintf(fid,'%i, ',iflag);
        for i = 1:nv
            if ( ( i == ih ) | ( i == ihf ) | ( i == iK0 )  | ( i == iK1 )  | ( i == iK2  ) | ( i == iKw )  | ( i == iKb ) | ...
                 ( i == iKS) | ( i == iKF ) | ( i == iK1p ) | ( i == iK2p ) | ( i == iK3p ) | ( i == iKSi ) )
                conv = 1;
                fprintf(fid, '%f, %f, ', conv * yobs(i), conv / sqrt( wobs(i) ) );
                fprintf(fid, '%f, %f, ', conv * y(i), conv / sqrt( w(i) ) );
            elseif ( (i == ifco2) | (i == ipco2) ) 
                % observed
                conv = 1e6;
                yo   = q( yobs(i) );                
                yu = q( yobs(i) - 1 / sqrt( wobs(i) ) );
                yl = q( yobs(i) + 1 / sqrt( wobs(i) ) );
                sigo  = 0.5 * ( yu + yl );
                % posterior
                yp = q( y(i) );
                yu = q( y(i) - 1 / sqrt( w(i) ) );
                yl = q( y(i) + 1 / sqrt( w(i) ) );
                sigp  = 0.5 * ( yu - yl );

                fprintf(fid,'%f, %f, ', conv * yo, conv * sigo);
                fprintf(fid,'%f, %f, ', conv * yp, conv * sigp);
            else
                % observed
                conv = 1e6;
                yo   = q( yobs(i) );
                yu = q( yobs(i) - 1 / sqrt( wobs(i) ) );
                yl = q( yobs(i) + 1 / sqrt( wobs(i) ) );
                sigo  = 0.5 * ( yu - yl );
                % posterior
                yp = q( y(i) );
                yu = q( y(i) - 1 / sqrt( w(i) ) );
                yl = q( y(i) + 1 / sqrt( w(i) ) );
                sigp  = 0.5 * ( yu - yl );

                fprintf(fid,'%e, %e, ', conv * yo, conv * sigo);
                fprintf(fid,'%e, %e, ', conv * yp, conv * sigp);
                
            end
        end
        fprintf(fid,'\n');
        out = [];
    end
end

