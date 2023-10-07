

% mksysV8
% in folder 'FP_QUODcarb/v7'

function sys = mksysV8(obs,abr,phscale)
%
% Private function for QUODcarb.m
% it creates the K amd M matrices and rph functions (residual ph)
%
% utility functions and constants
    LOG10 = log(10);
    p = @(x) -log10( x );  % inverse of q    
    q = @(x)  10.^( -x );  % inverse of p
    sys.p = p;
    sys.q = q;

    dqdx = @(x) - LOG10 * 10.^( -x );  % q'
    d2qdx2 = @(x) LOG10^2 * 10.^(-x ); % q"
    dpdx = @(x) -1 / (x .* LOG10);     % p'
    d2pdx2 = @(x) 1 / (x.^2 * LOG10);  % p"
    sys.dqdx = dqdx;
    sys.dpdx = dpdx;
    sys.d2pdx2 = d2pdx2;
    sys.d2qdx2 = d2qdx2;
    
    isgood = @(thing) (~sum(isnan(thing)) & ~isempty(thing));
    if (~isfield(obs,'m'))
        if opt.printmes ~= 0
            error('Need to provide temperature and pressure measurement.')
        end
    end
    field_names = union( fieldnames(obs), fieldnames(obs.m) );
    if ~ismember('sal',field_names)
        error('Error no salinity specified (obs.sal is missing).')
    end
    if ~ismember('T',field_names)
        error('Error no temperature specified (obs.m(:).T is missing).')
    end
    if ~ismember('P',field_names)
        error('Error no pressure specified (obs.m(:).P is missing).')
    end
    i = 0;
    % carbonate:
    i = i + 1; isal = i; sys.isal = isal;
    i = i + 1;  iTC = i;  sys.iTC  = iTC;
    i = i + 1;  iTA = i;  sys.iTA  = iTA;
    % borate: Kb = [h][boh4]/[boh3]
    i = i + 1;  iTB = i;  sys.iTB = iTB;
    % sulfate: Ks = [hf][so4]/[hso4]
    i = i + 1;  iTS = i;  sys.iTS = iTS;
    % fluoride: Kf = [h][F]/[HF]
    i = i + 1;  iTF = i;  sys.iTF = iTF;
    if ismember('phosphate',abr)
        % K1p = [h][h2po4]/[h3po4]
        % K2p = [h][hpo4]/[h2po4]
        % K3p = [h][po4]/[hpo4]
        i = i + 1;  iTP = i; sys.iTP = iTP;
    end
    if ismember('silicate',abr)
        % KSi = [h][siooh3]/[sioh4]
        i = i + 1;  iTSi = i; sys.iTSi = iTSi;
    end
    if ismember('ammonia',abr)
        % Knh4 = [h][nh3]/[nh4]
        i = i + 1;  iTNH3 = i; sys.iTNH3 = iTNH3;
    end
    if ismember('sulfide',abr)
        % Kh2s = [h][hs]/[h2s]
        i = i + 1;  iTH2S = i; sys.iTH2S = iTH2S;
    end
    if ismember('solubility',abr)
        % Kar = [co3][ca]/OmegaAr
        % Kca = [co3][ca]/OmegaCa
        i = i + 1;  iTCal = i; sys.iTCal = iTCal ;
    end

    nTP = length(obs.m); % number of different (T,P) systems

    for j = 1:nTP
        if (isgood(obs.m(j).T) && isgood(obs.m(j).P))
            i = i + 1;    iT = i;  m(j).iT = iT; 
            i = i + 1;    iP = i;  m(j).iP = iP;
        else
            error('Error must specify temperature and pressure.');
        end
        % K0 = [co2st]/fco2   
        % K1 = [h][hco3]/[co2st]
        % K2 = [h][co3]/[hco3]
        nrk = 3;            i = i + 1;
        m(j).iK0      = i;  i = i + 1;
        m(j).iK1      = i;  i = i + 1;
        m(j).iK2      = i;  i = i + 1;
        m(j).ifco2    = i;  i = i + 1;
        m(j).ico2st   = i;  i = i + 1;
        m(j).ihco3    = i;  i = i + 1;
        m(j).ico3     = i;  i = i + 1;
        m(j).iph      = i;  i = i + 1;
        m(j).iph_free = i;  i = i + 1;
        m(j).ip2f     = i;  i = i + 1;
        m(j).ipco2    = i;  i = i + 1; 
        m(j).ipfH     = i; % fH activity coefficient

        % Kw = [h][oh] water = {'Kw','oh' };
        nrk = nrk + 1;      i = i + 1;
        m(j).iKw = i;       i = i + 1;
        m(j).ioh = i;

        % Kb = [h][boh4]/[boh3]
        nrk = nrk + 1;      i = i + 1;
        m(j).iKb   = i;     i = i + 1;
        m(j).iboh4 = i;     i = i + 1;
        m(j).iboh3 = i;

        % Ks  = [hf][so4]/[hso4]
        nrk = nrk + 1;      i = i + 1;
        m(j).iKs     = i;   i = i + 1;
        m(j).iso4    = i;   i = i + 1;
        m(j).ihso4   = i;

        % Kf = [h][F]/[HF]
        nrk = nrk + 1;      i = i + 1;
        m(j).iKf = i;       i = i + 1;
        m(j).iF  = i;       i = i + 1;
        m(j).iHF = i;

        if ismember('phosphate',abr)
            % K1p = [h][h2po4]/[h3po4]
            % K2p = [h][hpo4]/[h2po4]
            % K3p = [h][po4]/[hpo4]
            nrk = nrk + 3;      i = i + 1;
            m(j).iK1p   = i;    i = i + 1;
            m(j).iK2p   = i;    i = i + 1;
            m(j).iK3p   = i;    i = i + 1;
            m(j).ih3po4 = i;    i = i + 1;
            m(j).ih2po4 = i;    i = i + 1;
            m(j).ihpo4  = i;    i = i + 1;
            m(j).ipo4   = i;
        end
        if ismember('silicate',abr)
            % KSi = [h][siooh3]/[sioh4]
            nrk = nrk + 1;      i = i + 1;
            m(j).iKsi    = i;   i = i + 1;
            m(j).isiooh3 = i;   i = i + 1;
            m(j).isioh4  = i;
        end
        if ismember('ammonia',abr)
            % Knh4 = [h][nh3]/[nh4]
            nrk = nrk + 1;      i = i + 1;
            m(j).iKnh4 = i;     i = i + 1;
            m(j).inh3  = i;     i = i + 1;
            m(j).inh4  = i;
        end
        if ismember('sulfide',abr)
            % Kh2s = [h][hs]/[h2s]
            nrk = nrk + 1;      i = i + 1;
            m(j).iKh2s = i;     i = i + 1;
            m(j).ihs   = i;     i = i + 1;
            m(j).ih2s  = i;
        end
        if ismember('solubility',abr)
            % Kar = [co3][ca]/OmegaAr
            % Kca = [co3][ca]/OmegaCa
            nrk = nrk + 1;      i = i + 1;
            m(j).iKar     = i;  i = i + 1;
            m(j).ica      = i;  i = i + 1;
            m(j).iOmegaAr = i;
            nrk = nrk + 1;      i = i + 1;
            m(j).iKca     = i;  i = i + 1;
            m(j).iOmegaCa = i;
        end
    end
    %
    nv = i;
    K = sparse(nTP*2*nrk,nv);
    row = 0;
    for j = 1:nTP
        % K0 = [CO2*]/fCO2 ==> -pK0 + pco2st - pfco2 = 0         
        row = row + 1;
        K(row,[ m(j).iK0, m(j).ico2st, m(j).ifco2]) = [-1, 1, -1];
        row = row + 1;
        m(j).kK0 = row;
        K(row, m(j).iK0) = 1; 
        
        % K1 = [HCO3][H]/[CO2*] ==> -pK1 + phco3 + ph - pco2st = 0
        row = row + 1;
        K(row,[ m(j).iK1, m(j).iph, m(j).ihco3, m(j).ico2st]) = ...
            [-1, 1, 1, -1];
        row = row + 1;
        m(j).kK1 = row;
        K(row, m(j).iK1) = 1;
        
        % K2 = [CO3][H]/[HCO3]
        row = row + 1;
        K(row,[ m(j).iK2, m(j).iph, m(j).ico3, m(j).ihco3 ]) = ...
            [-1, 1, 1, -1];
        row = row + 1;
        m(j).kK2 = row;
        K(row, m(j).iK2) = 1; 
        
        % fco2 = pco2 * p2f;
        row = row + 1;
        K(row,[ m(j).ifco2, m(j).ipco2, m(j).ip2f ]) = [-1 1 1];
        row = row + 1;
        m(j).kp2f = row;
        K(row, m(j).ip2f) = 1; 
        
        % Kw = [OH][H]
        row = row + 1;
        K(row,[ m(j).iKw, m(j).iph, m(j).ioh ]) = [-1, 1, 1];
        row = row+1;
        m(j).kKw = row;
        K(row, m(j).iKw) = 1; % Kw

        % Kb = [H][BOH4]/[BOH3]
        row = row+1;
        K(row,[ m(j).iKb, m(j).iph, m(j).iboh4, m(j).iboh3 ]) = ...
            [-1, 1, 1, -1];
        row = row+1;
        m(j).kKb = row;
        K(row, m(j).iKb) = 1; % Kb

        % KS  = [H]F[SO4]/[HSO4]
        row = row+1;
        K(row,[ m(j).iKs, m(j).iph_free, m(j).iso4, m(j).ihso4 ]) = ...
            [-1, 1, 1, -1];
        row = row+1;
        m(j).kKs = row;
        K(row, m(j).iKs) = 1; % Ks

        % Kf = [H]F[F]/[HF]     [H]Free like in CO2SYS
        row = row+1;
        K(row,[ m(j).iKf, m(j).iph_free, m(j).iF, m(j).iHF ]) = ...
            [-1, 1, 1, -1];
        row = row+1;
        m(j).kKf = row;
        K(row, m(j).iKf) = 1; % Kf

        if (ismember('phosphate',abr))
            % K1p = [H][H2PO4]/[H3PO4]
            row = row+1;
            K(row,[ m(j).iK1p, m(j).iph, m(j).ih2po4, m(j).ih3po4]) = ...
                [-1, 1, 1, -1];
            row = row + 1;
            m(j).kK1p = row;
            K(row, m(j).iK1p) = 1; % K1p

            % K2p = [H][HPO4]/[H2PO4]
            row = row + 1;
            K(row,[ m(j).iK2p, m(j).iph, m(j).ihpo4, m(j).ih2po4 ]) = ...
                [-1, 1, 1, -1];
            row = row + 1;
            m(j).kK2p = row;
            K(row, m(j).iK2p) = 1; % K2p

            % K3p = [H][PO4]/[HPO4]        
            row = row + 1;
            K(row,[ m(j).iK3p, m(j).iph, m(j).ipo4, m(j).ihpo4 ]) = ...
                [-1, 1, 1, -1];
            row = row+1;
            m(j).kK3p = row;
            K(row, m(j).iK3p) = 1; % K3p
        end
        if (ismember('silicate',abr))
            % KSi = [H][SiO(OH)3]/[Si(OH)4]
            row = row + 1;
            K(row,[ m(j).iKsi, m(j).iph, m(j).isiooh3, m(j).isioh4 ]) = ...
                [-1, 1, 1, -1];
            row = row + 1;
            m(j).kKsi = row;
            K(row, m(j).iKsi) = 1; % Ksi
        end
        if (ismember('ammonia',abr))
            % Knh4 = [H][NH3]/[NH4+]
            row = row + 1;
            K(row,[ m(j).iKnh4, m(j).iph, m(j).inh3, m(j).inh4]) = ...
                [-1, 1, 1, -1];
            row = row + 1;
            m(j).kKnh4 = row;
            K(row, m(j).iKnh4) = 1; % Knh4
        end    
        if (ismember('sulfide',abr))
            % Kh2s = [H][HS]/[H2S]
            row = row + 1;
            K(row,[ m(j).iKh2s, m(j).iph, m(j).ihs, m(j).ih2s]) = ...
                [-1, 1, 1, -1];
            row = row + 1;
            m(j).kKh2s = row;
            K(row, m(j).iKh2s) = 1; % Kh2s
        end
        if ismember('solubility',abr)
            % Kar = [co3][ca]/OmegaAr ==> -pKar + pco3 + pca - pOmegaAr = 0
            row = row + 1;
            K(row,[ m(j).iKar, m(j).ico3, m(j).ica, m(j).iOmegaAr]) = ...
                [-1, 1, 1, -1];
            row = row + 1;
            m(j).kKar = row;
            K(row, m(j).iKar) = 1; 

            % Kca = [co3][ca]/OmegaCa ==> -pKca + pco3 + pca - pOmegaCa = 0
            row = row + 1;
            K(row, [ m(j).iKca, m(j).ico3, m(j).ica, m(j).iOmegaCa]) = ...
                [-1, 1, 1, -1];
            row = row + 1;
            m(j).kKca = row;
            K(row, m(j).iKca) = 1;
        end
    end
    % "mass conservation" equations
    nr = 5; % TC and TA from the carbonate system, TB from borate system,
            % TS from sulfate system, TF from fluoride system
    if (ismember('phosphate',abr))
        nr = nr + 1;
    end
    if (ismember('silicate',abr))
        nr = nr + 1;
    end
    if (ismember('ammonia',abr))
        nr = nr + 1;
    end
    if (ismember('sulfide',abr))
        nr = nr + 1;
    end    
    if (ismember('solubility',abr))
        nr = nr + 1;
    end
    M = sparse(nTP*nr,nv);
    row = 0;
    for j = 1:nTP
        % Total alkalinity: TA - [HCO3] - 2[CO3] - [OH] ...
        row = row + 1;
        row_alk = row;
        m(j).jTA = row_alk;
        M(row_alk,[ iTA, m(j).ihco3, m(j).ico3, m(j).ioh ]) = ...
            [1, -1, -2, -1];

        % Total carbonate: TC - [CO2*] - [HCO3] - [CO3] = 0
        row = row + 1;
        M(row, [ iTC, m(j).ico2st, m(j).ihco3, m(j).ico3 ]) = ...
            [1, -1, -1, -1];
        m(j).jTC = row;
        
        % Total borate
        row = row + 1;
        M(row,[ iTB, m(j).iboh3, m(j).iboh4 ])   =  [1, -1, -1];
        M(row_alk, m(j).iboh4) = -1;
        m(j).jTB = row;
        
        % Total sulfate
        row = row + 1;
        M(row, [ iTS, m(j).ihso4, m(j).iso4 ])   =  [1, -1, -1];
        M(row_alk,[ m(j).iph_free, m(j).ihso4 ])  =  [1, 1]; 
        m(j).jTS = row;

        % Total fluoride
        row = row + 1;
        M(row, [ iTF, m(j).iHF, m(j).iF ]) =  [1, -1, -1];
        M(row_alk, m(j).iHF) =  1;
        m(j).jTF = row;
        
        if (ismember('phosphate',abr))
            row = row + 1;
            M(row, [ iTP,m(j).ih3po4, m(j).ih2po4, m(j).ihpo4, m(j).ipo4 ]) ...
                = [1, -1, -1, -1, -1];
            M(row_alk, [ m(j).ihpo4, m(j).ipo4, m(j).ih3po4 ]) ...
                = [-1, -2, 1]; 
            m(j).jTP = row;
        end
        if (ismember('silicate',abr))
            row = row + 1;
            M(row, [ iTSi, m(j).isioh4, m(j).isiooh3]) = [1, -1, -1];
            M(row_alk, m(j).isiooh3) = -1; 
            m(j).jTSi = row;
        end
        if (ismember('ammonia',abr))
            row = row+1;
            M(row, [ iTNH3, m(j).inh4, m(j).inh3 ]) = [1, -1, -1];
            M(row_alk,m(j).inh3) = -1; 
            m(j).jTNH3 = row;
        end
        if (ismember('sulfide',abr)) % H2S
            row = row+1;
            M(row,[ iTH2S, m(j).ih2s, m(j).ihs ]) = [1, -1, -1];
            M(row_alk,m(j).ihs) = -1; 
            m(j).jTH2S = row;
        end
        if (ismember('solubility',abr)) % total Calcium
            row = row + 1;
            M(row, [ iTCal, m(j).ica]) = [1, -1];
            m(j).jTCal = row;
        end

        % % RESIDUALS of ph_free equations based on input phscale 
        % %     (as functions solved in parse_zpK subfunction)
        
        % % ph_free = ph_tot - factor;
        % % rph_free = ph_tot - factor - ph_free; (r = RESIDUAL)
        
        % % convert to ph_tot with conversion factors FREE2tot & SWS2free
        % % where FREE2tot = 1 + TS/Ks; pFREE2tot = p(TS + Ks) - pKs;
        % % and SWS2free = 1/(1 + TS/Ks + TF/Kf); 
            % pSWS2free = -p( 1 + (TS/Ks) + (TF/Kf) ); 
        pFREE2tot = @(z) ( p( q( z(iTS)) + q( z(m(j).iKs)) ) ) ...       
                - (z(m(j).iKs)); % ( Ks/Ks + TS/Ks ) = (Ks + TS)/Ks
        pSWS2free = @(z) -p( 1 + ( q( z(iTS) ) / q( z(m(j).iKs) ) ) ...   
                + ( q(z(iTF) ) / q( z(m(j).iKf) ) ) );

        if phscale == 1 % total
            % ph_free = ph_tot - p(FREE2tot);
            m(j).rph_free = @(z) z(m(j).iph) - pFREE2tot(z) ...
                - z(m(j).iph_free);
        elseif phscale == 2 % sws
            % ph_free = ph_sws + (-pSWS2free);
            m(j).rph_free = @(z) z(m(j).iph) + pSWS2free(z) ...
                - z(m(j).iph_free);
        elseif phscale == 3 % free
            % ph_free = ph_free;
            m(j).rph_free = @(z) z(m(j).iph) - z(m(j).iph_free);
        elseif phscale == 4 % nbs
            % ph_free = ph_nbs + (-pSWS2free) + pfH;
            m(j).rph_free = @(z) z(m(j).iph) + pSWS2free(z) ...
                - z(m(j).ipfH) - z(m(j).iph_free);
        end

        % % calculate rph_free first derivatives grph_free
        if phscale == 2 || phscale == 4 % sws or nbs
            gpSWS2free_outer = @(z) -dpdx( 1 + ( q( z(iTS) ) / ... % outer of chain rule
                q( z(m(j).iKs) ) ) + ( q( z(iTF) ) / q( z(m(j).iKf) ) ) );
            gpSWS2free_pTS = @(z) gpSWS2free_outer(z) * ...
                ( dqdx(z(iTS) ) / q( z(m(j).iKs) ) );
            gpSWS2free_pKs = @(z) gpSWS2free_outer(z) * ...
                q( z(iTS) ) * -1 * ( q( z(m(j).iKs) )^(-2)) * ...
                dqdx( z(m(j).iKs) );
            gpSWS2free_pTF = @(z) gpSWS2free_outer(z) * ...
                ( dqdx( z(iTF)) / q( z(m(j).iKf) ) );
            gpSWS2free_pKf = @(z) gpSWS2free_outer(z) * ...
                q( z(iTF) ) * -1 * ( q( z(m(j).iKf) )^(-2)) * ...
                dqdx( z(m(j).iKf) );
            if phscale == 2 % sws
                g_pfH = @(z) 0;
                m(j).grph_free = @(z) [ gpSWS2free_pTS(z); ...
                    gpSWS2free_pKs(z);  gpSWS2free_pTF(z); ...
                    gpSWS2free_pKf(z);           g_pfH(z); ...
                      1              ;                -1 ];
                    % ^ grph_free_ph |grph_free_phfree ^
            elseif phscale == 4 % nbs
                g_pfH = @(z) -1;
                m(j).grph_free = @(z) [ gpSWS2free_pTS(z); gpSWS2free_pKs(z); ...
                    gpSWS2free_pTF(z); gpSWS2free_pKf(z); ...
                    g_pfH(z); 1              ; -1];
                    %         ^ grph_free_ph |  ^ grph_free_phfree
            end
        elseif phscale == 1 % tot
            gpF2t_outside = @(z) -dpdx( q( z(iTS)) + q( z(m(j).iKs)) );
            gpF2t_pTS = @(z) gpF2t_outside(z) * dqdx( z(iTS) );
            gpF2t_pKs = @(z) gpF2t_outside(z) * dqdx( z(m(j).iKs) ) + 1;
            m(j).grph_free = @(z) [ gpF2t_pTS(z); gpF2t_pKs(z); ...
                0       ; 0     ; ...   % g_pTF; g_pKf;
                0       ; 1     ; -1 ]; % g_pfH; g_ph; g_phfree;
        elseif phscale == 3 % free
            m(j).grph_free = @(z) [ 0; 0; 0; 0; 0; 1; -1 ];
            % g_pTS; g_pKs; g_pTF; g_pKf; g_pfH; g_ph; g_phfree;
        end


    end
    sys.M = M;
    sys.K = K;
    sys.m = m;   
end

