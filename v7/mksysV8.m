

% mksysV8
% in folder 'FP_QUODcarb/v7'

function sys = mksysV8(obs,abr,phscale)
%
% Private function for QUODcarb.m
% it creates the K and M matrices and rph functions (residual ph)
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
    if (~isfield(obs,'tp'))
        if opt.printmes ~= 0
            error('Need to provide temperature and pressure measurement.')
        end
    end
    field_names = union( fieldnames(obs), fieldnames(obs.tp) );
    if ~ismember('sal',field_names)
        error('Error no salinity specified (obs.sal is missing).')
    end
    if ~ismember('T',field_names)
        error('Error no temperature specified (obs.tp(:).T is missing).')
    end
    if ~ismember('P',field_names)
        error('Error no pressure specified (obs.tp(:).P is missing).')
    end
    i = 0;
    % carbonate:
    i = i + 1; isal = i;  sys.isal = isal;
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

    nTP = length(obs.tp); % number of different (T,P) systems

    for j = 1:nTP
        if (isgood(obs.tp(j).T) && isgood(obs.tp(j).P))
            i = i + 1;    iT = i;     tp(j).iT = iT; 
            i = i + 1;    iP = i;     tp(j).iP = iP;
        else
            error('Error must specify temperature and pressure.');
        end
        % K0 = [co2st]/fco2   
        % K1 = [h][hco3]/[co2st]
        % K2 = [h][co3]/[hco3]
        nrk = 3;             i = i + 1;
        tp(j).iK0      = i;  i = i + 1;
        tp(j).iK1      = i;  i = i + 1;
        tp(j).iK2      = i;  i = i + 1;
        tp(j).ifco2    = i;  i = i + 1;
        tp(j).ico2st   = i;  i = i + 1;
        tp(j).ihco3    = i;  i = i + 1;
        tp(j).ico3     = i;  i = i + 1;
        tp(j).iph      = i;  i = i + 1;
        tp(j).iph_free = i;  i = i + 1;
        tp(j).ip2f     = i;  i = i + 1;
        tp(j).ipco2    = i;  i = i + 1; 
        tp(j).ipfH     = i; % fH activity coefficient

        % Kw = [h][oh] water = {'Kw','oh' };
        nrk = nrk + 1;       i = i + 1;
        tp(j).iKw = i;       i = i + 1;
        tp(j).ioh = i;

        % Kb = [h][boh4]/[boh3]
        nrk = nrk + 1;       i = i + 1;
        tp(j).iKb   = i;     i = i + 1;
        tp(j).iboh4 = i;     i = i + 1;
        tp(j).iboh3 = i;

        % Ks  = [hf][so4]/[hso4]
        nrk = nrk + 1;       i = i + 1;
        tp(j).iKs     = i;   i = i + 1;
        tp(j).iso4    = i;   i = i + 1;
        tp(j).ihso4   = i;

        % Kf = [h][F]/[HF]
        nrk = nrk + 1;       i = i + 1;
        tp(j).iKf = i;       i = i + 1;
        tp(j).iF  = i;       i = i + 1;
        tp(j).iHF = i;

        if ismember('phosphate',abr)
            % K1p = [h][h2po4]/[h3po4]
            % K2p = [h][hpo4]/[h2po4]
            % K3p = [h][po4]/[hpo4]
            nrk = nrk + 3;       i = i + 1;
            tp(j).iK1p   = i;    i = i + 1;
            tp(j).iK2p   = i;    i = i + 1;
            tp(j).iK3p   = i;    i = i + 1;
            tp(j).ih3po4 = i;    i = i + 1;
            tp(j).ih2po4 = i;    i = i + 1;
            tp(j).ihpo4  = i;    i = i + 1;
            tp(j).ipo4   = i;
        end
        if ismember('silicate',abr)
            % KSi = [h][siooh3]/[sioh4]
            nrk = nrk + 1;       i = i + 1;
            tp(j).iKsi    = i;   i = i + 1;
            tp(j).isiooh3 = i;   i = i + 1;
            tp(j).isioh4  = i;
        end
        if ismember('ammonia',abr)
            % Knh4 = [h][nh3]/[nh4]
            nrk = nrk + 1;       i = i + 1;
            tp(j).iKnh4 = i;     i = i + 1;
            tp(j).inh3  = i;     i = i + 1;
            tp(j).inh4  = i;
        end
        if ismember('sulfide',abr)
            % Kh2s = [h][hs]/[h2s]
            nrk = nrk + 1;       i = i + 1;
            tp(j).iKh2s = i;     i = i + 1;
            tp(j).ihs   = i;     i = i + 1;
            tp(j).ih2s  = i;
        end
        if ismember('solubility',abr)
            % Kar = [co3][ca]/OmegaAr
            % Kca = [co3][ca]/OmegaCa
            nrk = nrk + 1;       i = i + 1;
            tp(j).iKar     = i;  i = i + 1;
            tp(j).ica      = i;  i = i + 1;
            tp(j).iOmegaAr = i;
            nrk = nrk + 1;       i = i + 1;
            tp(j).iKca     = i;  i = i + 1;
            tp(j).iOmegaCa = i;
        end
    end
    %
    nv = i;
    K = sparse(nTP*2*nrk,nv);
    row = 0;
    for j = 1:nTP
        % K0 = [CO2*]/fCO2 ==> -pK0 + pco2st - pfco2 = 0         
        row = row + 1;
        K(row,[ tp(j).iK0, tp(j).ico2st, tp(j).ifco2]) = [-1, 1, -1];
        row = row + 1;
        tp(j).kK0 = row;
        K(row, tp(j).iK0) = 1; 
        
        % K1 = [HCO3][H]/[CO2*] ==> -pK1 + phco3 + ph - pco2st = 0
        row = row + 1;
        K(row,[ tp(j).iK1, tp(j).iph, tp(j).ihco3, tp(j).ico2st]) = ...
            [-1, 1, 1, -1];
        row = row + 1;
        tp(j).kK1 = row;
        K(row, tp(j).iK1) = 1;
        
        % K2 = [CO3][H]/[HCO3]
        row = row + 1;
        K(row,[ tp(j).iK2, tp(j).iph, tp(j).ico3, tp(j).ihco3 ]) = ...
            [-1, 1, 1, -1];
        row = row + 1;
        tp(j).kK2 = row;
        K(row, tp(j).iK2) = 1; 
        
        % fco2 = pco2 * p2f;
        row = row + 1;
        K(row,[ tp(j).ifco2, tp(j).ipco2, tp(j).ip2f ]) = [-1 1 1];
        row = row + 1;
        tp(j).kp2f = row;
        K(row, tp(j).ip2f) = 1; 
        
        % Kw = [OH][H]
        row = row + 1;
        K(row,[ tp(j).iKw, tp(j).iph, tp(j).ioh ]) = [-1, 1, 1];
        row = row+1;
        tp(j).kKw = row;
        K(row, tp(j).iKw) = 1; % Kw

        % Kb = [H][BOH4]/[BOH3]
        row = row+1;
        K(row,[ tp(j).iKb, tp(j).iph, tp(j).iboh4, tp(j).iboh3 ]) = ...
            [-1, 1, 1, -1];
        row = row+1;
        tp(j).kKb = row;
        K(row, tp(j).iKb) = 1; % Kb

        % KS  = [H]F[SO4]/[HSO4]
        row = row+1;
        K(row,[ tp(j).iKs, tp(j).iph_free, tp(j).iso4, tp(j).ihso4 ]) = ...
            [-1, 1, 1, -1];
        row = row+1;
        tp(j).kKs = row;
        K(row, tp(j).iKs) = 1; % Ks

        % Kf = [H]F[F]/[HF]     [H]Free like in CO2SYS
        row = row+1;
        K(row,[ tp(j).iKf, tp(j).iph_free, tp(j).iF, tp(j).iHF ]) = ...
            [-1, 1, 1, -1];
        row = row+1;
        tp(j).kKf = row;
        K(row, tp(j).iKf) = 1; % Kf

        if (ismember('phosphate',abr))
            % K1p = [H][H2PO4]/[H3PO4]
            row = row+1;
            K(row,[ tp(j).iK1p, tp(j).iph, tp(j).ih2po4, ...
                tp(j).ih3po4]) = [-1, 1, 1, -1];
            row = row + 1;
            tp(j).kK1p = row;
            K(row, tp(j).iK1p) = 1; % K1p

            % K2p = [H][HPO4]/[H2PO4]
            row = row + 1;
            K(row,[ tp(j).iK2p, tp(j).iph, tp(j).ihpo4, ...
                tp(j).ih2po4 ]) = [-1, 1, 1, -1];
            row = row + 1;
            tp(j).kK2p = row;
            K(row, tp(j).iK2p) = 1; % K2p

            % K3p = [H][PO4]/[HPO4]        
            row = row + 1;
            K(row,[ tp(j).iK3p, tp(j).iph, tp(j).ipo4, ...
                tp(j).ihpo4 ]) = [-1, 1, 1, -1];
            row = row+1;
            tp(j).kK3p = row;
            K(row, tp(j).iK3p) = 1; % K3p
        end
        if (ismember('silicate',abr))
            % KSi = [H][SiO(OH)3]/[Si(OH)4]
            row = row + 1;
            K(row,[ tp(j).iKsi, tp(j).iph, tp(j).isiooh3, ...
                tp(j).isioh4 ]) = [-1, 1, 1, -1];
            row = row + 1;
            tp(j).kKsi = row;
            K(row, tp(j).iKsi) = 1; % Ksi
        end
        if (ismember('ammonia',abr))
            % Knh4 = [H][NH3]/[NH4+]
            row = row + 1;
            K(row,[ tp(j).iKnh4, tp(j).iph, tp(j).inh3, tp(j).inh4]) = ...
                [-1, 1, 1, -1];
            row = row + 1;
            tp(j).kKnh4 = row;
            K(row, tp(j).iKnh4) = 1; % Knh4
        end    
        if (ismember('sulfide',abr))
            % Kh2s = [H][HS]/[H2S]
            row = row + 1;
            K(row,[ tp(j).iKh2s, tp(j).iph, tp(j).ihs, tp(j).ih2s]) = ...
                [-1, 1, 1, -1];
            row = row + 1;
            tp(j).kKh2s = row;
            K(row, tp(j).iKh2s) = 1; % Kh2s
        end
        if ismember('solubility',abr)
            % Kar = [co3][ca]/OmegaAr ==> -pKar + pco3 + pca - pOmegaAr = 0
            row = row + 1;
            K(row,[ tp(j).iKar, tp(j).ico3, tp(j).ica, ...
                tp(j).iOmegaAr]) = [-1, 1, 1, -1];
            row = row + 1;
            tp(j).kKar = row;
            K(row, tp(j).iKar) = 1; 

            % Kca = [co3][ca]/OmegaCa ==> -pKca + pco3 + pca - pOmegaCa = 0
            row = row + 1;
            K(row, [ tp(j).iKca, tp(j).ico3, tp(j).ica, ...
                tp(j).iOmegaCa]) = [-1, 1, 1, -1];
            row = row + 1;
            tp(j).kKca = row;
            K(row, tp(j).iKca) = 1;
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
        tp(j).jTA = row_alk;
        M(row_alk,[ iTA, tp(j).ihco3, tp(j).ico3, tp(j).ioh ]) = ...
            [1, -1, -2, -1];

        % Total carbonate: TC - [CO2*] - [HCO3] - [CO3] = 0
        row = row + 1;
        M(row, [ iTC, tp(j).ico2st, tp(j).ihco3, tp(j).ico3 ]) = ...
            [1, -1, -1, -1];
        tp(j).jTC = row;
        
        % Total borate
        row = row + 1;
        M(row,[ iTB, tp(j).iboh3, tp(j).iboh4 ])   =  [1, -1, -1];
        M(row_alk, tp(j).iboh4) = -1;
        tp(j).jTB = row;
        
        % Total sulfate
        row = row + 1;
        M(row, [ iTS, tp(j).ihso4, tp(j).iso4 ])   =  [1, -1, -1];
        M(row_alk,[ tp(j).iph_free, tp(j).ihso4 ])  =  [1, 1]; 
        tp(j).jTS = row;

        % Total fluoride
        row = row + 1;
        M(row, [ iTF, tp(j).iHF, tp(j).iF ]) =  [1, -1, -1];
        M(row_alk, tp(j).iHF) =  1;
        tp(j).jTF = row;
        
        if (ismember('phosphate',abr))
            row = row + 1;
            M(row, [ iTP,tp(j).ih3po4, tp(j).ih2po4, tp(j).ihpo4, ...
                tp(j).ipo4 ]) = [1, -1, -1, -1, -1];
            M(row_alk, [ tp(j).ihpo4, tp(j).ipo4, tp(j).ih3po4 ]) ...
                = [-1, -2, 1]; 
            tp(j).jTP = row;
        end
        if (ismember('silicate',abr))
            row = row + 1;
            M(row, [ iTSi, tp(j).isioh4, tp(j).isiooh3]) = [1, -1, -1];
            M(row_alk, tp(j).isiooh3) = -1; 
            tp(j).jTSi = row;
        end
        if (ismember('ammonia',abr))
            row = row+1;
            M(row, [ iTNH3, tp(j).inh4, tp(j).inh3 ]) = [1, -1, -1];
            M(row_alk,tp(j).inh3) = -1; 
            tp(j).jTNH3 = row;
        end
        if (ismember('sulfide',abr)) % H2S
            row = row+1;
            M(row,[ iTH2S, tp(j).ih2s, tp(j).ihs ]) = [1, -1, -1];
            M(row_alk,tp(j).ihs) = -1; 
            tp(j).jTH2S = row;
        end
        if (ismember('solubility',abr)) % total Calcium
            row = row + 1;
            M(row, [ iTCal, tp(j).ica]) = [1, -1];
            tp(j).jTCal = row;
        end

        % % RESIDUALS of ph_free equations based on input phscale 
        % %     (as functions solved in parse_zpK subfunction)
        
        % % ph_free = ph_tot - factor;
        % % rph_free = ph_tot - factor - ph_free; (r = RESIDUAL)
        
        % % convert to ph_tot with conversion factors FREE2tot & SWS2free
        % % where FREE2tot = 1 + TS/Ks; pFREE2tot = p(TS + Ks) - pKs;
        % % and SWS2free = 1/(1 + TS/Ks + TF/Kf); 
            % pSWS2free = -p( 1 + (TS/Ks) + (TF/Kf) ); 
        pFREE2tot = @(z) ( p( q( z(iTS)) + q( z(tp(j).iKs)) ) ) ...       
                - (z(tp(j).iKs)); % ( Ks/Ks + TS/Ks ) = (Ks + TS)/Ks
        pSWS2free = @(z) -p( 1 + ( q( z(iTS) ) / q( z(tp(j).iKs) ) ) ...   
                + ( q(z(iTF) ) / q( z(tp(j).iKf) ) ) );

        if phscale == 1 % total
            % ph_free = ph_tot - p(FREE2tot);
            tp(j).rph_free = @(z) z(tp(j).iph) - pFREE2tot(z) ...
                - z(tp(j).iph_free);
        elseif phscale == 2 % sws
            % ph_free = ph_sws + (-pSWS2free);
            tp(j).rph_free = @(z) z(tp(j).iph) + pSWS2free(z) ...
                - z(tp(j).iph_free);
        elseif phscale == 3 % free
            % ph_free = ph_free;
            tp(j).rph_free = @(z) z(tp(j).iph) - z(tp(j).iph_free);
        elseif phscale == 4 % nbs
            % ph_free = ph_nbs + (-pSWS2free) + pfH;
            tp(j).rph_free = @(z) z(tp(j).iph) + pSWS2free(z) ...
                - z(tp(j).ipfH) - z(tp(j).iph_free);
        end

        % % calculate rph_free first derivatives grph_free
        if phscale == 2 || phscale == 4 % sws or nbs
            gpSWS2free_outer = @(z) -dpdx( 1 + ...  % outer of chain rule
                ( q( z(iTS) ) / q( z(tp(j).iKs) ) ) + ...
                ( q( z(iTF) ) / q( z(tp(j).iKf) ) ) );
            gpSWS2free_pTS = @(z) gpSWS2free_outer(z) * ...
                ( dqdx(z(iTS) ) / q( z(tp(j).iKs) ) );
            gpSWS2free_pKs = @(z) gpSWS2free_outer(z) * ...
                q( z(iTS) ) * -1 * ( q( z(tp(j).iKs) )^(-2)) * ...
                dqdx( z(tp(j).iKs) );
            gpSWS2free_pTF = @(z) gpSWS2free_outer(z) * ...
                ( dqdx( z(iTF)) / q( z(tp(j).iKf) ) );
            gpSWS2free_pKf = @(z) gpSWS2free_outer(z) * ...
                q( z(iTF) ) * -1 * ( q( z(tp(j).iKf) )^(-2)) * ...
                dqdx( z(tp(j).iKf) );
            if phscale == 2 % sws
                g_pfH = @(z) 0;
                tp(j).grph_free = @(z) [ gpSWS2free_pTS(z); ...
                    gpSWS2free_pKs(z);   gpSWS2free_pTF(z); ...
                    gpSWS2free_pKf(z);            g_pfH(z); ...
                      1              ;                 -1 ];
                    % ^ grph_free_ph | grph_free_phfree ^
            elseif phscale == 4 % nbs
                g_pfH = @(z) -1;
                tp(j).grph_free = @(z) [   gpSWS2free_pTS(z); ...
                    gpSWS2free_pKs(z);     gpSWS2free_pTF(z); ...
                    gpSWS2free_pKf(z);              g_pfH(z); ...
                       1             ;                  -1   ];
                    %  ^ grph_free_ph| grph_free_phfree ^
            end
        elseif phscale == 1 % tot
            gpF2t_outside = @(z) -dpdx( q( z(iTS)) + q( z(tp(j).iKs)) );
            gpF2t_pTS = @(z) gpF2t_outside(z) * dqdx( z(iTS) );
            gpF2t_pKs = @(z) gpF2t_outside(z) * dqdx( z(tp(j).iKs) ) + 1;
            tp(j).grph_free = @(z) [ gpF2t_pTS(z); gpF2t_pKs(z); ...
                0       ; 0     ; ...   % g_pTF; g_pKf;
                0       ; 1     ; -1 ]; % g_pfH; g_ph; g_phfree;
        elseif phscale == 3 % free
            tp(j).grph_free = @(z) [ 0; 0; 0; 0; 0; 1; -1 ];
            % g_pTS; g_pKs; g_pTF; g_pKf; g_pfH; g_ph; g_phfree;
        end


    end
    sys.M  = M;
    sys.K  = K;
    sys.tp = tp;   
end

