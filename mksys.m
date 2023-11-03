function sys = mksys(obs,phscale)
%
% Private function for QUODcarb.m
% it creates the K and M matrices
%
% utility functions and constants
    LOG10   = log(10);
    p       = @(x) -log10( x );  % inverse of q    
    q       = @(x)  10.^( -x );  % inverse of p
    sys.p   = p;
    sys.q   = q;
    dqdx    = @(x) - LOG10 * 10.^( -x );  % q'
    d2qdx2  = @(x) LOG10^2 * 10.^(-x ); % q"
    dpdx    = @(x) -1 / (x .* LOG10);     % p'
    d2pdx2  = @(x) 1 / (x.^2 * LOG10);  % p"
    sys.dqdx    = dqdx;
    sys.dpdx    = dpdx;
    sys.d2pdx2  = d2pdx2;
    sys.d2qdx2  = d2qdx2;
    isgood = @(thing) (~sum(isnan(thing)) & ~isempty(thing));
    if (~isfield(obs,'tp'))
        error('Need to provide temperature and pressure measurement.')
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
    % carbonate:
    isal = 1;  sys.isal  = isal;   % salinity
    ipTC = 2;  sys.ipTC  = ipTC;   % p(total carbonate)
    ipTA = 3;  sys.ipTA  = ipTA;   % p(total alkalinity)

    i = 3;

    % borate: Kb = [h][boh4]/[boh3]
    i = i + 1;
    ipTB = i;  sys.ipTB = ipTB; % p(total borate)
        
    % sulfate: Ks = [hf][so4]/[hso4]
    i = i + 1;
    ipTS = i;  sys.ipTS = ipTS; % p(total sulfate)
    
    % fluoride: Kf = [h][F]/[HF]
    i = i + 1;
    ipTF = i;  sys.ipTF = ipTF; % p(total fluoride)
    
    % phosphate: Kp1 = [h][h2po4]/[h3po4], Kp2 = [h][hpo4]/[h2po4], Kp3 = [h][po4]/[hpo4]
    i = i+1;
    ipTP = i; sys.ipTP = ipTP; % p(total phosphate)

    % KSi = [h][siooh3]/[sioh4]
    i = i+1;
    ipTSi = i; sys.ipTSi = ipTSi; % p(total silicate)
    
    % Knh4 = [h][nh3]/[nh4]
    i = i+1;
    ipTNH4 = i; sys.ipTNH4 = ipTNH4; % p(total amonia)
    
    % Kh2s = [h][hs]/[h2s]
    i = i+1;
    ipTH2S = i; sys.ipTH2S = ipTH2S; % p(total sulfide)
    
    % Kar = [co3][ca]/OmegaAr
    % Kca = [co3][ca]/OmegaCa
    i = i+1;
    ipTCa = i; sys.ipTCa = ipTCa ; % p(total calcium)

    nTP = length(obs.tp); % number of different (T,P) sub systems
    for j = 1:nTP % loop over (T,P) sub systems
        i = i + 1;    iT = i;     tp(j).iT = iT;  % temperature
        i = i + 1;    iP = i;     tp(j).iP = iP;  % pressure
            
        % K0 = [co2st]/fco2   
        % K1 = [h][hco3]/[co2st]
        % K2 = [h][co3]/[hco3]
        nrk = 3;              i = i + 1; 
        tp(j).ipK0      = i;  i = i + 1; 
        tp(j).ipK1      = i;  i = i + 1;
        tp(j).ipK2      = i;  i = i + 1;
        tp(j).ipfco2    = i;  i = i + 1;
        tp(j).ipco2st   = i;  i = i + 1;
        tp(j).iphco3    = i;  i = i + 1;
        tp(j).ipco3     = i;  i = i + 1;
        tp(j).iph       = i; 
        
        % Kb = [h][boh4]/[boh3]
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKb   = i;     i = i + 1;
        tp(j).ipboh4 = i;     i = i + 1;
        tp(j).ipboh3 = i;

        % Kw = [h][oh] 
        nrk = nrk + 1;       i = i + 1;
        tp(j).ipKw = i;      i = i + 1;
        tp(j).ipoh = i;      
     
        % Ks  = [hf][so4]/[hso4]
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKs     = i;   i = i + 1;
        tp(j).ipso4    = i;   i = i + 1;
        tp(j).iphso4   = i;
        
        % Kf = [h][F]/[HF]
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKf = i;       i = i + 1;
        tp(j).ipF  = i;       i = i + 1;
        tp(j).ipHF = i;
        % Kp1 = [h][h2po4]/[h3po4]
        % Kp2 = [h][hpo4]/[h2po4]
        % Kp3 = [h][po4]/[hpo4]
        nrk = nrk + 3;        i = i + 1;
        tp(j).ipKp1   = i;    i = i + 1;
        tp(j).ipKp2   = i;    i = i + 1;
        tp(j).ipKp3   = i;    i = i + 1;
        tp(j).iph3po4 = i;    i = i + 1;
        tp(j).iph2po4 = i;    i = i + 1;
        tp(j).iphpo4  = i;    i = i + 1;
        tp(j).ippo4   = i;
        
        % KSi = [h][siooh3]/[sioh4]
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKsi    = i;   i = i + 1;
        tp(j).ipsiooh3 = i;   i = i + 1;
        tp(j).ipsioh4  = i;
        
        % Knh4 = [h][nh3]/[nh4]
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKnh4 = i;     i = i + 1;
        tp(j).ipnh3  = i;     i = i + 1;
        tp(j).ipnh4  = i;
        
        % Kh2s = [h][hs]/[h2s]
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKh2s = i;     i = i + 1;
        tp(j).ipHS   = i;     i = i + 1;
        tp(j).ipH2S  = i;    

        % fco2 = pco2 * p2f;
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipp2f     = i;  i = i + 1;
        tp(j).ippco2    = i;  
        
        % Kar = [co3][ca]/OmegaAr
        % Kca = [co3][ca]/OmegaCa
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKar     = i;  i = i + 1;
        tp(j).ipca      = i;  i = i + 1;
        tp(j).ipOmegaAr = i;
        nrk = nrk + 1;        i = i + 1;
        tp(j).ipKca     = i;  i = i + 1;
        tp(j).ipOmegaCa = i;  i = i+1;
        tp(j).ipfH      = i; i = i + 1;
        
        %  ph scales 
        tp(j).iph_tot  = i;   i = i + 1;
        tp(j).iph_free = i;   i = i + 1;
        tp(j).iph_sws  = i;   i = i + 1;
        tp(j).iph_nbs  = i;
    end
    nv = i;
    K = sparse(nTP*(nrk+2),nv);
    row = 0;
    for j = 1:nTP
        kr = [];
        kc = [];
        % K0 = [CO2*]/fCO2 ==> -pK0 + pco2st - pfco2 = 0         
        row = row + 1;
        K(row,[ tp(j).ipK0, tp(j).ipco2st, tp(j).ipfco2]) = [ -1, 1, -1 ];
        kc = union(kc,[tp(j).ipK0,tp(j).ipco2st, tp(j).ipfco2]);
        kr = [kr, row];
        
        % K1 = [HCO3][H]/[CO2*] ==> -pK1 + phco3 + ph - pco2st = 0
        row = row + 1;
        K(row,[ tp(j).ipK1, tp(j).iph, tp(j).iphco3, tp(j).ipco2st]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[tp(j).ipK1, tp(j).iph, tp(j).iphco3, tp(j).ipco2st]);
        kr = [kr, row];
        
        % K2 = [CO3][H]/[HCO3]
        row = row + 1;
        K(row,[ tp(j).ipK2, tp(j).iph, tp(j).ipco3, tp(j).iphco3 ]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[tp(j).ipK2, tp(j).iph, tp(j).ipco3, tp(j).iphco3]);
        kr = [kr, row];
        
        % Kb = [H][BOH4]/[BOH3]
        row = row+1;
        K(row,[ tp(j).ipKb, tp(j).iph, tp(j).ipboh4, tp(j).ipboh3 ]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[tp(j).ipKb, tp(j).iph, tp(j).ipboh4, tp(j).ipboh3]);
        kr = [kr, row];

        % Kw = [OH][H]
        row = row + 1;
        K(row,[ tp(j).ipKw, tp(j).iph, tp(j).ipoh ]) = [ -1, 1, 1 ];
        kc = union(kc,[tp(j).ipKw, tp(j).iph, tp(j).ipoh]);
        kr = [kr, row];

        % Ks  = [H]free[SO4]/[HSO4] 
        row = row+1;
        K(row,[ tp(j).ipKs, tp(j).iph_free, tp(j).ipso4, tp(j).iphso4 ]) = [ -1, 1, 1, -1 ];
        kc = union(kc, [tp(j).ipKs, tp(j).iph_free, tp(j).ipso4, tp(j).iphso4]);
        kr = [kr, row];
        
        % Kf = [H]free[F]/[HF]     
        row = row+1;
        K(row,[ tp(j).ipKf, tp(j).iph_free, tp(j).ipF, tp(j).ipHF ]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKf, tp(j).iph_free, tp(j).ipF, tp(j).ipHF]);
        kr = [kr, row];
        
        % Kp1 = [H][H2PO4]/[H3PO4]
        row = row+1;
        K(row,[ tp(j).ipKp1, tp(j).iph, tp(j).iph2po4, tp(j).iph3po4]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKp1, tp(j).iph, tp(j).iph2po4, tp(j).iph3po4]);
        kr = [kr, row];
        
        % Kp2 = [H][HPO4]/[H2PO4]
        row = row + 1;
        K(row,[ tp(j).ipKp2, tp(j).iph, tp(j).iphpo4, tp(j).iph2po4 ]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKp2, tp(j).iph, tp(j).iphpo4, tp(j).iph2po4] );
        kr = [kr, row];
        
        % Kp3 = [H][PO4]/[HPO4]        
        row = row + 1;
        K(row,[ tp(j).ipKp3, tp(j).iph, tp(j).ippo4,tp(j).iphpo4 ]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKp3, tp(j).iph, tp(j).ippo4, tp(j).iphpo4]);
        kr = [kr, row];       
            
        % KSi = [H][SiO(OH)3]/[Si(OH)4]
        row = row + 1;
        K(row,[ tp(j).ipKsi, tp(j).iph, tp(j).ipsiooh3, tp(j).ipsioh4 ]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKsi, tp(j).iph, tp(j).ipsiooh3, tp(j).ipsioh4]);
        kr = [kr, row];
        
        % Knh4 = [H][NH3]/[NH4+]
        row = row + 1;
        K(row,[ tp(j).ipKnh4, tp(j).iph, tp(j).ipnh3, tp(j).ipnh4]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKnh4, tp(j).iph, tp(j).ipnh3, tp(j).ipnh4]);
        kr = [kr, row];
        
        % Kh2s = [H][HS]/[H2S]
        row = row + 1;
        K(row,[ tp(j).ipKh2s, tp(j).iph, tp(j).ipHS, tp(j).ipH2S]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKh2s, tp(j).iph, tp(j).ipHS, tp(j).ipH2S]);
        kr = [kr, row];
            
        % fco2 = pco2 * p2f;
        row = row + 1;
        K(row,[ tp(j).ipfco2, tp(j).ippco2, tp(j).ipp2f ]) = [ -1 1 1 ];
        kc = union(kc,[ tp(j).ipfco2, tp(j).ippco2, tp(j).ipp2f]);
        kr = [kr, row];
          
        % Kar = [co3][ca]/OmegaAr ==> -pKar + pco3 + pca - pOmegaAr = 0
        row = row + 1;
        K(row,[ tp(j).ipKar, tp(j).ipco3, tp(j).ipca, tp(j).ipOmegaAr]) = [ -1, 1, 1, -1 ];
        kc = union(kc, [ tp(j).ipKar, tp(j).ipco3, tp(j).ipca, tp(j).ipOmegaAr]);
        kr = [kr, row];
        
        % Kca = [co3][ca]/OmegaCa ==> -pKca + pco3 + pca - pOmegaCa = 0
        row = row + 1;
        K(row, [ tp(j).ipKca, tp(j).ipco3, tp(j).ipca, tp(j).ipOmegaCa]) = [ -1, 1, 1, -1 ];
        kc = union(kc,[ tp(j).ipKca, tp(j).ipco3, tp(j).ipca, tp(j).ipOmegaCa]);
        kr = [kr, row];
        
        row = row + 1;
        switch phscale % working ph to phscale
          case 1
            K(row,[ tp(j).iph  tp(j).iph_tot]) = [ -1 1 ];
            kc = union(kc,[ tp(j).iph, tp(j).iph_tot]);
          case 2
            K(row,[ tp(j).iph tp(j).iph_sws] ) = [ -1 1 ];
            kc = union(kc,[ tp(j).iph, tp(j).iph_sws]);
          case 3
            K(row,[ tp(j).iph tp(j).iph_free] ) = [ -1 1 ];
            kc = union(kc,[ tp(j).iph, tp(j).iph_free]);
          case 4
            K(row, [tp(j).iph tp(j).iph_nbs] ) = [ -1 1 ];
            kc = union(kc,[ tp(j).iph, tp(j).iph_nbs]);
        end        
        tp(j).kphscale = row;
        kr = [kr, row];

        % def'n for nbs
        row = row + 1;
        K(row,[ tp(j).iph_nbs, tp(j).iph_sws, tp(j).ipfH ]) = [ 1 -1 -1 ];
        kc = union(kc, [ tp(j).iph_nbs, tp(j).iph_sws, tp(j).ipfH]);            
        kr = [kr, row];
        tp(j).kr = kr;
        tp(j).kc = kc;
        nr = 12; % TA, TC, TB, TS, TF, TP, TSi TNH4 TH2S TCa ph_tot ph_sws ph_nbs
    end
    % "mass conservation" equations
    M = sparse(nTP*nr,nv);
    row = 0;

    for j = 1:nTP
        mr = [];
        mc = [];

        % Total alkalinity: TA - [HCO3] - 2[CO3] ( - [OH] )...
        row = row + 1;
        row_alk = row;
        mr = [mr, row];

        % carbonate 
        M(row_alk,[ ipTA, tp(j).iphco3, tp(j).ipco3]) = [ 1, -1, -2 ];
        mc = union(mc,[ipTA, tp(j).iphco3, tp(j).ipco3]);
        
        % Total carbonate: TC - [CO2*] - [HCO3] - [CO3] = 0
        row = row + 1;
        M(row, [ ipTC, tp(j).ipco2st, tp(j).iphco3, tp(j).ipco3 ]) = [ 1, -1, -1, -1 ];
        mc = union(mc,[ipTC, tp(j).ipco2st, tp(j).iphco3, tp(j).ipco3]);
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e2;
        
        % Total borate
        row = row + 1;
        M(row,[ ipTB, tp(j).ipboh3, tp(j).ipboh4 ])   =  [ 1, -1, -1 ];
        mc = union(mc,[ipTB, tp(j).ipboh3, tp(j).ipboh4]);
        M(row_alk, tp(j).ipboh4) = -1;
        M(row_alk, tp(j).ipoh) = -1;        
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e3;
        
        % Total sulfate
        row = row + 1;
        M(row, [ ipTS, tp(j).iphso4, tp(j).ipso4 ])   =  [ 1, -1, -1 ];
        M(row_alk,[ tp(j).iph_free, tp(j).iphso4 ])   =  [1, 1]; 
        mc = union(mc,[ipTS, tp(j).iphso4, tp(j).ipso4]);
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e1;
        
        % Total fluoride
        row = row + 1;
        M(row, [ ipTF, tp(j).ipHF, tp(j).ipF ]) =  [ 1, -1, -1 ];
        mc = union(mc,[ipTF, tp(j).ipHF, tp(j).ipF]);
        M(row_alk, tp(j).ipHF) =  1;
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e3;
        
        % Total phosphate
        row = row + 1;
        M(row, [ ipTP, tp(j).iph3po4, tp(j).iph2po4, tp(j).iphpo4, tp(j).ippo4 ]) = [ 1, -1, -1, -1, -1 ];
        mc = union(mc,[ipTP, tp(j).iph3po4, tp(j).iph2po4, tp(j).iphpo4, tp(j).ippo4]);
        M(row_alk, [ tp(j).iphpo4, tp(j).ippo4, tp(j).iph3po4 ]) = [ -1, -2, 1 ]; 
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e7;
        
        % Total silicate
        row = row + 1;
        M(row, [ ipTSi, tp(j).ipsioh4, tp(j).ipsiooh3]) = [ 1, -1, -1 ];
        mc = union(mc,[ipTSi, tp(j).ipsioh4, tp(j).ipsiooh3]);
        M(row_alk, tp(j).ipsiooh3) = -1; 
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e6;
        
        % Total amonia
        row = row+1;
        M(row, [ ipTNH4, tp(j).ipnh4, tp(j).ipnh3 ]) = [ 1, -1, -1 ];
        mc = union(mc,[ipTNH4, tp(j).ipnh4, tp(j).ipnh3]);
        M(row_alk,tp(j).ipnh3) = -1; 
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e8;
        
        % Total sulfide
        row = row+1;
        M(row,[ ipTH2S, tp(j).ipH2S, tp(j).ipHS ]) = [ 1, -1, -1 ];
        mc = union(mc,[ipTH2S, tp(j).ipH2S, tp(j).ipHS]);
        M(row_alk,tp(j).ipHS) = -1; 
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e8;
        M(row_alk,:) = M(row_alk,:)*1e2;
        
        % total Ca
        row = row + 1;
        M(row, [ ipTCa, tp(j).ipca]) = [ 1, -1 ];
        mc = union(mc,[ipTCa, tp(j).ipca]);
        mr = [mr,row];
        M(row,:) = M(row,:)*1e2;
        
        % ph_tot and ph_free relationship
        row = row + 1;
        M(row, [ tp(j).iph_tot, tp(j).iph_free, tp(j).iphso4 ] ) = [ 1 -1 -1 ];
        mc = union(mc,[tp(j).iph_tot, tp(j).iph_free, tp(j).iphso4]);
        mr = [mr,row];
        % rescale row
        M(row,:) = M(row,:)*1e8;

        % ph_sws and ph_free relationship
        row = row + 1;
        M(row,[ tp(j).iph_sws, tp(j).iph_free, tp(j).iphso4, tp(j).ipHF ]) = [ 1 -1 -1 -1 ];
        mc = union(mc,[tp(j).iph_sws, tp(j).iph_free, tp(j).iphso4, tp(j).ipHF]);
        mr = [mr,row];
        tp(j).mr = mr;
        tp(j).mc = mc.';
        % rescale row
        M(row,:) = M(row,:)*1e8;

    end
    for j = 1:nTP
        % STUFF needed to compute the Revelle buffer factor
        ifixed = [ ipTA, ipTB, ipTS, ipTF, ipTP, ipTSi, ipTNH4, ipTH2S, ipTCa, isal ];
        ifixed = [ifixed, tp(j).ipK0,  tp(j).ipK1,  tp(j).ipK2,  tp(j).ipKb,  tp(j).ipKw, tp(j).ipKs,  tp(j).ipKf, ...
                 tp(j).ipKp1, tp(j).ipKp2, tp(j).ipKp3, tp(j).ipKsi, tp(j).ipKnh4, tp(j).ipKh2s,...
                 tp(j).ipp2f, tp(j).ipKar, tp(j).ipKca, tp(j).ipfH,  tp(j).iP, tp(j).iT];
        tp(j).ifixed = ifixed;
        tp(j).ifree = setdiff(union(tp(j).mc,tp(j).kc),ifixed);

        ML = M(tp(j).mr,:); KL = K(tp(j).kr,:);
        ML = ML(:,tp(j).ifree); KL = KL(:,tp(j).ifree);
        
        %tp(j).c = @(xfree,zhat) [ML*q(xfree); KL*xfree] + [MR*q(zhat(ifixed)); KR*zhat(ifixed)];
        tp(j).dcdx_pTAfixed = @(xfree) [ML*diag(dqdx(xfree));KL];
                

        % STUFF needed to compute a dpfco2/dpTA holding pTC fixed 
        jfixed = [ ipTC, ipTB, ipTS, ipTF, ipTP, ipTSi, ipTNH4, ipTH2S, ipTCa, isal ];
        jfixed = [jfixed, tp(j).ipK0,  tp(j).ipK1,  tp(j).ipK2,  tp(j).ipKb,  tp(j).ipKw, tp(j).ipKs,  tp(j).ipKf, ...
                 tp(j).ipKp1, tp(j).ipKp2, tp(j).ipKp3, tp(j).ipKsi, tp(j).ipKnh4, tp(j).ipKh2s,...
                 tp(j).ipp2f, tp(j).ipKar, tp(j).ipKca, tp(j).ipfH,  tp(j).iP, tp(j).iT];
        tp(j).jfixed = jfixed;
        tp(j).jfree = setdiff(union(tp(j).mc,tp(j).kc),jfixed);

        ML = M(tp(j).mr,:); KL = K(tp(j).kr,:);
        ML = ML(:,tp(j).jfree); KL = KL(:,tp(j).jfree);
        tp(j).dcdx_pTCfixed = @(xfree) [ML*diag(dqdx(xfree));KL];
        

    end
    sys.M  = M;
    sys.K  = K;
    sys.tp = tp;
end

    

