import numpy as np

def init(opt, yobs, sys):
    q = sys['q']
    p = sys['p']
    
    y0 = yobs.copy()
    
    dic = q(yobs[sys['ipTC']])
    alk = q(yobs[sys['ipTA']])
    
    if np.isnan(dic):
        dic = 2200e-6
        y0[sys['ipTC']] = p(dic)
    
    if np.isnan(alk):
        alk = 2200e-6
        y0[sys['ipTA']] = p(alk)
    
    gam = dic / alk

    TB = q(yobs[sys['ipTB']])
    y0[sys['ipTB']] = p(TB)

    TS = q(yobs[sys['ipTS']])
    y0[sys['ipTS']] = p(TS)

    TF = q(yobs[sys['ipTF']])
    y0[sys['ipTF']] = p(TF)

    TP = q(yobs[sys['ipTP']])
    y0[sys['ipTP']] = p(TP)

    TSi = q(yobs[sys['ipTSi']])
    y0[sys['ipTSi']] = p(TSi)

    TNH4 = q(yobs[sys['ipTNH4']])
    y0[sys['ipTNH4']] = p(TNH4)

    TH2S = q(yobs[sys['ipTH2S']])
    y0[sys['ipTH2S']] = p(TH2S)

    TCa = q(yobs[sys['ipTCa']])
    y0[sys['ipTCa']] = p(TCa)
    
    nTP = len(sys['tp'])
    
    for i in range(nTP):
        K0 = q(y0[sys['tp'][i]['ipK0']])
        K1 = q(y0[sys['tp'][i]['ipK1']])
        K2 = q(y0[sys['tp'][i]['ipK2']])
        
        h = 0.5 * ((gam - 1) * K1 + np.sqrt((1 - gam) ** 2 * K1 ** 2 - 4 * K1 * K2 * (1 - 2 * gam)))
        
        hco3 = h * alk / (h + 2 * K2)
        co2st = h * hco3 / K1
        co3 = dic * K1 * K2 / (K1 * h + h ** 2 + K1 * K2)
        fco2 = co2st / K0
        
        y0[sys['tp'][i]['iph']] = p(h)
        y0[sys['tp'][i]['iphco3']] = p(hco3)
        y0[sys['tp'][i]['ipco2st']] = p(co2st)
        y0[sys['tp'][i]['ipco3']] = p(co3)
        y0[sys['tp'][i]['ipfco2']] = p(fco2)
        
        Kb = q(y0[sys['tp'][i]['ipKb']])
        boh4 = TB * Kb / (Kb + h)
        boh3 = TB - boh4
        y0[sys['tp'][i]['ipboh3']] = p(boh3)
        y0[sys['tp'][i]['ipboh4']] = p(boh4)
        
        Kw = q(y0[sys['tp'][i]['ipKw']])
        oh = Kw / h
        y0[sys['tp'][i]['ipoh']] = p(oh)
        
        Ks = q(y0[sys['tp'][i]['ipKs']])
        
        h_tot = h * (1 + TS / Ks)
        h_free = h
        hso4 = h_tot - h_free
        so4 = Ks * hso4 / h
        y0[sys['tp'][i]['iphso4']] = p(hso4)
        y0[sys['tp'][i]['ipso4']] = p(so4)
        
        Kf = q(y0[sys['tp'][i]['ipKf']])
        HF = TF / (1 + Kf / h_free)
        F = Kf * HF / h_free
        h_sws = h_tot + HF
        y0[sys['tp'][i]['ipF']] = p(F)
        y0[sys['tp'][i]['ipHF']] = p(HF)
        
        Kp1 = q(y0[sys['tp'][i]['ipKp1']])
        Kp2 = q(y0[sys['tp'][i]['ipKp2']])
        Kp3 = q(y0[sys['tp'][i]['ipKp3']])
        d = h ** 3 + Kp1 * h ** 2 + Kp1 * Kp2 * h + Kp1 * Kp2 * Kp3
        h3po4 = TP * h ** 3 / d
        h2po4 = TP * Kp1 * h ** 2 / d
        hpo4 = TP * Kp1 * Kp2 * h / d
        po4 = TP * Kp1 * Kp2 * Kp3 / d
        y0[sys['tp'][i]['iph3po4']] = p(h3po4)
        y0[sys['tp'][i]['iph2po4']] = p(h2po4)
        y0[sys['tp'][i]['iphpo4']] = p(hpo4)
        y0[sys['tp'][i]['ippo4']] = p(po4)
        
        Ksi = q(y0[sys['tp'][i]['ipKsi']])
        siooh3 = TSi / (1 + h / Ksi)
        sioh4 = TSi - siooh3
        y0[sys['tp'][i]['ipsiooh3']] = p(siooh3)
        y0[sys['tp'][i]['ipsioh4']] = p(sioh4)
        
        Knh4 = q(y0[sys['tp'][i]['ipKnh4']])
        nh3 = TNH4 / (1 + h / Knh4)
        nh4 = TNH4 - nh3
        y0[sys['tp'][i]['ipnh3']] = p(nh3)
        y0[sys['tp'][i]['ipnh4']] = p(nh4)
        
        Kh2s = q(y0[sys['tp'][i]['ipKh2s']])
        hs = TH2S / (1 + h / Kh2s)
        h2s = TH2S - hs
        y0[sys['tp'][i]['ipHS']] = p(hs)
        y0[sys['tp'][i]['ipH2S']] = p(h2s)
        
        p2f = q(y0[sys['tp'][i]['ipp2f']])
        pco2 = fco2 / p2f
        y0[sys['tp'][i]['ippco2']] = p(pco2)
        
        Kar = q(y0[sys['tp'][i]['ipKar']])
        OmegaAr = co3 * TCa / Kar
        Kca = q(y0[sys['tp'][i]['ipKca']])
        OmegaCa = co3 * TCa / Kca
        y0[sys['tp'][i]['ipca']] = p(TCa)
        y0[sys['tp'][i]['ipOmegaAr']] = p(OmegaAr)
        y0[sys['tp'][i]['ipOmegaCa']] = p(OmegaCa)
        
        y0[sys['tp'][i]['iph_tot']] = p(h_tot)
        y0[sys['tp'][i]['iph_sws']] = p(h_sws)
        y0[sys['tp'][i]['iph_free']] = p(h_free)
        fH = q(y0[sys['tp'][i]['ipfH']])
        h_nbs = h_free * fH
        y0[sys['tp'][i]['iph_nbs']] = p(h_nbs)
    
    nlam = sys['M'].shape[0] + sys['K'].shape[0]
    lam = np.zeros(nlam)
    z0 = np.concatenate([y0, lam])
    
    return z0
