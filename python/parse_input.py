import calc_pTOT as ct
import calc_pk as cp
import numpy as np

def parse_input(obs,sys,opt,nD):
    def isgood(thing):
        return (thing is not None) and not np.any(np.isnan(thing))
    p = sys['p']
    q = sys['q']
    # convert x+/-e into precision for p(x) (precision = 1/variance)
    def w(x, e):
        return np.abs(p(1 + e / x)) ** (-2)
    
    nv = sys['K'].shape[1]
    yobs = np.full((nD, nv), np.nan)
    wobs = np.full((nD, nv), np.nan)

    if 'tp' not in obs:
        if opt['printmes'] != 0:
            raise ValueError('Need to provide temperature and pressure measurement.')
    
    nTP = len(obs[0]['tp'])
    for i in range(nD): # loop over all the stations
        if 'sal' not in obs[i]:
            # make sure all the required fields in the obs struct exist
            if opt['printmes'] != 0:
                raise ValueError('Need to provide salinity measurement.')
        else:
            yobs[i, sys['isal']] = obs[i]['sal']

        if 'esal' not in obs[i]:
            obs[i]['esal'] = 0.002  # 1std = 0.002 PSU
            wobs[i, sys['isal']] = obs[i]['esal'] ** (-2)
            if opt['printmes'] != 0:
                print('Warning: Assuming salinity uncertainty is 0.002 PSU')
        else:
            wobs[i, sys['isal']] = obs[i]['esal'] ** (-2) # std e -> w
        if 'TC' not in obs[i] or not isgood(obs[i]['TC']):
            obs[i]['TC'] = np.nan
            yobs[i, sys['ipTC']] = np.nan
        else:
            yobs[i, sys['ipTC']] = p(obs[i]['TC'] * 1e-6) # convt to mol/kg

        if 'eTC' not in obs[i] or not isgood(obs[i]['eTC']):
            obs[i]['eTC'] = np.nan
            wobs[i, sys['ipTC']] = np.nan
        else:
            wobs[i, sys['ipTC']] = w(obs[i]['TC'], obs[i]['eTC']) # std e -> w
        if 'TA' not in obs[i] or not isgood(obs[i]['TA']):
            obs[i]['TA'] = np.nan
            yobs[i, sys['ipTA']] = np.nan
        else:
            yobs[i, sys['ipTA']] = p(obs[i]['TA'] * 1e-6) # convt to mol/kg
        if 'eTA' not in obs[i] or not isgood(obs[i]['eTA']):
            obs[i]['eTA'] = np.nan
            wobs[i, sys['ipTA']] = np.nan
        else:
            wobs[i, sys['ipTA']] = w(obs[i]['TA'], obs[i]['eTA']) # std e -> w
        # calculate totals that are a function of salinity
        pT, _, _, epT = ct.calc_pTOT(opt, obs[i]['sal'])
        # (see Ref's within calc_pTOT)
        pTB, TB, epTB = pT[0], q(pTB[0]) * 1e6, epT[0]
        pTS, TS, epTS = pT[1], q(pTS[1]) * 1e6, epT[1]
        pTF, TF, epTF = pT[2], q(pTF[2]) * 1e6, epT[2]
        pTCa, TCa, epTCa = pT[3], q(pTCa[3]) * 1e6, epT[3]
        # total borate
        if 'TB' not in obs[i] or not isgood(obs[i]['TB']):
            obs[i]['TB'] = np.nan
            yobs[i, sys['ipTB']] = pTB
        else:
            if obs[i]['TB'] == 0:
                obs[i]['TB'] = 1e-3  # umol/kg, reset minimum to 1 nanomolar
            yobs[i, sys['ipTB']] = p(obs[i]['TB'] * 1e-6) # convt µmol/kg to mol/kg
        if 'eTB' not in obs[i] or not isgood(obs[i]['eTB']):
            obs[i]['eTB'] = np.nan
            wobs[i, sys['ipTB']] = epTB ** (-2)  # convert to precision
        else:
            wobs[i, sys['ipTB']] = w(obs[i]['TB'], obs[i]['eTB']) # µmol/kg        

        # total sulfate
        if 'TS' not in obs[i] or not isgood(obs[i]['TS']):
            obs[i]['TS'] = np.nan
            yobs[i, sys['ipTS']] = pTS
        else:
            if obs[i]['TS'] == 0:
                obs[i]['TS'] = 1e-3 # µmol/kg reset minimum to 1 nanomolar
            yobs[i, sys['ipTS']] = p(obs[i]['TS'] * 1e-6) # mol/kg
        
        if 'eTS' not in obs[i] or not isgood(obs[i]['eTS']):
            obs[i]['eTS'] = np.nan
            wobs[i, sys['ipTS']] = epTS ** (-2) # convert to precision
        else:
            wobs[i, sys['ipTS']] = w(obs[i]['TS'], obs[i]['eTS'])

        # total fluoride
        if 'TF' not in obs[i] or not isgood(obs[i]['TF']):
            obs[i]['TF'] = np.nan
            yobs[i, sys['ipTF']] = pTF
        else:
            if obs[i]['TF'] == 0:
                obs[i]['TF'] = 1e-3 # µmol/kg reset minimum to 1 nanomolar
            yobs[i, sys['ipTF']] = p(obs[i]['TF'] * 1e-6) # convt µmol/kg to mol/kg
        if 'eTF' not in obs[i] or not isgood(obs[i]['eTF']):
            obs[i]['eTF'] = np.nan
            wobs[i, sys['ipTF']] = epTF ** (-2) 
        else:
            wobs[i, sys['ipTF']] = w(obs[i]['TF'], obs[i]['eTF'])
        
        # total phosphate
        if 'TP' not in obs[i] or not isgood(obs[i]['TP']):
            obs[i]['TP'] = np.nan
            yobs[i, sys['ipTP']] = p(1e-9) # reset minimum (mol/kg)
        else:
            if obs[i]['TP'] == 0: #zero po4 is very unlikely and breaks the code
                obs[i]['TP'] = 1e-3 # µmol/kg reset minimum to 1 nanomolar
            yobs[i, sys['ipTP']] = p(obs[i]['TP'] * 1e-6) # convt µmol/kg to mol/kg
        if 'eTP' not in obs[i] or not isgood(obs[i]['eTP']):
            obs[i]['eTP'] = np.nan
            wobs[i, sys['ipTP']] = w(1e-3, 1e-3)
        else:
            if ((obs[i]['eTP']) == 0):
                obs[i]['eTP'] = 1e-3 # umol/kg, reset minimum if zero
            wobs[i, sys['ipTP']] = w(obs[i]['TP'], obs[i]['eTP']) # mol/kg
        
        # total silicate
        if 'TSi' not in obs[i] or not isgood(obs[i]['TSi']):
            obs[i]['TSi'] = np.nan
            yobs[i, sys['ipTSi']] = p(1e-9) # reset minimum (mol/kg)
        else:
            if obs[i]['TSi'] == 0: #zero silicate very unlikely and breaks code
                obs[i]['TSi'] = 1e-3 # µmol/kg,  reset minimum to 1 nanomolar
            yobs[i, sys['ipTSi']] = p(obs[i]['TSi'] * 1e-6) 
        if 'eTSi' not in obs[i] or not isgood(obs[i]['eTSi']):
            obs[i]['eTSi'] = np.nan
            wobs[i, sys['ipTSi']] = w(1e-3, 1e-3) # mol/kg
            if isgood(obs[i]['TSi']) and opt['printmes'] != 0:
                print("Warning, no obs['eTSi'] input with obs['eTSi']. Assuming 1 nanomolar.")
        else:
            if (obs[i]['eTSi'] == 0):
                obs[i]['eTSi'] = 1e-3 #umol/kg, reset minimum to 1 nanomolar
            wobs[i, sys['ipTSi']] = w(obs[i]['TSi'], obs[i]['eTSi']) # mol/kg
        
        # total amonia
        if 'TNH4' not in obs[i] or not isgood(obs[i]['TNH4']):
            obs[i]['TNH4'] = np.nan
            yobs[i, sys['ipTNH4']] = p(1e-9) # reset minimum (mol/kg)
        else:
            if obs[i]['TNH4'] == 0: 
                obs[i]['TNH4'] = 1e-3; # umol/kg, reset minimum to 1 nanomolar
            yobs[i, sys['ipTNH4']] = p(obs[i]['TNH4'] * 1e-6)
        if 'eTNH4' not in obs[i] or not isgood(obs[i]['eTNH4']):
            obs[i]['eTNH4'] = np.nan
            wobs[i, sys['ipTP']] = w(1e-3, 5e-4) # 5e-4 umol/kg = 1 sigma
            if isgood(obs[i]['TNH4']) and opt['printmes'] != 0:
                print("Warning, no obs['eTNH4'] input with obs['eTNH4']. Assuming 5e-4 umol/kg.")
        else:
            wobs[i, sys['ipTNH4']] = w(obs[i]['TNH4'], obs[i]['eTNH4'])

        # total sulfide
        if 'TH2S' not in obs[i] or not isgood(obs[i]['TH2S']):
            obs[i]['TH2S'] = np.nan
            yobs[i, sys['ipTH2S']] = p(1e-9) # reset minimum (mol/kg)
        else:
            if obs[i]['TH2S'] == 0: 
                obs[i]['TH2S'] = 1e-3; # umol/kg, reset minimum to 1 nanomolar
            yobs[i, sys['ipTH2S']] = p(obs[i]['TH2S'] * 1e-6)
        if 'eTH2S' not in obs[i] or not isgood(obs[i]['eTH2S']):
            obs[i]['eTH2S'] = np.nan
            wobs[i, sys['ipTH2S']] = w(1e-3, 5e-4) # 5e-4 umol/kg = 1 sigma
            if isgood(obs[i]['TH2S']) and opt['printmes'] != 0:
                print("Warning, no obs['eTH2S'] input with obs['eTNH4']. Assuming 5e-4 umol/kg.")
        else:
            wobs[i, sys['ipTH2S']] = w(obs[i]['TH2S'], obs[i]['eTH2S'])
        
        # total calcium
        if 'TCa' not in obs[i] or not isgood(obs[i]['TCa']):
            obs[i]['TCa'] = np.nan
            yobs[i, sys['ipTCa']] = pTCa
        else:
            if obs[i]['TCa'] == 0: 
                obs[i]['TCa'] = 1e-9; # mol/kg, reset minimum to 1 nanomolar
            yobs[i, sys['ipTCa']] = p(obs[i]['TCa']) # assume user input of mol/kg
        if 'eTCa' not in obs[i] or not isgood(obs[i]['eTCa']):
            obs[i]['eTCa'] = np.nan
            wobs[i, sys['ipTCa']] = epTCa ** (-2)
            if isgood(obs[i]['TCa']) and opt['printmes'] != 0:
                print("Warning, no obs['eTCa'] input with obs['eTCa']. Assuming 6e-5 umol/kg.")
        else:
            wobs[i, sys['ipTCa']] = w(obs[i]['TCa'], obs[i]['eTCa'])

        for j in range(nTP):  # loop over (T,P) systems
            yobs[i,sys['tp'][j]['iT']]   = obs[i]['tp'][j]['T']
            yobs[i,sys['tp'][j]['iP']]   = obs[i]['tp'][j]['P']
            wobs[i,sys['tp'][j]['iT']]   = (obs[i]['tp'][j]['eT']) ** (-2)
            wobs[i,sys['tp'][j]['iP']]   = (obs[i]['tp'][j]['eP']) ** (-2)
            pK, _, epK = cp.calc_pK(opt, obs[i]['tp'][j]['T'], obs[i]['sal'], obs[i]['tp'][j]['P']) # T, S, P

            pK0, pK1, pK2 = pK[0], pK[1], pK[2]  
            pKb, pKw, pKs = pK[3], pK[4], pK[5]
            pKf, pKp1, pKp2 = pK[6], pK[7], pK[8]
            pKp3, pKsi, pKnh4 = pK[9], pK[10], pK[11] 
            pKh2s, pp2f, pKar = pK[12], pK[13], pK[14]
            pKca, pfH = pK[15], pK[16]

            epK0, epK1, epK2 = epK[0], epK[1], epK[2]
            epKb, epKw, epKs = epK[3], epK[4], epK[5]
            epKf, epKp1, epKp2 = epK[6], epK[7], epK[8]
            epKp3, epKsi, epKnh4 = epK[9], epK[10], epK[11]
            epKh2s, epp2f, epKar = epK[12], epK[13], epK[14]
            epKca, epfH = epK[15], epK[16]
            
            # add "observations" for the equilibrium constants
            # and transfer from obs struct to yobs and wobs
            
            # co2 solubility and fugacity
            if 'pK0' not in obs[i]['tp'][j] or isgood(obs[i]['tp']['pK0']):
                obs[i]['tp'][j]['pK0'] = np.nan
                yobs[[i], sys['tp'][j]['ipK0']] = pK0
            else:
                yobs[[i], sys['tp'][j]['ipK0']]  = obs[i]['tp'][j]['pK0']
            if 'epK0' not in obs[i]['tp'][j] or isgood(obs[i]['tp']['epK0']):
                obs[i]['tp'][j]['epK0'] = np.nan
                yobs[[i], sys['tp'][j]['ipK0']] = epK0 ** (-2)
            else:
                wobs[[i], sys['tp'][j]['ipK0']]  =(obs[i]['tp'][j]['pK0']) ** (-2)
            if 'co2st' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['co2st']):
                obs[i]['tp'][j]['co2st'] = np.nan
                yobs[[i], sys['tp'][j]['ipco2st']] = np.nan
            else:
                yobs[[i], sys['tp'][j]['ipco2st']] = p(obs[i]['tp'][j]['co2st']*1e-6) # convt µmol/kg to mol/kg
            if 'eco2st' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eco2st']):
                obs[i]['tp'][j]['eco2st'] = np.nan
                wobs[i, sys['tp'][j]['ipco2st']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipco2st']] = w[obs[i]['tp'][j]['co2st'], obs[i]['tp'][j]['eco2st']]
            if 'fco2' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['fco2']):
                obs[i]['tp'][j]['fco2'] = np.nan
                yobs[i, sys['tp'][j]['ipfco2']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipfco2']] = p(obs[i]['tp'][j]['fco2']*1e-6) # convt µatm to atm
            if 'efco2' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['efco2']):
                obs[i]['tp'][j]['efco2'] = np.nan
                wobs[i, sys['tp'][j]['ipfco2']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ifco2']] = w[obs[i]['tp'][j]['fco2'], obs[i]['tp'][j]['efco2']]
            if 'pp2f' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pp2f']):
                obs[i]['tp'][j]['pp2f'] = np.nan
                yobs[i, sys['tp'][j]['ipp2f']] = pp2f
            else:
                yobs[i, sys['tp'][j]['ipp2f']] = obs[i]['tp'][j]['pp2f']
            if 'epp2f' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epp2f']):
                obs[i]['tp'][j]['epp2f'] = np.nan
                wobs[i, sys['tp'][j]['ipp2f']] = (epp2f) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipp2f']] = (obs[i]['tp'][j]['epp2f']) ** (-2)
            if 'pco2' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pco2']):
                obs[i]['tp'][j]['pco2'] = np.nan
                yobs[i, sys['tp'][j]['ippco2']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ippco2']] = p(obs[i]['tp'][j]['pco2']*1e-6) # convt µatm to atm
            if 'epco2' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epco2']):
                obs[i]['tp'][j]['epco2'] = np.nan
                wobs[i, sys['tp'][j]['ippco2']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ippco2']] = w[obs[i]['tp'][j]['pco2'], obs[i]['tp'][j]['epco2']]

            # carbonate system 
            if 'pK1' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pK1']):
                obs[i]['tp'][j]['pK1'] = np.nan
                yobs[i, sys['tp'][j]['ipK1']] = pK1
            else:
                yobs[i, sys['tp'][j]['ipK1']] = obs[i]['tp'][j]['pK1']
            if 'epK1' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epK1']):
                obs[i]['tp'][j]['epK1'] = np.nan
                wobs[i, sys['tp'][j]['ipK1']] = (epK1) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipK1']] = (obs[i]['tp'][j]['epK1']) ** (-2)
            if 'pK2' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pK2']):
                obs[i]['tp'][j]['pK2'] = np.nan
                yobs[i, sys['tp'][j]['ipK2']] = pK2
            else:
                yobs[i, sys['tp'][j]['ipK2']] = obs[i]['tp'][j]['pK2']
            if 'epK2' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epK2']):
                obs[i]['tp'][j]['epK2'] = np.nan
                wobs[i, sys['tp'][j]['ipK2']] = (epK2) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipK2']] = (obs[i]['tp'][j]['epK2']) ** (-2)
            if 'hco3' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['hco3']):
                obs[i]['tp'][j]['hco3'] = np.nan
                yobs[i, sys['tp'][j]['iphco3']] = np.nan
            else:
                yobs[i, sys['tp'][j]['iphco3']] = p(obs[i]['tp'][j]['hco3']*1e-6) #  µmol/kg to mol/kg
            if 'ehco3' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['ehco3']):
                obs[i]['tp'][j]['ehco3'] = np.nan
                wobs[i, sys['tp'][j]['iphco3']] = np.nan
            else:
                wobs[i, sys['tp'][j]['iphco3']] = w[obs[i]['tp'][j]['hco3'], obs[i]['tp'][j]['ehco3']]
            if 'co3' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['co3']):
                obs[i]['tp'][j]['co3'] = np.nan
                yobs[i, sys['tp'][j]['ipco3']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipco3']] = p(obs[i]['tp'][j]['co3']*1e-6) #  µmol/kg to mol/kg
            if 'eco3' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eco3']):
                obs[i]['tp'][j]['eco3'] = np.nan
                wobs[i, sys['tp'][j]['ipco3']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipco3']] = w[obs[i]['tp'][j]['co3'], obs[i]['tp'][j]['eco3']]
            if 'ph' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['ph']):
                obs[i]['tp'][j]['ph'] = np.nan
                yobs[i, sys['tp'][j]['iph']] = np.nan
            else:
                yobs[i, sys['tp'][j]['iph']] = obs[i]['tp'][j]['ph']
            if 'eph' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eph']):
                obs[i]['tp'][j]['eph'] = np.nan
                wobs[i, sys['tp'][j]['iph']] = np.nan
            else:
                wobs[i, sys['tp'][j]['iph']] = (obs[i]['tp'][j]['eph']) ** (-2)
            if 'ph_tot' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['ph_tot']):
                obs[i]['tp'][j]['ph_tot'] = np.nan
                yobs[i, sys['tp'][j]['iph_tot']] = np.nan
            else:
                yobs[i, sys['tp'][j]['iph_tot']] = obs[i]['tp'][j]['ph_tot']
            if 'eph_tot' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eph_tot']):
                obs[i]['tp'][j]['eph_tot'] = np.nan
                wobs[i, sys['tp'][j]['iph_tot']] = np.nan
            else:
                wobs[i, sys['tp'][j]['iph_tot']] = (obs[i]['tp'][j]['eph_tot']) ** (-2)
            if 'ph_free' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['ph_free']):
                obs[i]['tp'][j]['ph_free'] = np.nan
                yobs[i, sys['tp'][j]['iph_free']] = np.nan
            else:
                yobs[i, sys['tp'][j]['iph_free']] = obs[i]['tp'][j]['ph_free']
            if 'eph_free' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eph_free']):
                obs[i]['tp'][j]['eph_free'] = np.nan
                wobs[i, sys['tp'][j]['iph_free']] = np.nan
            else:
                wobs[i, sys['tp'][j]['iph_free']] = (obs[i]['tp'][j]['eph_free']) ** (-2)
            if 'ph_sws' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['ph_sws']):
                obs[i]['tp'][j]['ph_sws'] = np.nan
                yobs[i, sys['tp'][j]['iph_sws']] = np.nan
            else:
                yobs[i, sys['tp'][j]['iph_sws']] = obs[i]['tp'][j]['ph_sws']
            if 'eph_sws' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eph_sws']):
                obs[i]['tp'][j]['eph_sws'] = np.nan
                wobs[i, sys['tp'][j]['iph_sws']] = np.nan
            else:
                wobs[i, sys['tp'][j]['iph_sws']] = (obs[i]['tp'][j]['eph_sws']) ** (-2)
            if 'ph_nbs' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['ph_nbs']):
                obs[i]['tp'][j]['ph_nbs'] = np.nan
                yobs[i, sys['tp'][j]['iph_nbs']] = np.nan
            else:
                yobs[i, sys['tp'][j]['iph_nbs']] = obs[i]['tp'][j]['ph_nbs']
            if 'eph_nbs' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eph_nbs']):
                obs[i]['tp'][j]['eph_nbs'] = np.nan
                wobs[i, sys['tp'][j]['iph_nbs']] = np.nan
            else:
                wobs[i, sys['tp'][j]['iph_nbs']] = (obs[i]['tp'][j]['eph_nbs']) ** (-2)
            if 'pfH' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pfH']):
                obs[i]['tp'][j]['pfH'] = np.nan
                yobs[i, sys['tp'][j]['ipfH']] = pfH
            else:
                yobs[i, sys['tp'][j]['ipFH']] = obs[i]['tp'][j]['pfH']
            if 'epfH' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epfH']):
                obs[i]['tp'][j]['epfH'] = np.nan
                wobs[i, sys['tp'][j]['ipfH']] = (epfH) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipfH']] = (obs[i]['tp'][j]['epfH']) ** (-2)
            
            # water dissociation
            if 'pKw' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pKw']):
                obs[i]['tp'][j]['pKw'] = np.nan
                yobs[i, sys['tp'][j]['ipKw']] = pKw
            else:
                yobs[i, sys['tp'][j]['ipKw']] = obs[i]['tp'][j]['pKw']
            if 'epKw' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epKw']):
                obs[i]['tp'][j]['epKw'] = np.nan
                wobs[i, sys['tp'][j]['ipKw']] = (epKw) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipKw']] = (obs[i]['tp'][j]['epKw']) ** (-2)
            if 'oh' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['oh']):
                obs[i]['tp'][j]['oh'] = np.nan
                yobs[i, sys['tp'][j]['ipoh']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipoh']] = p(obs[i]['tp'][j]['oh'] * 1e-6) # convt µmol/kg to mol/kg
            if 'eoh' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eoh']):
                obs[i]['tp'][j]['eoh'] = np.nan
                wobs[i, sys['tp'][j]['ipoh']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipoh']] = w(obs[i]['tp'][j]['oh'], obs[i]['tp'][j]['eoh']) 
            # borate system
            if 'pKb' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pKb']):
                obs[i]['tp'][j]['pKb'] = np.nan
                yobs[i, sys['tp'][j]['ipKb']] = pKb
            else:
                yobs[i, sys['tp'][j]['ipKb']] = obs[i]['tp'][j]['pKb']
            if 'epKb' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epKb']):
                obs[i]['tp'][j]['epKb'] = np.nan
                wobs[i, sys['tp'][j]['ipKb']] = (epKb) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipKb']] = (obs[i]['tp'][j]['epKb']) ** (-2) 
            if 'boh3' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['boh3']):
                obs[i]['tp'][j]['boh3'] = np.nan
                yobs[i, sys['tp'][j]['ipboh3']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipboh3']] = p(obs[i]['tp'][j]['boh3'] * 1e-6) # convt µmol/kg to mol/kg
            if 'eboh3' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eboh3']):
                obs[i]['tp'][j]['eboh3'] = np.nan
                wobs[i, sys['tp'][j]['ipboh3']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipboh3']] = w(obs[i]['tp'][j]['boh3'], obs[i]['tp'][j]['eboh3']) 
            if 'boh4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['boh4']):
                obs[i]['tp'][j]['boh4'] = np.nan
                yobs[i, sys['tp'][j]['ipboh4']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipboh4']] = p(obs[i]['tp'][j]['boh'] * 1e-6) # convt µmol/kg to mol/kg
            if 'eboh4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eboh4']):
                obs[i]['tp'][j]['eboh4'] = np.nan
                wobs[i, sys['tp'][j]['ipboh4']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipboh4']] = w(obs[i]['tp'][j]['boh4'], obs[i]['tp'][j]['eboh4']) 
            
            # sulfate system
            if 'pKs' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pKs']):
                obs[i]['tp'][j]['pKs'] = np.nan
                yobs[i, sys['tp'][j]['ipKs']] = pKs
            else:
                yobs[i, sys['tp'][j]['ipKs']] = obs[i]['tp'][j]['pKs'] 
            if 'epKs' not in obs[i]['epKs'][j] or isgood(obs[i]['tp'][j]['epKs']):
                obs[i]['tp'][j]['epKs'] = np.nan
                wobs[i, sys['tp'][j]['ipKs']] = (epKs) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipKs']] = (obs[i]['tp'][j]['epKs']) ** (-2)
            if 'hso4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['hso4']):
                obs[i]['tp'][j]['hso4'] = np.nan
                yobs[i, sys['tp'][j]['iphso4']] = np.nan
            else:
                yobs[i, sys['tp'][j]['iphso4']] = p(obs[i]['tp'][j]['hso4'] * 1e-6) # convt µmol/kg to mol/kg
            if 'ehso4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['ehso4']):
                obs[i]['tp'][j]['ehso4'] = np.nan
                wobs[i, sys['tp'][j]['iphso4']] = np.nan
            else:
                wobs[i, sys['tp'][j]['iphso4']] = w(obs[i]['tp'][j]['hso4'], obs[i]['tp'][j]['ehso4']) 
            if 'so4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['so4']):
                obs[i]['tp'][j]['so4'] = np.nan
                yobs[i, sys['tp'][j]['ipso4']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipso4']] = p(obs[i]['tp'][j]['so4'] * 1e-6) # convt µmol/kg to mol/kg
            if 'eso4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eso4']):
                obs[i]['tp'][j]['eso4'] = np.nan
                wobs[i, sys['tp'][j]['ipso4']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipso4']] = w(obs[i]['tp'][j]['so4'], obs[i]['tp'][j]['eso4']) 
            
            # fluoride system
            if 'pKf' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pKf']):
                obs[i]['tp'][j]['pKf'] = np.nan
                yobs[i, sys['tp'][j]['ipKf']] = pKf
            else:
                yobs[i, sys['tp'][j]['ipKf']] = obs[i]['tp'][j]['pKf']
            if 'epKf' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epKf']):
                obs[i]['tp'][j]['epKf'] = np.nan
                wobs[i, sys['tp'][j]['ipKf']] = (epKf) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipKf']] = (obs[i]['tp'][j]['epKf']) ** (-2) 
            if 'HF' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['HF']):
                obs[i]['tp'][j]['HF'] = np.nan
                yobs[i, sys['tp'][j]['ipHF']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipHF']] = p(obs[i]['tp'][j]['HF'] * 1e-6) # convt µmol/kg to mol/kg
            if 'eHF' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eHF']):
                obs[i]['tp'][j]['eHF'] = np.nan
                wobs[i, sys['tp'][j]['ipHF']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipHF']] = w(obs[i]['tp'][j]['HF'], obs[i]['tp'][j]['eHF'])
            if 'F' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['F']):
                obs[i]['tp'][j]['F'] = np.nan
                yobs[i, sys['tp'][j]['ipF']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipF']] = p(obs[i]['tp'][j]['F'] * 1e-6) # convt µmol/kg to mol/kg
            if 'eF' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eF']):
                obs[i]['tp'][j]['eF'] = np.nan
                wobs[i, sys['tp'][j]['ipF']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipF']] = w(obs[i]['tp'][j]['F'], obs[i]['tp'][j]['eF'])  
            
            # phosphate system
            if 'pKp1' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pKp1']):
                obs[i]['tp'][j]['pKp1'] = np.nan
                yobs[i, sys['tp'][j]['ipKp1']] = pKp1
            else:
                yobs[i, sys['tp'][j]['ipKp1']] = obs[i]['tp'][j]['pKp1']
            if 'epKp1' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epKp1']):
                obs[i]['tp'][j]['epKp1'] = np.nan
                wobs[i, sys['tp'][j]['ipKp1']] = (epKp1) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipKp1']] = (obs[i]['tp'][j]['epKp1']) ** (-2) 
            if 'pKp2' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pKp2']):
                obs[i]['tp'][j]['pKp2'] = np.nan
                yobs[i, sys['tp'][j]['ipKp2']] = pKp2
            else:
                yobs[i, sys['tp'][j]['ipKp2']] = obs[i]['tp'][j]['pKp2']
            if 'epKp2' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epKp2']):
                obs[i]['tp'][j]['epKp2'] = np.nan
                wobs[i, sys['tp'][j]['ipKp2']] = (epKp2) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipKp2']] = (obs[i]['tp'][j]['epKp2']) ** (-2) 
            if 'pKp3' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pKp3']):
                obs[i]['tp'][j]['pKp3'] = np.nan
                yobs[i, sys['tp'][j]['ipKp3']] = pKp3
            else:
                yobs[i, sys['tp'][j]['ipKp3']] = obs[i]['tp'][j]['pKp3'] 
            if 'epKp3' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epKp3']):
                obs[i]['tp'][j]['epKp3'] = np.nan
                wobs[i, sys['tp'][j]['ipKp3']] = (epKp3) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipKp3']] = (obs[i]['tp'][j]['epKp3']) ** (-2) 
            if 'h3po4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['h3po4']):
                obs[i]['tp'][j]['h3po4'] = np.nan
                yobs[i, sys['tp'][j]['iph3po4']] = np.nan
            else:
                yobs[i, sys['tp'][j]['iph3po4']] = p(obs[i]['tp'][j]['h3po4'] * 1e-6) # convt µmol/kg to mol/kg
            if 'eh3po4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eh3po4']):
                obs[i]['tp'][j]['eh3po4'] = np.nan
                wobs[i, sys['tp'][j]['iph3po4']] = np.nan
            else:
                wobs[i, sys['tp'][j]['iph3po4']] = w(obs[i]['tp'][j]['h3po4'], obs[i]['tp'][j]['eh3po4']) 
            if 'h2po4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['h2po4']):
                obs[i]['tp'][j]['h2po4'] = np.nan
                yobs[i, sys['tp'][j]['iph2po4']] = np.nan
            else:
                yobs[i, sys['tp'][j]['iph2po4']] = p(obs[i]['tp'][j]['h2po4'] * 1e-6) # convt µmol/kg to mol/kg
            if 'eh2po4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eh2po4']):
                obs[i]['tp'][j]['eh2po4'] = np.nan
                wobs[i, sys['tp'][j]['iph2po4']] = np.nan
            else:
                wobs[i, sys['tp'][j]['iph2po4']] = w(obs[i]['tp'][j]['h2po4'], obs[i]['tp'][j]['eh2po4']) 
            if 'hpo4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['hpo4']):
                obs[i]['tp'][j]['hpo4'] = np.nan
                yobs[i, sys['tp'][j]['iphpo4']] = np.nan
            else:
                yobs[i, sys['tp'][j]['iphpo4']] = p(obs[i]['tp'][j]['hpo4'] * 1e-6) # convt µmol/kg to mol/kg
            if 'ehpo4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['ehpo4']):
                obs[i]['tp'][j]['ehpo4'] = np.nan
                wobs[i, sys['tp'][j]['iphpo4']] = np.nan
            else:
                wobs[i, sys['tp'][j]['iphpo4']] = w(obs[i]['tp'][j]['hpo4'], obs[i]['tp'][j]['ehpo4']) 
            if 'po4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['po4']):
                obs[i]['tp'][j]['po4'] = np.nan
                yobs[i, sys['tp'][j]['ippo4']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ippo4']] = p(obs[i]['tp'][j]['po4'] * 1e-6) # convt µmol/kg to mol/kg
            if 'epo4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epo4']):
                obs[i]['tp'][j]['epo4'] = np.nan
                wobs[i, sys['tp'][j]['ippo4']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ippo4']] = w(obs[i]['tp'][j]['po4'], obs[i]['tp'][j]['epo4']) 
            
            # silicate system
            if 'pKsi' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pKsi']):
                obs[i]['tp'][j]['pKsi'] = np.nan
                yobs[i, sys['tp'][j]['ipKsi']] = pKsi
            else:
                yobs[i, sys['tp'][j]['ipKsi']] = obs[i]['tp'][j]['pKsi'] 
            if 'epKsi' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epKsi']):
                obs[i]['tp'][j]['epKsi'] = np.nan
                wobs[i, sys['tp'][j]['ipKsi']] = (epKsi) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipKsi']] = (obs[i]['tp'][j]['epKsi']) ** (-2)
            if 'siooh3' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['siooh3']):
                obs[i]['tp'][j]['siooh3'] = np.nan
                yobs[i, sys['tp'][j]['ipsiooh3']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipsiooh3']] = p(obs[i]['tp'][j]['siooh3'] * 1e-6) # convt µmol/kg to mol/kg
            if 'esiooh3' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['esiooh3']):
                obs[i]['tp'][j]['esiooh3'] = np.nan
                wobs[i, sys['tp'][j]['ipsiooh3']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipsiooh3']] = w(obs[i]['tp'][j]['siooh3'], obs[i]['tp'][j]['esiooh3'])     
            if 'sioh4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['sioh4']):
                obs[i]['tp'][j]['sioh4'] = np.nan
                yobs[i, sys['tp'][j]['ipsioh4']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipsioh4']] = p(obs[i]['tp'][j]['sioh4'] * 1e-6) # convt µmol/kg to mol/kg
            if 'eosioh4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['esioh4']):
                obs[i]['tp'][j]['esioh4'] = np.nan
                wobs[i, sys['tp'][j]['ipsioh4']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipsioh4']] = w(obs[i]['tp'][j]['sioh4'], obs[i]['tp'][j]['esioh4']) 
            
            # ammonia system
            if 'pKnh4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pKnh4']):
                obs[i]['tp'][j]['pKnh4'] = np.nan
                yobs[i, sys['tp'][j]['ipKnh4']] = pKnh4
            else:
                yobs[i, sys['tp'][j]['ipKnh4']] = obs[i]['tp'][j]['pKnh4']
            if 'epKnh4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epKnh4']):
                obs[i]['tp'][j]['epKnh4'] = np.nan
                wobs[i, sys['tp'][j]['ipKnh4']] = (epKnh4) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipKnh4']] = (obs[i]['tp'][j]['epKnh4']) ** (-2) 
            if 'nh4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['nh4']):
                obs[i]['tp'][j]['nh4'] = np.nan
                yobs[i, sys['tp'][j]['ipnh4']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipnh4']] = p(obs[i]['tp'][j]['nh4'] * 1e-6) # convt µmol/kg to mol/kg
            if 'enh4' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['enh4']):
                obs[i]['tp'][j]['enh4'] = np.nan
                wobs[i, sys['tp'][j]['ipnh4']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipnh4']] = w(obs[i]['tp'][j]['nh4'], obs[i]['tp'][j]['enh4']) 
            if 'nh3' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['nh3']):
                obs[i]['tp'][j]['nh3'] = np.nan
                yobs[i, sys['tp'][j]['ipnh3']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipnh3']] = p(obs[i]['tp'][j]['nh3'] * 1e-6) # convt µmol/kg to mol/kg
            if 'enh3' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['enh3']):
                obs[i]['tp'][j]['enh3'] = np.nan
                wobs[i, sys['tp'][j]['ipnh3']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipnh3']] = w(obs[i]['tp'][j]['nh3'], obs[i]['tp'][j]['enh3']) 
            
            # sulfide system
            if 'pKh2s' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pKh2s']):
                obs[i]['tp'][j]['pKh2s'] = np.nan
                yobs[i, sys['tp'][j]['ipKh2s']] = pKh2s
            else:
                yobs[i, sys['tp'][j]['ipKh2s']] = obs[i]['tp'][j]['pKh2s'] 
            if 'epKh2s' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epKh2s']):
                obs[i]['tp'][j]['epKh2s'] = np.nan
                wobs[i, sys['tp'][j]['ipKh2s']] = (epKh2s) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipKh2s']] = (obs[i]['tp'][j]['epKh2s']) ** (-2)
            if 'H2S' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['H2S']):
                obs[i]['tp'][j]['H2S'] = np.nan
                yobs[i, sys['tp'][j]['ipH2S']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipH2S']] = p(obs[i]['tp'][j]['H2S'] * 1e-6) # convt µmol/kg to mol/kg
            if 'eH2S' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eH2S']):
                obs[i]['tp'][j]['eH2S'] = np.nan
                wobs[i, sys['tp'][j]['ipH2S']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipH2S']] = w(obs[i]['tp'][j]['H2S'], obs[i]['tp'][j]['eH2S']) 
            if 'HS' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['HS']):
                obs[i]['tp'][j]['HS'] = np.nan
                yobs[i, sys['tp'][j]['ipHS']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipHS']] = p(obs[i]['tp'][j]['HS'] * 1e-6) # convt µmol/kg to mol/kg
            if 'eHS' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eHS']):
                obs[i]['tp'][j]['eHS'] = np.nan
                wobs[i, sys['tp'][j]['ipHS']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipHS']] = w(obs[i]['tp'][j]['HS'], obs[i]['tp'][j]['eHS']) 

            # calcium carbonate solubility system
            if 'pKar' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pKar']):
                obs[i]['tp'][j]['pKar'] = np.nan
                yobs[i, sys['tp'][j]['ipKar']] = pKar
            else:
                yobs[i, sys['tp'][j]['ipKar']] = obs[i]['tp'][j]['pKar']
            if 'epKar' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epKar']):
                obs[i]['tp'][j]['epKar'] = np.nan
                wobs[i, sys['tp'][j]['ipKar']] = (epKar) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipKar']] = (obs[i]['tp'][j]['epKar']) ** (-2)
            if 'pKca' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['pKca']):
                obs[i]['tp'][j]['pKca'] = np.nan
                yobs[i, sys['tp'][j]['ipKca']] = pKca
            else:
                yobs[i, sys['tp'][j]['ipKca']] = obs[i]['tp'][j]['pKca']
            if 'epKca' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['epKca']):
                obs[i]['tp'][j]['epKca'] = np.nan
                wobs[i, sys['tp'][j]['ipKca']] =(epKca) ** (-2)
            else:
                wobs[i, sys['tp'][j]['ipKca']] = (obs[i]['tp'][j]['epKca']) ** (-2)
            if 'OmegaAr' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['OmegaAr']):
                obs[i]['tp'][j]['OmegaAr'] = np.nan
                yobs[i, sys['tp'][j]['ipOmegaAr']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipOmegaAr']] = p(obs[i]['tp'][j]['OmegaAr']) # Omega is dimensionless
            if 'eOmegaAr' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eOmegaAr']):
                obs[i]['tp'][j]['eOmegaAr'] = np.nan
                wobs[i, sys['tp'][j]['ipOmegaAr']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipOmegaAr']] = w(obs[i]['tp'][j]['OmegaAr'], obs[i]['tp'][j]['eOmegaAr']) 
            if 'OmegaCa' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['OmegaCa']):
                obs[i]['tp'][j]['OmegaCa'] = np.nan
                yobs[i, sys['tp'][j]['ipOmegaCa']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipOmegaCa']] = p(obs[i]['tp'][j]['OmegaCa']) # Omega is dimensionless
            if 'eOmegaCa' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eOmegaCa']):
                obs[i]['tp'][j]['eOmegaCa'] = np.nan
                wobs[i, sys['tp'][j]['ipOmegaCa']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipOmegaCa']] = w(obs[i]['tp'][j]['OmegaCa'], obs[i]['tp'][j]['eOmegaCa']) 
            if 'ca' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['ca']):
                obs[i]['tp'][j]['ca'] = np.nan
                yobs[i, sys['tp'][j]['ipca']] = np.nan
            else:
                yobs[i, sys['tp'][j]['ipca']] = p(obs[i]['tp'][j]['ca'] * 1e-6) # convt µmol/kg to mol/kg
            if 'eca' not in obs[i]['tp'][j] or isgood(obs[i]['tp'][j]['eca']):
                obs[i]['tp'][j]['eca'] = np.nan
                wobs[i, sys['tp'][j]['ipca']] = np.nan
            else:
                wobs[i, sys['tp'][j]['ipca']] = w(obs[i]['tp'][j]['ca'], obs[i]['tp'][j]['eca']) 
