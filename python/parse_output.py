def parse_output(z, sigx, sys, f, C):
    # Populate est, output structure with best estimates
    
    # INPUT:
    #   z := system state including equilibrium constants and lagrange multipliers
    p = sys['p']
    q = sys['q']
    
    def ebar(j):
        return 0.5 * (q(z[j] - sigx[j]) - q(z[j] + sigx[j]))

    def ebar_l(j):
        return q(-sigx[j])

    def ebar_u(j):
        return q(sigx[j])
    
    # populate 'est' structure with best estimate:
    #   1. p(value) and p(error) where p(x) = -log10(x)
    #   2. value and average error about the value in 'q' 
    #           where q(x) = x^(-10)
    #   3. upper and lower bounds in 'q' space, not symmetric
    #           about the value in 'q' space

    # PLEASE ADD VARIABLES AT END
    # CHANGING THIS ORDER BREAKS PrintCSV
    est = {}
    est['sal'] = z[sys['isal']]
    est['esal'] = sigx[sys['isal']]

    # TC (DIC)
    est['pTC'] = z[sys['ipTC']]
    est['epTC'] = sigx[sys['ipTC']]
    est['TC'] = q(z[sys['ipTC']]) * 1e6
    est['eTC'] = ebar(sys['ipTC']) * 1e6
    est['eTC_l'] = ebar_l(sys['ipTC']) * 1e6
    est['eTC_u'] = ebar_u(sys['ipTC']) * 1e6

    # TA (Alkalinity)
    est['pTA'] = z[sys['ipTA']]
    est['epTA'] = sigx[sys['ipTA']]
    est['TA'] = q(z[sys['ipTA']]) * 1e6
    est['eTA'] = ebar(sys['ipTA']) * 1e6
    est['eTA_l'] = ebar_l(sys['ipTA']) * 1e6
    est['eTA_u'] = ebar_u(sys['ipTA']) * 1e6

    # TB (Borate)
    est['pTB'] = z[sys['ipTB']]
    est['epTB'] = sigx[['sys.ipTB']]
    est['TB'] = q(z[sys['ipTB']]) * 1e6
    est['eTB'] = ebar(sys['ipTB']) * 1e6
    est['eTB_l'] = ebar_l(sys['ipTB']) * 1e6
    est['eTB_u'] = ebar_u(sys['ipTB']) * 1e6

    # TS (Sulfate)
    est['pTS'] = z[sys['ipTS']]
    est['epTS'] = sigx[sys['ipTS']]
    est['TS'] = q(z[sys['ipTS']]) * 1e6
    est['eTS'] = ebar(sys['ipTS']) * 1e6
    est['eTS_l'] = ebar_l(sys['ipTS']) * 1e6
    est['eTS_u'] = ebar_u(sys['ipTS']) * 1e6

    # TF (Fluoride)
    est['pTF'] = z[sys['ipTF']]
    est['epTF'] = sigx[sys['ipTF']]
    est['TF'] = q(z[sys['ipTF']]) * 1e6
    est['eTF'] = ebar(sys['ipTF']) * 1e6
    est['eTF_l'] = ebar_l(sys['ipTF']) * 1e6
    est['eTF_u'] = ebar_u(sys['ipTF']) * 1e6

    # TP (Phosphate)
    est['pTP'] = z[sys['ipTP']]
    est['epTP'] = sigx[sys['ipTP']]
    est['TP'] = q(z[sys['ipTP']]) * 1e6
    est['eTP'] = ebar(sys['ipTP']) * 1e6
    est['eTP_l'] = ebar_l(sys['ipTP']) * 1e6
    est['eTP_u'] = ebar_u(sys['ipTP']) * 1e6

    # TSi (Silicate)
    est['pTSi'] = z[sys['ipTSi']]
    est['epTSi'] = sigx[sys['ipTSi']]
    est['TSi'] = q(z[sys['ipTSi']]) * 1e6
    est['eTSi'] = ebar(sys['ipTSi']) * 1e6
    est['eTSi_l'] = ebar_l(sys['ipTSi']) * 1e6
    est['eTSi_u'] = ebar_u(sys['ipTSi']) * 1e6

    # TNH4 (Ammonium)
    est['pTNH4'] = z[sys['ipTNH4']]
    est['epTNH4'] = sigx[sys['ipTNH4']]
    est['TNH4'] = q(z[sys['ipTNH4']]) * 1e6
    est['eTNH4'] = ebar(sys['ipTNH4']) * 1e6
    est['eTNH4_l'] = ebar_l(sys['ipTNH4']) * 1e6
    est['eTNH4_u'] = ebar_u(sys['ipTNH4']) * 1e6

    # TH2S (Sulfide)
    est['pTH2S'] = z[sys['ipTH2S']]
    est['epTH2S'] = sigx[sys['ipTH2S']]
    est['TH2S'] = q(z[sys['ipTH2S']]) * 1e6
    est['eTH2S'] = ebar(sys['ipTH2S']) * 1e6
    est['eTH2S_l'] = ebar_l(sys['ipTH2S']) * 1e6
    est['eTH2S_u'] = ebar_u(sys['ipTH2S']) * 1e6

    # TCa (Calcium)
    est['pTCa'] = z[sys['ipTCa']]
    est['epTCa'] = sigx[sys['ipTCa']]
    est['TCa'] = q(z[sys['ipTCa']]) * 1e6
    est['eTCa'] = ebar(sys['ipTCa']) * 1e6
    est['eTCa_l'] = ebar_l(sys['ipTCa']) * 1e6
    est['eTCa_u'] = ebar_u(sys['ipTCa']) * 1e6

    # Loop through temperature/pressure data points
    est['tp'] = []
    for i in range(len(sys['tp'])):
        # temp (deg C)
        est['tp'][i]['T'] = z[sys['tp'][i]['iT']]
        est['tp'][i]['eT'] = sigx[sys['tp'][i]['iT']]
        est['tp'][i]['eT_l'] = z[sys['tp'][i]['iT']] - sigx[sys['tp'][i]['iT']]
        est['tp'][i]['eT_u'] = z[sys['tp'][i]['iT']] + sigx[sys['tp'][i]['iT']]

        # pressure (dbar)
        est['tp'][i]['P'] = z[sys['tp'][i]['iP']]
        est['tp'][i]['eP'] = sigx[sys['tp'][i]['iP']]
        est['tp'][i]['eP_l'] = z[sys['tp'][i]['iP']] - sigx[sys['tp'][i]['iP']]
        est['tp'][i]['eP_u'] = z[sys['tp'][i]['iP']] + sigx[sys['tp'][i]['iP']]

        # fCO2
        est['tp'][i]['fco2'] = q(z[sys['tp'][i]['ipfco2']]) * 1e6
        est['tp'][i]['efco2'] = ebar(sys['tp'][i]['ipfco2']) * 1e6
        est['tp'][i]['efco2_l'] = ebar_l(sys['tp'][i]['ipfco2']) * 1e6
        est['tp'][i]['efco2_u'] = ebar_u(sys['tp'][i]['ipfco2']) * 1e6
        est['tp'][i]['pfco2'] = z[sys['tp'][i]['ipfco2']]
        est['tp'][i]['epfco2'] = sigx[sys['tp'][i]['ipfco2']]

        #pCO2
        est['tp'][i]['pco2'] = q(z[sys['tp'][i]['ippco2']]) * 1e6
        est['tp'][i]['epco2'] = ebar(sys['tp'][i]['ippco2']) * 1e6
        est['tp'][i]['epco2_l'] = ebar_l(sys['tp'][i]['ippco2']) * 1e6
        est['tp'][i]['epco2_u'] = ebar_u(sys['tp'][i]['ippco2']) * 1e6
        est['tp'][i]['ppco2'] = z[sys['tp'][i]['ippco2']]
        est['tp'][i]['eppco2'] = sigx[sys['tp'][i]['ippco2']]

        #HCO3
        est['tp'][i]['hco3'] = q(z[sys['tp'][i]['iphco3']]) * 1e6
        est['tp'][i]['ehco3'] = ebar(sys['tp'][i]['iphco3']) * 1e6
        est['tp'][i]['ehco3'] = ebar(sys['tp'][i]['iphco3']) * 1e6
        est['tp'][i]['ehco3_l']   = ebar_l(sys['tp'][i]['iphco3'])*1e6
        est['tp'][i]['ehco3_u']   = ebar_u(sys['tp'][i]['iphco3'])*1e6
        est['tp'][i]['phco3']     = z(sys['tp'][i]['iphco3'])
        est['tp'][i]['ephco3']    = sigx(sys['tp'][i]['iphco3'])

        # CO3
        est['tp'][i]['co3']       = q(z(sys['tp'][i]['ipco3']))*1e6 # convt mol/kg to µmol/kg
        est['tp'][i]['eco3']      = ebar(sys['tp'][i]['ipco3'])*1e6
        est['tp'][i]['eco3_l']    = ebar_l(sys['tp'][i]['ipco3'])*1e6
        est['tp'][i]['eco3_u']    = ebar_u(sys['tp'][i]['ipco3'])*1e6
        est['tp'][i]['pco3']      = z(sys['tp'][i]['ipco3'])
        est['tp'][i]['epco3']     = sigx(sys['tp'][i]['ipco3'])

        # CO2*
        est['tp'][i]['co2st']     = q(z(sys['tp'][i]['ipco2st'])) # convt mol/kg to µmol/kg
        est['tp'][i]['eco2st']    = ebar(sys['tp'][i]['ipco2st'])
        est['tp'][i]['eco2st_l']  = ebar_l(sys['tp'][i]['ipco2st'])
        est['tp'][i]['eco2st_u']  = ebar_u(sys['tp'][i]['ipco2st'])
        est['tp'][i]['pco2st']    = z(sys['tp'][i]['ipco2st'])
        est['tp'][i]['epco2st']   = sigx(sys['tp'][i]['ipco2st'])

        # pH on the scale opt.phscale used to compute the pK values
        est['tp'][i]['ph']        = z(sys['tp'][i]['iph']) # (9)
        est['tp'][i]['eph']       = sigx(sys['tp'][i]['iph'])
        est['tp'][i]['h']         = q(z(sys['tp'][i]['iph'])) * 1e6
        est['tp'][i]['eh']        = ebar(sys['tp'][i]['iph']) * 1e6
        est['tp'][i]['eh_l']      = ebar_l(sys['tp'][i]['iph']) * 1e6
        est['tp'][i]['eh_u']      = ebar_u(sys['tp'][i]['iph']) * 1e6

        # pH_free
        est['tp'][i]['ph_free']    = z(sys['tp'](i)['iph_free']) # (15)
        est['tp'][i]['eph_free']   = sigx(sys['tp'](i)['iph_free'])
        est['tp'][i]['h_free']     = q(z(sys['tp'][i]['iph_free'])) * 1e6
        est['tp'][i]['eh_free']    = ebar(sys['tp'][i]['iph_free']) * 1e6
        est['tp'][i]['eh_free_l']  = ebar_l(sys['tp'][i]['iph_free']) * 1e6
        est['tp'][i]['eh_free_u']  = ebar_u(sys['tp'][i]['iph_free']) * 1e6

        # pH_tot
        est['tp'][i]['ph_tot']    = z(sys['tp'][i]['iph_tot']) # (21)
        est['tp'][i]['eph_tot']   = sigx(sys['tp'][i]['iph_tot'])
        est['tp'][i]['h_tot']     = q(z(sys['tp'][i]['iph_tot'])) * 1e6
        est['tp'][i]['eh_tot']    = ebar(sys['tp'][i]['iph_tot']) * 1e6
        est['tp'][i]['eh_tot_l']  = ebar_l(sys['tp'][i]['iph_tot']) * 1e6
        est['tp'][i]['eh_tot_u']  = ebar_u(sys['tp'][i]['iph_tot']) * 1e6

        # pH_sws
        est['tp'][i]['ph_sws']    = z(sys['tp'][i]['iph_sws'])
        est['tp'][i]['eph_sws']   = sigx(sys['tp'][i]['iph_sws'])
        est['tp'][i]['h_sws']     = q(z(sys['tp'][i]['iph_sws'])) * 1e6
        est['tp'][i]['eh_sws']    = ebar(sys['tp'][i]['iph_sws']) * 1e6
        est['tp'][i]['eh_sws_l']  = ebar_l(sys['tp'][i]['iph_sws']) * 1e6
        est['tp'][i]['eh_sws_u']  = ebar_u(sys['tp'][i]['iph_sws']) * 1e6

        # pH_nbs
        est['tp'][i]['ph_nbs']    = z(sys['tp'][i]['iph_nbs'])
        est['tp'][i]['eph_nbs']   = sigx(sys['tp'][i]['iph_nbs'])
        est['tp'][i]['h_nbs']     = q(z(sys['tp'][i]['iph_nbs'])) * 1e6
        est['tp'][i]['eh_nbs']    = ebar(sys['tp'][i]['iph_nbs']) * 1e6
        est['tp'][i]['eh_nbs_l']  = ebar_l(sys['tp'][i]['iph_nbs']) * 1e6
        est['tp'][i]['eh_nbs_u']  = ebar_u(sys['tp'][i]['iph_nbs']) * 1e6

        # fH = activity coefficient
        est['tp'][i]['fH']        = q(z(sys['tp'][i]['ipfH'])) * 1e6
        est['tp'][i]['efH']       = ebar(sys['tp'][i]['ipfH']) * 1e6
        est['tp'][i]['efH_l']     = ebar_l(sys['tp'][i]['ipfH']) * 1e6
        est['tp'][i]['efH_u']     = ebar_u(sys['tp'][i]['ipfH']) * 1e6
        est['tp'][i]['pfH']       = z(sys['tp'][i]['ipfH'])
        est['tp'][i]['epfH']      = sigx(sys['tp'][i]['ipfH'])

        # p2f
        est['tp'][i]['p2f']       = q(z(sys['tp'][i]['ipp2f']))
        est['tp'][i]['ep2f']      = ebar(sys['tp'][i]['ipp2f'])
        est['tp'][i]['pp2f']      = z(sys['tp'][i]['ipp2f'])
        est['tp'][i]['epp2f']     = sigx(sys['tp'][i]['ipp2f'])

        # pK0 
        est['tp'][i]['pK0']       = z(sys['tp'][i]['ipK0']) # (79)
        est['tp'][i]['epK0']      = sigx(sys['tp'][i]['ipK0'])
        est['tp'][i]['K0']        = q(z(sys['tp'][i]['ipK0']))
        est['tp'][i]['eK0']       = ebar(sys['tp'][i]['ipK0'])
        est['tp'][i]['eK0_l']     = ebar_l(sys['tp'][i]['ipK0'])
        est['tp'][i]['eK0_u']     = ebar_u(sys['tp'][i]['ipK0'])

        # pK1
        est['tp'][i]['pK1']   = z(sys['tp'][i]['ipK1'])
        est['tp'][i]['epK1']  = sigx(sys['tp'][i]['ipK1'])
        est['tp'][i]['K1']    = q(z(sys['tp'][i]['ipK1']))
        est['tp'][i]['eK1']   = ebar(sys['tp'][i]['ipK1'])
        est['tp'][i]['eK1_l'] = ebar_l(sys['tp'][i]['ipK1']) 
        est['tp'][i]['eK1_u'] = ebar_u(sys['tp'][i]['ipK1'])

        # pK2
        est['tp'][i]['pK2']   = z(sys['tp'][i]['ipK2'])
        est['tp'][i]['epK2']  = sigx(sys['tp'][i]['ipK2'])
        est['tp'][i]['K2']    = q(z(sys['tp'][i]['ipK2']))
        est['tp'][i]['eK2']   = ebar(sys['tp'][i]['ipK2'])
        est['tp'][i]['eK2_l'] = ebar_l(sys['tp'][i]['ipK2'])
        est['tp'][i]['eK2_u'] = ebar_u(sys['tp'][i]['ipK2'])
        
        # OH 
        est['tp'][i]['oh']    = q(z(sys['tp'][i]['ipoh']))*1e6 # convt
        est['tp'][i]['eoh']   = ebar(sys['tp'][i]['ipoh'])*1e6
        est['tp'][i]['eoh_l'] = ebar_l(sys['tp'][i]['ipoh'])*1e6
        est['tp'][i]['eoh_u'] = ebar_u(sys['tp'][i]['ipoh'])*1e6
        est['tp'][i]['poh']   = z(sys['tp'][i]['ipoh'])
        est['tp'][i]['epoh']  = sigx(sys['tp'][i]['ipoh'])

        # pKw 
        est['tp'][i]['pKw']   = z(sys['tp'][i]['ipKw']) # (103)
        est['tp'][i]['epKw']  = sigx(sys['tp'][i]['ipKw'])
        est['tp'][i]['Kw']    = q(z(sys['tp'][i]['ipKw']))
        est['tp'][i]['eKw']   = ebar(sys['tp'][i]['ipKw'])
        est['tp'][i]['eKw_l'] = ebar_l(sys['tp'][i]['ipKw'])
        est['tp'][i]['eKw_u'] = ebar_u(sys['tp'][i]['ipKw'])

        # BOH4 borate
        est['tp'][i]['boh4']      = q(z(sys['tp'][i]['ipboh4']))*1e6 # convt mol/kg to µmol/kg
        est['tp'][i]['eboh4']     = ebar(sys['tp'][i]['ipboh4'])*1e6
        est['tp'][i]['eboh4_l']   = ebar_l(sys['tp'][i]['ipboh4'])*1e6
        est['tp'][i]['eboh4_u']   = ebar_u(sys['tp'][i]['ipboh4'])*1e6
        est['tp'][i]['pboh4']     = z(sys['tp'][i]['ipboh4'])
        est['tp'][i]['epboh4']    = sigx(sys['tp'][i]['ipboh4'])

        # BOH3
        est['tp'][i]['boh3']      = q(z(sys['tp'](i)['ipboh3']))*1e6
        est['tp'][i]['eboh3']     = ebar(sys['tp'][i]['ipboh3'])*1e6
        est['tp'][i]['eboh3_l']   = ebar_l(sys['tp'][i]['ipboh3'])*1e6
        est['tp'][i]['eboh3_u']   = ebar_u(sys['tp'][i]['ipboh3'])*1e6
        est['tp'][i]['pboh3']     = z(sys['tp'][i]['ipboh3'])
        est['tp'][i]['epboh3']    = sigx(sys['tp'][i]['ipboh3'])
            
        # pKb
        est['tp'][i]['pKb']       = z(sys['tp'][i]['ipKb']) # (121)
        est['tp'][i]['epKb']      = sigx(sys['tp'][i]['ipKb'])
        est['tp'][i]['Kb']        = q(z(sys['tp'][i]['ipKb']))
        est['tp'][i]['eKb']       = ebar(sys['tp'][i]['ipKb'])
        est['tp'][i]['eKb_l']     = ebar_l(sys['tp'][i]['ipKb'])
        est['tp'][i]['eKb_u']     = ebar_u(sys['tp'][i]['ipKb'])

        # SO4 sulfate
        est['tp'][i]['so4']       = q(z(sys['tp'][i]['ipso4'])) # mol/kg
        est['tp'][i]['eso4']      = ebar(sys['tp'][i]['ipso4'])
        est['tp'][i]['eso4_l']    = ebar_l(sys['tp'][i]['ipso4'])
        est['tp'][i]['eso4_u']    = ebar_u(sys['tp'][i]['ipso4'])
        est['tp'][i]['pso4']      = z(sys['tp'][i]['ipso4'])
        est['tp'][i]['epso4']     = sigx(sys['tp'][i]['ipso4'])

        # HSO4
        est['tp'][i]['hso4']      = q(z(sys['tp'][i]['iphso4']))*1e6
        est['tp'][i]['ehso4']     = ebar(sys['tp'][i]['iphso4'])*1e6
        est['tp'][i]['ehso4_l']   = ebar_l(sys['tp'][i]['iphso4'])*1e6
        est['tp'][i]['ehso4_u']   = ebar_u(sys['tp'][i]['iphso4'])*1e6
        est['tp'][i]['phso4']     = z(sys['tp'][i]['iphso4'])
        est['tp'][i]['ephso4']    = sigx(sys['tp'][i]['iphso4'])

        # pKs
        est['tp'][i]['pKs']       = z(sys['tp'][i]['ipKs']) # (145)
        est['tp'][i]['epKs']      = sigx(sys['tp'][i]['ipKs'])
        est['tp'][i]['Ks']        = q(z(sys['tp'][i]['ipKs']))
        est['tp'][i]['eKs']       = ebar(sys['tp'][i]['ipKs'])
        est['tp'][i]['eKs_l']     = ebar_l(sys['tp'][i]['ipKs'])
        est['tp'][i]['eKs_u']     = ebar_u(sys['tp'][i]['ipKs'])

        # F fluoride 
        est['tp'][i]['F']         = q(z(sys['tp'][i]['ipF']))*1e6 # convt
        est['tp'][i]['eF']        = ebar(sys['tp'][i]['ipF'])*1e6
        est['tp'][i]['eF_l']      = ebar_l(sys['tp'][i]['ipF'])*1e6
        est['tp'][i]['ef_u']      = ebar_u(sys['tp'][i]['ipF'])*1e6
        est['tp'][i]['pF']        = z(sys['tp'][i]['ipF'])
        est['tp'][i]['epF']       = sigx(sys['tp'][i]['ipF'])

        # HF 
        est['tp'][i]['HF']        = q(z(sys['tp'][i]['ipHF']))*1e6
        est['tp'][i]['eHF']       = ebar(sys['tp'][i]['ipHF'])*1e6
        est['tp'][i]['eHF_l']     = ebar_l(sys['tp'][i]['ipHF'])*1e6
        est['tp'][i]['eHF_u']     = ebar_u(sys['tp'][i]['ipHF'])*1e6
        est['tp'][i]['pHF']       = z(sys['tp'][i]['ipHF'])
        est['tp'][i]['epHF']      = sigx(sys['tp'][i]['ipHF'])

        # pKf
        est['tp'][i]['pKf']       = z(sys['tp'][i]['ipKf']) # (163)
        est['tp'][i]['epKf']      = sigx(sys['tp'][i]['ipKf'])
        est['tp'][i]['Kf']        = q(z(sys['tp'][i]['ipKf']))
        est['tp'][i]['eKf']       = ebar(sys['tp'][i]['ipKf'])
        est['tp'][i]['eKf_l']     = ebar_l(sys['tp'][i]['ipKf'])
        est['tp'][i]['eKf_u']     = ebar_u(sys['tp'][i]['ipKf'])

        # PO4
        est['tp'][i]['po4']       = q(z(sys['tp'][i]['ippo4']))*1e6 # convt
        est['tp'][i]['epo4']      = ebar(sys['tp'][i]['ippo4'])*1e6
        est['tp'][i]['epo4_l']    = ebar_l(sys['tp'][i]['ippo4'])*1e6
        est['tp'][i]['epo4_u']    = ebar_u(sys['tp'][i]['ippo4'])*1e6
        est['tp'][i]['ppo4']      = z(sys['tp'][i]['ippo4'])
        est['tp'][i]['eppo4']     = sigx(sys['tp'][i]['ippo4'])

        # HPO4
        est['tp'][i]['hpo4']      = q(z(sys['tp'][i]['iphpo4']))*1e6
        est['tp'][i]['ehpo4']     = ebar(sys['tp'][i]['iphpo4'])*1e6
        est['tp'][i]['ehpo4_l']   = ebar_l(sys['tp'][i]['iphpo4'])*1e6
        est['tp'][i]['ehpo4_u']   = ebar_u(sys['tp'][i]['iphpo4'])*1e6
        est['tp'][i]['phpo4']     = z(sys['tp'][i]['iphpo4'])
        est['tp'][i]['ephpo4']    = sigx(sys['tp'][i]['iphpo4'])

        # H2PO4
        est['tp'][i]['h2po4']     = q(z(sys['tp'][i]['iph2po4']))*1e6
        est['tp'][i]['eh2po4']    = ebar(sys['tp'][i]['iph2po4'])*1e6
        est['tp'][i]['eh2po4_l']  = ebar_l(sys['tp'][i]['iph2po4'])*1e6
        est['tp'][i]['eh2po4_u']  = ebar_u(sys['tp'][i]['iph2po4'])*1e6
        est['tp'][i]['ph2po4']    = z(sys['tp'][i]['iph2po4'])
        est['tp'][i]['eph2po4']   = sigx(sys['tp'][i]['iph2po4'])            

        # H3PO4
        est['tp'][i]['h3po4']     = q(z(sys['tp'][i]['iph3po4']))*1e6
        est['tp'][i]['eh3po4']    = ebar(sys['tp'][i]['iph3po4'])*1e6
        est['tp'][i]['eh3po4_l']  = ebar_l(sys['tp'][i]['iph3po4'])*1e6
        est['tp'][i]['eh3po4_u']  = ebar_u(sys['tp'][i]['iph3po4'])*1e6
        est['tp'][i]['ph3po4 ']   = z(sys['tp'][i]['iph3po4'])
        est['tp'][i]['eph3po4']   = sigx(sys['tp'][i]['iph3po4'])

        # pKp1
        est['tp'][i]['pKp1']      = z(sys['tp'][i]['ipKp1'])
        est['tp'][i]['epKp1']     = sigx(sys['tp'][i]['ipKp1'])
        est['tp'][i]['Kp1']       = q(z(sys['tp'][i]['ipKp1']))
        est['tp'][i]['eKp1']      = ebar(sys['tp'][i]['ipKp1'])
        est['tp'][i]['eKp1_l']    = ebar_l(sys['tp'][i]['ipKp1'])
        est['tp'][i]['eKp1_u']    = ebar_u(sys['tp'][i]['ipKp1'])

        # pKp2
        est['tp'][i]['pKp2']      = z(sys['tp'][i]['ipKp2'])
        est['tp'][i]['epKp2']     = sigx(sys['tp'][i]['ipKp2'])
        est['tp'][i]['Kp2']       = q(z(sys['tp'][i]['ipKp2']))
        est['tp'][i]['eKp2']      = ebar(sys['tp'][i]['ipKp2'])
        est['tp'][i]['pKp2_l']    = ebar_l(sys['tp'][i]['ipKp2'])
        est['tp'][i]['eKp2_u']    = ebar_u(sys['tp'][i]['ipKp2'])

        # pKp3
        est['tp'][i]['pKp3']      = z(sys['tp'][i]['ipKp3'])
        est['tp'][i]['epKp3']     = sigx(sys['tp'](i)['ipKp3'])
        est['tp'][i]['Kp3']       = q(z(sys['tp'][i]['ipKp3']))
        est['tp'][i]['eKp3']      = ebar(sys['tp'][i]['ipKp3'])
        est['tp'][i]['eKp3_l']    = ebar_l(sys['tp'][i]['ipKp3'])
        est['tp'][i]['eKp3_u']    = ebar_u(sys['tp'][i]['ipKp3'])
        
        # SiOH
        est['tp'][i]['sioh4']     = q(z(sys['tp'][i]['ipsioh4']))*1e6 # convt
        est['tp'][i]['esioh4']    = ebar(sys['tp'][i]['ipsioh4'])*1e6
        est['tp'][i]['esioh4_l']  = ebar_l(sys['tp'][i]['ipsioh4'])*1e6
        est['tp'][i]['esioh4_u']  = ebar_u(sys['tp'][i]['ipsioh4'])*1e6
        est['tp'][i]['psioh4']    = z(sys['tp'][i]['ipsioh4'])
        est['tp'][i]['epsioh4']   = sigx(sys['tp'][i]['ipsioh4'])

        # SiOH3
        est['tp'][i]['siooh3']    = q(z(sys['tp'][i]['ipsiooh3']))*1e6
        est['tp'][i]['esiooh3']   = ebar(sys['tp'][i]['ipsiooh3'])*1e6
        est['tp'][i]['esiooh3_l'] = ebar_l(sys['tp'][i]['ipsiooh3'])*1e6
        est['tp'][i]['esiooh3_'] = ebar_u(sys['tp'][i]['ipsiooh3'])*1e6
        est['tp'][i]['psiooh3']   = z(sys['tp'][i]['ipsiooh3'])
        est['tp'][i]['epsiooh3']  = sigx(sys['tp'][i]['ipsiooh3'])

        # pKsi
        est['tp'][i]['pKsi']      = z(sys['tp'][i]['ipKsi'])
        est['tp'][i]['epKsi']     = sigx(sys['tp'][i]['ipKsi'])
        est['tp'][i]['Ksi']       = q(z(sys['tp'][i]['ipKsi']))
        est['tp'][i]['eKsi']      = ebar(sys['tp'][i]['ipKsi'])
        est['tp'][i]['eKsi_l']    = ebar_l(sys['tp'][i]['ipKsi'])
        est['tp'][i]['eKsi_u']    = ebar_u(sys['tp'][i]['ipKsi'])

        # NH3
        est['tp'][i]['nh3']       = q(z(sys['tp'][i]['ipnh3']))*1e6 # convt
        est['tp'][i]['enh3']      = ebar(sys['tp'][i]['ipnh3'])*1e6
        est['tp'][i]['enh3_l']    = ebar_l(sys['tp'][i]['ipnh3'])*1e6
        est['tp'][i]['enh3_u']    = ebar_u(sys['tp'][i]['ipnh3'])*1e6
        est['tp'][i]['pnh3']      = z(sys['tp'][i]['ipnh3'])
        est['tp'][i]['epnh3']     = sigx(sys['tp'][i]['ipnh3'])

        # NH4
        est['tp'][i]['nh4']       = q(z(sys['tp'][i]['ipnh4']))*1e6
        est['tp'][i]['enh4']      = ebar(sys['tp'][i]['ipnh4'])*1e6
        est['tp'][i]['enh4_l']    = ebar_l(sys['tp'][i]['ipnh4'])*1e6
        est['tp'][i]['enh4_u']    = ebar_u(sys['tp'][i]['ipnh4'])*1e6
        est['tp'][i]['pnh4']      = z(sys['tp'][i]['ipnh4'])
        est['tp'][i]['epnh4']     = sigx(sys['tp'][i]['ipnh4'])

        # pKNH4
        est['tp'][i]['pKnh4']     = z(sys['tp'][i]['ipKnh4'])
        est['tp'][i]['epKnh4']    = sigx(sys['tp'][i]['ipKnh4'])
        est['tp'][i]['Knh4']      = q(z(sys['tp'][i]['ipKnh4']))
        est['tp'][i]['eKnh4']     = ebar(sys['tp'][i]['ipKnh4'])
        est['tp'][i]['eKnh4_l']   = ebar_l(sys['tp'][i]['ipKnh4'])
        est['tp'][i]['eKnh4_u']   = ebar_u(sys['tp'][i]['ipKnh4'])

        # HS
        est['tp'][i]['HS']        = q(z(sys['tp'][i]['ipHS']))*1e6 # convt
        est['tp'][i]['eHS']       = ebar(sys['tp'][i]['ipHS'])*1e6
        est['tp'][i]['eHS_l']     = ebar_l(sys['tp'][i]['ipHS'])*1e6
        est['tp'][i]['eHS_u']     = ebar_u(sys['tp'][i]['ipHS'])*1e6
        est['tp'][i]['pHS']       = z(sys['tp'][i]['ipHS'])
        est['tp'][i]['epHS']      = sigx(sys['tp'][i]['ipHS'])

        # H2S
        est['tp'][i]['H2S']       = q(z(sys['tp'][i]['ipH2S']))*1e6
        est['tp'][i]['eH2S']      = ebar(sys['tp'][i]['ipH2S'])*1e6
        est['tp'][i]['eH2S_l']    = ebar_l(sys['tp'][i]['ipH2S'])*1e6
        est['tp'][i]['eHS2_u']    = ebar_u(sys['tp'][i]['ipH2S'])*1e6
        est['tp'][i]['pH2S']      = z(sys['tp'][i]['ipH2S'])
        est['tp'][i]['epH2S']     = sigx(sys['tp'][i]['ipH2S'])

        # pKh2s
        est['tp'][i]['pKh2s']     = z(sys['tp'][i]['ipKh2s'])
        est['tp'][i]['epKh2s']    = sigx(sys['tp'][i]['ipKh2s'])
        est['tp'][i]['Kh2s']      = q(z(sys['tp'][i]['ipKh2s']))
        est['tp'][i]['eKh2s']     = ebar(sys['tp'][i]['ipKh2s'])
        est['tp'][i]['eKh2s_l']   = ebar_l(sys['tp'][i]['ipKh2s'])
        est['tp'][i]['eKh2s_u']   = ebar_u(sys['tp'][i]['ipKh2s'])

        # Ca
        est['tp'][i]['ca']        = q(z(sys['tp'][i]['ipca']))*1e6
        est['tp'][i]['eca']       = ebar(sys['tp'][i]['ipca'])*1e6
        est['tp'][i]['eca_l']     = ebar_l(sys['tp'][i]['ipca'])*1e6
        est['tp'][i]['eca_u']     = ebar_u(sys['tp'][i]['ipca'])*1e6
        est['tp'][i]['pca']       = z(sys['tp'][i]['ipca'])
        est['tp'][i]['epca']      = sigx(sys['tp'][i]['ipca'])

        # Omega_Ar
        est['tp'][i]['OmegaAr']    = q(z(sys['tp'][i]['ipOmegaAr'])) # unitless
        est['tp'][i]['eOmegaAr']   = ebar(sys['tp'][i]['ipOmegaAr'])
        est['tp'][i]['eOmegaAr_l'] = ebar_l(sys['tp'][i]['ipOmegaAr'])
        est['tp'][i]['eOmegaAr_u ']= ebar_u(sys['tp'][i]['ipOmegaAr'])
        est['tp'][i]['pOmegaAr']   = z(sys['tp'][i]['ipOmegaAr'])
        est['tp'][i]['epOmegaAr']  = sigx(sys['tp'][i]['ipOmegaAr'])

        # pKar
        est['tp'][i]['pKar']      = z(sys['tp'][i]['ipKar'])
        est['tp'][i]['epKar']     = sigx(sys['tp'][i]['ipKar'])
        est['tp'][i]['Kar']       = q(z(sys['tp'][i]['ipKar']))
        est['tp'][i]['eKar']      = ebar(sys['tp'][i]['ipKar'])
        est['tp'][i]['eKar_l']    = ebar_l(sys['tp'][i]['ipKar'])
        est['tp'][i]['eKar_u']    = ebar_u(sys['tp'][i]['ipKar'])

        # Omega_Ca
        est['tp'][i]['OmegaCa']    = q(z(sys['tp'][i]['ipOmegaCa']))
        est['tp'][i]['eOmegaCa']   = ebar(sys['tp'][i]['ipOmegaCa'])
        est['tp'][i]['eOmegaCa_l'] = ebar_l(sys['tp'][i]['ipOmegaCa'])
        est['tp'][i]['eOmegaCa_u'] = ebar_u(sys['tp'][i]['ipOmegaCa'])
        est['tp'][i]['pOmegaCa']   = z(sys['tp'][i]['ipOmegaCa'])
        est['tp'][i]['epOmegaCa']  = sigx(sys['tp'][i]['ipOmegaCa'])

        # pKca
        est['tp'][i]['pKca']      = z(sys['tp'][i]['ipKca'])
        est['tp'][i]['epKca']     = sigx(sys['tp'][i]['ipKca'])
        est['tp'][i]['Kca']       = q(z(sys['tp'][i]['ipKca']))
        est['tp'][i]['eKca']      = ebar(sys['tp'][i]['ipKca'])
        est['tp'][i]['eKca_l']    = ebar_l(sys['tp'][i]['ipKca'])
        est['tp'][i]['eKca_u']    = ebar_u(sys['tp'][i]['ipKca'])

    # PLEASE ADD NEW VARIABLES HERE AT END
    est['f']       = f # residual f value, from limp
    est['C']       = C
    return est