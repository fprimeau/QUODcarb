import csv

def PrintCSV(est,obs,iflag,opt):
# subfunction of QUODcarb to print results to a CSV
# uses filename opt['fname'] 
#               opt['fname'] = 'QUODcarb_output.csv'(default)
    if opt['printcsv'] == 1:
        nD = len(obs)
        for i in range(nD):
            if i == 0:
                with open(opt['fname'], 'w') as fname:
                # make column headers
                    make_headers(est(i),opt,fname)
            # fill one row with data
            parse_CSV(est(i),obs(i),iflag(i),opt,fname)


def make_headers(est,opt,fid):
# Print the column headers  
    nTP = len(est['tp'])
    # row 1
    fid.write('  , ')
    fid.write('est = output, ')
    fid.write('e = 1sigma, ')

    fn = list(est.keys())
    fnl = len(fn) - 1
    for i in range(4, fnl, 6):  # MATLAB 5:6:fnl is equivalent to Python range(4, fnl, 6)
        fid.write('  , ')  # est
        fid.write('  , ')  # est.e

    for j in range(nTP):
        fnj = list(est['tp'][0].keys())
        for i in range(0, 8, 4):  # MATLAB 1:4:8 is equivalent to Python range(0, 8, 4)
            fid.write(f'tp({j+1}), ')  # est
            fid.write('  , ')  # est.e

        for i in range(8, 39, 6):  # MATLAB 9:6:39 is equivalent to Python range(8, 39, 6)
            fid.write(f'tp({j+1}), ')  # est
            fid.write('  , ')  # est.e

        for i in range(44, 75, 6):  # MATLAB 45:6:75 is equivalent to Python range(44, 75, 6)
            fid.write(f'tp({j+1}), ')  # est
            fid.write('  , ')  # est.e

        for i in range(78, len(fnj), 6):  # MATLAB 79:6:length(fnj) is equivalent to Python range(78, len(fnj), 6)
            fid.write(f'tp({j+1}), ')  # est
            fid.write('  , ')  # est.e

    fid.write('\n')  # finish first row

    # row 2
    fid.write('iflag, ')

    fid.write('est.sal, ')  # sal is first
    fid.write('est.esal, ')  # esal

    for i in range(4, fnl, 6):  # temperature independent totals
        fid.write(f'est.{fn[i]}, ')
        fid.write(f'est.e{fn[i]}, ')  # e

    for j in range(nTP):
        fnj = list(est['tp'][j].keys())
        for i in range(0, 8, 4):  # T, P
            fid.write(f'est.{fnj[i]}, ')
            fid.write(f'est.e{fnj[i]}, ')  # e = error

        for i in range(8, 39, 6):  # fco2, pco2, hco3, co3, co2*, ph
            fid.write(f'est.{fnj[i]}, ')
            fid.write(f'est.e{fnj[i]}, ')  # e = error

        for i in range(44, 75, 6):  # ph_free, ph_tot, ph_sws, ph_nbs, fH, p2f
            fid.write(f'est.{fnj[i]}, ')  
            fid.write(f'est.e{fnj[i]}, ')  # e = error

        for i in range(78, 288, 6):  # all the rest (except Revelle if on)
            fid.write(f'est.{fnj[i]}, ')  
            fid.write(f'est.e{fnj[i]}, ')  # e = error

        if opt['Revelle'] == 1:
            fid.write(f'est.{fnj[-2]}, ')  # Revelle
            fid.write(f'est.{fnj[-1]}, ')  # dpfco2dpTA
    
    fid.write('\n')  # finish second row

    # row 3
    fid.write('(0=good), ')

    fid.write('(PSU), ')  # est.sal
    fid.write('(PSU), ')  # est.esal

    for i in range(2, fnl, 6):  # temperature independent totals
        fid.write('(umol/kg), ')
        fid.write('(umol/kg), ')

    for j in range(nTP):
        fid.write('deg C, ')  # est['T']
        fid.write('deg C, ')
        fid.write('dbar, ')  # est['P']
        fid.write('dbar, ')

        for i in range(8, 16, 6):  # fco2, pco2
            fid.write('(uatm), ')
            fid.write('(uatm), ')

        for i in range(20, 34, 6):  # hco3, co3, co2*
            fid.write('(umol/kg), ')
            fid.write('(umol/kg), ')

        fid.write('(p units), ')  # ph, log10 unitless
        fid.write('  , ')

        for i in range(44, 64, 6):  # ph_free, ph_tot, ph_sws, ph_nbs
            fid.write('(p units), ')  # log10 unitless
            fid.write('  , ')

        for i in range(68, 76, 6):  # fH, p2f aka FugFac
            fid.write('(unitless), ')
            fid.write('  , ')

        for i in range(78, 92, 6):  # pK0, pK1, pK2
            fid.write('(p units), ')  # log10 unitless
            fid.write('  , ')

        fid.write('(umol/kg), ')  # oh (97)
        fid.write('(umol/kg), ')
        fid.write('(p units), ')  # pKw (103)
        fid.write('  , ')

        for i in range(3):  # (boh3, boh4, pKb) (so4, hso4, pKs)(F, HF, pKf)
            fid.write('(umol/kg), ')  # est 1st
            fid.write('(umol/kg), ')
            fid.write('(umol/kg), ')  # est 2nd
            fid.write('(umol/kg), ')
            fid.write('(p units), ')  # est['pK']
            fid.write('  , ')

        for i in range(4):  # po4, hpo4, h2po4, h3po4
            fid.write('(umol/kg), ')
            fid.write('(umol/kg), ')

        for i in range(3):  # pKp1, pKp2, pKp3
            fid.write('(p units), ')
            fid.write('  , ')

        for i in range(3):  # (sioh4, siooh3, pKsi)(nh3, nh4, pKnh4)(HS, H2S, pKh2s)
            fid.write('(umol/kg), ')  # est 1st
            fid.write('(umol/kg), ')
            fid.write('(umol/kg), ')  # est 2nd
            fid.write('(umol/kg), ')
            fid.write('(p units), ')  # est['pK']
            fid.write('  , ')

        fid.write('(umol/kg), ')  # ca
        fid.write('(umol/kg), ')

        for i in range(2):  # (Omega_Ar, pKar)(Omega_Ca, pKca)
            fid.write('(umol/kg), ')  # est 1st
            fid.write('(umol/kg), ')
            fid.write('(p units), ')  # est['pK']
            fid.write('  , ')

        if opt['Revelle'] == 1:
            fid.write('(unitless), ')  # Revelle
            fid.write('  , ')  # dpfco2dpTA

    fid.write('\n')  # finish third row


def parse_CSV(varargin):
# Print the output
    est     = varargin[0]; # posterior marginal precision
    obs     = varargin[1]; # measurement
    iflag   = varargin[2]; # Newton solver convergence flag: 0 converged, 1 not converged
    opt     = varargin[3]; # options
    fid   = varargin[4]; # filename with fopen command

    nTP     = len(est['tp'])

    with open(fid, 'w', newline='') as file:
        writer = csv.writer(file)
        row = [iflag]

        #salinity
        row.extend([est['sal'], est['esal']])
    # TC (DIC)
        row.extend([est['TC'], est['eTC']])
        # TA
        row.extend([est['TA'], est['eTA']])
        # TB borate
        row.extend([est['TB'], est['eTB']])
        # TS sulfate
        row.extend([est['TS'], est['eTS']])
        # TF fluoride
        row.extend([est['TF'], est['eTF']])
        # TP phosphate
        row.extend([est['TP'], est['eTP']])
        # TSi silicate
        row.extend([est['TSi'], est['eTSi']])
        # TNH4 nitrate
        row.extend([est['TNH4'], est['eTNH4']])
        # TH2S sulfide
        row.extend([est['TH2S'], est['eTH2S']])
        # TCa calcium solubility
        row.extend([est['TCa'], est['eTCa']])
        
        for j in range(nTP):
            tp = est['tp'][j]
            # Temp
            row.extend([tp['T'], tp['eT']])
            # Pres
            row.extend([tp['P'], tp['eP']])
            # fco2
            row.extend([tp['fco2'], tp['efco2']])
            # pco2
            row.extend([tp['pco2'], tp['epco2']])
            # hco3
            row.extend([tp['hco3'], tp['ehco3']])
            # co3
            row.extend([tp['co3'], tp['eco3']])
            # co2* (co2st)
            row.extend([tp['co2st'], tp['eco2st']])
            # ph
            row.extend([tp['ph'], tp['eph']])
            # ph_free
            row.extend([tp['ph_free'], tp['eph_free']])
            # ph_tot
            row.extend([tp['ph_tot'], tp['eph']])
            # ph_sws
            row.extend([tp['ph_sws'], tp['eph']])
            # ph_nbs
            row.extend([tp['ph_nbs'], tp['eph']])
            # fH = activity coefficient
            row.extend([tp['fH'], tp['efH']])
            # pp2f
            row.extend([tp['pp2f'], tp['epp2f']])
            # pK0
            row.extend([tp['pK0'], tp['epK0']])
            # pK1
            row.extend([tp['pK1'], tp['epK1']])
            # pK2
            row.extend([tp['pK2'], tp['epK2']])
            # oh
            row.extend([tp['oh'], tp['eoh']])
            # pKw = [h][oh]
            row.extend([tp['pKw'], tp['epKw']])
            # boh4
            row.extend([tp['boh4'], tp['eboh4']])
            # boh3
            row.extend([tp['boh3'], tp['eboh3']])
            # pKb = [h][boh4]/[boh3]
            row.extend([tp['pKb'], tp['epKb']])
            # so4
            row.extend([tp['so4'], tp['eso4']])
            # hso4
            row.extend([tp['hso4'], tp['ehso4']])
            # pKs  = [hf][so4]/[hso4]
            row.extend([tp['pKs'], tp['epKs']])
            # [F]
            row.extend([tp['F'], tp['eF']])
            # [HF] hydrogen fluoride
            row.extend([tp['HF'], tp['eHF']])
            # pKf = [h][F]/[HF]
            row.extend([tp['pKf'], tp['epKf']])
            # po4
            row.extend([tp['po4'], tp['epo4']])
            # hpo4
            row.extend([tp['hpo4'], tp['ehpo4']])
            # h2po4
            row.extend([tp['h2po4'], tp['eh2po4']])
            # h3po4
            row.extend([tp['h3po4'], tp['eh3po4']])
            # pKp1 = [h][h2po4]/[h3po4]
            row.extend([tp['pKp1'], tp['epKp1']])
            # pKp2 = [h][hpo4]/[h2po4]
            row.extend([tp['pKp2'], tp['epKp2']])
            # pKp3 = [h][po4]/[hpo4]
            row.extend([tp['pKp3'], tp['epKp3']])
            # sioh4
            row.extend([tp['sioh4'], tp['esioh4']])
            # siooh3
            row.extend([tp['siooh3'], tp['esiooh3']])
            # pKSi = [h][siooh3]/[sioh4]
            row.extend([tp['pKsi'], tp['epKsi']])
            # nh3
            row.extend([tp['nh3'], tp['enh3']])
            # nh4
            row.extend([tp['nh4'], tp['enh4']])
            # pKnh4 = [h][nh3]/[nh4]
            row.extend([tp['pKnh4'], tp['epKnh4']])
            # hs
            row.extend([tp['HS'], tp['eHS']])
            # h2s
            row.extend([tp['H2S'], tp['eH2S']])
            # pKh2s = [h][hs]/[h2s]
            row.extend([tp['pKh2s'], tp['epKh2s']])
            # ca
            row.extend([tp['ca'], tp['eca']])
            # OmegaAr
            row.extend([tp['OmegaAr'], tp['eOmegaAr']])
            # OmegaCa
            row.extend([tp['OmegaCa'], tp['eOmegaCa']])
            # pKar = [ca][co3]/[omegaAr]
            row.extend([tp['pKar'], tp['epKar']])
            # pKca = [ca][co3]/[omegaCa]
            row.extend([tp['pKca'], tp['epKca']])
            
            if opt['Revelle'] == 1:
                row.extend([tp['Revelle'], tp['dpfco2dpTA']])
