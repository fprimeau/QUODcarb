import numpy as np
import calc_pTOT as cpt
import calc_pk as cp

def update_y(y, x, obs, sys, opt):
    nTP = len(sys['tp'])
    M = sys['M']
    K = sys['K']
    nv = M.shape[1]
    gy = np.zeros((nv, len(y)))
    ggy = np.zeros((nv, len(y), len(y)))
    sal = x[sys['isal']]
    e3 = np.finfo(float).eps ** 3
    ie3 = 1j * e3

    # Calculate totals (see Ref's in calc_pTOT)
    pT, gpT, ggpT, _ = cpt.calc_pTOT(opt, sal)
    pTB, gpTB, ggpTB = pT[0], gpT[0], ggpT[0]
    pTS, gpTS, ggpTS = pT[1], gpT[1], ggpT[1]
    pTF, gpTF, ggpTF = pT[2], gpT[2], ggpT[2]
    pTCa, gpTCa, ggpTCa = pT[3], gpT[3], ggpT[3]

    # Update y
    if np.isnan(obs.TB):
        y[sys['ipTB']] = pTB
        gy[sys['ipTB'], sys['isal']] = gpTB
        ggy[sys['ipTB'], sys['isal'], sys['isal']] = ggpTB

    if np.isnan(obs.TS):
        y[sys['ipTS']] = pTS
        gy[sys['ipTS'], sys['isal']] = gpTS
        ggy[sys['ipTS'], sys['isal'], sys['isal']] = ggpTS

    if np.isnan(obs.TF):
        y[sys['ipTF']] = pTF
        gy[sys['ipTF'], sys['isal']] = gpTF
        ggy[sys['ipTF'], sys['isal'], sys['isal']] = ggpTF

    if np.isnan(obs.TCa):
        y[sys['ipTCa']] = pTCa
        gy[sys['ipTCa'], sys['isal']] = gpTCa
        ggy[sys['ipTCa'], sys['isal'], sys['isal']] = ggpTCa

    for i in range(nTP):
        # Use complex step method to get ∂T, ∂S, ∂P
        pK, gpK = cp.calc_pK(opt, x[sys['tp'][i]['iT']], x[sys['isal']], x[sys['tp'][i]['iP']])
        
        pK_T, gpK_T = cp.calc_pK(opt, x[sys['tp'][i]['iT']] + ie3, x[sys['isal']], x[sys['tp'][i]['iP']])
        pK_S, gpK_S = cp.calc_pK(opt, x[sys['tp'][i]['iT']], x[sys['isal']] + ie3, x[sys['tp'][i]['iP']])
        pK_P, gpK_P = cp.calc_pK(opt, x[sys['tp'][i]['iT']], x[sys['isal']], x[sys['tp'][i]['iP']] + ie3)
        
        ggpK = np.zeros((len(pK), 3, 3))
        ggpK[:, 0, :] = np.imag(gpK_T) / e3
        ggpK[:, 1, :] = np.imag(gpK_S) / e3
        ggpK[:, 2, :] = np.imag(gpK_P) / e3
        
        iTSP = [sys['tp'][i]['iT'], sys['isal'], sys['tp'][i]['iP']]
        if np.isnan(obs.tp[i].pK0):
            y[sys['tp'][i]['ipK0']] = pK[0]
            gy[sys['tp'][i]['ipK0'], iTSP] = gpK[0, :]
            ggy[sys['tp'][i]['ipK0'], iTSP, iTSP] = ggpK[0, :, :]
        
        if np.isnan(obs.tp[i].pK1):
            y[sys['tp'][i]['ipK1']] = pK[1]
            gy[sys['tp'][i]['ipK1'], iTSP] = gpK[1, :]
            ggy[sys['tp'][i]['ipK1'], iTSP, iTSP] = ggpK[1, :, :]
        
        if np.isnan(obs.tp[i].pK2):
            y[sys['tp'][i]['ipK2']] = pK[2]
            gy[sys['tp'][i]['ipK2'], iTSP] = gpK[2, :]
            ggy[sys['tp'][i]['ipK2'], iTSP, iTSP] = ggpK[2, :, :]
        
        if np.isnan(obs.tp[i].pKb):
            y[sys['tp'][i]['ipKb']] = pK[3]
            gy[sys['tp'][i]['ipKb'], iTSP] = gpK[3, :]
            ggy[sys['tp'][i]['ipKb'], iTSP, iTSP] = ggpK[3, :, :]
        
        if np.isnan(obs.tp[i].pKw):
            y[sys['tp'][i]['ipKw']] = pK[4]
            gy[sys['tp'][i]['ipKw'], iTSP] = gpK[4, :]
            ggy[sys['tp'][i]['ipKw'], iTSP, iTSP] = ggpK[4, :, :]
        
        if np.isnan(obs.tp[i].pKs):
            y[sys['tp'][i]['ipKs']] = pK[5]
            gy[sys['tp'][i]['ipKs'], iTSP] = gpK[5, :]
            ggy[sys['tp'][i]['ipKs'], iTSP, iTSP] = ggpK[5, :, :]
        
        if np.isnan(obs.tp[i].pKf):
            y[sys['tp'][i]['ipKf']] = pK[6]
            gy[sys['tp'][i]['ipKf'], iTSP] = gpK[6, :]
            ggy[sys['tp'][i]['ipKf'], iTSP, iTSP] = ggpK[6, :, :]
        
        if np.isnan(obs.tp[i].pKp1):
            y[sys['tp'][i]['ipKp1']] = pK[7]
            gy[sys['tp'][i]['ipKp1'], iTSP] = gpK[7, :]
            ggy[sys['tp'][i]['ipKp1'], iTSP, iTSP] = ggpK[7, :, :]
        
        if np.isnan(obs.tp[i].pKp2):
            y[sys['tp'][i]['ipKp2']] = pK[8]
            gy[sys['tp'][i]['ipKp2'], iTSP] = gpK[8, :]
            ggy[sys['tp'][i]['ipKp2'], iTSP, iTSP] = ggpK[8, :, :]
        
        if np.isnan(obs.tp[i].pKp3):
            y[sys['tp'][i]['ipKp3']] = pK[9]
            gy[sys['tp'][i]['ipKp3'], iTSP] = gpK[9, :]
            ggy[sys['tp'][i]['ipKp3'], iTSP, iTSP] = ggpK[9, :, :]
        
        if np.isnan(obs.tp[i].pKsi):
            y[sys['tp'][i]['ipKsi']] = pK[10]
            gy[sys['tp'][i]['ipKsi'], iTSP] = gpK[10, :]
            ggy[sys['tp'][i]['ipKsi'], iTSP, iTSP] = ggpK[10, :, :]
        
        if np.isnan(obs.tp[i].pKnh4):
            y[sys['tp'][i]['ipKnh4']] = pK[11]
            gy[sys['tp'][i]['ipKnh4'], iTSP] = gpK[11, :]
            ggy[sys['tp'][i]['ipKnh4'], iTSP, iTSP] = ggpK[11, :, :]
        
        if np.isnan(obs.tp[i].pKh2s):
            y[sys['tp'][i]['ipKh2s']] = pK[12]
            gy[sys['tp'][i]['ipKh2s'], iTSP] = gpK[12, :]
            ggy[sys['tp'][i]['ipKh2s'], iTSP, iTSP] = ggpK[12, :, :]
        
        if np.isnan(obs.tp[i].pp2f):
            y[sys['tp'][i]['ipp2f']] = pK[13]
            gy[sys['tp'][i]['ipp2f'], iTSP] = gpK[13, :]
            ggy[sys['tp'][i]['ipp2f'], iTSP, iTSP] = ggpK[14, :, :]
        
        if np.isnan(obs.tp[i].pKar):
            y[sys['tp'][i]['ipKar']] = pK[15]
            gy[sys['tp'][i]]['ipKar']['iTSP'] = gpK[15, :]
            ggy[sys['tp'][i]['ipKar'], iTSP, iTSP] = ggpK[15, :, :]
            y[sys['tp'][i]['ipKca']] = pK[16]
            gy[sys['tp'][i]['ipKca'], iTSP] = gpK[16, :]
            ggy[sys['tp'][i]['ipKca'], iTSP, iTSP] = ggpK[16, :, :]
    return [y,gy,ggy]
