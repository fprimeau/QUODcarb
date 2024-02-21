# driver to go with QUODcarb

# g = GOMECC data

# format long

from compare.m import compare
import numpy as np
import copy
import QUODcarb

data = np.load('datag.npy')

nD = 5  # 1186; until UW start, or before: #length(in) 

# choose options for opt structure
opt = {
    'K1K2': 10,     # option for K1K2 formulation
    'KSO4': 1,      # option for KSO4 formulation
    'KF': 2,        # option for KF formulation
    'TB': 2,        # option for TB formulation
    'phscale': 1,   # 1 = tot, 2 = sws, 3 = free, 4 = NBS
    'printcsv': 0,  # 1print est to CSV? 1 = on, 0 = off
    #'fname': 'QUODcarb_output.csv',  # Uncomment and set the filename if printcsv is on
    #'fname': 'output_csv/Q5_Nov3.csv',
    'co2press': 0,  # 1 = on, 0 = off
    'Revelle': 0,   # 1 = on, 0 = off
    'printmes': 0   # 1 = on, 0 = off
}

# read in GOMECC data and put into obs structure

obs = []

for i in range(1, nD + 1):
    obs_i = {}

    # measurements independent of (T,P)
    obs_i['TC'] = data[i-1, 4]  # (umol/kg)
    obs_i['eTC'] = 2.01  # TC error ±2.01 umol/kg
    obs_i['TA'] = data[i-1, 6]
    obs_i['eTA'] = 1.78  # TA error ±2.01 umol/kg
    obs_i['sal'] = data[i-1, 2]  # PSU
    obs_i['esal'] = 0.002  # new = 0.001

    # nutrients P and Si also independent of (T,P)
    obs_i['TP'] = data[i-1, 23]
    obs_i['eTP'] = data[i-1, 23] * 0.001  # 0.01% meas uncertainty
    obs_i['TSi'] = data[i-1, 17]
    obs_i['eTSi'] = data[i-1, 17] * 0.001  # 0.01% meas uncertainty

    # first (T,P)-dependent measurement
    '''
    'T': deg C, CTD temp
    'eT': ±0.02 degC
    'P': dbar
    'eP': (max) ± 0.63 dbar    
    '''
    
    obs_i['tp'] = [{'T': data[i-1, 1], 'eT': 0.02, 'P': data[i-1, 0], 'eP': 0.63}]

    # second(T,P)-dependent measurement
    '''
    'T':  degC
    'eT': from cruise report
    TODO:  Not sure about the P option
    'P':  NOT in situ
    'ph': total scale
    'co3':(µmol/kg)
    'eco3':std ±2µmol/kg

    '''

    obs_i['tp'].append({'T': 25, 'eT': 0.05, 'P': 0.0, 'eP': 0.63, 'ph': data[i-1, 8], 'eph': 0.0004,
                        'co3': data[i-1, 14], 'eco3': 2.0})

    # third (T,P)-dependent measurement
    '''
    'T':       degC
    'eT':      from cruise report
    'P':       dbar (surface pressure for pco2)
    'pco2':    (µatm)
    'epco2':   0.21% relative std error (avg)
    '''
    obs_i['tp'].append({'T': 20, 'eT': 0.03, 'P': 0.0, 'eP': 0.63,
                        'pco2': data[i-1, 11], 'epco2': data[i-1, 11] * 0.0021})

    obs.append(obs_i)

obs_backup = copy.deepcopy(obs)


# Q5: All five input
# CT AT pH pCO2 CO3 (Q5) (fid5)

obs = copy.deepcopy(obs_backup)
#future TODO:make sure the return from the function is right
newtuple = QUODcarb(obs,opt)
est = newtuple[0]
obs = newtuple[1]
sys = newtuple[2]
iflag = newtuple[3]

est05 = est


# Q2: Input pairs
# TA TC (Q2) (fid2)
for i in range(nD):
    obs[i]['tp'][1]['ph'] = np.nan
    obs[i]['tp'][1]['eph'] = np.nan
    obs[i]['tp'][2]['pco2'] = np.nan
    obs[i]['tp'][2]['epco2'] = np.nan
    obs[i]['tp'][1]['co3'] = np.nan
    obs[i]['tp'][1]['eco3'] = np.nan

est,obs, _, _ = QUODcarb(obs,opt); # [est, obs, sys, iflag]
est02  = est
fid2   = 'compare_outs/compare_TC_TA.csv'; 
tp     = 2 # second tp system for ph in there
#future TODO: make sure compare function works
A      = compare(obs,est,opt,tp,1,fid2); # 1 for input pair TA TC

# TA ph (Q2) (fid3)

obs = obs_backup.copy()

for i in range(nD):
    obs[i]['TC'] = np.nan
    obs[i]['eTC'] = np.nan
    obs[i]['tp'][2]['pco2'] = np.nan
    obs[i]['tp'][2]['epco2'] = np.nan
    obs[i]['tp'][1]['co3'] = np.nan
    obs[i]['tp'][1]['eco3'] = np.nan

est, obs, _, _ = QUODcarb(obs, opt)
est03 = est
# tp      = 2
# fid3    = 'compare_outs/compare_TA_ph.csv'
# [A]     = compare(obs,est,opt,tp,3,fid3)
