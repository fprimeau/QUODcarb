# driver to go with QUODcarb

# g = GOMECC data
# format long

#from compare.m import compare
import scipy
import numpy as np
import copy
import mksys
import parse_input as pi
#import QUODcarb

data = scipy.io.loadmat('data.mat')
nD = len(data)

#   choose options for opt structure
opt = {
    'K1K2': 4,     # option for K1K2 formulation
    'KSO4': 1,      # option for KSO4 formulation
    'KF': 2,        # option for KF formulation
    'TB': 2,        # option for TB formulation
    'phscale': 1,   # 1 = tot, 2 = sws, 3 = free, 4 = NBS
    'printcsv': 0,  # 1print est to CSV? 1 = on, 0 = off
    #'fname': 'QUODcarb_output.csv',  # Uncomment and set the filename if printcsv is on
    #'fname': 'output_csv/Q5_Nov3.csv',
    'co2press': 1,  # 1 = on, 0 = off
    'Revelle': 0,   # 1 = on, 0 = off
    'printmes': 0   # 1 = on, 0 = off
}

obs = []

# read in GOMECC data and put into obs structure
for i in range(nD):
    obs.append({})
    # measurements independent of (T,P)
    obs[i]['TC'] = data['data'][i][4]     # (umol/kg)
    obs[i]['eTC'] =  2.01                   # TC error ±2.01 umol/kg
    obs[i]['TA'] = data['data'][i][6]
    obs[i]['TC'] = 1.78                    # TA error ±2.01 umol/kg
    obs[i]['sal'] = data['data'][i][2]    # PSU
    obs[i]['esal'] = 0.002                 # new = 0.001

    # nutrients P and Si also independent of (T,P)
    obs[i]['TP'] = data['data'][i][23]
    obs[i]['eTP'] = data['data'][i][23] * 0.001    # 0.01% meas uncertainty
    obs[i]['TSi'] = data['data'][i][17]
    obs[i]['eTSi'] = data['data'][i][17] * 0.001   # 0.01% meas uncertainty
    obs[i]['tp'] = []
    # first(T,P)-dependent measurement
    """
    'T': deg C, CTD temp
    'eT': ±0.02 degC
    'P': dbar
    'eP': (max) ± 0.63 dbar 
    """
    obs[i]['tp'].append({})
    obs[i]['tp'][0]['T'] = data['data'][i][1]
    obs[i]['tp'][0]['eT'] = 0.02
    obs[i]['tp'][0]['P'] = data['data'][i][0]
    obs[i]['tp'][0]['eP'] = 0.63
    # second(T,P)-dependent measurement
    """
    'T':  degC
    'eT': from cruise report
    'P':  in(i+ad,1);  NOT in situ
    'ph': total scale
    'co3':(µmol/kg)
    'eco3':std ±2µmol/kg
    """
    obs[i]['tp'].append({})
    obs[i]['tp'][1]['T'] = 25
    obs[i]['tp'][1]['eT'] = 0.05
    obs[i]['tp'][1]['P'] = 0.0
    obs[i]['tp'][1]['eP'] = 0.63
    obs[i]['tp'][1]['ph'] = data['data'][i][8]
    obs[i]['tp'][1]['eph'] = 0.0004, 
    obs[i]['tp'][1]['co3'] = data['data'][i][13]
    obs[i]['tp'][1]['eco3'] = 2.0
    # third (T,P)-dependent measurement
    """
    'T':       degC
    'eT':      from cruise report
    'P':       dbar (surface pressure for pco2)
    'pco2':    (µatm)
    'epco2':   0.21% relative std error (avg)
    """
    obs[i]['tp'].append({})
    obs[i]['tp'][2]['T'] = 20
    obs[i]['tp'][2]['eT'] = 0.03
    obs[i]['tp'][2]['P'] = 0.0
    obs[i]['tp'][2]['eP'] = 0.63
    obs[i]['tp'][2]['pco2'] = data['data'][i, 10]
    obs[i]['tp'][2]['epco2'] = data['data'][i][10] * 0.0021

obs = np.array(obs)
obs_backup = copy.deepcopy(obs)

# Q5: All five input
# CT AT pH pCO2 CO3 (Q5) (fid5)

# obs = copy.deepcopy(obs_backup)

# newtuple = QUODcarb(obs,opt)
# est = newtuple[0]
# obs = newtuple[1]
# sys = newtuple[2]
# iflag = newtuple[3]

# est05 = est
# Q2: Input pairs
# TA TC (Q2) (fid2)
# for i in range(nD):
#     obs[i]['tp']['ph'] = np.nan
#     obs[i]['tp']['eph'] = np.nan
#     obs[i]['tp']['pco2'] = np.nan
#     obs[i]['tp']['epco2'] = np.nan
#     obs[i]['tp']['co3'] = np.nan
#     obs[i]['tp']['eco3'] = np.nan

#est,obs, _, _ = QUODcarb(obs,opt); # [est, obs, sys, iflag]
# est02  = est
# fid2   = 'compare_outs/compare_TC_TA.csv'; 
# tp     = 2 # second tp system for ph in there
#A      = compare(obs,est,opt,tp,1,fid2); # 1 for input pair TA TC

# TA ph (Q2) (fid3)

# obs = obs_backup.copy()

# for i in range(nD):
#     obs[i]['TC'] = np.nan
#     obs[i]['eTC'] = np.nan
#     obs[i]['tp']['pco2'] = np.nan
#     obs[i]['tp']['epco2'] = np.nan
#     obs[i]['tp']['co3'] = np.nan
#     obs[i]['tp']['eco3'] = np.nan

# est, obs, _, _ = QUODcarb(obs, opt)
#est03 = est
# tp      = 2
# fid3    = 'compare_outs/compare_TA_ph.csv'
# [A]     = compare(obs,est,opt,tp,3,fid3)

sys = mksys.mksys(obs[0], opt['phscale'])
pi.parse_input(obs,sys,opt,nD)    