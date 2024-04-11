# driver to go with QUODcarb

# g = GOMECC data
# format long

#from compare.m import compare
import scipy.io
import numpy as np
import copy
#import QUODcarb

data = scipy.io.loadmat('datag.mat')
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
    obs.append([])
    # measurements independent of (T,P)
    obs[i].append({'TC': data['datag'][i][4]})  # (umol/kg)
    obs[i].append({'eTC': 2.01})                # TC error ±2.01 umol/kg
    obs[i].append({'TA': data['datag'][i][6]})
    obs[i].append({'TC': 1.78})                # TA error ±2.01 umol/kg
    obs[i].append({'sal': data['datag'][i][2]}) # PSU
    obs[i].append({'esal': 0.002}) # new = 0.001

    # nutrients P and Si also independent of (T,P)
    obs[i].append({'TP': data['datag'][i][23]})
    obs[i].append({'eTP': data['datag'][i][23] * 0.001})    # 0.01% meas uncertainty
    obs[i].append({'TSi': data['datag'][i][17]})
    obs[i].append({'eTSi': data['datag'][i][17] * 0.001})   # 0.01% meas uncertainty
    obs[i].append({'tp': []})

    # first(T,P)-dependent measurement
    """
    'T': deg C, CTD temp
    'eT': ±0.02 degC
    'P': dbar
    'eP': (max) ± 0.63 dbar 
    """

    first_dep_meas = {'T': data['datag'][i][1], 'eT': 0.02, 'P': data['datag'][i][0], 'eP': 0.63}
    obs[i][10]['tp'].append(first_dep_meas)

    # second(T,P)-dependent measurement
    """
    'T':  degC
    'eT': from cruise report
    'P':  in(i+ad,1);  NOT in situ
    'ph': total scale
    'co3':(µmol/kg)
    'eco3':std ±2µmol/kg
    """
    second_dep_meas = {'T': 25, 'eT': 0.05, 'P': 0.0, 'eP': 0.63, 'ph': data['datag'][i][8], 'eph': 0.0004, 'co3': data['datag'][i][13], 'eco3': 2.0}
    obs[i][10]['tp'].append(second_dep_meas)

    # third (T,P)-dependent measurement
    """
    'T':       degC
    'eT':      from cruise report
    'P':       dbar (surface pressure for pco2)
    'pco2':    (µatm)
    'epco2':   0.21% relative std error (avg)
    """

    third_dep_meas = {'T': 20, 'eT': 0.03, 'P': 0.0, 'eP': 0.63, 'pco2': data['datag'][i, 10], 'epco2': data['datag'][i][10] * 0.0021}
    obs[i][10]['tp'].append(third_dep_meas)

obs = np.array(obs)
obs_backup = copy.deepcopy(obs)

# Q5: All five input
# CT AT pH pCO2 CO3 (Q5) (fid5)

obs = copy.deepcopy(obs_backup)

# newtuple = QUODcarb(obs,opt)
# est = newtuple[0]
# obs = newtuple[1]
# sys = newtuple[2]
# iflag = newtuple[3]

# est05 = est
# Q2: Input pairs
# TA TC (Q2) (fid2)
for i in range(nD):
    obs[i][10]['tp'][1]['ph'] = np.nan
    obs[i][10]['tp'][1]['eph'] = np.nan
    obs[i][10]['tp'][2]['pco2'] = np.nan
    obs[i][10]['tp'][2]['epco2'] = np.nan
    obs[i][10]['tp'][1]['co3'] = np.nan
    obs[i][10]['tp'][1]['eco3'] = np.nan

#est,obs, _, _ = QUODcarb(obs,opt); # [est, obs, sys, iflag]
# est02  = est
# fid2   = 'compare_outs/compare_TC_TA.csv'; 
# tp     = 2 # second tp system for ph in there
#A      = compare(obs,est,opt,tp,1,fid2); # 1 for input pair TA TC

# TA ph (Q2) (fid3)

obs = obs_backup.copy()

for i in range(nD):
    obs[i][0]['TC'] = np.nan
    obs[i][1]['eTC'] = np.nan
    obs[i][10]['tp'][2]['pco2'] = np.nan
    obs[i][10]['tp'][2]['epco2'] = np.nan
    obs[i][10]['tp'][1]['co3'] = np.nan
    obs[i][10]['tp'][1]['eco3'] = np.nan

# est, obs, _, _ = QUODcarb(obs, opt)
#est03 = est
tp      = 2
# fid3    = 'compare_outs/compare_TA_ph.csv'
# [A]     = compare(obs,est,opt,tp,3,fid3)
