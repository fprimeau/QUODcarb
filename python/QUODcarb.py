import check_opt as co
import mksys as ms
import parse_input as pi
import parse_output as po
import PrintCSV as pc
import init as i
import limp as l
import newtn as n
import numpy as np


def QUODcarb(obs,opt):
#
# OUTPUT:
#   est := posterior estimates of co2-system variables, equilibrium constants, and precisions
#   obs := same as input except that the pK's have been added
# iflag := 0 if solver converged to specified accuracy 
#          1 after reaching maximum number of iterations without converging
#          2 if it was one order of magnitude away from converging at
#               maximum iteration
#   CSV := optional CSV output, if turned on in options
#
# INPUT:
#   obs  := co2-system measured quantities with precisions 
#   opt  := solver options input by user 
#
# -------------------------------------------------------------------------
#
# SYNTAX example:
#   obs.sal         = salinity;     (PSU)           
#   obs.esal        = sal_error;    (±sigma)        
#   obs.TC          = total_c;      (umol/kg-SW)    
#   obs.eTC         = TC_error;     (±sigma)        
#   obs.TA          = total_alk;    (umol/kg-SW)    
#   obs.eTA         = alk_error;    (±sigma)        
#   obs.tp(1).T     = temp;         (deg C)         
#   obs.tp(1).eT    = temp_error;   (±sigma)        
#   obs.tp(1).P     = pressure;     (dbar, 0 = surface)          
#   obs.tp(1).eP    = pres_error;   (±sigma)     
#   obs.tp(1).ph    = ph_meas;      
#   obs.tp(1).eph   = ph_error;     (±sigma)
#
#   opt['K1K2']        = 10;           % (Lueker et al 2000)
#   opt['KSO4']        = 1;            % (Dickson et al 1990a) 
#   opt['KF']          = 2;            % (Perez and Fraga 1987
#   opt['TB']          = 2;            % (Lee et al. 2010)
#   opt['phscale']     = 1;            % (1=tot, 2=free, 3=sws, 4=nbs)
#   opt['printcsv']    = 1;            % (1=on, 0=off)
#   opt['fname']       = 'output.csv'; % (CSV filename)
#   opt['printmes']    = 1;            % (1=on, 0=off)
#   opt['co2press']    = 1;            % (1=on, 0=off)
#   opt['Revelle']     = 1;            % (1=on, 0=off)
#
# --------------------------------------------------------------------------
# 
# INPUT OPTIONS:
#   opt.K1K2  -> choice of K1 and K2 formulation    (T-range)   (S-range)
#           1 = Roy et al, 1993                     (T:0-45)    (S:5-45)
#           2 = Goyet & Poisson, 1989               (T:-1-40)   (S:10-50)
#           3 = Hansson, 1973          REFIT by Dickson & Millero, 1987
#                                                   (T:2-35)    (S:20-40)
#           4 = Mehrbach et al., 1973  REFIT by Dickson & Millero, 1987
#                                                   (T:2-35)    (S:20-40)
#           5 = Hansson, 1973 and Mehrbach, 1973 
#                                      REFIT by Dickson & Millero, 1987
#                                                   (T:2-35)    (S:20-40)
#           x = x(GEOSECS)            ~NOT AVAILABLE IN QUODcarb~
#           x = x(Peng)               ~NOT AVAILABLE IN QUODcarb~
#           x = x(Millero, 1979)      ~NOT AVAILABLE IN QUODcarb~
#           9 = Cai and Wang, 1998                  (T:2-35)    (S:0-49)
#          10 = Lueker et al., 2000    (DEFAULT)    (T:2-35)    (S:19-43)
#          11 = Mojica Prieto and Millero, 2002     (T:0-45)    (S:5-42)
#          12 = Millero et al., 2002                (T:-1.6-35) (S:34-37)
#          13 = Millero et al., 2006                (T:0-50)    (S:1-50)
#          14 = Millero et al., 2010                (T:0-50)    (S:1-50)
#          15 = Waters, Millero, and Woosley, 2014  (T:0-50)    (S:1-50)
#          16 = Sulpis et al., 2020                 (T:-1.7-32) (S:31-38)
#          17 = Schockman and Byrne, 2021           (T:15-35)   (S:19-41)
#
#   opt['KSO4']  -> choice of KSO4 formulation
#           1 = Dickson (1990a)         (DEFAULT)
#           2 = Khoo et al., 1977
#           3 = Waters and Millero, 2013
#
#   opt['KF']    -> choice of KF formulation
#           1 = Dickson and Riley, 1979
#           2 = Perez and Fraga, 1987   (DEFAULT)
#
#   opt['TB']    -> choice of total borate formulation
#           1 = Uppstrom, 1979
#           2 = Lee et al., 2010        (DEFAULT)
#
#   opt['co2press'] -> turn on or off the pressure dependencies for K0 and
#           pCO2 to fCO2 fugacity factor (p2f)
#           1 = on                      (DEFAULT)
#           2 = off
#
# --------------------------------------------------------------------------
#
# OUTPUT:
#    est ->  'est' structure with best estimate contains:
#                1. p(value) and p(error) where p(x) = -log10(x)
#                2. value and average error about the value in 'q'
#                    where q(x) = x^(-10)
#                3. upper and lower bounds in 'q' space, not symmetric
#                    about the value in 'q' space
#    csv ->  csv file with most of est populated in a spreadsheet, 
#                  contains column headers with labels and units                 
#                     -does not include upper and lower errors
# 
# --------------------------------------------------------------------------
# 
#  Changes? -> the only things you may want to change are: 
#                1. tolerance level of Newton solver -> line 111
#                2. Max Iteration number -> MAXIT in newtn.m
# 
# --------------------------------------------------------------------------


    opt = co.check_opt(opt);                   # check opt structure

    sys = ms.mksys(obs(1),opt.phscale);        # create indexing
    nD  = len(obs)  
    nv  = sys['K'].shape[1]

    # populate obs, yobs, wobs at each datapoint
    [obs,yobs,wobs] = pi.parse_input(obs,sys,opt,nD)

    est = []    # list of dictionaries where each dictionary has posterior estimates of CO2 system
                # variables, equilibrium constants, and precisions
    iflag = []  # list of integers indicating solver's status for each dataset point

    for i in range(nD): # loop over the full data set

        z0 = i.init(opt, yobs[i,:],sys)             # initialize
        if 'tol' not in opt:
            tol = 1e-6;                          # tolerance
        else:
            tol = opt['tol']
        # negative of the log of the posterior 
        # aka the log improbability (limp for short)
        fun = lambda z: l.limp(z, yobs[i, :], wobs[i, :], obs[i], sys, opt)

        # find the maximum of the posterior probability 
        [zhat,J,iflag[i]] = n.newtn(z0,fun,tol)
        if iflag[i] != 0 and opt['printmes'] != 0:
            print(f"Newton\'s method iflag = {iflag[i]} at i = {i}")
        [_,_,f] = l.limp(zhat,yobs[i,:],wobs[i,:],obs[i],sys,opt)
        # residual f value, to tack onto est

        # calculate the marginalized posterior uncertainty using Laplace's approximation
        C = np.linalg.inv(J)
        C = C[:nv, :nv]
        sigx = np.sqrt(np.diag(C))
        if opt['printmes'] != 0:
            if np.isnan(sigx).sum() > 0 or np.isinf(sigx).sum() > 0: 
                print(f'NaN found in output means faulty run. i = {i}')

        # populate est
        [est[i]] = po.parse_output(zhat,sigx,sys,f,C);    
        
        # calculate the Revelle factor if opt.Revelle = 1 ('on')
        if opt['Revelle'] == 1:
            for j in range(len(sys['tp'])):
                # Revelle
                ifree   = sys['tp'][j]['ifree']

                ei      = np.zeros(len(ifree))
                ei[0]   = 1
                jac     = sys['tp'][j]['dcdx_pTAfixed'](zhat[ifree])
                z = ei - np.dot(jac.T, np.linalg.solve(np.dot(jac, jac.T), np.dot(jac, ei)))
                est[i]['tp'][j]['Revelle'] = z[1] / z[0]
            
                # dpfCO2dpTA (similar to Revelle but TC held fixed)
                jfree = sys['tp'][j]['jfree']
                ej = np.zeros(len(jfree))
                ej[0] = 1
                jac = sys['tp'][j]['dcdx_pTCfixed'](zhat[jfree])
                z = ej - np.dot(jac.T, np.linalg.solve(np.dot(jac, jac.T), np.dot(jac, ej)))
                est[i]['tp'][j]['dpfco2dpTA'] = z[1] / z[0]

    # PrintCSV if opt.printcsv = 1 using filename opt.fname
    pc.PrintCSV(est,obs,iflag,opt)

    return est,obs,sys,iflag,opt