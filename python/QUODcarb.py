import mksys
import numpy as np
from scipy.linalg import inv, sqrtm
import PrintCSV

def QUODcarb(obs, opt):         #returns [est,obs,sys,iflag]
    """
    OUTPUT:
    est := posterior estimates of co2-system variables, equilibrium constants, and precisions
    obs := same as input except that the pK's have been added
    iflag := 0 if solver converged to specified accuracy 
            1 after reaching maximum number of iterations without converging
            2 if it was one order of magnitude away from converging at
                maximum iteration
    CSV := optional CSV output, if turned on in options

    INPUT:
    obs  := co2-system measured quantities with precisions 
    opt  := solver options input by user 

    -------------------------------------------------------------------------

    SYNTAX example:
    obs.sal         = salinity;     (PSU)           
    obs.esal        = sal_error;    (±sigma)        
    obs.TC          = total_c;      (umol/kg-SW)    
    obs.eTC         = TC_error;     (±sigma)        
    obs.TA          = total_alk;    (umol/kg-SW)    
    obs.eTA         = alk_error;    (±sigma)        
    obs.tp(1).T     = temp;         (deg C)         
    obs.tp(1).eT    = temp_error;   (±sigma)        
    obs.tp(1).P     = pressure;     (dbar)          
    obs.tp(1).eP    = pres_error;   (±sigma)     
    obs.tp(1).ph    = ph_meas;      
    obs.tp(1).eph   = ph_error;     (±sigma)

    opt.K1K2        = 10;           % (Lueker et al 2000)
    opt.KSO4        = 1;            % (Dickson et al 1990a) 
    opt.KF          = 2;            % (Perez and Fraga 1987
    opt.TB          = 2;            % (Lee et al. 2010)
    opt.phscale     = 1;            % (1=tot, 2=free, 3=sws, 4=nbs)
    opt.printcsv    = 1;            % (1=on, 0=off)
    opt.fname       = 'output.csv'; % (CSV filename)
    opt.printmes    = 1;            % (1=on, 0=off)
    opt.co2press    = 1;            % (1=on, 0=off)
    opt.Revelle     = 1;            % (1=on, 0=off)

    --------------------------------------------------------------------------

    INPUT OPTIONS:
    opt.K1K2  -> choice of K1 and K2 formulation
            1 = Roy et al, 1993
            2 = Goyet & Poisson, 1989
            3 = Hansson, 1973          REFIT by Dickson and Millero, 1987
            4 = Mehrbach et al., 1973  REFIT by Dickson and Millero, 1987
            5 = Hansson, 1973 and Mehrbach, 1973 
                                        REFIT by Dickson and Millero, 1987
        x(6) = x(GEOSECS)            ~NOT AVAILABLE IN QUODcarb~
        x(7) = x(Peng)               ~NOT AVAILABLE IN QUODcarb~
        x(8) = x(Millero, 1979)      ~NOT AVAILABLE IN QUODcarb~
            9 = Cai and Wang, 1998
            10 = Lueker et al., 2000    (DEFAULT)
            11 = Mojica Prieto and Millero, 2002
            12 = Millero et al., 2000
            13 = Millero et al., 2002
            14 = Millero et al., 2006
            15 = Waters, Millero, and Woosley, 2014
            16 = Sulpis et al., 2020
            17 = Schockman and Byrne, 2021

    opt.KSO4  -> choice of KSO4 formulation
            1 = Dickson (1990a)         (DEFAULT)
            2 = Khoo et al., 1977
            3 = Waters and Millero, 2013

    opt.KF    -> choice of KF formulation
            1 = Dickson and Riley, 1979
            2 = Perez and Fraga, 1987   (DEFAULT)

    opt.TB    -> choice of total borate formulation
            1 = Uppstrom, 1979
            2 = Lee et al., 2010        (DEFAULT)

    opt.co2press -> turn on or off the pressure dependencies for K0 and
            pCO2 to fCO2 fugacity factor (p2f)
            1 = on                      (DEFAULT)
            2 = off

    --------------------------------------------------------------------------

    OUTPUT:
    est ->  'est' structure with best estimate contains:
                1. p(value) and p(error) where p(x) = -log10(x)
                2. value and average error about the value in 'q'
                    where q(x) = x^(-10)
                3. upper and lower bounds in 'q' space, not symmetric
                    about the value in 'q' space
    csv ->  csv file with most of est populated in a spreadsheet, 
                    contains column headers with labels and units                 
                    -does not include upper and lower errors

    --------------------------------------------------------------------------

    Changes? -> the only things you may want to change are: 
                1. tolerance level of Newton solver -> line 111
                2. Max Iteration number -> MAXIT in newtn.m

    --------------------------------------------------------------------------
    """

    opt = check_opt(opt)                   # check opt structure

    sys = mksys(obs[0], opt['phscale'])        # create indexing
    nD  = len(obs)
    nv  = sys.K.shape[1]

    # populate obs, yobs, wobs at each datapoint
    [obs,yobs, wobs] = parse_input(obs,sys,opt,nD)

    for i in range(len(nD)): # loop over the full data set

        z0      = init(yobs[i,:],sys)           # initialize
        tol     = 2e-6                          # tolerance
        # negative of the log of the posterior 
        # aka the log improbability (limp for short)
        def fun(z):
            return limp(z, yobs[i, :], wobs[i, :], obs[i], sys, opt)

        # find the maximum of the posterior probability 
        [zhat, J, iflag[i]] = newtn(z0,fun,tol)

        if (iflag[i] != 0) and (opt.printmes != 0):
            print('Newton\'s method iflag = {iflag[i]} at i = {i}\n')
        
        # if (0):
        #     % check derivatives
        #     n = length(zhat)
        #     I = eye(n)
        #     e = sqrt(eps)
        #     Htest = zeros(n,n)
        #     for ii = 1:n
        #         gp = fun(zhat+I(:,ii)*e)
        #         gm = fun(shat-I(:,ii)*e)
        #         H(:,ii) = (gp-gm)/(2*e)
        #     [ii,jj,iv] = find(H)
        #     [iii,jjj,iiv] = find(J)
        #     keyboard
        #     re = abs(iv - iiv)./iiv; % 1e-3

        # calculate the marginalized posterior uncertainty using Laplace's approximation
        C = inv(J)
        C = C[:nv, :nv]
        sigx = np.sqrt(np.diag(C))
        if opt.printmes != 0:
            if np.sum(np.isnan(sigx)) > 0 or np.sum(np.isinf(sigx)) > 0: 
                print('NaN found in output means faulty run. i = {I}\n')

        # populate est
        est[i] = parse_output(zhat, sigx, sys)   

        # calculate the Revelle factor if opt.Revelle = 1
        if opt.Revelle == 1:
            for j in range(len(sys.tp)):
                # Revelle
                ifree   = sys.tp[j].ifree
                I_TA    = speye(len(ifree))
                ei      = np.ones(len(ifree))
                jfree   = sys.tp(j).jfree
                I_TC    = speye(len(jfree))
                ej      = np.ones(len(jfree))
                # dpfCO2dpTA (similar to Revelle but TC held fixed)
                jac     = sys.tp[j].dcdx_pTAfixed(zhat[ifree])
                z       = ei - ( jac.' ) * ( ( jac * jac.' ) \ ( jac*ei ) )
                est[i].tp[j].Revelle = z(2)/z(1)
                jac     = sys.tp(j).dcdx_pTCfixed(zhat(jfree))
                z       = ej - ( jac.' ) * ( ( jac * jac.' ) \ ( jac*ej ) )
                est[i].tp[j].dpfco2dpTA = z(2)/z(1)

    # PrintCSV if opt.printcsv = 1 using filename opt.fname
    PrintCSV(est, obs, iflag, opt)

# -------------------------------------------------------------------------
# Subfunctions
# -------------------------------------------------------------------------

def limp(z,y,w,obs,sys,opt):    # returns [g,H,f]
    """ 
    [g,H,f] = limp(z,y,w,obs,sys,opt)

    Negative log probability for the co2 system  a.k.a. log improbability i.e., limp!

    INPUT:

        y  := measured components, non-measured components set to NaN
        w  := measurement precisions, same size and order as y, non-measured components set to anything, they are ignored
        gpK  := first order derivatives of pK wrt T,S,P
        ggpK := second order derivatives of pK wrt T,S,P

    OUTPUT:
    
        f := limp
        g := grad f w.r.t. x
        h := hessian of f w.r.t. x
    """

    p = sys['p']
    q = sys['q']
    M = sys['M']
    K = sys['K']
    nrk = K.shape[0]
    nTP = len(sys['tp'])
    nv = M.shape[1]
    x = z[:nv]  # measured variables

    # fill ypK, gypK, ggypK, and with associated calculated pK, gpK, and ggpK values
    # update the pK values based on the new estimate of (T,P)

    y, gy, ggy = update_y(y, x, obs, sys, opt)
    # Make a vector of measured quantities
    id = np.where(~np.isnan(y))[0]
    y = y[id]
    gy = gy[id, :]
    ggy = ggy[id, :, :]

    # Make a precision matrix
    W = np.diag(w[id])

    # Build a matrix that Picks out the measured components of x
    I = np.eye(nv)  # for chain rule
    PP = I[id, :]  # picking/pick out the measured ones
    e = np.dot(PP, x) - y  # calculated - measured (minus)
    ge = PP - gy

    nlam = M.shape[0] + K.shape[0]
    lam = z[nv:]  # Lagrange multipliers

    # constraint equations
    c = np.concatenate((np.dot(M, q(x)), np.dot(K, x)))
    f = 0.5 * np.dot(np.dot(e.T, W), e) + np.dot(lam.T, c)  # limp, method of lagrange multipliers
    # -(-1/2 sum of squares) + constraint eqns, minimize f => grad(f) = 0
    dcdx = np.concatenate((np.dot(M, np.diag(sys['dqdx'](x))), K))
    #iTSP = [[sys.tp(:).iT], sys.isal, [sys.tp(:).iP]]
    g = np.concatenate((np.dot(np.dot(e.T, W), ge) + np.dot(lam.T, dcdx), c))

    ddq = np.diag(sys['d2qdx2'](x))  # q"

    nr, nc = M.shape
    gg = np.zeros((nc, nc))
    for row in range(nr):
        gg += lam[row] * np.diag(M[row, :]) * ddq

    tmp = np.zeros((len(x), len(x)))
    eW = np.dot(e.T, W)
    for jj in range(ggy.shape[0]):
        tmp += eW[jj] * ggy[jj]

    H = np.block([[np.dot(np.dot(ge.T, W), ge) - tmp + gg, dcdx.T],  # derivatives wrt lambdas
                  [dcdx, np.zeros((nlam, nlam))]])  # derivatives wrt var's
    g = g.reshape(-1, 1)  # make sure g is returned as a column vector

    return f, g, H

# -----------------------------------------------------------------------------------

def check_opt():
    # check opt input
    isbad = lambda thing: (thing is None) or (sum(np.isnan(thing)) > 0)

    # opt.printmes
    if 'printmes' not in opt or isbad(opt['printmes']):
        opt['printmes'] = 1  # default on

    # opt.K1K2
    if 'K1K2' not in opt or isbad(opt['K1K2']):
        if opt['printmes'] != 0:
            print('No K1K2 formulation chosen. Assuming opt.K1K2 = 4')
        opt['K1K2'] = 4  # default K1K2 setting
    elif opt['K1K2'] > 18 or opt['K1K2'] < 1 or opt['K1K2'] in [6, 7, 8]:
        if opt['printmes'] != 0:
            print('Invalid K1K2 formulation chosen. Assuming opt.K1K2 = 4')
        opt['K1K2'] = 4  # default K1K2 setting
    
    # opt.TB
    if 'TB' not in opt or isbad(opt['TB']):
        if opt['printmes'] != 0:
            print('No K1K2 formulation chosen. Assuming opt.K1K2 = 4')
        opt['K1K2'] = 4  # default K1K2 setting
    elif opt['K1K2'] > 18 or opt['K1K2'] < 1 or opt['K1K2'] in [6, 7, 8]:
        if opt['printmes'] != 0:
            print('Invalid K1K2 formulation chosen. Assuming opt.K1K2 = 4')
        opt['K1K2'] = 4  # default K1K2 setting

    # opt.KSO4
    if 'KSO4' not in opt or isbad(opt['KSO4']):
        if opt.get('printmes', 0) != 0:
            print('No KSO4 formulation chosen. Assuming opt.KSO4 = 1')
        opt['KSO4'] = 1  # default opt.KSO4 setting
    elif opt['KSO4'] > 3 or opt['KSO4'] < 1:
        if opt.get('printmes', 0) != 0:
            print('Invalid KSO4 formulation chosen. Assuming opt.KSO4 = 1')
        opt['KSO4'] = 1  # default opt.KSO4 setting

    # opt.KF
    if 'KF' not in opt or isbad(opt['KF']):
        if opt.get('printmes', 0) != 0:
            print('No KF formulation chosen. Assuming opt.KF = 2')
        opt['KF'] = 2  # default KF
    elif opt['KF'] > 2 or opt['KF'] < 1:
        if opt.get('printmes', 0) != 0:
            print('Invalid KF formulation chosen. Assuming opt.KF = 2')
        opt['KF'] = 2

    # opt.phscale
    if 'phscale' not in opt or isbad(opt['phscale']):
        raise ValueError('No opt.phscale chosen, must choose 1 = tot, 2 = sws, 3 = free, 4 = NBS')
    elif opt['phscale'] > 4 or opt['phscale'] < 1:
        raise ValueError('Invalid opt.phscale chosen, must choose 1 = tot, 2 = sws, 3 = free, 4 = NBS')

    # opt.printcsv and opt.fname
    if 'printcsv' not in opt or isbad(opt['printcsv']):
        opt['printcsv'] = 0  # default off
    elif opt['printcsv'] > 1 or opt['printcsv'] < 0:
        if opt.get('printmes', 0) != 0:
            print('Invalid CSV opt chosen. Assuming opt.csv = 1')
    else:
        if 'fname' not in opt or isbad(opt['fname']):
            opt['fname'] = 'QUODcarb_output.csv'
            if opt.get('printmes', 0) != 0:
                print("Invalid CSV filename. Assuming opt.fname = 'QUODcarb_output.csv'")
    
    # opt.co2press
    if 'co2press' not in opt or isbad(opt['co2press']):
        opt['co2press'] = 1  # on
        if opt.get('printmes', 0) != 0:
            print('No opt.co2press chosen. Assuming opt.co2press = 1 (on).')

    # opt.Revelle
    if 'Revelle' not in opt or isbad(opt['Revelle']):
        opt['Revelle'] = 0
        if opt['printmes'] != 0:
            print('No opt.Revelle chosen. Assuming opt.Revelle = 0 (off).')

    # opt.turnoff
    return opt

# ------------------------------------------------------------------------

def parse_input(obs, sys, opt, nD):
    isgood = lambda thing: (thing is not None) and (not np.isnan(thing))
    p = sys["p"]
    q = sys["q"]
    # convert x+/-e into precision for p(x) (precision = 1/variance)
    w = lambda x, e: np.abs(p(1 + e/x))**(-2)

    nv = sys["K"].shape[1]
    yobs = np.full((nD, nv), np.nan)
    wobs = np.full((nD, nv), np.nan)

    if "tp" not in obs:
        if opt["printmes"] != 0:
            raise ValueError("Need to provide temperature and pressure measurement.")
    
    nTP = len(obs[0]["tp"])
    for i in range(nD):  # loop over all the stations
        # make sure all the required fields in the obs struct exist
        if "sal" not in obs[i]:
            if opt["printmes"] != 0:
                raise ValueError("Need to provide salinity measurement.")
        else:
            yobs[i, sys["isal"]] = obs[i]["sal"]

        if "esal" not in obs[i]:
            obs[i]["esal"] = 0.002  # std = 0.002 PSU
            wobs[i, sys["isal"]] = obs[i]["esal"]**(-2)
            if opt["printmes"] != 0:
                print("Warning: Assuming salinity uncertainty is 0.002 PSU")
        else:
            wobs[i, sys["isal"]] = obs[i]["esal"]**(-2)  # std e -> w

        if 'TC' not in obs[i] or not obs[i]['TC']:
            obs[i]['TC'] = np.nan
            yobs[i, sys['ipTC']] = np.nan
        else:
            yobs[i, sys['ipTC']] = p(obs[i]['TC'] * 1e-6)  # Assuming p is a function
            
        if "TC" not in obs[i] or not obs[i]['TC']:
            obs[i]['TC'] = np.nan #[]
            yobs[i]['ipTC'] = np.nan
        else:
            yobs[i, sys["ipTC"]] = p((obs(i).TC*1e-6)) #convt to mol/kg

        if 'eTC' not in obs[i] or not obs[i]['eTC']:
            obs[i]['TC'] = np.nan
            wobs[i, sys['ipTC']] = np.nan
        else:
            wobs[i, sys['ipTC']] = w(obs[i]['TC'], obs[i]['eTC']) #std e -> w

        if 'TA' not in obs[i] or not isgood(obs[i]['TA']):
            obs[i]['TA'] = np.nan #[]
            yobs[i, sys['ipTA']] = np.nan
        else:
            yobs[i, sys['ipTA']] = p(obs[i]['TA'] * 1e-6) #convt to mol/kg

        if 'eTA' not in obs[i] or not isgood(obs[i]['eTA']):
            obs[i]['eTA'] = np.nan
            wobs[i, sys['ipTA']] = np.nan
        else:
            wobs[i, sys['ipTA']] = w(obs[i]['TA'], obs[i]['eTA']) #std e -> w

        # calculate totals that are a function of salinity
        pT, _, _, epT = calc_pTOT(opt, obs[i]['sal'])
        pTB, epTB = pT[0], epT[0]
        pTS, epTS = pT[1], epT[1]
        pTF, epTF = pT[2], epT[2]
        pTCa, epTCa = pT[3], epT[3]  # (see Ref's within calc_pTOT)
        #total borate
        if 'TB' not in obs[i] or not isgood(obs[i]['TB']):
            obs[i]['TB'] = np.nan
            yobs[i, sys['ipTB']] = pTB
        else:
            if obs[i]['TB'] == 0:
                obs[i]['TB'] = 1e-3 # umol/kg, reset minimum to 1 nanomolar
            yobs[i, sys['ipTB']] = p(obs[i]['TB'] * 1e-6) # convt µmol/kg to mol/kg
        
        if 'eTB' not in obs[i] or not isgood(obs[i]['eTB']):
            obs[i]['eTB'] = np.nan
            wobs[i, sys['ipTB']] = (epTB)^(-2) # convert to precision
        else:
            wobs[i, sys['ipTB']] = w(obs[i]['TB'], obs[i]['eTB']) #mol/kg

        # total sulfate
        
                        
# stopped at line 455
