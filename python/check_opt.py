import numpy as np
# check opt input
def isbad(thing):
    return thing is None or (np.isnan(thing).sum() if isinstance(thing, np.ndarray) else False)

def check_opt(opt):
    # opt['printmes']
    if 'printmes' not in opt or isbad(opt['printmes']):
        opt['printmes'] = 1  # default on
    
    # opt['K1K2']
    if 'K1K2' not in opt or isbad(opt['K1K2']):
        if opt['printmes'] != 0:
            print('No K1K2 formulation chosen. Assuming opt.K1K2 = 4')
        opt['K1K2'] = 4  # default K1K2 setting
    elif opt['K1K2'] > 18 or opt['K1K2'] < 1 or opt['K1K2'] in {6, 7, 8}:
        if opt['printmes'] != 0:
            print('Invalid K1K2 formulation chosen. Assuming opt.K1K2 = 4')
        opt['K1K2'] = 4  # default K1K2 setting

    # opt['TB']
    if 'TB' not in opt or isbad(opt['TB']):
        if opt['printmes'] != 0:
            print('No TB formulation chosen. Assuming opt.TB = 2')
        opt['TB'] = 2
    elif opt['TB'] > 2 or opt['TB'] < 1:
        if opt['printmes'] != 0:
            print('Invalid TB formulation chosen. Assuming opt.TB = 2')
        opt['TB'] = 2

    # opt['KSO4']
    if 'KSO4' not in opt or isbad(opt['KSO4']):
        if opt['printmes'] != 0:
            print('No KSO4 formulation chosen. Assuming opt.KSO4 = 1')
        opt['KSO4'] = 1  # default opt.KSO4 setting
    elif opt['KSO4'] > 3 or opt['KSO4'] < 1:
        if opt['printmes'] != 0:
            print('Invalid KSO4 formulation chosen. Assuming opt.KSO4 = 1')
        opt['KSO4'] = 1  # default opt.KSO4 setting

    # opt['KF']
    if 'KF' not in opt or isbad(opt['KF']):
        if opt['printmes'] != 0:
            print('No KF formulation chosen. Assuming opt.KF = 2')
        opt['KF'] = 2  # default KF
    elif opt['KF'] > 2 or opt['KF'] < 1:
        if opt['printmes'] != 0:
            print('Invalid KF formulation chosen. Assuming opt.KF = 2')
        opt['KF'] = 2

    # opt['phscale']
    if 'phscale' not in opt or isbad(opt['phscale']):
        raise ValueError('No opt.phscale chosen, must choose 1 = tot, 2 = sws, 3 = free, 4 = NBS')
    elif opt['phscale'] > 4 or opt['phscale'] < 1:
        raise ValueError('Invalid opt.phscale chosen, must choose 1 = tot, 2 = sws, 3 = free, 4 = NBS')

    # opt['printcsv'] and opt['fname']
    if 'printcsv' not in opt or isbad(opt['printcsv']):
        opt['printcsv'] = 0  # default off
    elif opt['printcsv'] > 1 or opt['printcsv'] < 0:
        if opt['printmes'] != 0:
            print('Invalid CSV opt chosen. Assuming opt.csv = 1')
    else:
        if 'fname' not in opt or isbad(opt['fname']):
            opt['fname'] = 'QUODcarb_output.csv'
            if opt['printmes'] != 0:
                print('Invalid CSV filename. Assuming opt.fname = "QUODcarb_output.csv"')

    # opt['co2press']
    if 'co2press' not in opt or isbad(opt['co2press']):
        opt['co2press'] = 1  # on
        if opt['printmes'] != 0:
            print('No opt.co2press chosen. Assuming opt.co2press = 1 (on).')

    # opt['Revelle']
    if 'Revelle' not in opt or isbad(opt['Revelle']):
        opt['Revelle'] = 0
        if opt['printmes'] != 0:
            print('No opt.Revelle chosen. Assuming opt.Revelle = 0 (off).')

    if 'mpk' not in opt:
        mpk0 = lambda T, S, P: 1
        mpk1 = lambda T, S, P: 1
        mpk2 = lambda T, S, P: 1
        mpk = lambda T, S, P: [mpk0(T, S, P), mpk1(T, S, P), mpk2(T, S, P), 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        opt['mpk'] = mpk

    if 'gmpk' not in opt:
        gmpk0 = lambda T, S, P: [0, 0, 0]
        gmpk1 = lambda T, S, P: [0, 0, 0]
        gmpk2 = lambda T, S, P: [0, 0, 0]
        gmpk = lambda T, S, P: [
            gmpk0(T, S, P),
            gmpk1(T, S, P),
            gmpk2(T, S, P),
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0]
        ]
        opt['gmpk'] = gmpk

    if 'empk' not in opt:
        empk0 = 1
        empk1 = 1
        empk2 = 1
        empk = [empk0, empk1, empk2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        opt['empk'] = empk

    # opt['turnoff']
    return opt
