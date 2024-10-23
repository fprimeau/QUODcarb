## QUODcarb.m
`function [est,obs,sys,iflag,opt] = QUODcarb(obs,opt)`

OUTPUT:
- `est :=` posterior estimates of co2-system variables, equilibrium constants, and precisions
- `obs :=` same as input except that the variable names are all there
- `sys :=` indexing scheme as created in mksys, useful for debugging
- `iflag := `
  - `0` if solver converged to specified accuracy 
  - `1` after reaching maximum number of iterations without converging
  - `2` if it was one order of magnitude away from converging at maximum iteration
- `opt :=` same as input except that all fields are now filled in 
- `CSV :=` optional CSV output, if turned on in options

INPUT:
- `obs  :=` co2-system measured quantities with precisions 
- `opt  :=` solver options input by user 

SYNTAX example:
- `obs` structure
  - `obs.sal         =` salinity;     (PSU)           
  - `obs.usal        =` sal_uncert;    (±sigma) (default 0.002)     
  - `obs.TC          =` total_c;      (umol/kg-SW)    
  - `obs.uTC         =` TC_uncert;     (±sigma) (default 2 umol/kg)
  - `obs.TA          =` total_alk;    (umol/kg-SW)    
  - `obs.uTA         =` alk_uncert;    (±sigma) (default 2 umol/kg)
  - `obs.tp(1).T     =` temp;         (deg C)         
  - `obs.tp(1).uT    =` temp_uncert;   (±sigma) (default 0.1 degC)
  - `obs.tp(1).P     =` pressure;     (dbar, 0 = surface)          
  - `obs.tp(1).uP    =` pres_uncert;   (±sigma) (default 0.1 dbar)
  - `obs.tp(1).ph    =` ph_meas;      
  - `obs.tp(1).uph   =` ph_uncert;     (±sigma) (default 0.010)
- `opt` structure
  - `opt.K1K2        = 10;`           % (Lueker et al 2000)
  - `opt.KSO4        = 1;`            % (Dickson et al 1990a) 
  - `opt.KF          = 2;`            % (Perez and Fraga 1987
  - `opt.TB          = 2;`           % (Lee et al. 2010)
  - `opt.phscale     = 1;`           % (1=tot, 2=free, 3=sws, 4=nbs)
  - `opt.printcsv    = 1;`            % (1=on, 0=off)
  - `opt.fname       = 'output.csv';` % (CSV filename)
  - `opt.printmes    = 1;`            % (1=on, 0=off)
  - `opt.co2press    = 1;`            % (1=on, 0=off)
  - `opt.Revelle     = 1;`            % (1=on, 0=off)

INPUT OPTIONS:
  - `opt.K1K2`  -> choice of K1 and K2 formulation
      - (T-range)(S-range)
    - `1 =` Roy et al, 1993
      - (T:0-45)    (S:5-45)
    - `2 =` Goyet & Poisson, 1989
      - (T:-1-40)   (S:10-50)
    - `3 =` Hansson, 1973 REFIT by Dickson & Millero, 1987
      - (T:2-35)    (S:20-40)
    - `4 =` Mehrbach et al., 1973 REFIT by Dickson & Millero, 1987
      - (T:2-35)    (S:20-40)
    - `5 =` Hansson, 1973 and Mehrbach, 1973 REFIT by Dickson & Millero, 1987
      - (T:2-35)    (S:20-40)
    - x = x(GEOSECS)  ~NOT AVAILABLE IN QUODcarb~
    - x = x(Peng)  ~NOT AVAILABLE IN QUODcarb~
    - `8 =` Millero, 1979 PURE WATER ONLY
      - (T:0-50)    (S: 0.)
    - `9 =` Cai and Wang, 1998
      - (T:2-35)    (S:0-49)
    - `10 =` Lueker et al., 2000 (DEFAULT)
      - (T:2-35)    (S:19-43)
    - `11 =` Mojica Prieto and Millero, 2002
      - (T:0-45)    (S:5-42)
    - `12 =` Millero et al., 2002
      - (T:-1.6-35) (S:34-37)
    - `13 =` Millero et al., 2006
      - (T:0-50)    (S:1-50)
    - `14 =` Millero et al., 2010
      - (T:0-50)    (S:1-50)
    - `15 =` Waters, Millero, and Woosley, 2014
      - (T:0-50)    (S:1-50)
    - `16 =` Sulpis et al., 2020
      - (T:-1.7-32) (S:31-38)
    - `17 =` Schockman and Byrne, 2021
      - (T:15-35)   (S:19-41)
  - `opt.KSO4` -> choice of KSO4 formulation
    - `1 =` Dickson (1990a) (DEFAULT)
    - `2 =` Khoo et al., 1977
    - `3 =` Waters and Millero, 2013
  - `opt.KF` -> choice of KF formulation
    - `1 =` Dickson and Riley, 1979
    - `2 =` Perez and Fraga, 1987 (DEFAULT)
  - `opt.TB` -> choice of total borate formulation
    - `1 =` Uppstrom, 1979
    - `2 =` Lee et al., 2010 (DEFAULT)
  - `opt.phscale` -> pH scale of choice for calculations, the equilibrium constants are returned on this scale, REQUIRED
    - `1 =` total scale
    - `2 =` free scale
    - `3 =` seawater scale (SWS)
    - `4 =` National Bureau of Standards scale (NBS) 
  - `opt.printcsv` -> option to print the output est structure to a csv file
    - `1 =` on
    - `0 =` off 
  - `opt.fname` -> chosen filename for CSV, if opt.printcsv = 1, must be a .csv file with no spaces
    - `opt.fname = 'QUODcarb_output.csv'` (DEFAULT)
  - `opt.printmes` -> should QUODcarb print messages to your command line? recommended 'on' for new users
    - `1 =` on (DEFAULT)
    - `0 =` off 
  - `opt.co2press` -> turn on or off the pressure dependencies for K0 and pCO2 to fCO2 fugacity factor (p2f)
    - `1 =` on (DEFAULT)
    - `0 =` off
  - `opt.Revelle` -> calculate the Revelle factor? reminder: QUODcarb cannot calculate uRevelle
    - `1 =` on
    - `0 =` off (DEFAULT)
  - `opt.tol` -> change default tolerance to some other value, recommended to stay within 1e-5 and 1e-8
    - `opt.tol = 1e-6;` (DEFAULT) 

OUTPUT:
  - `est` -> 'est' structure with best estimate, as given in 'parse_input', contains:
    1. p(value) and p(error) where p(x) = -log10(x)
    2. value and average error about the value in 'q', where q(x) = x^(-10)
    3. upper and lower bounds in 'q' space, not symmetric about the value in 'q' space
  - `obs` -> 'obs' structure with updated values as given in 'parse_input' so all variable names have a value, nan or otherwise
    - restates what user input and adds omitted fields or default values
  - `sys` -> 'sys' structure with the system's variable names and associated indexing scheme, as created in 'mksys'
  - `iflag` ->
    - 0 if solver converged to specified tolerance (opt.tol)
    - 1 after reaching maximum number of iterations without converging, see MAXIT in newtn.m
    - 2 if it was one order of magnitude away from converging at maximum iteration, consider keeping or rerunning with a larger tolerance
  - `opt` -> 'opt' structure with all variable names as given in 'check_opt'
    - restates what user input and adds omitted fields  
  - `csv` ->  if opt.printcsv = 1, csv file with most of `est` populated in a spreadsheet, contains column headers with labels and units
    - does not contain upper and lower errors


Changes? -> the only things you may want to change are:
- Max Iteration number -> `MAXIT` in newtn.m
- Print Iteration number and F0 value to screen -> `iprint` in newtn.m







