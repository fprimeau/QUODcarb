## About
This repository includes software compatible with MATLAB for calculating 
 marine CO2-system variables in the exactly-determined and over-determined 
 cases. The calculation also estimates the posterior uncertainties and 
 requires input measurement uncertainties.

## History
CO2SYS was initially developed by Lewis and Wallace (1998) for MS DOS, 
 later adapted for MS Excel and MATLAB by Pierrot (2006). The code was 
 vectorized, refined, and optimized for computational speed by van Heuven 
 et al. (2011). Options for error propagation were added by Orr et al. 
 (2018). The CO2SYSv3 update by Sharp (2022, 2023) added numerous options 
 and extended capabilities of the MATLAB software. The QUODcarb software
 builds upon these previous versions. 

## Installation and Use
Download the files in this repository and put them together in a directory
 that is accessible by MATLAB. 
To perform CO2-system calculations, use QUODcarb.m as directed in this 
documentation and in the comments at the top of the file. Download the
 compare.m file only if you would like to compare the QUODcarb output to
 CO2SYSv3’s output.

## Quick Start
See example scripts (example1.m, example2.m, example3.m) showing possible
 inputs, including correct formulation and capitalization. The function
 QUODcarb requires an input observations structure (obs) and an input
 options structure (opt). All possible inputs and outputs can be found
 in the document ‘QUODcarb_inputs_outputs’.

## Default Settings
There are a multitude of options, here are all the default settings:
### `opt` structure
- `opt.K1K2 = 4`	K1K2 formulation
 	- Mehrbach et al. (1973) refit by Dickson and Millero (1987)
- `opt.KSO4 = 1`		KSO4 formulaation
 	- Dickson et al. (1990a)
 - `opt.KF = 2`			KF formulation
	- Perez and Fraga (1987)
 - `opt.TB = 2`			Total Borate formulation
	- Lee et al. (2010)
 - `opt.printcsv = 0`		Print output to CSV file?
	- Default setting is off/not printing to csv
 - `opt.printmes = 1`		Print messages to screen?
	- Default setting is on/do print to screen
 - `opt.co2press = 0`		Pressure correction on P2F (FugFac) and K0
	- Default setting is off/not correcting these for pressure
 - `opt.Revelle = 0`		Calculate the Revelle factor?
	- Default setting is off/do not calculate Revelle factor
###	`obs` structure
 - `obs.esal` = 0.002 (1 sigma, PSU)
 - `obs.eTC` = 5 (1 sigma, µmol/kg-SW)
 - `obs.eTA` = 5 (1 sigma, µmol/kg-SW)
 - `obs.tp(i).eph` = 0.02	(1 sigma, -log10 scale)
 - `obs.tp(i).epco2` = 5 (1 sigma, µatm)
 - `obs.tp(i).eco3` = 5 (1 sigma, µmol/kg-SW)

## Tips for Success
- Capitalization DOES matter, please be careful to capitalize or use
     lower case exactly as in the examples. 
	- Example: obs.m(1).ph ≠ obs.m(1).pH
- As with CO2SYS, QUODcarb requires two separate CO2-system measurements. The five possible inputs include: (input at least two)
	- `obs.TC` (± `obs.eTC`)
	- `obs.TA` (± `obs.eTA`)
	- `obs.tp(1).ph` (± `obs.tp(1).eph`)
	- `obs.tp(2).pco2` (± `obs.tp(2).epco2`)
		- (or `obs.tp(2).fco2` ± `obs.tp(2).efco2`)
	- `obs.tp(3).co3` (± `obs.tp(3).eco3`)
        - pCO2 and CO3 can also go into tp(1) if they are at the same
          temperature
        - This does not show the salinity, nor the tp(1) through tp(3)
          temperatures and pressures, please do not forget them.
- Even when working with just TC and TA, tp(1) still needs at least a
     temperature and pressure.
        - This is the same as CO2SYS, just formatted differently.
- The input to QUODcarb must always start at `obs(1)` or else it does not
 initialize properly, do not try to give it `obs(2)` with an empty `obs(1)`.

## Citation
Not published yet, documentation will be updated when available.

## References
Lewis and Wallace (1998)
Sharp et al (2023)
Pierrot (2006)
Van Heuven (2011)
Orr et al (2018)





