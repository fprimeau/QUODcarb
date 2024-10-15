## About
This repository includes QUODcarb software compatible with MATLAB for 
 calculating marine CO2-system variables in the exactly-determined and
 over-determined cases. QUODcarb, as presented in Fennell & Primeau 
 (2024), stands for Quantifying Uncertainty in an Over-Determined marine
 CARBonate dataset. The calculation also estimates the output 
 uncertainties and requires input measurement uncertainties. See 
 Fennell & Primeau (2024) for methods and introduction to the solver. 
 
 We intend for QUODcarb to have the
 same capabilities as CO2SYSv3, in addition to the following improvements:
 - Allows over-determined input
 	- i.e. two, three, four, or five (and beyond) measurements of a single water parcel may be input into a single QUODcarb calculation
 - Allows multiple input temperatures and pressures
 	-  i.e. pH at 25 Celsius and pCO2 at 20 Celsius input in tp(1) and tp(2) for a single QUODcarb calculation
 - Uncertainty calculation inherent to the solver
 	- allows for self-consistent uncertainty estimates in the over-determined cases
  - Possibility for the output values to change based on the measured system for the equilibrium constants and mass balance totals
  	-  pK1, pK2, TB etc., in the over-determined cases only
  - Possibility to turn 'off' the pK1, pK2, or TB formulations
  	- when three or more measurements are available AT A SINGLE TEMPERTURE AND PRESSURE (tp(1) only, no tp(2))
   	- see top of `src/QUODcarb.m` for full list of input options

## History
CO2SYS was initially developed by Lewis and Wallace (1998) for MS DOS, 
 later adapted for MS Excel and MATLAB by Pierrot (2006). The code was 
 vectorized, refined, and optimized for computational speed by van Heuven 
 et al. (2011). Options for error propagation were added by Orr et al. 
 (2018). The CO2SYSv3 update by Sharp (2022, 2023) added numerous options 
 and extended capabilities of the MATLAB software. The QUODcarb software
 builds upon these previous versions. See `resources/CO2SYS_code_notes.pdf` for a more in-depth description of the updates to the code by the various previous authors.

## Installation and Use
Download the files of interest in the `src` repository and put them together in a directory
 that is accessible by MATLAB. 
To perform CO2-system calculations, use `src/QUODcarb.m` as directed in this 
documentation and in the comments at the top of the file. The `src/QUODcarb.m` file is a standalone function, it does not need the other files to function. Download the example files only if you are interested in trying them yourself or interested in seeing the example syntax. Download the
 `src/compare.m` file only if you would like to compare the QUODcarb output to
 CO2SYSv3’s output (see `src/example4.m` for example use of `src/compare.m`).

## Syntax Example
QUODcarb requires on input an `obs` structure and an `opt` structure with, respectively, the input conditions and measurements and the options for equilibrium constants etc. 
- `obs` structure with temperature, pressure, etc. may be added to in a for loop at each `obs(i).TC` etc.
	- `obs.sal = salinity;`
	- `obs.usal = sal_uncertainty;` % (± sigma)
	- `obs.TC = total_carbon;` % DIC (umol/kg)
	- `obs.uTC = TC_uncertainty;` % (± sigma)
	- `obs.TA = total_alk;`	% (umol/kg)
	- `obs.uTA = TA_uncertainty;` % (± sigma)
	- `obs.tp(1).T = temp;` % deg Celcius
	- `obs.tp(1).uT = temp_uncertainty;` % (± sigma)
	- `obs.tp(1).P = pressure;` % (dbar)
	- `obs.tp(1).uP = press_uncertainty;` % (± sigma)
	- `obs.tp(1).ph = ph_meas;` % any ph scale
	- `obs.tp(1).uph = ph_uncertainty;` % (± sigma)

- `opt` structure with the input options, the equilibrium constant options are the same as for CO2SYSv3
	- `opt.K1K2 = 10;` % Lueker et al., 2000
	- `opt.KSO4 = 1;` % Dickson et al., 1990a
	- `opt.KF = 2;`	% Perez and Fraga, 1987
	- `opt.TB = 1;`	% Lee et al., 2010
	- `opt.phscale = 1;` % (1=tot, 2=free, 3=sws, 4=nbs)
	- `opt.printcsv	= 1;` % (1=on, 0=off)
	- `opt.fname = 'output.csv';` % CSV filename
	- `opt.printmes = 1;` % (1=on, 0=off)
	- `opt.co2press = 1;` % (1=on, 0=off)
	- `opt.Revelle = 1;` % (1=on, 0=off) Calculate Revelle factor?
 - see the top of the `src/QUODcarb.m` file for all possible input options
 - see `resources/QUODcarb_inputs_outputs.pdf` for a list of the possible input and outut parameters and variable names

## Quick Start
See example scripts (`src/example1.m`, `src/example2.m`, `src/example3.m`) showing possible
 inputs, including correct formulation and capitalization. The function
 QUODcarb requires an input observations structure (`obs`) and an input
 options structure (`opt`). All possible inputs and outputs can be found
 in the document `resources/QUODcarb_inputs_outputs`. To see how we parsed the data
 for the paper, take a look at `paper_scripts/driver.m` (can download the 
 GOMECC-3 dataset via `paper_scripts/parse_gomecc3_data.m` to try it out for yourself).

## Recreate Fennell & Primeau 2024 Figures
To recreate the figures from the first QUODcarb paper, see the directory `paper_scripts`. Start by running `paper_scripts/parse_gomecc_data.m` to load the GOMECC-3 Dataset to your local computer. Then run `paper_scripts/driver.m` to run QUODcarb through all 1119 datapoints 26 times, one for each possible calculation. This requires making new directories on your device to store the output .mat files. The code takes ~4hrs to run on a remote device and 20GB of storage. Once the 26 .mat files are loaded you may choose any of the plotting scripts to plot the desired figures.

## Tips for Success
- Capitalization DOES matter, please be careful to capitalize and/or use
     lower case exactly as in the examples. 
	- Example: obs.m(1).ph ≠ obs.m(1).pH
 	- See the full list in `resources/QUODcarb_inputs_outputs.pdf`
- As with CO2SYS, QUODcarb requires two separate CO2-system measurements. The possible inputs include: (input at least two)
	- `obs.TC` (± `obs.uTC`)
	- `obs.TA` (± `obs.uTA`)
	- `obs.tp(1).ph` (± `obs.tp(1).uph`)
	- `obs.tp(2).pco2` (± `obs.tp(2).upco2`)
		- or `obs.tp(2).fco2` ± `obs.tp(2).ufco2`
  		- or `obs.tp(2).co2st` ± `obs.tp(2).uco2st`
	- `obs.tp(3).co3` (± `obs.tp(3).uco3`)
 		- or `obs.tp(3).hco3` ± `obs.tp(3).uhco3`
        - pCO2 and CO3 can also go into tp(1) if they are at the same
          temperature
        - This does not show the salinity, nor the tp(1) through tp(3)
          temperatures and pressures, please do not forget them.
- Even when working with just TC and TA, tp(1) still needs at least a
     temperature and pressure.
        - This is the same as CO2SYS, just formatted differently.
- The input to QUODcarb must always start at `obs(1)` or else it does not
 initialize properly, do not try to give it `obs(2)` with an empty `obs(1)`.
- The input uncertainty must be a non-zero number as -log10(0) is undefined and breaks the code.

## Citation
Not published yet, documentation will be updated when available.

## References
- For a full list of references, please see `resources/Works_Cited.pdf`. This includes all the possible formulations and other miscellaneous works that are important to the history and functionality of the calculations of the CO2 system in seawater. 

Lewis, E., & Wallace, D. W. R. (1998). Program developed for CO2 system calculations. ORNL/CDIAC-105. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.

Sharp, J. D., Pierrot, D., Humphreys, M. P., Epitalon, J. M., Orr, J. C., Lewis, E. R., & Wallace, D. W. R. (2023). CO2SYSv3 for MATLAB (Version v3.2.1). http://doi.org/10.5281/zenodo.3950562

Pierrot, D., Lewis, E., & Wallace, D. W. R. (2006). CO2SYS DOS Program developed for CO2 system calculations. ORNL/CDIAC-105. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, US Department of Energy, Oak Ridge, TN.

Van Heuven, S., Pierrot, D., Rae, J. W. N., Lewis, E., & Wallace, D. W. R. (2011). MATLAB Program developed for CO2 system calculations. ORNL/CDIAC-105b. Carbon Dioxide Information Analysis Center, Oak Ridge National Laboratory, Oak Ridge, TN.

Orr. J. C., Epitalon, J. M., Dickson, A. G., & Gattuso, J. P. (2018). Routine uncertainty propagation for the marine carbon dioxide system. Marine Chemistry 207, 84-107.





