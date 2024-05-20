# Alternative raw processing for MRR-Pro - RaProM-Pro.py

RaProM-Pro is a novel MRR processing methodology developed for MRR-Pro (Micro Rain Radar Doppler profiler manufactured by Metek GmbH) with enhanced spectra processing and Doppler dealiasing. RaProM-Pro can work with the <i>spectrum_raw</i> or with <i>spectrum_reflectivity</i> from the Manufacturer's netcdf. RaProM-Pro produces a number of output fields which include equivalent reflectivity (Ze), Doppler fall speed and derived parameters such as spectral width, skewness, kurtosis, a simplified precipitation type classification (drizzle, rain, mixed, snow, graupel and hail) and additional variables depending on the precipitation type.<br/><br/>
**Note1:  More information about the processing of RaProM-Pro, with several examples and results of the bright band features can be found on the article:  Garcia-Benadí A, Bech J, Gonzalez S, Udina M, Codina B. A New Methodology to Characterise the Radar Bright Band Using Doppler Spectral Moments from Vertically Pointing Radar Observations. Remote Sens. 2021, 13, 4323. https://doi.org/10.3390/rs13214323**<br/><br/>
**Note2: The scripts works for MRR-Pro data. There is another version of this program for MRR-2 data called RaProM.py (https://github.com/AlbertGBena/RaProM). More information available at: Garcia-Benadi A, Bech J, Gonzalez S, Udina M, Codina B, Georgis JF (2020). Precipitation Type Classification of Micro Rain Radar Data Using an Improved Doppler Spectral Processing Methodology. Remote Sensing, 12(24), 4113 https://doi.org/10.3390/rs12244113** <br />

More information at: Garcia-Benadí A, Bech J, Gonzalez S, Udina M, Codina B. A New Methodology to Characterise the Radar Bright Band Using Doppler Spectral Moments from Vertically Pointing Radar Observations. Remote Sens. 2021, 13, 4323. https://doi.org/10.3390/rs13214323

## Versions and dependences

The main script is called RaProM-Pro.py and it is available in Python 3.8. The following libraries are necessary::

	numpy , version 1.14.5 or later
	miepython, version 1.3.0 or later (matplotlib is necessary for this library works)
	netCDF4, version 1.2.7 or later (cftime is necessary for this library works)

	Other necessary libaries are datetime, math, os, glob and sys
	
The libraries can be installed with pip, using these sentences:

	pip install numpy
	pip install miepython
	pip install netCDF4
	pip install matplotlib
	pip install cftime
	

**The script works with the MRR-Pro files.**

## How to cite

If you use this script for your publication, please cite as:<br/>
Garcia-Benadí A, Bech J, Gonzalez S, Udina M, Codina B. A New Methodology to Characterise the Radar Bright Band Using Doppler Spectral Moments from Vertically Pointing Radar Observations. Remote Sens. 2021, 13, 4323. https://doi.org/10.3390/rs13214323


## Outputs
The script produces the following outputs from MRR-Pro netcdf data:<br />
**W:** Fall speed with aliasing correction (m s<sup>-1</sup>)<br />
**spectral width:** Spectral width of the dealiased velocity distribution (m s<sup>-1</sup>)<br />
**skewness:** Skewness of the dealiased velocity distribution<br />
**kurtosis:** Kurtosis of the dealiased velocity distribution<br />
**DBPIA:** Path Integrated Attenuation (dB) calculated assuming all hydrometeors are in liquid phase regardless of hydrometeor type classification<br />
**Type:** Predominant hydrometeor type numerical value where possible values are: -20 (hail), -15 (graupel), -10 (snow), 0 (mixed), 5 (drizzle), 10 (rain) and 20 (unknown precipitation)<br />
**LWC:** Liquid Water Content (g m<sup>-3</sup>) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**RR:** Rain Rate (mm h<sup>-1</sup>) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**LWC_all:** Liquid Water Content (g m<sup>-3</sup>) calculated assuming all hydrometeors are in liquid phase regardless of hydrometeor type classification<br />
**RR_all:** Rain Rate (mm h<sup>-1</sup>) calculated assuming all hydrometeors are in liquid phase regardless of hydrometeor type classification<br />
**SR:** Snow Rate (mm h<sup>-1</sup>)<br />
**Z:** Radar reflectivity (dBZ) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**Z_all:** Radar reflectivity (dBZ) calculated assuming all hydrometeors are in liquid phase regardless of hydrometeor type classification <br />
**Za:** Attenuated radar reflectivity (dBZ) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**Ze:** Equivalent radar reflectivity (dBZ)<br />
**Zea:** Attenuated equivalent radar reflectivity (dBZ)<br />
**N(D):** Drop Size Distribution (log10(m<sup>-3</sup> mm<sup>-1</sup>)) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**N(D)_all:** Drop Size Distribution (log10(m<sup>-3</sup> mm<sup>-1</sup>)) calculated assuming all hydrometeors are in liquid phase regardless of hydrometeor type classification<br />
**SNR:** Signal to noise ratio from signal without dealiasing (dB)<br />
**Noise:** Noise from spectra reflectivity (m<sup>-1</sup>)<br />
**N<sub>w</sub>:** Intercept of the gamma distribution normalized to the Liquid Water Content (log10(m<sup>-3</sup> mm<sup>-1</sup>)) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**D<sub>m</sub>:** Mean mass-weighted raindrop diameter (mm) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**N<sub>w</sub>_all:** Intercept of the gamma distribution normalized (log10(m<sup>-3</sup> mm<sup>-1</sup>)) calculated assuming all hydrometeors are in liquid phase regardless of  hydrometeor type classification<br />
**D<sub>m</sub>_all:** Mean mass-weighted raindrop diameter (mm) calculated assuming all hydrometeors are in liquid phase regardless of hydrometeor type classification <br />
**BB<sub>bottom</sub>:** Bright Band bottom height  (m) (above sea level)<br />
**BB<sub>top</sub>:** Bright Band top height (m) (above sea level)<br />
**BB<sub>peak</sub>:** Bright Band peak height (m) (above sea level)<br />
**TyPrecipi:** Precipitation regime numerical value where possible values are: 5 (convective), 0 (transition) and -5 (stratiform) calculated using only liquid hydrometeors according to hydrometeor type classification<br />
**TyPrecipi_all:** Precipitation regime numerical value where possible values are: 5 (convective), 0 (transition) and -5 (stratiform) calculated assuming all hydrometeors are in liquid phase regardless of hydrometeor type classification<br />
<br />


## How to execute the script
The script can be executed from a command line at the system prompt (see MS-Windows example):<br />
<br />
![commandWindow](https://user-images.githubusercontent.com/35369817/67784656-64703d00-fa6c-11e9-94fa-0e616d703168.JPG)
<br />
at the directory where RaProM-Pro.py has been copied:
```
python RaProM-Pro.py

```
The script has some additional command line execution options. Please note that their use may imply a substantial increase of the netcdf output file (see below). Command line options are (more than one is possible, in any order):<br /> 
<br /> 
**<i>-spe3D</i>**: this option saves spectral reflectivity raw values after noise and dealising in the spe3D variable. Which contains the spectral reflectivity in function of time, height and speed dealiased. Output file size increases about a factor of 8.<br />
<br /> 
**<i>-dsd3D</i>**: this option saves the Drop Size Distribution vales in the dsd3D variable which contains the Drop Size Distribution in function of time, height and drop diameters.  Output file size increases about a factor of 4.<br />
<br /> 
**<i>-hxxx</i>**: this option forces the MRR antenna height to be at xxx meters above sea level which is important if the height was not correctly configured in the original raw data file. xxx can be a float or an integer value.<br />
<br /> 
**<i>-Myyy</i>**: this option modifies the MRR radar constant so it will affect Z, RR and other variables. M is the multiplicative bias calculated by comparing the MRR rainfall (R_MRR) with a reference rainfall value such as a rain gauge (R_REF), M=R_MRR/R_REF. M can be a float or an integer, typically close to 1.<br />
<br /> 


The syntax of these options are:

```
python RaProM-Pro.py -spe3D

```
```
python RaProM-Pro.py -dsd3D

```
```
python RaProM-Pro.py -h100.3

```
This example forces the antenna height to be at 100.3 m above sea level.<br />
```
python RaProM-Pro.py -M0.78

```
This example assumes a multiplicative bias of 0.87 between MRR and reference rainfall.<br />


The script asks the directory where the netcdf files to be processed are located (it will process all the MRR-Pro netcdf files of the selected folder), for example:
```
C:\mrrPro\test\
```
**NOTE 1: The path must end with \\ in Windows or a / in Linux**<br />
**NOTE 2:  Please avoid blank spaces and special characters in your file path**<br />
**NOTE 3: In macOS, depending on your path environment configuration, it may not be necessary to indicate the complete path so "./" may be enough**<br />

The script indicates the number of netcdf files in the folder and starts the process.

The result is stored in a netcdf file with the same name but finished "-processed"


## Contact
If you have any question, please contact with Albert at albert.garcia@meteo.ub.edu  or   albert.garcia-benadi@upc.edu
