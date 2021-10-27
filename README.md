# Alternative raw processing by MRR-Pro - RaProM-Pro.py

RaProM-Pro is a novel MRR processing methodology eveloped for MRR-Pro (Micro Rain Radar Doppler profiler manufactured by Metek GmbH) with enhanced spectra processing and Doppler dealiasing. RaProM-Pro produces a number of output fields which include equivalent reflectivity (Ze), Doppler fall speed and derived parameters such as spectral width, skewness, kurtosis, a simplified precipitation type classification (drizzle, rain, mixed, snow, graupel and hail) and additional variables depending on the precipitation type.<br/><br/>
**Note1:  More information about the processing of RaProM-Pro, with several examples and results of the bright band features can be found on the article:  Garcia-Benadí A, Bech J, Gonzalez S, Udina M, Codina B. A New Methodology to Characterise the Radar Bright Band Using Doppler Spectral Moments from Vertically Pointing Radar Observations. Remote Sens. 2021, 13, 4323. https://doi.org/10.3390/rs13214323**<br/>
**Note2: The scripts works for MRR-Pro data. There is another version of this program for MRR-2 data called RaProM.py (https://github.com/AlbertGBena/RaProM). More information available at: Garcia-Benadi A, Bech J, Gonzalez S, Udina M, Codina B, Georgis JF (2020). Precipitation Type Classification of Micro Rain Radar Data Using an Improved Doppler Spectral Processing Methodology. Remote Sensing, 12(24), 4113 https://doi.org/10.3390/rs12244113** <br />

More information at: Garcia-Benadí A, Bech J, Gonzalez S, Udina M, Codina B. A New Methodology to Characterise the Radar Bright Band Using Doppler Spectral Moments from Vertically Pointing Radar Observations. Remote Sens. 2021, 13, 4323. https://doi.org/10.3390/rs13214323

## Versions and dependences

The main script is called RaProM-Pro.py and it is available in Python 3.8. The following libraries are necessary::

	numpy , version 1.14.5 or later

	miepython, version 1.3.0 or later

	netCDF4, version 1.2.7 or later

The libraries can be installed with pip, using these sentences:

	pip install numpy
	pip install miepython
	pip install netCDF4

**The script works with the MRR-Pro files.**

## How to cite

If you use this script for your publication, please cite as:<br/>
Garcia-Benadí A, Bech J, Gonzalez S, Udina M, Codina B. A New Methodology to Characterise the Radar Bright Band Using Doppler Spectral Moments from Vertically Pointing Radar Observations. Remote Sens. 2021, 13, 4323. https://doi.org/10.3390/rs13214323


## Outputs
The script produces the following outputs from MRR-Pro netcdf data:<br />
**W:** fall speed with aliasing correction (m s<sup>-1</sup>)<br />
**spectral width:** spectral width of the dealiased velocity distribution (m s<sup>-1</sup>)<br />
**skewness:** skewness of the dealiased velocity distribution<br />
**kurtosis:** kurtosis of the dealiased velocity distribution<br />
**DBPIA:** Path Integrated Attenuation without the hydrometeor type consideration<br />
**Type:** Hydrometeor type (unknown[20], rain [10], drizzle [5], mixed [0], snow [-10], graupel [-15] and hail [-20])<br />
**LWC:** Liquid water content (g m<sup>-3</sup>)<br />
**RR:** Rain rate (mm h<sup>-1</sup>)<br />
**LWC_all:** Liquid water content (g m<sup>-3</sup>) supposing that all hydrometeors are in liquid phase<br />
**RR_all:** Rain rate (mm h<sup>-1</sup>) supposing that all hydrometeors are in liquid phase<br />
**SR:** Snow rate (mm h<sup>-1</sup>)<br />
**Z:** Reflectivity considering only liquid drops (dBZ)<br />
**Z_all:** Reflectivity (dBZ) supposing that all hydrometeors are in liquid phase <br />
**Za:** Attenuated Reflectivity considering only liquid drops (dBZ)<br />
**Ze:** Equivalent Reflectivity (dBZ)<br />
**Zea:** Attenuated Equivalent Reflectivity (dBZ)<br />
**N(D):** Drop Size Distribution (log10(m<sup>-3</sup> mm<sup>-1</sup>)) only for liquid type<br />
**N(D)_all:** Drop Size Distribution (log10(m<sup>-3</sup> mm<sup>-1</sup>)) suposing that all hydrometeor are in liquid phase<br />
**SNR:** Signal noise relation from signal without dealiasing (dB)<br />
**Noise:** Noise from spectra reflectivity (m<sup>-1</sup>)<br />
**N<sub>w</sub>:** Intercept of the gamma distribution normalized to the liquid water content (log10(m<sup>-3</sup> mm<sup>-1</sup>))<br />
**D<sub>m</sub>:** Mean mass-weighted raindrop diameter (mm)<br />
**N<sub>w</sub>_all:** Intercept of the gamma distribution normalized (log10(m<sup>-3</sup> mm<sup>-1</sup>)) supposing that all hydrometeors are in liquid phase<br />
**D<sub>m</sub>_all:** Mean mass-weighted raindrop diameter (mm) supposing that all hydrometeors are in liquid phase <br />
**BB<sub>bottom</sub>:** Bright Band bottom height  (m) (above sea level)<br />
**BB<sub>top</sub>:** Bright Band top height (m) (above sea level)<br />
**BB<sub>peak</sub>:** Bright Band peak height (m) (above sea level)<br />
**TyPrecipi:** Rainfall type where the value 5 is convective, 0 is transition and -5 is stratiform<br />
**TyPrecipi_all:** Rainfall type supposing that all hydrometeors are in liquid phase where the value 5 is convective, 0 is transition and -5 is stratiform<br />
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
