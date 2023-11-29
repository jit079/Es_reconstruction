# Es_reconstruction

The scripts are written in Python and use the libraries: numpy, matplotlib, netCDF4, and h5py. The versions I used on my computer are: Python - 3.10.8, numpy - 1.22.1, matplotlib - 3.7.0, netCDF4 - 1.6.2, h5py - 3.6.0.

The code reconstructs hyperspectral downwelling solar irradiance from 310.25 to 899.75 nm using only multispectral Ed (10 nm resolution) at 412, 489, 555, and 705 nm, using SeaBird OCR sensor. Users need to define the output wavelength that they would like to reconstruct Ed at. The finest reconstructed resolution can be 0.5 nm.

## Files needed to run the program

**Software Files (in the "code" folder):**

_reconstruct_Es.py_: includes the functions to read linear coefficients, read Tg, read extraterrestrial solar irradiance E0, read ancillary data from MERRA2, and the function to reconstruct Ed.

_luts.py_: includes the functions to do look-up table interpolation

_example.py_: demonstrates how to run the code using an example

**Data files (in the "auxdata" folder):**

Data files that are necessary to successfully run the scripts mainly consist of (1) one Look-Up Table (LUT), (2) a txt file containing linear coefficients derived for 310.25 - 899.75 nm (every 0.5 nm), and (3) a txt file containing the extraterrestrial solar irradiance (1st column - wavelength (nm), 2nd column - E0 (W/m2/nm). They are placed in the “auxdata” folder. 

## How to run the program

1. Download the MERRA2 files corresponding to the date when the multispectral Ed data is collected. (The MERRA2 data for the example code has been provided in the 'merra2/' folder.)

2. To run the example script, _“python pace_example.py”_.

## Output

The reconstructed Ed data is saved in 'output/' as a txt file (1st column - wavelength (nm), 2nd columm - reconstructed Ed (W/m2/nm).

