
# Hinode/EIS Observations

Here we have collected some observed intensities for the five Fe XIII lines of interest. All of
the observations are taken from a single raster of an active region (eis_l1_20130708_002042).

* eis_l1_20130708_002042.fe_density.h5 -> HDF5 file containing the intensities
* eis_l1_20130708_002042.fe_density.a.jpg and eis_l1_20130708_002042.fe_density.b.jpg -> plots of the intensities ratios
* test_read_fe_13_density.py -> A python routine illustrating how the intensities can be read
* fe_13.monte_carlo.no_error.h5 -> CHIANTI atomic data for these lines with no errors
* fe_density.txt -> A text file containing the intensities and errors, for cross-checking with the HDF5 file
* fe_13_density.pro -> The IDL routine used to create the HDF5 file; requires special routines to run and is provided for informational purposes.