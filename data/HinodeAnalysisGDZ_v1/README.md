
Files related to the analysis of EIS spectra using GDZ's perturbed CHIANTI atomic data are located
here. 

* fe_13_fit_intensities.pro: read the EIS intensities and the perturbed emissivities, compute a
  simple isothermal fit to the intensities, save the results to a file.

  - eis_l1_20130708_002042.fe_density.h5: a file containing a set of 1000 observed EIS Fe XIII
    intensities.

  - fe_13.monte_carlo_normal_nsim=1000.h5: Perturbed emissivities

  - fe_13_fit_intensities.normal.h5: output file with saved fits

* fe_13_fit_intensities_plot.pro: display the calculated fits.

* chianti_default_hpw.pro: Compare the perturbed atomic data with the default CHIANTI
  calculation. Requires other routines to run.