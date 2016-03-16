# ISSIAtomicData

The goal here is to compute different realizations of the CHIANTI atomic data for use in
interpreting observed density sensitive line ratios. The calculations so far are simply proof of
concept to illustrate the procedure and to provide example data files.

See chianti_ratio_mc.md for information.

atomic_data_v1_hpw: First version of perturbed atomic data

* chianti_ratio_mc.pro -> driver routine
* setup_ion.pro -> modified CHIANTI routine for perturbing collision strengths
* read_wgfa2.pro -> modified CHIANTI routine for perturbing A values
* nrl_save_hdf.pro -> IDL routine for writing HDF5
* nrl_restore_hdf.pro -> IDL routine for reading HDF5
* test_chianti_mc.py -> test python script for reading HDF5 files
* test_chianti_mc.R -> test R script for reading HDF5 files
* fe_13.monte_carlo.h5 -> 100 realizations for Fe XIII (selected lines)
* o_7.monte_carlo.h5 -> 100 realizations for O VII (selected lines)
* o_8.monte_carlo.h5 -> 100 realizations for O VIII (selected lines)

atomic_data_v2_gdz: Second version of perturbed atomic data
