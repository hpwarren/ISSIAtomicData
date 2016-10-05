
#### ISSIAtomicData/phase2_20161006/01_chianti_errors

Here we generate many realizations of emissivities computed with CHIANTI. Only the emsissivites for
lines of interest are saved. The primary file of interest is fe_13.monte_carlo_normal_nsim=1000.h5

* RUN_CALC_EMISS_MC: This routine runs everthing and takes no inputs

	- calc_emiss_mc: Loop over realizations, collect emissivities for selected lines, save to HDF.
	
	- emiss_calc: Modified version of standard CHIANTI routine, adds keyword for
		adding uncertainties (either normal or uniform)
		
	- setup_ion: Called by emiss_calc, perturbs A-values and collision rates. Some nuances

		+ Looks for the file 'setup_ion.seed.txt' to get the initial seed, if this file doesn't exist a
      seed is generated and saved to the file. This should allow for subsequent runs to generate
      the same sequence of random perturbations to the atomic data. We're not saving all of the
      emissivity calculations, but we may want to go back and do this for a specific realization of
      interest.

		+ The magnitude of the perturbation is set by 'ups_sigma.txt' and 'a-values_sigma.txt' where
      uncertainties have been assigned based on recent calculations. Note that the minimum
      uncertainty on the collision rates is 5\% and not 2\%, as in the previous version of the
      file.

* PLOT_EMISS_MC: plots the emissivities

![A sampling of the randomly perturbed emissivities for the lines of interest. The red line is the
 default CHIANTI calculation. ](plot_emiss_mc.jpg)

* Output file

The output file is saved using

```


nrl_save_hdf,  ioneq_file=ioneq_file, $
               logn=logn, logt_max=logt_max, nsim=nsim, $
               emissivity=emissivity, transition=transition, wavelength=wavelength, $
               time_stamp=time_stamp, $
               file=out_file+'_nsim='+trim(nsim)+'.h5'
```								 