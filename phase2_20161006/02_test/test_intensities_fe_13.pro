
;;
;; Here we generate some fake intensities from an assumed density and a theoretical relationship
;; between the density and the path length. 
;;

pro test_intensities_fe_13

  ;; --- read in the atomic data

  chianti_file = '../01_chianti_errors/fe_13.monte_carlo_normal_nsim=1000.h5'

  nrl_restore_hdf, logn=logn, logt_max=logt_max, emissivity=emissivity, wavelength=wavelength, $
                   file=chianti_file

  ;; --- select some random densities and estimate the intensities for all of the lines

  n_intensities = 1000
  seed = 42L
  logn_obs = 8.5 + 2.0*randomu(seed, n_intensities)
  p_obs = 2*!boltzmann*10.0^(logn_obs + logt_max)

  exposure = 60.0 ;; EIS exposure time in seconds, for eis_throughput
  slit = 2 ;; EIS slit width in arcsecs, for eis_throughput

  n_lines = n_elements(wavelength) 
  intensities = fltarr(n_intensities, n_lines)
  intensities_error = fltarr(n_intensities, n_lines)
  ds_obs = fltarr(n_intensities)
  for i=0, n_intensities-1 do begin
    ds = (2560.*1.0D+5/p_obs[i]) ;; from Martens et al. 2000, Eq 24.
    dens = 10.d0^logn_obs[i]
    em = dens*dens*ds
    err = fltarr(n_lines)
    ints = fltarr(n_lines)
    for j=0, n_lines-1 do begin
      emiss = interpol(emissivity[0, *, j], logn, logn_obs[i]) ;; 0 is the default CHIANTI
      ints[j] = emiss*em 
      eis = eis_throughput(wavelength[j], ints_erg=ints[j], exposure=exposure, slit=slit)
      err[j] = eis.ints_erg_err
    endfor
    intensities[i, *] = ints
    intensities_error[i, *] = err
    ds_obs[i] = ds
    
    print, logn_obs[i], logt_max, p_obs[i], ds, format='(3f10.2,e12.2)'
    print, ints, format='(7f10.2)'
    print, err, format='(7f10.2)'
    print, 100*err/ints, format='(7f10.2)'
    print
  endfor

  intensity_file = 'test_intensities_fe_13.h5'
  time_stamp = systime(0)

  nrl_save_hdf, logn=logn, logt_max=logt_max, emissivity=emissivity, wavelength=wavelength, $
                intensities=intensities, intensities_error=intensities_error, $
                n_intensities=n_intensities, $
                time_stamp=time_stamp, logn_obs=logn_obs, ds_obs=ds_obs, $
                file=intensity_file

end
