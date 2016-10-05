
;;
;;
;;
;;

pro test_intensities_fe_13

  ;; --- read in the atomic data

  chianti_file = '../01_chianti_errors/fe_13.monte_carlo_normal_nsim=1000.h5'

  nrl_restore_hdf, logn=logn, logt_max=logt_max, emissivity=emissivity, wavelength=wavelength, $
                   file=chianti_file

  ;; --- multiply by the factors needed to get the units correct for the emissivity

  sngl_ion = 'fe_13'
  chianti_path = !xuvtop
  chianti_path_ioneq = concat_dir(chianti_path,'ioneq')
  ioneq_file = concat_dir(chianti_path_ioneq,'chianti.ioneq')
  convertname, sngl_ion, iz, ion
  read_ioneq, ioneq_file, logt, ioneq, ioneq_ref
  this_ioneq = ioneq[*, iz-1, ion-1]
  logt_max = ch_tmax(sngl_ion, /log)
  diff = min(abs(logt_max - logt), p)
  this_ioneq = interpol(this_ioneq, logt, logt_max)
  nH_ne = (proton_dens(logt_max, /hydrogen))[0]
  abund_fe = 10.0^(8.10 - 12.0) ;; coronal abundance for Fe

  n_lines = (size(emissivity, /dim))[2]
  n_sim = (size(emissivity, /dim))[0]
  density = 10.0^logn
  chianti_factor = abund_fe*this_ioneq*nH_ne/(4*!pi*density)
  for i=0, n_sim-1 do begin
    for j=0, n_lines-1 do begin
      emissivity[i, *, j] *= chianti_factor
    endfor
  endfor

  ;; --- select some random densities and estimate the intensities for all of the lines

  n_intensities = 1000
  logn_obs = 8.5 + 2.5*randomu(seed, n_intensities)
  logn_obs = logn_obs[sort(logn_obs)]
  p_obs = 2*!boltzmann*10.0^(logn_obs + logt_max)

  exposure = 60.0
  slit = 2

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
      emiss = interpol(emissivity[0, *, j], logn, logn_obs[i])
      ints[j] = emiss*em ;; 0 is the default CHIANTI calculation
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
