
pro chianti_default

  ;; --- read GDZ file

  perturb_file = 'fe_13.monte_carlo_normal.h5'
  nrl_restore_hdf, logn=logn, emissivity=emissivity_perturb, wavelength=wavelength_perturb, $
                   file=perturb_file

  ;; --- compute CHIANTI

  sngl_ion = 'fe_13'

  logt_max = ch_tmax(sngl_ion, /log)

  chianti_path = !xuvtop
  chianti_path_ioneq = concat_dir(chianti_path,'ioneq')
  ioneq_file = concat_dir(chianti_path_ioneq,'chianti.ioneq')

  convertname, sngl_ion, iz, ion

  emiss = emiss_calc(iz, ion, temp=logt_max, dens=logn, ioneq_file=ioneq_file, /quiet)
  emiss = emiss[where(emiss.flag eq 0)]

  ;; --- compare

  wvl_list = 202.044
  diff = min(abs(emiss.lambda-wvl_list), p)
  help, /str, emiss[p]

  plot, logn, reform(emiss[p].em), /ylog, xtitle='Log Density (cm !a-3!n)', $
        ytitle='Emissivity'

  m = where(wavelength_perturb eq wvl_list)
  print, wavelength_perturb[m[0]], wvl_list
  for k=0, 99 do begin
    oplot, logn, reform(emissivity_perturb[k, *, m[0]]), linestyle=1
  endfor
    
end
