
pro fe_13_fit_intensities_spec, ps=ps

  file = 'eis_l1_20130708_002042.003.055.spec.txt'
  src  = rd_tfile(file, /auto)
  wave = float(reform(src[0, *]))
  spec = float(reform(src[3, *]))
  err  = float(reform(src[4, *]))

  sngl_ion = 'fe_13'
  log_t = ch_tmax(sngl_ion, /log)
  log_n = 9.48
  log_ds = 9.14
  log_em = 2*log_n + log_ds
  ioneq_path = concat_dir(!xuvtop, 'ioneq')
  ioneq_name = concat_dir(ioneq_path, 'chianti.ioneq')
  abund_path = concat_dir(!xuvtop, 'abundance') 
  abund_name = concat_dir(abund_path, 'sun_coronal_1992_feldman.abund')

  wmin  = 190.0
  wmax  = 210.0
  dwave = (wave[1]-wave[0])/4.

  file = 'fe_13_fit_intensities_spec.idl'
  if not(file_exist(file)) then begin
    ch_synthetic, wmin, wmax, output=output,$
                 density=10.0^log_n, $
                 logt_isothermal=log_t, $
                 logem_isothermal=log_em,$
                 ioneq_name=ioneq_name
    
    make_chianti_spec, output, lambda, o,$
                       bin_size=dwave,$
                       abund_name=abund_name,$
                       continuum=1

    save, output, o, file=file
  endif
  restore, file
  
  ch_wave = o.lambda
  ch_spec = o.spectrum*0.85
  fwhm = 4*3*dwave
  ch_spec = eve_gaussfold(ch_wave, ch_spec, fwhm)

  hpw_setup_ps, ps=ps, w=10.0, h=5.0, /land, /color, file='fe_13_fit_intensities_spec'
  hpw_setup_xwindow, 2000, 500
  hpw_thicken_lines
  linecolors
  !p.position = hpw_plot_pos(1, 1, spacing=[11, 4, 2, 2])
  @hpw_setup_symbols

  plot, wave, spec, psym=10, xrange=[199, 205], xstyle=1, $
        ytitle='Intensity (erg cm!a-2!n s!a-1!n sr!a-1!n '+angstrom+'!a-1!n)', $
        xtitle='Wavelength ('+angstrom+')'
  oplot_err, wave, spec, yerr=err, psym=3

  oplot, ch_wave, ch_spec, color=2, thick=3

  w = [200.021, 201.121, 203.152, 203.826, 202.044]  
  for i=0, n_elements(w)-1 do begin
    plots, w[i]+0.2*[-1,1], !y.crange[1], thick=6, color=2
    new = convert_coord(w[i], !y.crange[1], /data, /to_normal)
    xyouts, new[0], new[1]*1.01, trim(w[i], '(f10.3)'), align=0.5, /normal
  endfor

  hpw_setup_ps, ps=ps, /close
  hpw_clean_display

end
