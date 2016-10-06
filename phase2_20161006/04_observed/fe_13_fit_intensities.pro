
;; #################################################################################################

function fe_13_compute_intensities, logn, emissivity, input_logn, input_logds

  n_lines = (size(emissivity, /dim))[1]
  ints_model = fltarr(n_lines)

  for i=0, n_lines-1 do begin
    e = interpol(emissivity[*, i], logn, input_logn)
    EM = 10.0^(2.*input_logn)*10.0^input_logds
    ints_model[i] = e*EM
  endfor

  return, ints_model
end

;; #################################################################################################

function fe_13_compute_deviates, p, $
                                 ints=ints, $
                                 err=err, $
                                 logn=logn, $
                                 emissivity=emissivity

  model_logn = p[0]
  model_logds = p[1]

  model_ints = fe_13_compute_intensities(logn, emissivity, model_logn, model_logds)

  deviates = (model_ints - ints)/err

  return, deviates
end

;; #################################################################################################

pro fe_13_fit_intensities, ps=ps

  file = 'eis_l1_20130708_002042.fe_density.h5'

  nrl_restore_hdf, intensities=intensities, intensities_error=intensities_error, $
                   index=index, eis_files=eis_files, eis_nx=eis_nx, eis_ny=eis_ny, $
                   logn=logn, emissivity=emissivity, wavelength=wavelength, $
                   n_lines=n_lines, ints_204_min=ints_204_min, $
                   file=file

  ;; -----------------------------------------------------------------------------------------------
  ;; --- multiply emissivity by ioneq and other factors

  sngl_ion = 'fe_13'
  chianti_path = !xuvtop
  chianti_path_ioneq = concat_dir(chianti_path,'ioneq')
  ioneq_file = concat_dir(chianti_path_ioneq,'chianti.ioneq')
  convertname, sngl_ion, iz, ion
  read_ioneq, ioneq_file, logt, ioneq, ioneq_ref
  this_ioneq = ioneq[*,iz-1,ion-1]
  logt_max = ch_tmax(sngl_ion, /log)
  diff = min(abs(logt_max - logt), p)
  this_ioneq = interpol(this_ioneq, logt, logt_max)
  nH_ne = (proton_dens(logt_max, /hydrogen))[0]
  abund_fe = 10.0^(8.10 - 12.0)

  n_lines = (size(emissivity, /dim))[2]
  n_sim = (size(emissivity, /dim))[0]
  density = 10.0^logn
  chianti_factor = abund_fe*this_ioneq*nH_ne/(4*!pi*density)
  for i=0, n_sim-1 do begin
    for j=0, n_lines-1 do begin
      emissivity[i, *, j] *= chianti_factor
    endfor
  endfor

  ;; -----------------------------------------------------------------------------------------------
  ;; --- Fit with the perturbed atomic data

  n_chianti = (size(emissivity, /dim))[0]
  res = fltarr(n_chianti, 3)
  for n=0, n_chianti-1 do begin

    this_emissivity = reform(emissivity[n, *, *])
    ints = intensities[n, *]
    err = intensities_error[n, *]

    guess = [9.0, 9.0]
    fa = {ints: ints, err: err, emissivity: this_emissivity, logn: logn}
    fit = mpfit('fe_13_compute_deviates', guess, functargs=fa, /quiet, perror=perr, $
                bestnorm=chi2)

    model = fe_13_compute_intensities( logn, this_emissivity, fit[0], fit[1])

    print
    print, '  model log_n = '+trim(fit[0], '(f10.2)') + ' +- '+trim(perr[0], '(f10.3)')
    print, ' model log_ds = '+trim(fit[1], '(f10.2)') + ' +- '+trim(perr[1], '(f10.3)')
    print, '         chi2 = '+trim(chi2, '(f10.1)')
    print, 'Line', 'Imodel', 'Iobs', 'SigmaI', 'dI/Sigma', 'dI/I', format='(6a10)'
    for i=0, n_elements(model)-1 do begin
      var1 = abs(model[i]-ints[i])/err[i]
      var2 = 100*abs(model[i]-ints[i])/ints[i]
      print, wavelength[i], model[i], ints[i], err[i], var1, var2, format='(f10.3, 4f10.2, f10.1)'
    endfor

    res[n, *] = [fit[0], fit[1], chi2]
    if n eq 0 then pause
  endfor

  ;; ------------------------------------------------------------------------------------------

  hpw_setup_ps, ps=ps, w=10.0, h=6.0, /land, file='fe_13_fit_intensities'
  hpw_setup_xwindow, 1000, 600
  hpw_thicken_lines
  !p.multi = [0, 2, 1]

  bs = 0.05
  hist = histogram(res[*,0], binsize=bs, locations=xhist)
  xhist += bs/2.
  plot, xhist, hist, psym=10, $
        xtitle='Inferred Density (log cm!a-3!n)'
  plots, median(res[*,0]), !y.crange, thick=3

  ssw_legend, ['Input', 'Median'], linestyle=[0,2], box=0, /right, $
              spacing=1.5, pspacing=1.5

  bs = 0.05
  hist = histogram(res[*,1], binsize=bs, locations=xhist)
  xhist += bs/2.
  plot, xhist, hist, psym=10, $
        xtitle='Inferred Path Length (log cm)'
  plots, median(res[*,1]), !y.crange, thick=3

  hpw_setup_ps, ps=ps, /close
  hpw_clean_display

end
