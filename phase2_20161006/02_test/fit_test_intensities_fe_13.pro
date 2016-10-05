
;; #################################################################################################

function fe_13_compute_intensities, logn, emissivity, input_logn, input_logds

  n_lines = (size(emissivity, /dim))[1]
  ints_model = fltarr(n_lines)

  for i=0, n_lines-1 do begin
    e = interpol(emissivity[*, i], logn, input_logn)
    dens = 10.d0^input_logn
    ds = 10.d0^input_logds
    EM = dens*dens*ds
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

pro fit_test_intensities_fe_13, ps=ps

  intensity_file = 'test_intensities_fe_13.h5'

  nrl_restore_hdf, logn=logn, logt_max=logt_max, emissivity=emissivity, wavelength=wavelength, $
                   intensities=intensities, intensities_error=intensities_error, $
                   n_intensities=n_intensities, logn_obs=logn_obs, ds_obs=ds_obs, $
                   time_stamp=time_stamp, $
                   file=intensity_file

  ;; -----------------------------------------------------------------------------------------------
  ;; --- fit with the default CHIANTI data

  p = n_intensities/2
  ints = reform(intensities[p, *])
  err = reform(intensities_error[p, *])
  emissivity0 = reform(emissivity[0, *, *])
  logn_obs = logn_obs[p]
  ds_obs = alog10(ds_obs[p])

  ;; --- perturb the intensities
  n_lines = n_elements(ints)
  ints = ints + err*randomn(seed, n_lines)

  ;; -----------------------------------------------------------------------------------------------
  ;; --- Fit with the perturbed atomic data

  n_chianti = (size(emissivity, /dim))[0]
  res = fltarr(n_chianti, 3)
  for n=0, n_chianti-1 do begin

    this_emissivity = reform(emissivity[n, *, *])

    guess = [9.0, 9.0]
    fa = {ints: ints, err: err, emissivity: this_emissivity, logn: logn}
    fit = mpfit('fe_13_compute_deviates', guess, functargs=fa, /quiet, perror=perr, $
                bestnorm=chi2)

    model = fe_13_compute_intensities( logn, this_emissivity, fit[0], fit[1])

    model = fe_13_compute_intensities( logn, this_emissivity, logn_obs, ds_obs)

    print
    print, '  model log_n = '+trim(fit[0], '(f10.2)') + ' +- '+trim(perr[0], '(f10.3)') + $
           '  ['+trim(logn_obs,'(f10.2)')+']'
    print, ' model log_ds = '+trim(fit[1], '(f10.2)') + ' +- '+trim(perr[1], '(f10.3)') + $
           '  ['+trim(ds_obs,'(f10.2)')+']'           
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

  opf = str_replace(intensity_file, '.h5', '.fits.h5')
  nrl_save_hdf, res=res, file=opf

  hpw_setup_ps, ps=ps, w=10.0, h=6.0, /land, file='fit_test_intensities_fe_13'
  hpw_setup_xwindow, 1000, 600
  hpw_thicken_lines
  !p.multi = [0, 2, 1]

  bs = 0.03
  hist = histogram(res[*,0], binsize=bs, locations=xhist)
  xhist += bs/2.
  plot, xhist, hist, psym=10, $
        xtitle='Inferred Density (log cm!a-3!n)'
  plots, logn_obs, !y.crange, linestyle=2
  plots, median(res[*,0]), !y.crange, thick=3

  ssw_legend, ['Input', 'Median'], linestyle=[0,2], box=0, /right, $
              spacing=1.5, pspacing=1.5

  bs = 0.03
  hist = histogram(res[*,1], binsize=bs, locations=xhist)
  xhist += bs/2.
  plot, xhist, hist, psym=10, $
        xtitle='Inferred Path Length (log cm)'
  plots, ds_obs, !y.crange, linestyle=2
  plots, median(res[*,1]), !y.crange, thick=3

  hpw_setup_ps, ps=ps, /close
  hpw_clean_display

end
