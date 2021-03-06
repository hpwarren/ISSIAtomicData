
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

;;
;; Fit a *single* set of intensities with all of the realizations of the atomic data, plot the
;; distribution of n and ds. Note that the single set of intensities is perturbed by a set of
;; normally distributed random deviations.
;;

pro fit_test_intensities_fe_13, ps=ps, noerror=noerror

  intensity_file = 'test_intensities_fe_13.h5'

  nrl_restore_hdf, logn=logn, logt_max=logt_max, emissivity=emissivity, wavelength=wavelength, $
                   intensities=intensities, intensities_error=intensities_error, $
                   n_intensities=n_intensities, logn_obs=logn_obs, ds_obs=ds_obs, $
                   time_stamp=time_stamp, $
                   file=intensity_file

  ;; -----------------------------------------------------------------------------------------------
  ;; --- select a set of intensities

  diff = min(abs(logn_obs - 9.9), p)
    
  ints = reform(intensities[p, *])
  err = reform(intensities_error[p, *])
  logn_obs = logn_obs[p]
  ds_obs = alog10(ds_obs[p])

  ;; --- perturb the intensities
  ints_saved = ints
  if not(keyword_set(noerror)) then begin
    seed = 42
    n_lines = n_elements(ints)
    ints = ints + err*randomn(seed, n_lines)
  endif

  ;; -----------------------------------------------------------------------------------------------
  ;; --- Fit with the perturbed atomic data

  n_chianti = (size(emissivity, /dim))[0]
  res = fltarr(n_chianti, 3)
  for n=0, n_chianti-1 do begin

    this_emissivity = reform(emissivity[n, *, *])

    guess = [9.5, 9.5]
    fa = {ints: ints, err: err, emissivity: this_emissivity, logn: logn}
    fit = mpfit('fe_13_compute_deviates', guess, functargs=fa, /quiet, perror=perr, $
                bestnorm=chi2, dof=dof)

    model = fe_13_compute_intensities( logn, this_emissivity, fit[0], fit[1])

    print
    print, '  model log_n = '+trim(fit[0], '(f10.2)') + ' +- '+trim(perr[0], '(f10.3)') + $
           '  ['+trim(logn_obs,'(f10.2)')+']'
    print, ' model log_ds = '+trim(fit[1], '(f10.2)') + ' +- '+trim(perr[1], '(f10.3)') + $
           '  ['+trim(ds_obs,'(f10.2)')+']'           
    print, '         chi2 = '+trim(chi2, '(f10.1)')
    print, ' normalized chi2 = '+trim(chi2/dof, '(f10.1)')        
    print, 'Line', 'Iobs', 'SigmaI', 'Imodel', 'dI/I', 'dI/Sigma', format='(6a10)'    
    for i=0, n_elements(model)-1 do begin
      var1 = abs(model[i]-ints[i])/err[i]
      var2 = 100*abs(model[i]-ints[i])/ints[i]
      ;; print, wavelength[i], model[i], ints[i], err[i], var1, var2, format='(f10.3, 4f10.2, f10.1)'
      print, wavelength[i],  ints[i], err[i], model[i], var2, var1, format='(f10.3, 4f10.1, f10.1)'
    endfor

    res[n, *] = [fit[0], fit[1], chi2]
    if n eq 0 then pause
  endfor

  opf = str_replace(intensity_file, '.h5', '.fits.h5')
  nrl_save_hdf, res=res, file=opf

  ;; ------------------------------------------------------------------------------------------

  hpw_setup_ps, ps=ps, w=6.0, h=7.0, file='fit_test_intensities_fe_13'
  hpw_setup_xwindow, 500, 1000
  hpw_thicken_lines
  @hpw_setup_symbols
  !p.multi = [0, 1, 2]

  bs = 0.03
  hist = histogram(res[*,0], binsize=bs, locations=xhist)
  xhist += bs/2.
  plot, xhist, hist, psym=10, $
        xtitle='Inferred Density (log cm!a-3!n)'

  ssw_legend, ['mean = ' + trim(mean(res[*,0]), '(f10.2)'), $
               s_sigma + ' = ' +trim(stddev(res[*,0]), '(f10.2)')], box=0, spacing=1.5

  bs = 0.03
  hist = histogram(res[*,1], binsize=bs, locations=xhist)
  xhist += bs/2.
  plot, xhist, hist, psym=10, $
        xtitle='Inferred Path Length (log cm)'

  ssw_legend, ['mean = ' + trim(mean(res[*,1]), '(f10.2)'), $
               s_sigma + ' = ' +trim(stddev(res[*,1]), '(f10.2)')], box=0, spacing=1.5  

  hpw_setup_ps, ps=ps, /close
  hpw_clean_display

  print, res[*, 2]

end
