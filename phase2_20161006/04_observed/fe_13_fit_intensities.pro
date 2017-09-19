
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
  ;; --- Fit with the perturbed atomic data

  diff = min(abs(intensities[*,0] - 1473.1), p)

  n_pixel = 216
  ints = intensities[n_pixel, *]
  err = intensities_error[n_pixel, *]  

  n_chianti = (size(emissivity, /dim))[0]
  res = fltarr(n_chianti, 3)
  perr = fltarr(n_chianti, 2)
  for n=0, n_chianti-1 do begin

    this_emissivity = reform(emissivity[n, *, *])

    guess = [9.5, 9.0]
    fa = {ints: ints, err: err, emissivity: this_emissivity, logn: logn}
    fit = mpfit('fe_13_compute_deviates', guess, functargs=fa, /quiet, perror=fit_perr, $
                bestnorm=chi2, dof=dof)

    model = fe_13_compute_intensities( logn, this_emissivity, fit[0], fit[1])

    print
    print, '       n chianti = ' + trim(n)
    print, '         n pixel = ' + trim(n_pixel)
    print, '     model log_n = '+trim(fit[0], '(f10.2)') + ' +- '+trim(fit_perr[0], '(f10.3)')
    print, '    model log_ds = '+trim(fit[1], '(f10.2)') + ' +- '+trim(fit_perr[1], '(f10.3)')
    print, '            chi2 = '+trim(chi2, '(f10.1)')
    print, ' normalized chi2 = '+trim(chi2/dof, '(f10.1)')    
    print, 'Line', 'Iobs', 'SigmaI', 'Imodel', 'dI/I', 'dI/Sigma', format='(6a10)'
    for i=0, n_elements(model)-1 do begin
      var1 = abs(model[i]-ints[i])/err[i]
      var2 = 100*abs(model[i]-ints[i])/ints[i]
      print, wavelength[i],  ints[i], err[i], model[i], var2, var1, format='(f10.3, 4f10.1, f10.1)'
    endfor

    res[n, *] = [fit[0], fit[1], chi2]
    perr[n, *] = [fit_perr[0], fit_perr[1]]

    if n eq 0 then pause
    if n eq 471 then pause
    if n eq 216 then pause
    
  endfor

  med_log_n = median(res[*,0])
  std_log_n = stddev(res[*,0])
  med_log_ds = median(res[*,1])
  std_log_ds = stddev(res[*,1])  

  print, res[0,0], perr[0,0], format='(2f10.3)'
  print, med_log_n, std_log_n, format='(2f10.3)'
  print
  print, res[0,1], perr[0,1], format='(2f10.3)'
  print, med_log_ds, std_log_ds, format='(2f10.3)'

  ;; ------------------------------------------------------------------------------------------

  hpw_setup_ps, ps=ps, w=10.0, h=6.0, /land, file='fe_13_fit_intensities'
  hpw_setup_xwindow, 1000, 600
  hpw_thicken_lines
  !p.multi = [0, 2, 1]

  bs = 0.05
  hist = histogram(res[*, 0], binsize=bs, locations=xhist)
  xhist += bs/2.
  plot, xhist, hist, psym=10, $
        xtitle='Inferred Density (log cm!a-3!n)'

  ssw_legend, ['Input', 'Median'], linestyle=[0,2], box=0, /right, $
              spacing=1.5, pspacing=1.5

  bs = 0.05
  hist = histogram(res[*,1], binsize=bs, locations=xhist)
  xhist += bs/2.
  plot, xhist, hist, psym=10, $
        xtitle='Inferred Path Length (log cm)'

  hpw_setup_ps, ps=ps, /close
  hpw_clean_display

end
