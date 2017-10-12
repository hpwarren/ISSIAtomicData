
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

pro fe_13_fit_intensities_all

  file = 'eis_l1_20130708_002042.fe_density.h5'

  nrl_restore_hdf, intensities=intensities, intensities_error=intensities_error, $
                   index=index, eis_files=eis_files, eis_nx=eis_nx, eis_ny=eis_ny, $
                   logn=logn, emissivity=emissivity, wavelength=wavelength, $
                   n_lines=n_lines, ints_204_min=ints_204_min, $
                   file=file

  ;; -----------------------------------------------------------------------------------------------
  ;; --- Fit with the perturbed atomic data

  n_intensities = (size(intensities))[1]

  fit_ne = fltarr(n_intensities)
  fit_ne_err = fltarr(n_intensities)
  fit_ds = fltarr(n_intensities)
  fit_ds_err = fltarr(n_intensities)
  fit_chi2 = fltarr(n_intensities)
  ints_model = fltarr(n_intensities, n_lines)
  text = list()
  text.add, '# least-squares fits to the observed intensities'
  text.add, '# using the default CHIANTI atomic data'
  text.add, '# '
  for n=0, n_intensities-1 do begin

    ints = intensities[n, *]
    err = intensities_error[n, *]

    this_emissivity = reform(emissivity[0, *, *]) ;; use default

    guess = [9.5, 9.0]
    fa = {ints: ints, err: err, emissivity: this_emissivity, logn: logn}
    fit = mpfit('fe_13_compute_deviates', guess, functargs=fa, /quiet, perror=fit_perr, $
                bestnorm=chi2, dof=dof)

    model = fe_13_compute_intensities( logn, this_emissivity, fit[0], fit[1])

    text.add, ''
    text.add, '         n pixel = ' + trim(n)
    text.add, '     model log_n = '+trim(fit[0], '(f10.2)') + ' +- '+trim(fit_perr[0], '(f10.3)')
    text.add, '    model log_ds = '+trim(fit[1], '(f10.2)') + ' +- '+trim(fit_perr[1], '(f10.3)')
    text.add, '            chi2 = '+trim(chi2, '(f10.1)')
    text.add, ' normalized chi2 = '+trim(chi2/dof, '(f10.1)')    
    text.add, string('Line', 'Iobs', 'SigmaI', 'Imodel', 'dI/I', 'dI/Sigma', format='(6a10)')
    for i=0, n_elements(model)-1 do begin
      var1 = abs(model[i]-ints[i])/err[i]
      var2 = 100*abs(model[i]-ints[i])/ints[i]
      text.add, string(wavelength[i],  ints[i], err[i], model[i], var2, var1, $
                       format='(f10.3, 4f10.1, f10.1)')
    endfor

    fit_ne[n] = fit[0]
    fit_ne_err[n] = fit_perr[0]
    fit_ds[n] = fit[1]
    fit_ds_err[n] = fit_perr[1]
    fit_chi2[n] = chi2
    ints_model[n, *] = model
    
  endfor

  opf = str_replace(file, '.h5', '.default_fits.txt')
  openw, unit, opf, /get
  for i=0, n_elements(text)-1 do printf, unit, text[i]
  free_lun, unit
  print
  print, 'saved to ' + opf  
  
  opf = str_replace(file, '.h5', '.default_fits.h5')
  print
  print, 'saved to ' + opf
  nrl_save_hdf, fit_ne=fit_ne, fit_ds=fit_ds, fit_chi2=fit_chi2, $
                fit_ne_err=fit_ne_err, fit_ds_err=fit_ds_err, ints_model=ints_model, $
                file=opf

end
