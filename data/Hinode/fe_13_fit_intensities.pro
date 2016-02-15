
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

  n_lines = (size(emissivity, /dim))[1]
  density = 10.0^logn
  chianti_factor = abund_fe*this_ioneq*nH_ne/(4*!pi*density)  
  for i=0, n_lines-1 do begin
    emissivity[*, i] *= chianti_factor
  endfor
  
  ;; -----------------------------------------------------------------------------------------------
  ;; --- fit with the default CHIANTI data

  mx = max(intensities[*, 2], p)
  ints = reform(intensities[p, *])
  err = reform(intensities_error[p, *])

  guess = [9.0, 9.0]
  fa = {ints: ints, err: err, emissivity: emissivity, logn: logn}
  fit = mpfit('fe_13_compute_deviates', guess, functargs=fa, /quiet, perror=perr, $
              bestnorm=chi2)

  model = fe_13_compute_intensities( logn, emissivity, fit[0], fit[1])

  print, '  model log_n = '+trim(fit[0], '(f10.2)') + ' +- '+trim(perr[0], '(f10.3)')
  print, ' model log_ds = '+trim(fit[1], '(f10.2)') + ' +- '+trim(perr[1], '(f10.3)')
  print, '         chi2 = '+trim(chi2, '(f10.1)')  
  print, 'Line', 'Imodel', 'Iobs', 'SigmaI', 'dI/Sigma', 'dI/I', format='(6a10)'
  for i=0, n_elements(model)-1 do begin
    var1 = abs(model[i]-ints[i])/err[i]
    var2 = 100*abs(model[i]-ints[i])/ints[i]
    print, wavelength[i], model[i], ints[i], err[i], var1, var2, format='(f10.3, 4f10.2, f10.1)'
  endfor

  res0 = [fit[0], fit[1], chi2]

  pause
  
  ;; -----------------------------------------------------------------------------------------------
  ;; --- Fit with the perturbed atomic data

  all_perturb_files = file_search('../../fe_13.monte_carlo.*.h5', count=nfiles)
  for i=0, nfiles-1 do begin
    print, i+1, '  ', all_perturb_files[i]
  endfor
  read, ii, prompt=' select a file > '
  perturb_file = all_perturb_files[ii-1]
  nrl_restore_hdf, emissivity=emissivity_perturb, $
                   CHIANTI_PERTURB_SPL1=CHIANTI_PERTURB_SPL1, $
                   file=perturb_file
  opf = 'fe_13_fit_intensities.'+trim(fix(CHIANTI_PERTURB_SPL1*100),'(i2.2)')+'.h5'

  n_perturb = (size(emissivity_perturb, /dim))[0]
  res = fltarr(n_perturb, 3)
  for n=0, n_perturb-1 do begin

    this_emissivity = reform(emissivity_perturb[n, *, *])
    for k=0, n_lines-1 do begin
      this_emissivity[*, k] = this_emissivity[*, k]*chianti_factor
    endfor

    guess = [9.0, 9.0]
    fa = {ints: ints, err: err, emissivity: this_emissivity, logn: logn}
    fit = mpfit('fe_13_compute_deviates', guess, functargs=fa, /quiet, perror=perr, $
                bestnorm=chi2)                

    model = fe_13_compute_intensities( logn, this_emissivity, fit[0], fit[1])
    chi2_est = total(((model - ints)/err)^2)
    
    print, '  model log_n = '+trim(fit[0], '(f10.2)') + ' +- '+trim(perr[0], '(f10.3)')
    print, ' model log_ds = '+trim(fit[1], '(f10.2)') + ' +- '+trim(perr[1], '(f10.3)')
    print, '         chi2 = '+trim(chi2, '(f10.1)')
    print, '         chi2 = '+trim(chi2_est, '(f10.1)')    
    print, 'Line', 'Imodel', 'Iobs', 'SigmaI', 'dI/Sigma', 'dI/I', format='(6a10)'
    for i=0, n_elements(model)-1 do begin
      var1 = abs(model[i]-ints[i])/err[i]
      var2 = 100*abs(model[i]-ints[i])/ints[i]
      print, wavelength[i], model[i], ints[i], err[i], var1, var2, format='(f10.3, 4f10.2, f10.1)'
    endfor

    res[n, *] = [fit[0], fit[1], chi2]
    
  endfor

  nrl_save_hdf, res0=res0, res=res, file=opf

end
