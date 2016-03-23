
;; #################################################################################################
;;
;; Given a density and path length use the atomic data to compute the expected intensities for
;; these Fe XIII lines
;;

function fe_13_compute_intensities, logn, emissivity, input_logn, input_logds

  n_lines = (size(emissivity, /dim))[1]
  ints_model = fltarr(n_lines)

  for i=0, n_lines-1 do begin
    e = interpol(emissivity[*, i], logn, input_logn)
    EM = 10.0^(2.*input_logn + input_logds)
    ints_model[i] = e*EM
  endfor

  return, ints_model
end

;; #################################################################################################
;;
;; Compute the deviation of the model and observed intensities. For use with MPFIT.

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
;; The abundance, ionization fraction, and other factors have been left off of the emissivity
;; calculation. This makes the resulting emission measures have very funny units. Let's put
;; them into the emissivity so that we can compare with other calculations.
;;

pro compute_chianti_factor, logn, emissivity, chianti_factor, logt_max

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

end

;; #################################################################################################
;;
;; A routine to analyze observed EIS Fe XIII intensities with perturbed CHIANTI atomic data
;;

pro fe_13_fit_intensities

  ;; Read the observed EIS intensities. This is a file of 1000 randomly selected sets of Fe XIII
  ;; line intensities. Note that emissivity0 is the unperturbed atomic data for these lines as a
  ;; function of density. Computed with fe_13_density.pro in data/Hinode.

  file = 'eis_l1_20130708_002042.fe_density.h5'

  nrl_restore_hdf, intensities=intensities, intensities_error=intensities_error, $
                   index=index, eis_files=eis_files, eis_nx=eis_nx, eis_ny=eis_ny, $
                   logn=logn, emissivity=emissivity0, wavelength=wavelength0, $
                   n_lines=n_lines, ints_204_min=ints_204_min, $
                   file=file

  compute_chianti_factor, logn, emissivity0, chianti_factor, logt_max

  nrl_save_hdf, logn=logn, logt_max=logt_max, chianti_factor=chianti_factor, $
                file='chianti_factor.fe13.h5'

  ;; -----------------------------------------------------------------------------------------------
  ;; fit with the default CHIANTI data

  mx = max(intensities[*, 2], p)
  print, p
  ints = reform(intensities[p, *])
  err = reform(intensities_error[p, *])

  guess = [9.0, 9.0]
  fa = {ints: ints, err: err, emissivity: emissivity0, logn: logn}
  fit = mpfit('fe_13_compute_deviates', guess, functargs=fa, /quiet, perror=perr, $
              bestnorm=chi2)

  model = fe_13_compute_intensities( logn, emissivity0, fit[0], fit[1])
  log_em_model = (2*fit[0] + fit[1])

  print, '     logt_max = '+trim(logt_max, '(f10.2)')
  print, '  model log_n = '+trim(fit[0], '(f10.2)') + ' +- '+trim(perr[0], '(f10.3)')
  print, ' model log_ds = '+trim(fit[1], '(f10.2)') + ' +- '+trim(perr[1], '(f10.3)')
  print, ' model log_em = '+trim(log_em_model, '(f10.2)')
  print, '         chi2 = '+trim(chi2, '(f10.1)')
  print, 'Line', 'Imodel', 'Iobs', 'SigmaI', 'dI/Sigma', 'dI/I', format='(6a10)'
  for i=0, n_elements(model)-1 do begin
    var1 = abs(model[i]-ints[i])/err[i]
    var2 = 100*abs(model[i]-ints[i])/ints[i]
    print, wavelength0[i], model[i], ints[i], err[i], var1, var2, format='(f10.3, 4f10.2, f10.1)'
  endfor

  res0 = [fit[0], fit[1], chi2]

  pause

  ;; -----------------------------------------------------------------------------------------------
  ;;  now read GDZ perturbed CHIANTI atomic data and fit intensities. Assume that the density scale
  ;;  is the same

  perturb_file = 'fe_13.monte_carlo_normal_nsim=1000.h5'
  nrl_restore_hdf, emissivity=emissivity_perturb, wavelength=wavelength, file=perturb_file

  nwave = n_elements(wavelength) 
  wave_index = lonarr(nwave)
  for i=0, nwave-1 do begin
    m = where(wavelength eq wavelength0[i])
    wave_index[i] = m[0]
  endfor
  
 
  n_perturb = (size(emissivity_perturb, /dim))[0]
  res = fltarr(n_perturb, 3)
  for n=0, n_perturb-1 do begin

    this_emissivity_gdz = reform(emissivity_perturb[n, *, *])
    this_emissivity = this_emissivity_gdz
    for k=0, n_lines-1 do begin
      this_emissivity[*, k] = this_emissivity_gdz[*, wave_index[k]]*chianti_factor
      ;; print, wavelength[wave_index[k]], wavelength0[k]
    endfor

    guess = [9.0, 9.0]
    fa = {ints: ints, err: err, emissivity: this_emissivity, logn: logn}
    fit = mpfit('fe_13_compute_deviates', guess, functargs=fa, /quiet, perror=perr, $
                bestnorm=chi2)

    model = fe_13_compute_intensities( logn, this_emissivity, fit[0], fit[1])
    log_em_model = (2*fit[0] + fit[1])
    chi2_est = total(((model - ints)/err)^2)

    print
    print, ' Run # '+trim(n+1)
    print, '     logt_max = '+trim(logt_max, '(f10.2)')
    print, '  model log_n = '+trim(fit[0], '(f10.2)') + ' +- '+trim(perr[0], '(f10.3)')
    print, ' model log_ds = '+trim(fit[1], '(f10.2)') + ' +- '+trim(perr[1], '(f10.3)')
    print, ' model log_em = '+trim(log_em_model, '(f10.2)')
    print, '         chi2 = '+trim(chi2, '(f10.1)')
    print, 'Line', 'Imodel', 'Iobs', 'SigmaI', 'dI/Sigma', 'dI/I', format='(6a10)'
    for i=0, n_elements(model)-1 do begin
      var1 = abs(model[i]-ints[i])/err[i]
      var2 = 100*abs(model[i]-ints[i])/ints[i]
      print, wavelength0[i], model[i], ints[i], err[i], var1, var2, format='(f10.3, 4f10.2, f10.1)'
    endfor

    res[n, *] = [fit[0], fit[1], chi2]

    if n le 5 then pause

  endfor

  opf = 'fe_13_fit_intensities.normal.h5'  
  nrl_save_hdf, res0=res0, res=res, file=opf

end
