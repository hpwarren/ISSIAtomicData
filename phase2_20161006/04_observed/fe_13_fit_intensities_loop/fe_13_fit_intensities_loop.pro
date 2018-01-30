
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

pro fe_13_fit_set, n_pixel, intensities, intensities_error, wavelength, logn, emissivity, $
                   res, perr, text

  ints = intensities[n_pixel-1, *]
  err = intensities_error[n_pixel-1, *]

  n_chianti = (size(emissivity, /dim))[0]
  res = fltarr(n_chianti, 3)
  perr = fltarr(n_chianti, 2)
  text = list()
  for n=0, n_chianti-1 do begin

    this_emissivity = reform(emissivity[n, *, *])

    guess = [9.5, 9.0]
    fa = {ints: ints, err: err, emissivity: this_emissivity, logn: logn}
    fit = mpfit('fe_13_compute_deviates', guess, functargs=fa, /quiet, perror=fit_perr, $
                bestnorm=chi2, dof=dof)

    model = fe_13_compute_intensities( logn, this_emissivity, fit[0], fit[1])

    text.add, ''
    text.add, '       n chianti = ' + trim(n)
    text.add, '         n pixel = ' + trim(n_pixel)
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

    res[n, *] = [fit[0], fit[1], chi2]
    perr[n, *] = [fit_perr[0], fit_perr[1]]

  endfor

end

;; #################################################################################################

pro fe_13_fit_intensities_loop

  file = '../eis_l1_20130708_002042.fe_density.h5'

  nrl_restore_hdf, intensities=intensities, intensities_error=intensities_error, $
                   index=index, eis_files=eis_files, eis_nx=eis_nx, eis_ny=eis_ny, $
                   logn=logn, emissivity=emissivity, wavelength=wavelength, $
                   n_lines=n_lines, ints_204_min=ints_204_min, $
                   file=file

  break_file, file, disk, dir, opf_base, ext, /last

  n_pixel_max = 1000
  for n_pixel=1, n_pixel_max do begin

    ;; --- compute solution
    fe_13_fit_set, n_pixel, intensities, intensities_error, wavelength, $
                   logn, emissivity, res, perr, output

    ;; --- compute statistics
    med_log_n = median(res[*,0])
    std_log_n = stddev(res[*,0])
    med_log_ds = median(res[*,1])
    std_log_ds = stddev(res[*,1])

    bs = 0.0125
    hist_log_n = histogram(res[*, 0], binsize=bs, locations=xhist_log_n)
    xhist_log_n += bs/2.
    yfit_log_n = mpfitpeak(xhist_log_n, hist_log_n, fit_log_n, /positive, nterms=3)

    bs = 0.025
    hist_log_ds = histogram(res[*,1], binsize=bs, locations=xhist_log_ds)
    xhist_log_ds += bs/2.
    yfit_log_ds = mpfitpeak(xhist_log_ds, hist_log_ds, fit_log_ds, /positive, nterms=3)

    ;; --- output file names

    opd = concat_dir('data', trim(n_pixel, '(i4.4)'))
    hpw_output_dir, opd

    opf_hist = opf_base + '.hist.'+trim(n_pixel, '(i4.4)')+'.h5'
    opf_text = opf_base + '.hist.'+trim(n_pixel, '(i4.4)')+'.txt'
    opf_plot = opf_base + '.hist.'+trim(n_pixel, '(i4.4)')+'.jpg'

    opf_hist = concat_dir(opd, opf_hist)
    opf_text = concat_dir(opd, opf_text)
    opf_plot = concat_dir(opd, opf_plot)

    ;; --- save the histograms
    nrl_save_hdf, xhist_log_n=xhist_log_n, hist_log_n=hist_log_n, $
                  xhist_log_ds=xhist_log_ds, hist_log_ds=hist_log_ds, $
                  med_log_n=med_log_n, std_log_n=std_log_n, $
                  med_log_ds=med_log_ds, std_log_ds=std_log_ds, $
                  fit_log_n=fit_log_n, fit_log_ds=fit_log_ds, $
                  fits=res, file=opf_hist

    ;; --- save information on the fits
    openw, unit, opf_text, /get
    foreach t, output do printf, unit, t
    free_lun, unit

    ;; --- plot the histograms
    if n_pixel eq 1 then w = window(dimension=[800, 500], /buffer) else w.erase

    title = ' med log_n = ' + trim(med_log_n, '(f10.2)') + ' $\sigma$ = ' + $
            trim(std_log_n, '(f10.3)')
    p = plot(xhist_log_n, hist_log_n, layout=[2,1,1], /stairstep, /current, $
             xtitle='Log n', title=title)
    p = plot(xhist_log_n, yfit_log_n, thick=2, linestyle=0, /overplot, /current)

    title =' med log_n = ' + trim(med_log_ds, '(f10.2)') + ' $\sigma$ = ' + $
           trim(std_log_ds, '(f10.3)')
    p = plot(xhist_log_ds, hist_log_ds, layout=[2,1,2], /stairstep, /current, $
             xtitle='Log ds', title=title)
    p = plot(xhist_log_ds, yfit_log_ds, thick=2, linestyle=0, /overplot, /current)

    label = text(10, 10, 'Pixel = ' + trim(n_pixel), /device, target=p)

    p.save, opf_plot, resolution=96

    hpw_progress_hash, n_pixel-1, n_pixel_max-1, message='Working on ' + trim(n_pixel, '(i4.4)')

  endfor

end
