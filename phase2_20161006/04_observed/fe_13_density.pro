
pro fe_13_density, ps=ps, check=check

  n_intensities_out = 1000
  ints_204_min = 1000.

  ;; ------------------------------------------------------------------------------------------
  ;; --- get inputs (emissivitis and intensities)

  chianti_file = '../01_chianti_errors/fe_13.monte_carlo_normal_nsim=1000.h5'

  nrl_restore_hdf, logn=logn, logt_max=logt_max, emissivity=emissivity, wavelength=wavelength, $
                   file=chianti_file  

  ;; --- read the line profile fits
  files = ['eis_l1_20130708_002042.fe_13_196_525.2c.fit.genx', $
           'eis_l1_20130708_002042.fe_13_200_021.1c.fit.genx', $
           'eis_l1_20130708_002042.fe_13_201_121.2c.fit.genx', $
           'eis_l1_20130708_002042.fe_13_202_044.1c.fit.genx', $           
           'eis_l1_20130708_002042.fe_13_203_165.3c.fit.genx', $
           'eis_l1_20130708_002042.fe_13_203_826.2c.fit.genx', $
           'eis_l1_20130708_002042.fe_13_209_916.3c.fit.genx']
  nfiles = n_elements(files) 

  ;; --- group into multi-dimensional arrays
  for i=0, n_elements(files)-1 do begin
    restgen, fit, file=files[i]
    if i eq 0 then begin
      intensities = fltarr(nfiles, fit.nx, fit.ny)
      intensities_error = fltarr(nfiles, fit.nx, fit.ny)
      wavelength_obs = fltarr(nfiles)
      line_ids = strarr(nfiles)
    endif
    intensities[i, *, *] = fit.int[*, *, fit.component-1]
    intensities_error[i, *, *] = fit.e_int[*, *, fit.component-1]
    line_ids[i] = fit.line_id[fit.component-1]
    wavelength_obs[i] = float((str2arr(line_ids[i], ' '))[2])
  endfor

  intensities = reform(intensities, [nfiles, fit.nx*fit.ny])
  intensities_error = reform(intensities_error, [nfiles, fit.nx*fit.ny])

  print, 'These should match if the emissivities match the observations'
  print, wavelength, format='(7f10.3)'
  print, wavelength_obs, format='(7f10.3)'

  ;; ------------------------------------------------------------------------------------------
  ;; --- we don't need all of the intensities, randomly select some

  intensities_out = fltarr(n_intensities_out, nfiles)
  intensities_error_out = fltarr(n_intensities_out, nfiles)
  index_out = lonarr(n_intensities_out)

  ;; --- look for all sets of intensities with fe_13_203_826 > ints_min
  diff = min(abs(wavelength_obs - 203.826), pw)
  match = where(intensities[pw, *] ge ints_204_min, nmatch)
  if nmatch lt n_intensities_out then stop
  
  ;; --- randomize the matches
  match = match[ permute(nmatch) ]

  ;; --- make sure the old example set of intensities is on the list
  d = min(abs(intensities[pw, match] - 10620.57), p_index)
  saved = match[0]  
  match[0] = match[p_index]
  match[p_index] = saved

  ;; --- loop over intensities but weed out any bad sets
  cntr = 0
  index = 0
  openw, unit, 'fe_13_density.txt', /get
  printf, unit, '#'
  for i=0, nfiles-1 do printf, unit, '# '+files[i]
  printf, unit, '#'  
  printf, unit, 'cntr', 'index', wavelength, format='(2a10, 7f10.3)'
  while (cntr lt n_intensities_out) do begin
    m = match[index]
    this_set = intensities[*, m]
    bad = where(this_set le 0 or finite(this_set) ne 1, n_bad)
    if n_bad eq 0 then begin
      intensities_out[cntr, *] = this_set
      intensities_error_out[cntr, *] = intensities_error[*, m]
      index_out[cntr] = m
      printf, unit, cntr, m, this_set, format='(i10, i10, 7f10.2)'
      printf, unit, intensities_error[*, m], format='(f30.2, 6f10.2)'
      ++cntr      
    endif
    ++index
  endwhile
  free_lun, unit

  ;; ------------------------------------------------------------------------------------------
  ;; --- a quick check that we can recover the intensities correctly from the original files

  if keyword_set(check) then begin
    print
    print, intensities_out[0, *], format='(7f10.2)'
    print, intensities_error_out[0, *], format='(7f10.2)'
    m = index_out[0]
    i_list = list()
    i_err_list = list()
    i = m mod fit.nx
    j = m / fit.nx
    print, m, i, j
    for n=0, n_elements(files)-1 do begin
      restgen, fit, file=files[n]
      i_list.add, fit.int[i, j, fit.component-1]
      i_err_list.add, fit.e_int[i, j, fit.component-1]
    endfor
    print, i_list, format='(7f10.2)'
    print, i_err_list, format='(7f10.2)'
    pause
  endif

  ;; ------------------------------------------------------------------------------------------
  ;; --- save to output file

  break_file, files[0], disk, dir, name, ext, /last
  opf = (str2arr(name, '.'))[0] + '.fe_density'
  
  time_stamp = systime(0)

  nrl_save_hdf, intensities=intensities_out, intensities_error=intensities_error_out, $
                index=index_out, eis_files=files, eis_nx=fit.nx, eis_ny=fit.ny, $
                n_intensities = n_intensities, n_lines=nfiles, ints_204_min=ints_204_min, $
                logn=logn, emissivity=emissivity, wavelength=wavelength, $               
                time_stamp=time_stamp, $
                file=opf+'.h5'

  ;; ------------------------------------------------------------------------------------------
  ;; --- display densities

  hpw_setup_ps, ps=ps, w=10.0, h=5.0, file=opf+'.a', /land, /color
  hpw_setup_xwindow, 1000, 500, 0
  hpw_thicken_lines
  hpw_setup_dots, size=0.65
  !p.multi = [0, 3, 2]
  !p.charsize=2
  linecolors

  diff = min(abs(wavelength_obs - 202.044), pw)
  idx = indgen(nfiles)
  good = where(idx ne pw, n_good)
  idx = idx[good]

  dens = list()
  title = list()
  for k=0, n_good-1 do begin

    this_title = trim(wavelength[idx[k]],'(f10.3)')+'/'+trim(wavelength[pw],'(f10.3)')
    ratio_theory = reform(emissivity[0, *, idx[k]])/reform(emissivity[0, *, pw])    

    ;; --- plot the densities inferred from the data
    if k ne n_good-1 then begin
      plot, logn, ratio_theory, /xstyle, title=this_title, $
            xtitle='Log Density (cm!a-3!n)', ytitle='Ratio (energy units)'
      logn_obs = fltarr(n_intensities_out)
      for i=0, n_intensities_out-1 do begin
        ratio_obs = intensities_out[i, idx[k]]/intensities_out[i, pw]
        logn_obs[i] = spline(ratio_theory, logn, ratio_obs)
        plots, logn_obs[i], ratio_obs, psym=8
      endfor
      dens.add, logn_obs
    endif else begin
      plot, intensities_out[*, idx[k]]/intensities_out[*, pw], psym=8, /noerase, $
            title=this_title
      plots, !x.crange, mean(ratio_theory), color=2
      print, mean(ratio_theory), mean(intensities_out[*, idx[k]]/intensities_out[*, pw])
    endelse
    title.add, this_title

  endfor

  break_file, files[0], disk, dir, name, ext, /last
  ps_title = (str2arr(name, '.'))[0] + '.fe_density'
  xyouts, 0.5, 0.5, ps_title, align=0.5, charsize=0.75, /normal

  hpw_setup_ps, ps=ps, /close
  hpw_clean_display

  ;; ------------------------------------------------------------------------------------------
  ;; --- display correlations between densities

  hpw_setup_ps, ps=ps, w=10.0, h=5.0, file=opf+'.b', /land, /color
  hpw_thicken_lines  
  hpw_setup_xwindow, 1000, 500, 1
  hpw_setup_dots, size=0.75
  !p.multi = [0, 2, 2]

  range = minmax(logn)

  idx1 = [0, 1, 2, 3]
  idx2 = [4, 4, 4, 4]

  for i=0, 3 do begin
    plot, dens[idx1[i]], dens[idx2[i]], psym=8, xrange=range, yrange=range, $
        xtitle=title[idx1[i]], ytitle=title[idx2[i]]
    plots, !x.crange, !y.crange
    r = correlate(dens[idx1[i]], dens[idx2[i]])
    ssw_legend, 'r = '+trim(r, '(f10.2)'), box=0
  endfor

  break_file, files[0], disk, dir, name, ext, /last
  ps_title = (str2arr(name, '.'))[0] + '.fe_density'
  xyouts, 0.5, 0.5, ps_title, align=0.5, charsize=0.75, /normal

  hpw_setup_ps, ps=ps, /close
  hpw_clean_display
  
end
