
pro fe_13_density, ps=ps, check=check

  n_intensities_out = 1000
  ints_204_min = 1000.

  ;; ------------------------------------------------------------------------------------------
  ;; --- get inputs (emissivitis and intensities)

  ;; --- this file has the errors set to zero and only one set of realizations
  nrl_restore_hdf, logn=logn, emissivity=emissivity, wavelength=wavelength, $                   
                   file='fe_13.monte_carlo.no_error.h5'
  emissivity = reform(emissivity[0, *, *])

  ;; --- read the line profile fits
  files = ['eis_l1_20130708_002042.fe_13_200_021.1c.fit.genx', $
           'eis_l1_20130708_002042.fe_13_201_121.2c.fit.genx', $
           'eis_l1_20130708_002042.fe_13_202_044.1c.fit.genx', $           
           'eis_l1_20130708_002042.fe_13_203_152.2c.fit.genx', $           
           'eis_l1_20130708_002042.fe_13_203_826.2c.fit.genx']
  nfiles = n_elements(files) 

  ;; --- group into multi-dimensional arrays
  for i=0, n_elements(files)-1 do begin
    restgen, fit, file=files[i]
    if i eq 0 then begin
      intensities = fltarr(nfiles, fit.nx, fit.ny)
      intensities_error = fltarr(nfiles, fit.nx, fit.ny)
    endif
    intensities[i, *, *] = fit.int[*, *, fit.component-1]
    intensities_error[i, *, *] = fit.e_int[*, *, fit.component-1]
  endfor

  intensities = reform(intensities, [nfiles, fit.nx*fit.ny])
  intensities_error = reform(intensities_error, [nfiles, fit.nx*fit.ny])

  ;; ------------------------------------------------------------------------------------------
  ;; --- we don't need all of the intensities, randomly select some

  intensities_out = fltarr(n_intensities_out, nfiles)
  intensities_error_out = fltarr(n_intensities_out, nfiles)
  index_out = lonarr(n_intensities_out)

  ;; --- look for all sets of intensities with fe_13_203_826 > ints_min
  match = where(intensities[4, *] ge ints_204_min, nmatch)
  if nmatch lt n_intensities_out then stop
  
  ;; --- randomize the matches
  match = match[ permute(nmatch) ]

  ;; --- loop over intensities but weed out any bad sets
  cntr = 0
  index = 0
  openw, unit, 'fe_density.txt', /get
  printf, unit, '#'
  for i=0, nfiles-1 do printf, unit, '# '+files[i]
  printf, unit, '#'  
  printf, unit, 'cntr', 'index', wavelength, format='(2a10, 5f10.3)'
  while (cntr lt n_intensities_out) do begin
    m = match[index]
    this_set = intensities[*, m]
    bad = where(this_set le 0 or finite(this_set) ne 1, n_bad)
    if n_bad eq 0 then begin
      intensities_out[cntr, *] = this_set
      intensities_error_out[cntr, *] = intensities_error[*, m]
      index_out[cntr] = m
      printf, unit, cntr, m, this_set, format='(i10, i10, 5f10.2)'
      printf, unit, intensities_error[*, m], format='(f30.2, 4f10.2)'
      ++cntr      
    endif
    ++index
  endwhile
  free_lun, unit

  ;; ------------------------------------------------------------------------------------------
  ;; --- a quick check that we can recover the intensities correctly from the original files

  if keyword_set(check) then begin
    print, intensities_out[0, *], format='(5f10.2)'
    print, intensities_error_out[0, *], format='(5f10.2)'
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
    print, i_list, format='(5f10.2)'
    print, i_err_list, format='(5f10.2)'
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

  hpw_setup_ps, ps=ps, w=10.0, h=5.0, file=opf+'.a', /land
  hpw_setup_xwindow, 1000, 500, 0
  hpw_thicken_lines
  hpw_setup_dots, size=0.75
  !p.multi = [0, 2, 2]

  idx = [0, 1, 3, 4]

  dens = list()
  title = list()
  for k=0, 3 do begin

    ;; --- plot CHIANTI temperature ratio
    this_title = trim(wavelength[idx[k]],'(f10.3)')+'/'+trim(wavelength[2],'(f10.3)')
    ratio_theory = reform(emissivity[*, idx[k]])/reform(emissivity[*, 2])
    plot, logn, ratio_theory, /xstyle, title=this_title, $
          xtitle='Log Density (cm!a-3!n)', ytitle='Ratio (energy units)'

    ;; --- plot the densities inferred from the data
    logn_obs = fltarr(n_intensities_out)
    for i=0, n_intensities_out-1 do begin
      ratio_obs = intensities_out[i, idx[k]]/intensities_out[i, 2]
      logn_obs[i] = spline(ratio_theory, logn, ratio_obs)
      plots, logn_obs[i], ratio_obs, psym=8
    endfor
    dens.add, logn_obs
    title.add, this_title

  endfor

  break_file, files[0], disk, dir, name, ext, /last
  ps_title = (str2arr(name, '.'))[0] + '.fe_density'
  xyouts, 0.5, 0.5, ps_title, align=0.5, charsize=0.75, /normal

  hpw_setup_ps, ps=ps, /close
  hpw_clean_display

  ;; ------------------------------------------------------------------------------------------
  ;; --- display correlations between densities

  hpw_setup_ps, ps=ps, w=10.0, h=5.0, file=opf+'.b', /land
  hpw_thicken_lines  
  hpw_setup_xwindow, 1000, 500, 1
  hpw_setup_dots, size=0.75
  !p.multi = [0, 2, 2]

  range = minmax(logn)

  idx1 = [0, 1, 2, 0]
  idx2 = [3, 3, 3, 1]

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
