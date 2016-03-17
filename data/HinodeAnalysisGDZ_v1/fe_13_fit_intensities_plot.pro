
pro fe_13_fit_intensities_plot, ps=ps

  files = ['fe_13_fit_intensities.normal.h5']
  results = list()
  for i=0, n_elements(files)-1 do begin
    nrl_restore_hdf, res=res, res0=res0, file=files[i]
    print, res0
    results.add, res
  endfor

  hpw_setup_ps, ps=ps, w=9.0, h=9.0, /color, file='fe_13_fit_intensities_plot'
  hpw_setup_xwindow, 1000, 1000
  hpw_setup_dots, size=0.6
  hpw_thicken_lines
  !p.multi = [0, 2, 2, 0, 0]
  linecolors

  xr = [9.0, 10.0]
  yr = [8.0, 10.0]
  
  for i=0, n_elements(files)-1 do begin
    res = results[i]
    if i eq 0 then begin
      plot, res[*, 0], res[*, 1], psym=8, xrange=xr, yrange=yr, $
            xtitle='Log Electron Density', ytitle='Log ds', /nodata
    endif
    color = i + 2
    oplot, res[*, 0], res[*, 1], psym=8, color=color
  endfor

  plots, res0[0], res0[1], psym=8, color=7, symsize=2.5

  for i=0, n_elements(files)-1 do begin
    res = results[i]
    if i eq 0 then begin
      plot, res[*, 2], res[*, 0], psym=8, xrange=[0, 500], yrange=[9, 10], $
            ytitle='Log Electron Density', xtitle='Chi2', /nodata
    endif
    color = i + 2
    oplot, res[*, 2], res[*, 0], psym=8, color=color
  endfor

  plots, res0[2], res0[0], psym=8, color=7, symsize=2.5


  r = results.ToArray(dim=1)
  c = r[*, 2]
  m = where(c lt 250)
  
  bs = 0.05
  hist = histogram(r[m, 0], binsize=bs, locations=xhist)
  xhist += bs/2.
  plot, xhist, hist, psym=10, xrange=xr, xtitle='Log Electron Density', $
        ytitle='Number of Samples'
  print, median(r[m, 0]), res0[0]
  print, stddev(r[m, 0])

  bs = 0.05
  hist = histogram(r[m, 1], binsize=bs, locations=xhist)
  xhist += bs/2.
  plot, xhist, hist, psym=10, xrange=yr, xtitle='Log ds', $
        ytitle='Number of Samples'        
  print, median(r[m, 1]), res0[1]
  print, stddev(r[m, 1])

  hpw_setup_ps, ps=ps, /close
  hpw_clean_display

end
