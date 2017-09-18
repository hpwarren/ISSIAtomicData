
FUNCTION logticks_exp, axis, index, value
   exponent   = LONG( ALOG10( value ) )
   tickmark = '10!a' + STRTRIM( STRING( exponent ), 2 ) + '!N'
   RETURN, tickmark
END

pro plot_title, pos, title

  x = mean(pos[[0,2]])
  y = pos[3]
  xyouts, x, y+0.015, title, align=0.5, /normal

end

pro plot_emiss_mc, ps=ps

  files = file_search('fe_13.monte_carlo_*.h5', count=n_files)
  case n_files of
    0: return
    1: file = files[0]
    else: begin
      for i=0, n_files-1 do print, i+1, ' ', files[i]
      read, ii, prompt=' Select a file > '
      file = files[ii-1]
    end
  endcase

  restore, file=str_replace(file, '.h5', '.save')

  nrl_restore_hdf, emissivity=emissivity, logn=logn, wavelength=wavelength, $
                   file=file

  sz = size(emissivity, /dim)
  n_sim = sz[0]
  n_logn = sz[1]
  n_lines = sz[2]

  idx_max = 471

  ;; ------------------------------------------------------------------------------------------

  hpw_setup_ps, ps=ps, w=11.0, h=5.0, /land, /color, file='plot_emiss_mc'
  hpw_setup_xwindow, 1100, 700
  hpw_thicken_lines
  pos = hpw_plot_pos(4, 2, spacing=[9,4,6,3])
  linecolors
  erase

  for i=0, n_lines-1 do begin
    ii = i mod 4
    jj = i / 4
    !p.position = pos[ii, jj, *]
    for j=0, n_sim-1 do begin
      y = reform(emissivity[j,*,i])
      if j eq 0 then begin
        com = "plot, logn, y, /ylog, /noerase, /nodata"
        if jj eq 1 then com += ",xtitle='Log N (cm!a-3!n)'"
        if jj ne 1 then com += ",xtickname=replicate(' ',60)"
        res = execute(com)
        plot_title, pos[ii,jj,*], 'Fe XIII '+trim(wavelength[i],'(f10.3)')
      endif
      oplot, logn, y, linestyle=0, color=180
    endfor
    y = reform(emissivity[0,*,i]) 
    oplot, logn, y, thick=3, color=2
    y = reform(emissivity[idx_max,*,i]) 
    oplot, logn, y, thick=3, color=10
    y = reform(emissivity[0,*,i])
    res = execute(com)
  endfor

  dx = pos[0,0,2] - pos[0,0,0]
  xyouts, 0.025, 0.5, 'Emissivity (erg cm!a3!n s!a-1!n sr!a-1!n)', align=0.5, /normal,$
          orientation=90, charsize=1.25

  hpw_setup_ps, ps=ps, /close
  hpw_clean_display

  ;; ------------------------------------------------------------------------------------------

  hpw_setup_ps, ps=ps, w=11.0, h=5.0, /land, /color, file='plot_emiss_mc_rat'
  hpw_setup_xwindow, 1100, 700
  hpw_thicken_lines
  pos = hpw_plot_pos(3, 2, spacing=[9,4,6,3])
  linecolors
  erase

  cntr = 0
  for i=0, n_lines-1 do begin
    if i eq 3 then continue
    ii = i mod 4
    jj = i / 4
    !p.position = pos[ii, jj, *]
    for j=0, n_sim-1 do begin
      y = reform(emissivity[j,*,i])/reform(emissivity[j,*,3])
      if j eq 0 then begin
        com = "plot, logn, y, /ylog, /noerase, /nodata, ytickformat='logticks_exp'"
        if jj eq 1 then com += ",xtitle='Log N (cm!a-3!n)'"
        if jj ne 1 then com += ",xtickname=replicate(' ',60)"
        res = execute(com)
        plot_title, pos[ii,jj,*], 'Fe XIII '+trim(wavelength[i],'(f10.3)')+'/202.044'
      endif
      oplot, logn, y, linestyle=0, color=180
    endfor
    y = reform(emissivity[0,*,i])/reform(emissivity[0,*,3]) 
    oplot, logn, y, thick=3, color=2
    y = reform(emissivity[idx_max,*,i])/reform(emissivity[idx_max,*,3]) 
    oplot, logn, y, thick=3, color=10
    y = reform(emissivity[0,*,i])/reform(emissivity[0,*,3]) 
    res = execute(com)
  endfor

  dx = pos[0,0,2] - pos[0,0,0]
  xyouts, 0.0225, 0.5, 'Ratio (ergs)', align=0.5, /normal,$
          orientation=90, charsize=1.25

  hpw_setup_ps, ps=ps, /close
  hpw_clean_display

end
