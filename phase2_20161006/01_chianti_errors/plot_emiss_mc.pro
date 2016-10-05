
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

  nrl_restore_hdf, emissivity=emissivity, logn=logn, wavelength=wavelength, $
                   file=file

  sz = size(emissivity, /dim)
  n_sim = sz[0]
  n_logn = sz[1]
  n_lines = sz[2]

  hpw_setup_ps, ps=ps, w=11.0, h=7.0, /land, /color, file='plot_emiss_mc'
  hpw_setup_xwindow, 1100, 700
  hpw_thicken_lines
  pos = hpw_plot_pos(4, 2, /ytitles, /xtitles, spacing=[7,4,2,2])
  linecolors
  erase

  for i=0, n_lines-1 do begin
    !p.position = pos[i mod 4, i / 4, *]
    for j=0, n_sim-1, 10 do begin
      y = reform(emissivity[j,*,i])/10.0^logn
      if j eq 0 then begin
        plot, logn, y, /ylog, /noerase, yrange=[1.0E-21, 1.0E-18], $
              title=trim(wavelength[i],'(f10.3)'), xtitle='Log T (K)', /nodata
      endif
      oplot, logn, y, linestyle=1, color=128
    endfor
    y = reform(emissivity[0,*,i])/10.0^logn    
    oplot, logn, y, thick=3, color=2
  endfor

  hpw_setup_ps, ps=ps, /close
  hpw_clean_display

end
