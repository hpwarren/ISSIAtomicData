
pro mhd2ints, ps=ps

  file = 'data_0180.idl'
  restore, file

  density = nn
  temperature = te

  nx = n_elements(xco)
  ny = n_elements(yco)
  nz = n_elements(zco)

  ipf = 'mhd2ints.save'
  if not(file_exist(ipf)) then begin

    ints1_sum = fltarr(nx, ny)
    ints2_sum = fltarr(nx, ny)
    nz = n_elements(zco)

    for iz=1, nz/2-1 do begin
      hpw_progress_hash, iz-1, nz/2-1

      logn = alog10(density[*,*,iz])
      logt = alog10(temperature[*,*,iz])

      w1 = 202.044
      o1 = nrl_atodat(w1, logt=logt, logn=logn, emiss=emiss1, /ion_only)
      w2 = 203.826
      o2 = nrl_atodat(w2, logt=logt, logn=logn, emiss=emiss2, /ion_only)

      dz = zco[iz]-zco[iz-1]
      em = density[*,*,iz]*density[*,*,iz]*dz
      ints1 = emiss1*em
      ints2 = emiss2*em

      ints1_sum += ints1
      ints2_sum += ints2

      if iz mod 10 eq 0 then begin
        hpw_setup_xwindow, 700, 700
        pos = hpw_plot_pos(2, 2, /nolabels)
        loadct, 5, /silent

        img = bytscl(temperature[*,*,iz])
        hpw_display_image, [0,1], [0, 1], img, /nox, /noy, /black, pos=pos[0,0,*], /noerase
        ssw_legend, 'Temperature', box=0, textcolor=color
        ssw_legend, 'iz = '+trim(iz), box=0, /right, textcolor=color

        img = bytscl(density[*,*,iz])
        hpw_display_image, [0,1], [0, 1], img, /nox, /noy, /black, pos=pos[1,0,*], /noerase
        ssw_legend, 'Density', box=0, textcolor=color

        loadct, 3, /silent

        img = bytscl(ints1)
        hpw_display_image, [0,1], [0, 1], img, /nox, /noy, /black, pos=pos[0,1,*], /noerase
        ssw_legend, 'Fe XIII 202.044', box=0, textcolor=color

        img = bytscl(ints2)
        hpw_display_image, [0,1], [0, 1], img, /nox, /noy, /black, pos=pos[1,1,*], /noerase
        ssw_legend, 'Fe XIII 203.826', box=0, textcolor=color
      endif

    endfor

    save, ints1_sum, ints2_sum, file=ipf

  endif else restore, ipf

  x0 = 20
  y0 = 16
  ds = 10

  hpw_setup_ps, ps=ps, w=10.0, h=5.0, /land, /color, file='mhd2ints.1'
  hpw_setup_xwindow, 1000, 500, 1
  hpw_thicken_lines
  linecolors
  if keyword_set(ps) then color = 254 else color = 255
  pos = hpw_plot_pos(2, 1, /nolabels)

  loadct, 3, /silent

  img = bytscl(ints1_sum)
  hpw_display_image, [0,nx-1], [0, ny-1], img, /nox, /noy, /black, pos=pos[0,*], /noerase
  ssw_legend, 'Fe XIII 202.044', box=0, textcolor=color
  hpw_drawbox, x0, y0, ds, ds

  img = bytscl(ints2_sum)
  hpw_display_image, [0, nx-1], [0, ny-1], img, /nox, /noy, /black, pos=pos[1,*], /noerase
  ssw_legend, 'Fe XIII 203.826', box=0, textcolor=color
  hpw_drawbox, x0, y0, ds, ds

  hpw_setup_ps, ps=ps, /close

  x = findgen(nx)
  y = findgen(ny)
  x1 = x0 - ds/2.
  x2 = x0 + ds/2.
  y1 = y0 - ds/2.
  y2 = y0 + ds/2.

  x = x[x1:x2]
  y = y[y1:y2]
  temp = temperature[x1:x2, y1:y2, 0:nz/2-1]
  dens = density[x1:x2, y1:y2, 0:nz/2-1]
  i1 = ints1_sum[x1:x2, y1:y2]
  i2 = ints2_sum[x1:x2, y1:y2]

  mx = max(ints2_sum[x1:x2, y1:y2], p)
  px = p mod n_elements(x)
  py = p / n_elements(x)

  t = temperature[px, py, 0:nz/2-1]
  d = density[px, py, 0:nz/2-1]
  z = zco[0:nz/2-1]/1.0E+8

  diff = min(abs(alog10(t) - 6.25), p)

  hpw_setup_ps, ps=ps, w=10.0, h=4.0, /land, /color, file='mhd2ints.2'
  hpw_setup_xwindow, 1000, 600, 2
  hpw_thicken_lines
  linecolors
  !p.position = hpw_plot_pos(1, 1, spacing=[8,4,8,2])

  plot, z, d, /ylog, xstyle=1, ystyle=8, ytitle='Density (cm!a-3!n)', $
        xtitle='Distance along the Loop (Mm)'
  plots, z[p], d[p], psym=2

  axis, yaxis=1, yrange=[1.0E+3, 1.0E+7], ystyle=1, /save, /ylog, color=10, $
        ytitle='Temperature (K)'
  oplot, z, t, color=10
  plots, z[p], t[p], psym=2

  hpw_setup_ps, ps=ps, /close

  peak_i1 = i1[px, py]
  peak_i2 = i2[px, py]
  print, 'peak = ', peak_i2, '/', peak_i1, '=', peak_i2/peak_i1, $
         format='(a10, f8.1, a4, f8.1, a4, f8.2)
  mean_i1 = mean(ints1_sum[x1:x2, y1:y2])
  mean_i2 = mean(ints2_sum[x1:x2, y1:y2])
  print, 'mean = ', mean_i2, '/', mean_i1, '=', mean_i2/mean_i1, $
         format='(a10, f8.1, a4, f8.1, a4, f8.2)


end
