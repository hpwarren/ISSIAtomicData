    title = ' med log_n = ' + trim(med_log_n, '(f10.2)') + ' $\sigma$ = ' + $
          trim(std_log_n, '(f10.3)')

      
  p = plot(xhist_log_n, yfit_log_n, thick=2, linestyle=0, /overplot, /current)

  title =' med log_n = ' + trim(med_log_ds, '(f10.2)') + ' $\sigma$ = ' + $
         trim(std_log_ds, '(f10.3)')
  p = plot(xhist_log_ds, hist_log_ds, layout=[2,1,2], /stairstep, /current, $
           xtitle='Log ds', title=title)
  p = plot(xhist_log_ds, yfit_log_ds, thick=2, linestyle=0, /overplot, /current)

  label = text(10, 10, 'Pixel = ' + trim(n_pixel), /device, target=p)
