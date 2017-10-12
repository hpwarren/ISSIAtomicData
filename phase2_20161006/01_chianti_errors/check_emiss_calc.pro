
pro check_emiss_calc

  sngl_ion = 'fe_13'

  logn1 =  7.0
  logn2 = 12.0
  dlogn =  0.2
  nlogn = fix((logn2-logn1)/dlogn + 1)
  logn  = logn1 + dlogn*findgen(nlogn)

  chianti_path = !xuvtop
  chianti_path_ioneq = concat_dir(chianti_path,'ioneq')
  ioneq_file = concat_dir(chianti_path_ioneq,'chianti.ioneq')

  convertname, sngl_ion, iz, ion

  logt_max = ch_tmax(sngl_ion, /log)

  emiss = emiss_calc(iz, ion, temp=logt_max, dens=logn, ioneq_file=ioneq_file, /quiet)

  wvl_list = [202.044, 209.619, 209.916]
  dl = 0.1
  for i=0, n_elements(wvl_list)-1 do begin
    diff = abs(emiss.lambda-wvl_list[i])
    match = where(diff lt dl, n_nearby)
    mx = max(emiss[match].em)
    for j=0, n_nearby-1 do begin
      jj = match[j]
      if max(emiss[jj].em) gt 1.0E-2*mx then begin
        print,emiss[jj].lambda, emiss[jj].level1, emiss[jj].level2, max(emiss[jj].em)
      endif
    endfor
    print
  endfor

end
