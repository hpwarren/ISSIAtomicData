
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

  density = 10.0^logn
  chianti_factor = abund_fe*this_ioneq*nH_ne/(4*!pi*density)
  emissivity *= chianti_factor

end

pro chianti_default_hpw, ps=ps

  ;; --- read GDZ file

  perturb_file = 'fe_13.monte_carlo_normal_nsim=1000.h5'
  nrl_restore_hdf, logn=logn, emissivity=emissivity_perturb, wavelength=wavelength_perturb, $
                   file=perturb_file

  ;; --- compute default versions of CHIANTI emissivities

  sngl_ion = 'fe_13'

  logt_max = ch_tmax(sngl_ion, /log)

  chianti_path = !xuvtop
  chianti_path_ioneq = concat_dir(chianti_path,'ioneq')
  ioneq_file = concat_dir(chianti_path_ioneq,'chianti.ioneq')

  convertname, sngl_ion, iz, ion

  emiss = emiss_calc(iz, ion, temp=logt_max, dens=logn, ioneq_file=ioneq_file, /quiet)
  emiss = emiss[where(emiss.flag eq 0)]

  ;; --- plot default and perturbed atomic data

  hpw_setup_ps, ps=ps, w=11.0, h=6.0, /color, file='chianti_default_hpw'
  hpw_setup_xwindow, 1200, 600, 1  
  hpw_thicken_lines
  linecolors
  pos = plot_pos(3, 2, /xtitle, /ytitle, spacing=[8, 5, 2, 2])
  erase

  for i=0, n_elements(wavelength_perturb)-1 do begin

    wvl_list = wavelength_perturb[i]
    dlambda = 0.1
    match = where(abs(emiss.lambda-wvl_list) lt dlambda, n_nearby)
    
    ;; --- sum the relevant transitions
    nlogn = n_elements(logn) 
    em = fltarr(nlogn, n_nearby)
    for n=0, n_nearby-1 do begin
      em[*, n] = reform(emiss[match[n]].em)
    endfor
    if n_nearby gt 1 then begin
      em = total(em, 2)
    endif else begin
      em = reform(em)
    endelse

    ii = i mod 3 & jj = i / 3
    !p.position = pos[ii, jj, *]
    title = trim(wvl_list, '(f10.3)')
    plot, logn, em, /ylog, xtitle='Log Density (cm !a-3!n)', $
          ytitle='Emissivity', /noerase, xrange=[8, 10], title=title
    
    m = where(wavelength_perturb eq wvl_list)
    dim = size(emissivity_perturb, /dim)
    n_perturb = dim[0]
    for k=0, n_perturb-1 do begin
      oplot, logn, reform(emissivity_perturb[k, *, m[0]]), linestyle=1
    endfor
    
    oplot, logn, em, color=2

  endfor
    
  hpw_setup_ps, ps=ps, /close
    
end
