
;;
;; This rountine generates realizations of randomly perturbed atomic data and saves the
;; emissivieis for select lines to an output file. See the CHIANTI routines setup_ion.pro and
;; read_wgfa2.pro for details.
;;
;; See chianti_ratio_mc.input for how to specify the ion and the lines of interest.
;;
;; Calculations are performed for a range of densities at a single temperature, the peak
;; temperature of formation. These calculations are stored in the EMISSIVITY array.
;;
;; Optionally, calculations for a range of temperatures can also be computed if the calc_t keyword
;; is set. These calculations are stored in a different array, EMISSIVITY_T. If the keyword isn't
;; used the array exists but is filled with zeros.
;;
;; Note that transitions within dw_max are included in the emissivity, accounting for self-blends
;;

pro chianti_ratio_mc, calc_t=calc_t

  ;; --- select an input from chianti_ratio_mc.input
  src = rd_tfile('chianti_ratio_mc.input', 1, /nocomment)
  nsrc = n_elements(src)
  for i=0, nsrc-1 do print, i+1, ' ', src[i]
  read, ii, prompt=' select an input > '
  if ii le 0 then return

  ;; --- parse the input
  line = src[ii-1]
  sngl_ion = strtrim((str2arr(line, '|'))[0], 2)
  wvl_list = (str2arr(line, '|'))[1]
  wvl_list = float(str2arr(wvl_list, ','))
  wvl_list = wvl_list[sort(wvl_list)]

  ;; ------------------------------------------------------------------------------------------
  ;; --- setup the calculation

  logt_max = ch_tmax(sngl_ion, /log)  

  logt1 = 5.5
  logt2 = 7.5
  dlogt = 0.1
  nlogt = fix((logt2 - logt1)/dlogt + 1)
  logt  = logt1 + dlogt*findgen(nlogt)

  logn1 =  7.0
  logn2 = 12.0
  dlogn =  0.2
  nlogn = fix((logn2-logn1)/dlogn + 1)
  logn  = logn1 + dlogn*findgen(nlogn)  

  ;; --- these variables are read by setup_ion.pro and read_wgfa2.pro
  PERTURB = 0.10
  setenv, 'CHIANTI_PERTURB_SPL1='+trim(PERTURB)
  setenv, 'CHIANTI_PERTURB_SPL2='+trim(PERTURB)
  setenv, 'CHIANTI_PERTURB_AVAL='+trim(PERTURB)

  ;; --- number of realizations to be run, 100 is a lot for Fe XII and Fe XIII!
  nsim = 1000

  ;; --- transitions within dw_max are included in the emissivity, accounting for blends
  dw_max = 0.1

  ;; --- select the default CHIANTI ionization balance calculation
  chianti_path = !xuvtop
  chianti_path_ioneq = concat_dir(chianti_path,'ioneq')
  ioneq_file = concat_dir(chianti_path_ioneq,'chianti.ioneq')

  convertname, sngl_ion, iz, ion

  ;; ------------------------------------------------------------------------------------------
  ;; --- loop over the realizations of the atomic data
  
  n_wvl = n_elements(wvl_list)
  data = {wvl: 0.0, trans: '', em: dblarr(nsim, nlogn), em_t: dblarr(nsim, nlogt, nlogn)}
  data = replicate(data, n_wvl)

  for k=0, nsim-1 do begin

    print, k+1, ' of ', nsim, format='(i3, a4, i4.4)'

    tic
    ;; --- compute for a range of temperatures
    if keyword_set(calc_t) then begin
      emiss_t = emiss_calc(iz, ion, temp=logt, dens=logn, ioneq_file=ioneq_file, /quiet)
      emiss_t = emiss_t[where(emiss_t.flag eq 0)]
    endif

    ;; --- compute only for logt_max
    emiss = emiss_calc(iz, ion, temp=logt_max, dens=logn, ioneq_file=ioneq_file, /quiet)
    emiss = emiss[where(emiss.flag eq 0)]
    toc

    ;; --- loop over the wavelengths
    if k eq 0 then text = list()
    for i=0, n_wvl-1 do begin
      ;; --- find nearby transtions 
      diff = abs(emiss.lambda-wvl_list[i])
      match = where(diff lt dw_max, n_nearby)
      if k eq 0 then begin
        text.add, ' ion = '+sngl_ion
        text.add, ' this wave = '+trim(wvl_list[i], '(f10.3)')
        text.add, ' nearby transiions = '+trim(n_nearby)
        for mm=0, n_nearby-1 do begin
          m = match[mm]
          if mm eq 0 then text.add, $
              string( 'index', 'wave', 'dwave', 'emiss', format='(a7, a10, a10, a12)')
          text.add, string(m, emiss[m].lambda, diff[m], max(emiss[m].em), $
                 format='(i7, f10.4, f10.4, e12.2)')
        endfor
      endif

      ;; --- sum the relevant transitions
      em = fltarr(nlogn, n_nearby)
      em_t = fltarr(nlogt, nlogn, n_nearby)
      for n=0, n_nearby-1 do begin
        em[*, n] = reform(emiss[match[n]].em)
        if keyword_set(calc_t) then em_t[*, *, n] = reform(emiss_t[match[n]].em)
      endfor
      if n_nearby gt 1 then begin
        em = total(em, 2)
        em_t = total(em_t, 3)
      endif else begin
        em = reform(em)
        em_t = reform(em_t)
      endelse

      ;; --- keep a record of what was added
      trans = trim(emiss[match].level1)+'-'+trim(emiss[match].level2)
      trans = arr2str(trans, ' / ')
      if k eq 0 then begin
        text.add, ' transitions = '+trans
        text.add, ''
      endif
      
      data[i].trans = trans
      data[i].wvl = wvl_list[i]
      data[i].em[k, *] = em
      data[i].em_t[k, *, *] = em_t
    endfor

  endfor

  ;; --- pad trans to avoid variable length strings, which get truncated :(
  s = data.trans
  len = max(strlen(s))
  s = strpad(s, len, /after, fill=' ')
  data.trans = s

  text = text.ToArray()
  foreach t, text do print, t

  emissivity_t = data.em_t  
  emissivity   = data.em
  transition   = data.trans
  wavelength   = data.wvl

  ;; ------------------------------------------------------------------------------------------
  ;; --- save results

  CHIANTI_PERTURB_SPL1 = float(getenv('CHIANTI_PERTURB_SPL1'))
  CHIANTI_PERTURB_SPL2 = float(getenv('CHIANTI_PERTURB_SPL2'))
  CHIANTI_PERTURB_AVAL = float(getenv('CHIANTI_PERTURB_AVAL'))

  opf = sngl_ion+'.monte_carlo.'+trim(fix(100*PERTURB), '(i2.2)')+'.h5'
  time_stamp = systime(0)

  nrl_save_hdf, CHIANTI_PERTURB_SPL1=CHIANTI_PERTURB_SPL1, $
                CHIANTI_PERTURB_SPL2=CHIANTI_PERTURB_SPL2, $
                CHIANTI_PERTURB_AVAL=CHIANTI_PERTURB_AVAL, $
                ioneq_file=ioneq_file, $
                logt=logt, logn=logn, logt_max=logt_max, nsim=nsim, $
                emissivity=emissivity, transition=transition, wavelength=wavelength, $
                emissivity_t=emissivity_t, text=text, $
                time_stamp=time_stamp, $
                file=opf

end
