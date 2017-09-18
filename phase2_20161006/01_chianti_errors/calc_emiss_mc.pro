;+
;
; PROJECT:  CHIANTI
;
; NAME:
; CALC_EMISS_MC
;
; PURPOSE:
;
; This rountine calculates the emissivities for an ion using the
; CHIANTI routine EMISS_CALC, then selects only the lines that are
; within a dl  range of the input wavelengths wvl_list. The
; emissivities of the lines are then summed up and saved.
; This is to account for self-blends.
;
; The routine is run a nsim number of times.  It assumes that the
; modification of setup_ion.pro is used. In this modification, the
; (scaled) effective collision strengths and the A-values are
; randomly varied within a set of 'uncertainties' that are defined in
; two input  files.
; The routine also assumes that the modified EMISS_CALC is used.
;
; Calculations are performed for a fixed range of densities at a
; single temperature, the peak temperature of formation. These
; calculations are stored in the EMISSIVITY array.
;
;
; MODIFICATION HISTORY:
;
; routine modified from an original one written by Harry Warren.
;
;       v.1, 3 March 2016, Giulio Del Zanna
;       v.2, 7 March 2016, GDZ, added keyword to calculate random
;               numbers using a normal distribtion.
;       v.3, 10 March 2016, GDZ, fixed a major BUG, the default
;            delta-lambda was 1 Angstrom, instead of 0.1 Angstrom !
;       V.4, 16 Mar 2016, GDZ, added the unperturbed case, stored in
;            the first array
;       V.5, 11 Oct 2016, HPW, added multiplication by
;             abund_fe*this_ioneq*nH_ne/(4*!pi*density)
;
; VERSION     :   5
;
;-

pro calc_emiss_mc, sngl_ion, wvl_list, dl=dl, nsim=nsim, normal=normal, $
                   out_file=out_file

; --- transitions within dl are included in the emissivity,
; accounting for self-blends
  if n_elements(dl) eq 0 then dl=0.1

; --- number of realizations to be run,
  if n_elements(nsim) eq 0 then nsim=100

; --- set the density range
  logn1 =  7.0
  logn2 = 12.0
  dlogn =  0.2
  nlogn = fix((logn2-logn1)/dlogn + 1)
  logn  = logn1 + dlogn*findgen(nlogn)

; --- select the default CHIANTI ionization balance calculation
  chianti_path = !xuvtop
  chianti_path_ioneq = concat_dir(chianti_path,'ioneq')
  ioneq_file = concat_dir(chianti_path_ioneq,'chianti.ioneq')

  convertname, sngl_ion, iz, ion

; ------------------------------------------------------------------------------------------

; setup the output

  n_wvl = n_elements(wvl_list)
  data = {wvl: 0.0, trans: '', em: dblarr(nsim, nlogn)}
  data = replicate(data, n_wvl)

  logt_max = ch_tmax(sngl_ion, /log)

; --- loop over the realizations of the atomic data

  for k=0, nsim-1 do begin

    print, k+1, ' of ', nsim, format='(i4.4, a4, i4.4)'

; --- compute only for logt_max

    ;; --- store the unperturbed as the  first one
    if k eq 0 then begin
      emiss = emiss_calc(iz, ion, temp=logt_max, dens=logn, ioneq_file=ioneq_file, /quiet)
    endif else begin
      emiss = emiss_calc(iz, ion, temp=logt_max, dens=logn, ioneq_file=ioneq_file, /quiet, $
                         /add_uncertainties, normal=normal)
    endelse
    emiss = emiss[where(emiss.flag eq 0)]

    ;; --- loop over the input wavelengths
    for i=0, n_wvl-1 do begin
      ;; --- find nearby transtions
      diff = abs(emiss.lambda-wvl_list[i])
      match = where(diff lt dl, n_nearby)

      ;; --- sum the relevant transitions
      if n_nearby eq 1 then em=reform(emiss[match].em) else if n_nearby gt 1 then $
          em=total(reform(emiss[match[*]].em),2) else if  n_nearby eq 0 then $
              message,'Error ! no lines found !'

      ;; --- keep a record of what was added
      trans = trim(emiss[match].level1)+'-'+trim(emiss[match].level2)
      trans = arr2str(trans, ' / ')

      data[i].trans = trans
      data[i].wvl = wvl_list[i]
      data[i].em[k, *] = em

    endfor
  endfor

  ;; --- pad trans to avoid variable length strings, which get truncated :(
  
  s = data.trans
  len = max(strlen(s))
  s = strpad(s, len, /after, fill=' ')
  data.trans = s

  ;; --- Multiply by the factors needed to get the units correct for the emissivity. In earlier
  ;;     versions of the software this was done downstream but it is best done here.
  
  emissivity = data.em

  read_ioneq, ioneq_file, logt, ioneq, ioneq_ref
  this_ioneq = ioneq[*, iz-1, ion-1]
  this_ioneq = interpol(this_ioneq, logt, logt_max)
  nH_ne = (proton_dens(logt_max, /hydrogen))[0]
  abund_fe = 10.0^(8.10 - 12.0) ;; coronal abundance for Fe

  n_lines = (size(emissivity, /dim))[2]
  n_sim = (size(emissivity, /dim))[0]
  density = 10.0^logn
  chianti_factor = abund_fe*this_ioneq*nH_ne/(4*!pi*density)
  for i=0, n_sim-1 do begin
    for j=0, n_lines-1 do begin
      emissivity[i, *, j] *= chianti_factor
    endfor
  endfor

  data.em = emissivity

  ;; --- save output
  if n_elements(out_file) eq 0 then out_file=sngl_ion+'.monte_carlo'
  if keyword_set(normal) then out_file=out_file+'_normal'

  ;; --- save as an IDL savefle:  
  save, data, logn, nsim, file=out_file+'_nsim='+trim(nsim)+'.save', /compress

  ;; --- save as an HDF file for use with R, no structures
  transition = data.trans
  wavelength = data.wvl
  time_stamp = systime(0)
  
  nrl_save_hdf,  ioneq_file=ioneq_file, $
                 logn=logn, logt_max=logt_max, nsim=nsim, $
                 emissivity=emissivity, transition=transition, wavelength=wavelength, $
                 time_stamp=time_stamp, $
                 file=out_file+'_nsim='+trim(nsim)+'.h5'

end
