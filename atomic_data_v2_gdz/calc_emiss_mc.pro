;+
;
; PROJECT:  CHIANTI
;
; NAME:
;	CALC_EMISS_MC
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
;
;
; VERSION     :   4
;
;-


pro calc_emiss_mc, sngl_ion, wvl_list, dl=dl, nsim=nsim, normal=normal, $
                   out_file=out_file
  
  if n_params() eq 0 then begin
    ;; --- a default case
    sngl_ion = 'fe_13'
    wvl_list = [200.021, 201.121, 203.152, 203.826, 202.044]
    nsim = 1
    normal = 1
  endif

; --- transitions within dl are included in the emissivity, 
; accounting for self-blends
if n_elements(dl) eq 0 then dl=0.1 


; --- number of realizations to be run,

if n_elements(nsim) eq 0 then nsim=100


; set the density range 

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

  hpw_progress_hash, k, nsim-1
   
; --- compute only for logt_max
   
; store the unperturbed as the  first one    
   if k eq 0 then emiss = emiss_calc(iz, ion, temp=logt_max, dens=logn, ioneq_file=ioneq_file, /quiet) $
                          else $   
   emiss = emiss_calc(iz, ion, temp=logt_max, dens=logn, ioneq_file=ioneq_file, /quiet,/add_uncertainties, normal=normal)
   emiss = emiss[where(emiss.flag eq 0)]
   
   ; --- loop over the input  wavelengths

   for i=0, n_wvl-1 do begin
      ; --- find nearby transtions 
      diff = abs(emiss.lambda-wvl_list[i])
      match = where(diff lt dl, n_nearby)

      ; --- sum the relevant transitions
      
      if n_nearby eq 1 then em=reform(emiss[match].em) else if n_nearby gt 1 then $
         em=total(reform(emiss[match[*]].em),2) else if  n_nearby eq 0 then $
            message,'Error ! no lines found !'
      
      
      ; --- keep a record of what was added
      trans = trim(emiss[match].level1)+'-'+trim(emiss[match].level2)
      trans = arr2str(trans, ' / ')
      

      data[i].trans = trans
      data[i].wvl = wvl_list[i]
      data[i].em[k, *] = em

   endfor
endfor 

; --- pad trans to avoid variable length strings, which get truncated :(
s = data.trans
len = max(strlen(s))
s = strpad(s, len, /after, fill=' ')
data.trans = s  


; save as an IDL savefle: 

if n_elements(out_file) eq 0 then out_file=sngl_ion+'.monte_carlo'

if keyword_set(normal) then out_file=out_file+'_normal'

save, data,logn, file=out_file+'_nsim='+trim(nsim)+'.save',/compress


 emissivity   = data.em
  transition   = data.trans
  wavelength   = data.wvl
  
  time_stamp = systime(0)
  
 nrl_save_hdf,  ioneq_file=ioneq_file, $
                logn=logn, logt_max=logt_max, nsim=nsim, $
                emissivity=emissivity, transition=transition, wavelength=wavelength, $                
                time_stamp=time_stamp, $
                file=out_file+'_nsim='+trim(nsim, '(i4.4)')+'.h5'
 

end 
