;+
;
; PROJECT:  CHIANTI
;
;       CHIANTI is an Atomic Database Package for Spectroscopic Diagnostics of
;       Astrophysical Plasmas. It is a collaborative project involving the Naval
;       Research Laboratory (USA), the University of Florence (Italy), the
;       University of Cambridge and the Rutherford Appleton Laboratory (UK). 
;
; NAME:
;	SETUP_ION
;
; PURPOSE:
;
;	Several of the CHIANTI routines require the atomic data for ions 
;       to be read into common blocks, and so this routine performs this 
;       task for individual ions.
;
;       Note that the inputs WVLMIN and WVLMAX merely serve to check 
;       whether the ion has lines lying in this wavelength range. If yes, 
;       then the common blocks are filled with the data; if no, then the 
;       common blocks are left empty. 
;
;       An important point to note is that the 
;       wgfa and upsilon common blocks contain their data in 2-D array 
;       form, suitable for solving the level balance equations. When 
;       dealing with dielectronic files, there will be some transitions 
;       which will have two A-values, one representing autoionisation 
;       rather than radiative decay. For the level balance equations, 
;       these two A-values need to be summed together. However, for 
;       computing the line emissivity, one needs to only use the radiative 
;       decay A-value. For this purpose, the routine outputs the 1-D 
;       arrays read by READ_WGFA2 which contain the two separate A-values, 
;       allowing other routines to pick out just the radiative decay 
;       A-value. These outputs are listed below (LVL1, LVL2, etc.).
;
;
; CALLING SEQUENCE:
;
;       SETUP_ION,Ion,Wvlmin,Wvlmax,Wvltst,Lvl1,Lvl2,Wvl1,Gf1,A_value1
;
;
; INPUTS
;
;	ION     String specifying ion, e.g., 'c_2'
;       
;       WVLMIN  Minimum wavelength of interest (angstroms). If there are 
;               no lines with wavelengths between WVLMIN and WVLMAX, then 
;               no data is read into the common blocks. If both WVLMIN 
;               and WVLMAX are set to -1, then no checks are made on the 
;               wavelengths (i.e., the common blocks will be filled).
;
;       WVLMAX  Maximum wavelength of interest (angstroms). See above.
;	
; OPTIONAL INPUT
;
;       PATH    Can be used to specify a different directory to !xuvtop 
;               for the location of the data files. E.g., 
;               path='/data/user/chianti/o_4'
;
; KEYWORDS
;
;       NOPROT  Switches off the inclusion of proton rates.
;
;       ALL     By default the routine only looks for wavelengths which 
;               are positive when checking if they lie in the wavelength 
;               range. By setting /ALL then negative wavelengths are also 
;               included.
;
; OUTPUTS:
;
;	WVLTST  The number of transitions lying between WVLMIN and WVLMAX. 
;               If WVLMIN and WVLMAX are both -1, then WVLTST will 
;               be the total number of transitions in the .wgfa file 
;               (including any with zero wavelength).
;
;       LVL1    Number of final level of transition (1-D array)
;
;       LVL2    Number of initial level of transition (1-D array)
;
;       WVL1    Wavelenths (in angstroms) of spectral lines formed by this 
;               ion (1-D array)
;
;       GF1     Oscillator strength (1-D array)
;
;       A_VALUE1 Radiative transition probability (1-D array)
;
; OPTIONAL OUTPUTS
;
;       ANYLINES Contains the indices of the elements of lvl1, lvl2, wvl1, 
;                gf1, a_value1 arrays for which the wavelengths lie in 
;                specified wavelength range.
;
; EXAMPLE:
;
;    IDL> setup_ion,'c_2',1000.,1500.,wvltst,Lvl1,Lvl2,Wvl1,Gf1,A_value1
;    IDL> print,wvltst
;          17
;    IDL> setup_ion,'c_2',-1,-1,wvltst,Lvl1,Lvl2,Wvl1,Gf1,A_value1
;    IDL> print,wvltst
;          43
;    IDL> setup_ion,'c_2',-1,0,wvltst,Lvl1,Lvl2,Wvl1,Gf1,A_value1
;    IDL> print,wvltst
;           0
;
; COMMON BLOCKS:
;
;	common elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
;       common wgfa, wvl,gf,a_value
;       common upsilon,t_type,deu,c_ups,splups
;
; CALLS
;
;    READ_WGFA2, READ_SPLUPS, READ_ELVLC, ION2FILENAME
;
;
; MODIFICATION HISTORY:
;
;       Ver.1, Sep-99, Ken Dere
;       Ver.2, 10-Oct-00, Peter Young
;                Allowed wvlmin and wvlmax to take -1 values.
;                Removed elements common block.
;                Removed for-loop.
;                Added PATH.
;       Ver.3, 12-Nov-01, Peter Young
;                Added proton common block and /noprot keyword
;       Ver.4, 11-Dec-01, Peter Young
;                Added keyword /all.
;
;       V.5, Giulio Del Zanna (GDZ)
;                   generalized directory concatenation to work for
;                   Unix, Windows  and VMS. 
;
;       V.6, 25-Feb-2004, Enrico Landi (EL)
;                   Added call to read_ionrec.pro to read ionization 
;                   and recombination files.
;
;       V.7, 8-Jan-2014, Peter Young
;                   Modified to read the new .scups files instead of
;                   the old .splups files.
;
;       v.8, 1-May-2014 Giulio Del Zanna
;            fixed the reading of the new .scups files.
;
;       v.9, 27-Jun-2014, Peter Young
;            Modified how splstr is used (retain the data, info tags).
;
;       v.10, 3 March 2016, Giulio Del Zanna (GDZ)
;             added the option to read uncertainty files and perform a 
;             random variation of the rates, the A-values and the
;             collision strengths. We modify the scaled Upsilons since
;             they are linearly dependent on the Upsilons.
;             Note: the IDL function randomu is used to create
;             uniformly-distributed values within +/- sigma, where the
;             sigma values are read by two input files.
;
;             The program searches two files in the working directory,
;             a-values_sigma.txt   ups_sigma.txt 
;             with the lower,upper indices and the uncertainty (%).
;             If no values are found, a default_uncertainty value is
;             used (by default it is set to 30%). 
;
;       v.11, 7 March 2016, GDZ, added keyword normal
;       v.12, 16 March 2016, GDZ, fixed a major bug.
;
;
; VERSION     :   12
;
;-


PRO setup_ion,ions,wvlmin,wvlmax,wvltst,lvl1,lvl2,wvl1,gf1,a_value1, $
              path=path, noprot=noprot, all=all, anylines=anylines, $
              add_uncertainties=add_uncertainties,$
              default_uncertainty=default_uncertainty, normal=normal


COMMON elvlc,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
COMMON wgfa, wvl,gf,a_value
COMMON upsilon,splstr
COMMON proton, pstr, pe_ratio
COMMON ionrec,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status


IF n_params() LT 1 THEN BEGIN
   print,'  IDL>  setup_ion,ions,wvlmin,wvlmax,wvltst,lvl1,lvl2,wvl,gf,a_value'
   print,'               [ path= , noprot= ]'
   return
ENDIF


IF N_ELEMENTS(path) NE 0 THEN fname = concat_dir(path, ions) $
ELSE ion2filename,ions,fname

if n_elements(default_uncertainty) eq 0 then $
   default_uncertainty=30


wname=fname+'.wgfa'
elvlcname=fname+'.elvlc'
;upsname=fname+'.splups'
upsname=fname+'.scups'
pname=fname+'.psplups'

;
; Read .ci and .rec files
;

read_ionrec,fname,rec_rate,ci_rate,temp_ionrec,luprec,lupci,status

;
; Read .wgfa file
;
read_wgfa2,wname,lvl1,lvl2,wvl1,gf1,a_value1,wgfaref


IF (wvlmin EQ -1) AND (wvlmax EQ -1) THEN BEGIN
   wvltst = n_elements(lvl1)
   anylines=findgen(wvltst)
ENDIF ELSE BEGIN
   IF keyword_set(all) THEN BEGIN
      anylines=where((abs(wvl1) GE wvlmin) AND (abs(wvl1) LE wvlmax) AND $
                     (a_value1 NE 0.),wvltst)
   ENDIF ELSE BEGIN
      anylines=where((wvl1 GE wvlmin) AND (wvl1 LE wvlmax) AND $
                     (a_value1 NE 0.),wvltst)
   ENDELSE
ENDELSE


IF wvltst GT 0 THEN BEGIN
   
;------- GDZ addition here ------------------------------------------
   
  if keyword_set(add_uncertainties) then begin

    ;; --- 
    common get_mc_seed, init_seed, seed
    if n_elements(init_seed) eq 0 then begin
      ;; --- save the seed in case we need to recompute the emissivities 
      seed_file = 'setup_ion.seed.txt'
      if not file_exist(seed_file) then begin
        print, ' *** no pre-existing seed file, creating one *** '
        seed = randomu(s, /long)
        openw, unit, seed_file, /get
        printf, unit, seed
        free_lun, unit
      endif
      seed = long((rd_tfile(seed_file))[0])
      print, ' *** using seed = '+trim(seed)
      init_seed = 0
    endif
      
; read the uncertainty file. search in the working directory.

      if not file_exist('a-values_sigma.txt') then $
         message,'Error, a-values_sigma.txt file not found !' else begin 
         
         da=read_ascii('a-values_sigma.txt')
         lower=fix(reform(da.field1[0,*]))
         upper=fix(reform(da.field1[1,*]))
         perc=reform(da.field1[2,*])
         
         ori_a_value1=a_value1
         
         ind2 = where(wvl1 NE 0., nnl)
         if nnl eq 0 then message, 'Error ! no lines ??! '

         if keyword_set(normal) then rnd=randomn(seed, nnl) else $
             rnd = -1.+2*randomu(seed, nnl)
         
         for ii=0L,nnl-1 do begin 
            
; check if the uncertainty is present, otherwise use the  default
            ind=where(lower eq lvl1[ind2[ii]] and upper eq lvl2[ind2[ii]] ,nn)
            if nn eq 1 then  unc=perc[ind[0]]  else if nn eq 0 then unc=default_uncertainty
            
            a_value1[ind2[ii]]=a_value1[ind2[ii]]*((1. + unc/100.*rnd[ii]) > 0.01) 
            
         endfor 
         
      endelse       
   endif 
   
; plot_oo, ori_a_value1, a_value1, psym=6, syms=0.3, xr=[1, 1e12], yr=  [1, 1e12] 
; plot, a_value1/ ori_a_value1,psym=6, yr=[0,2] 
   
   
;------- end of GDZ addition here ------------------------------------------
   
                                ;
                                ; read proton rates
                                ;
   if keyword_set(noprot) then begin
      pstr=-1
   endif else BEGIN
      read_splups, pname, pstr, pref, /prot
   endelse

   nlvls=max([lvl1,lvl2])
   wvl=fltarr(nlvls,nlvls)
   gf=fltarr(nlvls,nlvls)
   a_value=fltarr(nlvls,nlvls)

   ind1 = where(wvl1 EQ 0.)
   ind2 = where(wvl1 NE 0.)

   wvl[lvl1-1,lvl2-1]= abs(wvl1)
   gf[lvl1-1,lvl2-1]=gf1
   
   IF ind1[0] NE -1 THEN BEGIN
      a_value[lvl1[ind1]-1,lvl2[ind1]-1] = a_value1[ind1]
      a_value[lvl1[ind2]-1,lvl2[ind2]-1] = $
         a_value[lvl1[ind2]-1,lvl2[ind2]-1] + a_value1[ind2]
   ENDIF ELSE BEGIN
      a_value[lvl1[ind2]-1,lvl2[ind2]-1] = a_value1[ind2]
   ENDELSE

                                ;
                                ; Read .elvlc file
                                ;
   read_elvlc,elvlcname,l1a,term,conf,ss,ll,jj,ecm,eryd,ecmth,erydth,eref
   mult=2.*jj+1.
                                ;
   g=where(ecm EQ 0.)
   IF max(g) GT 0 THEN ecm[g]=ecmth[g]

                                ;
                                ; Read .splups file (updated to new scups format, 8-Jan-2014)
                                ;

   read_scups, upsname, splstr
   
   
;------- GDZ addition here ------------------------------------------
   
   if keyword_set(add_uncertainties) then begin      
      
      
      ori_splups=splstr.data
      
; read the uncertainty file. search in the working directory.
      
      if not file_exist('ups_sigma.txt') then $
         message,'Error, ups_sigma.txt file not found !' else begin 
         
         da=read_ascii('ups_sigma.txt')
         lower=fix(reform(da.field1[0,*]))
         upper=fix(reform(da.field1[1,*]))
         perc=reform(da.field1[2,*])
         
         if keyword_set(normal) then rnd=randomn(seed, nnl) else $
; this is a uniform  distribution between -1 and +1 
         rnd = -1.+2*randomu(seed,splstr.info.ntrans) 
         
; this is a different distribution 
;   randomn(seed, splstr.info.ntrans)
         
         
         for ii=0L,splstr.info.ntrans-1 do begin &$
            
; check if the uncertainty is present, otherwise use the  default
            ind=where(lower eq splstr.data[ii].lvl1 and upper eq splstr.data[ii].lvl2 ,nn) &$
            if nn eq 1 then  unc=perc[ind[0]]  else if nn eq 0 then unc=default_uncertainty &$
            
            good = where(splstr.data[ii].spl gt 0, n_good)
            
            if n_good gt 0 then splstr.data[ii].spl[good] = $
               splstr.data[ii].spl[good]*((1. + unc/100.*rnd[ii]) > 0.01 ) &           
            
            endfor 
                        
         endelse        
      

endif   
   
   
   
ENDIF 


END
