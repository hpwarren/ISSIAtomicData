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
; VERSION     :   9, 27-Jun-2014
;
;-


PRO setup_ion,ions,wvlmin,wvlmax,wvltst,lvl1,lvl2,wvl1,gf1,a_value1, $
                 path=path, noprot=noprot, all=all, anylines=anylines

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
  
  ;; ##########################################################################################
  ;; Modifications by HPW to introduce white noise perturbations to the collison
  ;; strengths. Transitions to the ground state and other transitions are scaled differently.

  ;; --- use environment variable to control the behavior of the rountine
  ;; in the driver routine use:
  ;;   setenv, 'CHIANTI_PERTURB_SPL1=0.1'
  ;;   setenv, 'CHIANTI_PERTURB_SPL2=0.2'  
  CHIANTI_PERTURB_SPL1 = float(getenv('CHIANTI_PERTURB_SPL1')) ;; ground state
  CHIANTI_PERTURB_SPL2 = float(getenv('CHIANTI_PERTURB_SPL2')) ;; other states
  
  if keyword_set(CHIANTI_PERTURB_SPL1) OR keyword_set(CHIANTI_PERTURB_SPL2) then begin
    print, ' | Perturbing the effective collision strengths '
    print, ' |   Ground config sigma = '+trim(CHIANTI_PERTURB_SPL1,'(f10.2)')
    print, ' |   Other config sigma  = '+trim(CHIANTI_PERTURB_SPL2,'(f10.2)')

    ;; --- option to output information for each level, delete opf to generate a new file
    opf = ions+'.monte_carlo.txt'
    if not(file_exist(opf)) then begin
      write_to_file = 1
      openw, unit, opf, /get
      print, ' | writing example to file'
      printf, unit, '# This is an example, not all realizations of the atomic data are saved'
      printf, unit, '# Time Stamp = '+systime(0)
      printf, unit, '# '
      printf, unit, '# Perturbing the effective collision strengths '
      printf, unit, '#   Ground config sigma = '+trim(CHIANTI_PERTURB_SPL1,'(f10.2)')
      printf, unit, '#   Other config sigma  = '+trim(CHIANTI_PERTURB_SPL2,'(f10.2)')
      printf, unit, '#'
      printf, unit, '#', 'CONF', 'LEVEL1', 'LEVEL2', 'SIGMA', 'RND', format='(a1,a12,2a13,2a13)'
    endif else write_to_file = 0

    ;; --- find the levels for the ground configuration
    index_ground_conf_levels = where(conf eq 1, n_ground_conf_levels)
    ground_conf_levels = l1a[index_ground_conf_levels]
    
    rnd = randomn(seed, splstr.info.ntrans)
    
    for i=0,splstr.info.ntrans-1 do begin

      ;; --- does this level match a ground configuration level?
      match_ground_levels = where(ground_conf_levels eq splstr.data[i].lvl1, n_match_ground_levels)
      if n_match_ground_levels gt 0 then begin
        CHIANTI_PERTURB = CHIANTI_PERTURB_SPL1
      endif else CHIANTI_PERTURB = CHIANTI_PERTURB_SPL2
      perturb = 1 + CHIANTI_PERTURB*rnd[i]      

      ;; --- perturb; spl is -1 for missing data, so don't overwrite
      good = where(splstr.data[i].spl ge 0, n_good)
      if n_good gt 0 then splstr.data[i].spl[good] = splstr.data[i].spl[good]*perturb
      
      if write_to_file then begin
        match_level = where(l1a eq splstr.data[i].lvl1, n_match)
        if n_match eq 0 then stop
        printf, unit, conf[match_level[0]], splstr.data[i].lvl1, splstr.data[i].lvl2, $
                CHIANTI_PERTURB, perturb, format='(3i13, f13.1, f13.2)'
      endif
    endfor

    if write_to_file then free_lun, unit

  endif

  ;; End modifications by HPW
  ;; ##########################################################################################

; PRY 27-Jun-2014: I commented out the two lines below.
;  splref=splstr.info.comments
;  splstr=splstr.data
; read_splups,upsname,splstr,splref

ENDIF


END
