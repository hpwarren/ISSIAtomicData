;+
;summary_ea
;	generates effective areas separately for +ve and -ve ACIS/MEG, ACIS/LEG, and HRC/LEG segments for full Capella dataset over 18-23 Ang region
;	computes exposure*EA for each ObsID for each arm, adds them up, and writes out to rdb file
;
;vinay k (2016-oct-30)
;-

;	initialize
evtrdb='Capella_18t23.rdb'
outrdb='CapellaChandra_garf_18t23.rdb'
sep=string(byte(9B))

if n_elements(evt) eq 0 then evt=rdb(evtrdb)
obsid=evt.OBSID & cm2=evt.CM2 & wvl=evt.WVL & obsexp=evt.OBSEXP & bkg=evt.BKGFACTOR & garm=evt.GARM

o1=where(garm eq 1,n1) & o2=where(garm eq 2,n2) & o3=where(garm eq 3,n3) & o4=where(garm eq 4,n4)
obs1=obsid[o1] & obs2=obsid[o2] & obs3=obsid[o3] & obs4=obsid[o4]
uobs1=obs1[uniq(obs1,sort(obs1))] & nu1=n_elements(uobs1)
uobs2=obs2[uniq(obs2,sort(obs2))] & nu2=n_elements(uobs2)
uobs3=obs3[uniq(obs3,sort(obs3))] & nu3=n_elements(uobs3)
uobs4=obs4[uniq(obs4,sort(obs4))] & nu4=n_elements(uobs4)

wmin=18. & wmax=23. & wbin=0.005 & numw=long((wmax-wmin)/wbin+0.5)+1L
wea=findgen(numw)*wbin+wmin & wmid=findgen(numw-1L)*wbin+wmin+0.5*wbin
ea1=dblarr(numw-1L) & ea1m=ea1 & ea1p=ea1
ea2=dblarr(numw-1L) & ea2m=ea2 & ea2p=ea2
ea3=dblarr(numw-1L) & ea3m=ea3 & ea3p=ea3
ea4=dblarr(numw-1L) & ea4m=ea4 & ea4p=ea4

for i=0,nu1-1 do begin
  ok=where(obsid eq uobs1[i] and cm2 gt 1,mok)
  if mok gt 0 then begin
    xp=median(obsexp[ok])
    xx=wvl[ok] & yy=cm2[ok]
    op=where(xx gt 0,mop)
    om=where(xx lt 0,mom)
    if mop gt 0 then begin & xxp=xx[op] & yyp=yy[op] & os=sort(xxp) & ea1p=ea1p+xp*(interpol(yyp[os],xxp[os],wmid)>0) & endif
    if mom gt 0 then begin & xxm=-xx[om] & yym=yy[om] & os=sort(xxm) & ea1m=ea1m+xp*(interpol(yym[os],xxm[os],wmid)>0) & endif
  endif
endfor

for i=0,nu2-1 do begin
  ok=where(obsid eq uobs2[i] and cm2 gt 1,mok)
  if mok gt 0 then begin
    xp=median(obsexp[ok])
    xx=wvl[ok] & yy=cm2[ok]
    op=where(xx gt 0,mop)
    om=where(xx lt 0,mom)
    if mop gt 0 then begin & xxp=xx[op] & yyp=yy[op] & os=sort(xxp) & ea2p=ea2p+xp*(interpol(yyp[os],xxp[os],wmid)>0) & endif
    if mom gt 0 then begin & xxm=-xx[om] & yym=yy[om] & os=sort(xxm) & ea2m=ea2m+xp*(interpol(yym[os],xxm[os],wmid)>0) & endif
  endif
endfor

for i=0,nu3-1 do begin
  ok=where(obsid eq uobs3[i] and cm2 gt 1,mok)
  if mok gt 0 then begin
    xp=median(obsexp[ok])
    xx=wvl[ok] & yy=cm2[ok]
    op=where(xx gt 0,mop)
    om=where(xx lt 0,mom)
    if mop gt 0 then begin & xxp=xx[op] & yyp=yy[op] & os=sort(xxp) & ea3p=ea3p+xp*(interpol(yyp[os],xxp[os],wmid)>0) & endif
    if mom gt 0 then begin & xxm=-xx[om] & yym=yy[om] & os=sort(xxm) & ea3m=ea3m+xp*(interpol(yym[os],xxm[os],wmid)>0) & endif
  endif
endfor

for i=0,nu4-1 do begin
  ok=where(obsid eq uobs4[i] and cm2 gt 1,mok)
  if mok gt 0 then begin
    xp=median(obsexp[ok])
    xx=wvl[ok] & yy=cm2[ok]
    op=where(xx gt 0,mop)
    om=where(xx lt 0,mom)
    if mop gt 0 then begin & xxp=xx[op] & yyp=yy[op] & os=sort(xxp) & ea4p=ea4p+xp*(interpol(yyp[os],xxp[os],wmid)>0) & endif
    if mom gt 0 then begin & xxm=-xx[om] & yym=yy[om] & os=sort(xxm) & ea4m=ea4m+xp*(interpol(yym[os],xxm[os],wmid)>0) & endif
  endif
endfor

ea1=ea1p+ea1m
ea2=ea2p+ea2m
ea3=ea3p+ea3m
ea4=ea4p+ea4m

!p.multi=[0,2,2]

plot,wmid,ea1/1e6,thick=2,charsize=1.4,xtitle='[wvl]',ytitle='[cm!u2!n Msec]',title='ACIS/HEG',/xs
oplot,wmid,ea1p/1e6,thick=2,col=200
oplot,wmid,ea1m/1e6,thick=2,line=2
plot,wmid,ea2/1e6,thick=2,charsize=1.4,xtitle='[wvl]',ytitle='[cm!u2!n Msec]',title='ACIS/MEG',/xs
oplot,wmid,ea2p/1e6,thick=2,col=200
oplot,wmid,ea2m/1e6,thick=2,line=2
plot,wmid,ea3/1e6,thick=2,charsize=1.4,xtitle='[wvl]',ytitle='[cm!u2!n Msec]',title='ACIS/LEG',/xs
oplot,wmid,ea3p/1e6,thick=2,col=200
oplot,wmid,ea3m/1e6,thick=2,line=2
plot,wmid,ea4/1e6,thick=2,charsize=1.4,xtitle='[wvl]',ytitle='[cm!u2!n Msec]',title='HRC/LEG',/xs
oplot,wmid,ea4p/1e6,thick=2,col=200
oplot,wmid,ea4m/1e6,thick=2,line=2
xyouts,!x.crange[0],!y.crange[0]-0.2*(!y.crange[1]-!y.crange[0]),'red +ve, dash -ve',align=1,charsize=1.4,charthick=2

!p.multi=0

openw,uo,outrdb,/get_lun
printf,uo,'# Combined Chandra effective areas for different detector/grating combinations for Capella observations'
printf,uo,'# made on '+systime()+' with summary_ea.pro'
printf,uo,'# '
printf,uo,'# Columns:
printf,uo,'# WAVE -- wavelength [Ang]'
printf,uo,'# EA_AH -- exposure weighted effective area for ACIS/HEG [cm^2*sec]'
printf,uo,'# EA_AH_P -- exposure weighted effective area for +1 order ACIS/HEG [cm^2*sec]'
printf,uo,'# EA_AH_M -- exposure weighted effective area for -1 order ACIS/HEG [cm^2*sec]'
printf,uo,'# EA_AM -- exposure weighted effective area for ACIS/MEG [cm^2*sec]'
printf,uo,'# EA_AM_P -- exposure weighted effective area for +1 order ACIS/MEG [cm^2*sec]'
printf,uo,'# EA_AM_M -- exposure weighted effective area for -1 order ACIS/MEG [cm^2*sec]'
printf,uo,'# EA_AL -- exposure weighted effective area for ACIS/LEG [cm^2*sec]'
printf,uo,'# EA_AL_P -- exposure weighted effective area for +1 order ACIS/LEG [cm^2*sec]'
printf,uo,'# EA_AL_M -- exposure weighted effective area for -1 order ACIS/LEG [cm^2*sec]'
printf,uo,'# EA_HL -- exposure weighted effective area for HRC/LEG [cm^2*sec]'
printf,uo,'# EA_HL_P -- exposure weighted effective area for +1 order HRC/LEG [cm^2*sec]'
printf,uo,'# EA_HL_M -- exposure weighted effective area for -1 order HRC/LEG [cm^2*sec]'
printf,uo,'# '
printf,uo,'WAVE'+sep+'EA_AH'+sep+'EA_AH_P'+sep+'EA_AH_M'+sep+'EA_AM'+sep+'EA_AM_P'+sep+'EA_AM_M'+sep+'EA_AL'+sep+'EA_AL_P'+sep+'EA_AL_M'+sep+'EA_HL'+sep+'EA_HL_P'+sep+'EA_HL_M'
printf,uo,'N'+sep+'N'+sep+'N'+sep+'N'+sep+'N'+sep+'N'+sep+'N'+sep+'N'+sep+'N'+sep+'N'+sep+'N'+sep+'N'+sep+'N'
for i=0L,numw-2L do begin
  cc=strtrim(string(wmid[i],'(f8.5)'),2)+sep+$
     strtrim(string(ea1[i],'(f12.1)'),2)+sep+strtrim(string(ea1p[i],'(f12.1)'),2)+sep+strtrim(string(ea1m[i],'(f12.1)'),2)+sep+$
     strtrim(string(ea2[i],'(f12.1)'),2)+sep+strtrim(string(ea2p[i],'(f12.1)'),2)+sep+strtrim(string(ea2m[i],'(f12.1)'),2)+sep+$
     strtrim(string(ea3[i],'(f12.1)'),2)+sep+strtrim(string(ea3p[i],'(f12.1)'),2)+sep+strtrim(string(ea3m[i],'(f12.1)'),2)+sep+$
     strtrim(string(ea4[i],'(f12.1)'),2)+sep+strtrim(string(ea4p[i],'(f12.1)'),2)+sep+strtrim(string(ea4m[i],'(f12.1)'),2)
  printf,uo,cc
endfor
close,uo & free_lun,uo

end
