;+
;Capella_OVIII_flux_variations
;	make a plot of the Capella line flux variations for O VIII to put into response to referee
;
;vinay k (2023-nov-28)
;-

peasecolr & loadct,3 & blindcolr,/white & thk=3 & csz=3 & !p.background=255 & window,0,xsize=1500,ysize=1000

if n_tags(ctaln) eq 0 then ctaln=rdb('aln_ct_17t23.rdb')
if n_tags(ctalp) eq 0 then ctalp=rdb('alp_ct_17t23.rdb')
if n_tags(ctamn) eq 0 then ctamn=rdb('amn_ct_17t23.rdb')
if n_tags(ctamp) eq 0 then ctamp=rdb('amp_ct_17t23.rdb')
if n_tags(eaaln) eq 0 then eaaln=rdb('aln_ea_17t23.rdb')
if n_tags(eaalp) eq 0 then eaalp=rdb('alp_ea_17t23.rdb')
if n_tags(eaamn) eq 0 then eaamn=rdb('amn_ea_17t23.rdb')
if n_tags(eaamp) eq 0 then eaamp=rdb('amp_ea_17t23.rdb')
if n_tags(obs) eq 0 then obs=rdb('Capella_obslog.rdb')

iam=where(obs.CONFIG eq 'am',miam) & ial=where(obs.CONFIG eq 'al',mial)

obsid=obs.OBSID
obsam=obsid[iam] & obsal=obsid[ial]
cobsam=strtrim(obsam,2) & cobsal=strtrim(obsal,2)
sobsam=string(obsam,'(i5.5)') & sobsal=string(obsal,'(i5.5)')
yram=(obs.yrobs)[iam] & yral=(obs.yrobs)[ial]

tmp=[obsam,obsal] & os=sort(tmp) & print,strtrim(tmp[os],2)

o8camp=fltarr(miam) & o8camn=o8camp & o8calp=fltarr(mial) & o8caln=o8calp
o8campe=fltarr(miam) & o8camne=o8campe & o8calpe=fltarr(mial) & o8calne=o8calpe
o8Aamp=fltarr(miam) & o8Aamn=o8Aamp & o8Aalp=fltarr(mial) & o8Aaln=o8Aalp
o8famp=fltarr(miam) & o8famn=o8famp & o8falp=fltarr(mial) & o8faln=o8falp
o8fampe=fltarr(miam) & o8famne=o8fampe & o8falpe=fltarr(mial) & o8falne=o8falpe

namamp=tag_names(ctamp) & namamn=tag_names(ctamn) & namalp=tag_names(ctalp) & namaln=tag_names(ctaln)
print,long(strmid(namamp,8,5))-long(strmid(namamn,8,5))
print,long(strmid(namalp,8,5))-long(strmid(namaln,8,5))
namAam=tag_names(eaamp) & namAal=tag_names(eaalp)

wlo=ctamp.WLOW & whi=ctamp.WHIGH & wmid=(wlo+whi)/2.

backscal=1.0
for i=0,miam-1 do begin
  o1=where(namamp eq 'SCT_POS_'+sobsam[i],mo1)
  o2=where(namamp eq 'BCT_POS_'+sobsam[i],mo2)
  o3=where(namAam eq 'EAP_'+sobsam[i],mo3)
  cts=ctamp.(o1[0])
  ctb=ctamp.(o2[0])
  ea=eaamp.(o3[0])
  ok=where(wlo ge 18.5 and whi le 19.5,mok)
  o8camp[i]=total(cts[ok]-ctb[ok]/backscal) & o8campe[i]=sqrt(total(cts[ok]+ctb[ok]/backscal^2))
  o8Aamp[i]=mean(ea[ok])
  o8famp[i]=o8camp[i]/o8Aamp[i] & o8fampe[i]=o8campe[i]/o8Aamp[i]
endfor

for i=0,miam-1 do begin
  o1=where(namamn eq 'SCT_NEG_'+sobsam[i],mo1)
  o2=where(namamn eq 'BCT_NEG_'+sobsam[i],mo2)
  o3=where(namAam eq 'EAP_'+sobsam[i],mo3)
  cts=ctamn.(o1[0])
  ctb=ctamn.(o2[0])
  ea=eaamn.(o3[0])
  ok=where(wlo ge 18.5 and whi le 19.5,mok)
  o8camn[i]=total(cts[ok]-ctb[ok]/backscal) & o8camne[i]=sqrt(total(cts[ok]+ctb[ok]/backscal^2))
  o8Aamn[i]=mean(ea[ok])
  o8famn[i]=o8camn[i]/o8Aamn[i] & o8famne[i]=o8camne[i]/o8Aamn[i]
endfor

for i=0,mial-1 do begin
  o1=where(namalp eq 'SCT_POS_'+sobsal[i],mo1)
  o2=where(namalp eq 'BCT_POS_'+sobsal[i],mo2)
  o3=where(namAal eq 'EAP_'+sobsal[i],mo3)
  cts=ctalp.(o1[0])
  ctb=ctalp.(o2[0])
  ea=eaalp.(o3[0])
  ok=where(wlo ge 18.5 and whi le 19.5,mok)
  o8calp[i]=total(cts[ok]-ctb[ok]/backscal) & o8calpe[i]=sqrt(total(cts[ok]+ctb[ok]/backscal^2))
  o8Aalp[i]=mean(ea[ok])
  o8falp[i]=o8calp[i]/o8Aalp[i] & o8falpe[i]=o8calpe[i]/o8Aalp[i]
endfor

for i=0,mial-1 do begin
  o1=where(namaln eq 'SCT_NEG_'+sobsal[i],mo1)
  o2=where(namaln eq 'BCT_NEG_'+sobsal[i],mo2)
  o3=where(namAal eq 'EAP_'+sobsal[i],mo3)
  cts=ctaln.(o1[0])
  ctb=ctaln.(o2[0])
  ea=eaaln.(o3[0])
  ok=where(wlo ge 18.5 and whi le 19.5,mok)
  o8caln[i]=total(cts[ok]-ctb[ok]/backscal) & o8calne[i]=sqrt(total(cts[ok]+ctb[ok]/backscal^2))
  o8Aaln[i]=mean(ea[ok])
  o8faln[i]=o8caln[i]/o8Aaln[i] & o8falne[i]=o8calne[i]/o8Aaln[i]
endfor

yr=[yram,yram,yral,yral] & o8f=[o8famp,o8famn,o8falp,o8faln] & o8fe=[o8fampe,o8famne,o8falpe,o8falne]
os=sort(yr) & yr=yr[os] & o8f=o8f[os] & o8fe=o8fe[os]

plot,yr,o8f,psym=8,xr=[1999,2020],yr=[2,8],/xs,/ys,$
	xtitle='Year',ytitle='O VIII flux [ph/ks/cm^2]',title='Capella variations',$
	xticklen=1,yticklen=1,xgridstyle=1,ygridstyle=1,$
	thick=thk,xthick=thk,ythick=thk,charthick=thk,charsize=csz,col=0
poaintsym,'circle',/pfil,psiz=1
oplot,yram,o8famp,psym=8,col=1,symsize=csz & for i=0,miam-1 do oplot,yram[i]*[1,1],o8famp[i]+o8fampe[i]*[-1,1],thick=thk,col=1
oplot,yral,o8falp,psym=8,col=2,symsize=csz & for i=0,mial-1 do oplot,yral[i]*[1,1],o8falp[i]+o8falpe[i]*[-1,1],thick=thk,col=2
plots,!x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]),!y.crange[1]-0.09*(!y.crange[1]-!y.crange[0]),psym=8,symsize=csz,col=0
xyouts,!x.crange[0]+0.08*(!x.crange[1]-!x.crange[0]),!y.crange[1]-0.1*(!y.crange[1]-!y.crange[0]),'amp',charsize=csz,charthick=thk,col=1
xyouts,!x.crange[0]+0.15*(!x.crange[1]-!x.crange[0]),!y.crange[1]-0.1*(!y.crange[1]-!y.crange[0]),'alp',charsize=csz,charthick=thk,col=2
poaintsym,'box',/pfil,psiz=1
oplot,yram,o8famn,psym=8,col=6,symsize=csz & for i=0,miam-1 do oplot,yram[i]*[1,1],o8famn[i]+o8famne[i]*[-1,1],thick=thk,col=6
oplot,yral,o8faln,psym=8,col=5,symsize=csz & for i=0,mial-1 do oplot,yral[i]*[1,1],o8faln[i]+o8falne[i]*[-1,1],thick=thk,col=5
plots,!x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]),!y.crange[1]-0.15*(!y.crange[1]-!y.crange[0]),psym=8,symsize=csz,col=0
xyouts,!x.crange[0]+0.08*(!x.crange[1]-!x.crange[0]),!y.crange[1]-0.16*(!y.crange[1]-!y.crange[0]),'amn',charsize=csz,charthick=thk,col=6
xyouts,!x.crange[0]+0.15*(!x.crange[1]-!x.crange[0]),!y.crange[1]-0.16*(!y.crange[1]-!y.crange[0]),'aln',charsize=csz,charthick=thk,col=5
write_png,'Capella_OVIII_flux_variations.png',tvrd(/true)

end
