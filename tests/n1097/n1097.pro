; Compares my model with Feng's model for NGC 1097

; Set output file
set_plot, 'ps'
device,/encapsulated,filename='n1097.eps'

; Uses astro routine plotsym to define filled circles
plotsym, 0, 0.4, /fill

; Defines aspect ratio, look for aspect.pro in Google
; http://www.dfanning.com/tips/ps_aspect.html
;plotPosition = ASPECT(1.)

; Plot
readcol,'/users/grads/nemmen/Work/doutorado/adaf_code/tests/templates/n1097_adaf.dat',anu,anulum
readcol,'/users/grads/nemmen/Work/doutorado/adaf_code/tests/templates/n1097_thin.dat',dnu,dnulum
plot, anu, anulum, xrange=[17,19],yrange=[39.5,41], xstyle=1, ystyle=1, linestyle=0
oplot, dnu, dnulum, linestyle=0
;plot, anu, anulum, xrange=[8.5,20],yrange=[36,45], xstyle=1, ystyle=1, linestyle=0

; Box with description of different lines
legend, ['Nemmen+06 model','mdot=1E-3, p=0.2, delta=0.1','mdot=1E-3, p=0.2, delta=0.3','Best fit: mdot=1.1E-3, p=0.2, delta=0.1'], lines =[0,1,2,3], /right, /bottom

;plot, lognu, lognulum, psym=8, xrange=[-1.9,1.4],yrange=[-1.4,2.1], xstyle=1, ystyle=1, xtitle="log(P!djet!n/10!u43!n erg s!u-1!n)", ytitle="log(P!dBondi!n/10!u43!n erg s!u-1!n)", position=plotPosition

readcol,'spec_a01', anu,anulum
readcol,'spec_a01_ssd', dnu,dnulum
oplot, anu, anulum, linestyle=1
oplot, dnu, dnulum, linestyle=1

;readcol,'spec_a02',anu,anulum
;readcol,'spec_a02_ssd',dnu,dnulum
;oplot, anu, anulum, linestyle=2
;oplot, dnu, dnulum, linestyle=2

;readcol,'spec_std',anu,anulum
;readcol,'spec_std_ssd',dnu,dnulum
;oplot, anu, anulum, linestyle=3
;oplot, dnu, dnulum, linestyle=3

readcol,'spec_b01', anu,anulum
readcol,'spec_b01_ssd', dnu,dnulum
oplot, anu, anulum, linestyle=2
oplot, dnu, dnulum, linestyle=2

;readcol,'spec_b02',anu,anulum
;readcol,'spec_b02_ssd',dnu,dnulum
;oplot, anu, anulum, linestyle=5
;oplot, dnu, dnulum, linestyle=5

readcol,'spec_01',anu,anulum
readcol,'spec_01_ssd',dnu,dnulum
oplot, anu, anulum, linestyle=3
oplot, dnu, dnulum, linestyle=3 

;readcol,'/users/grads/nemmen/tmp.dat',anu,anulum
;oplot, anu, anulum, psym=7

; X-ray spectrum
; ===============

gamma=1.64 ; best-fit value of Gamma    *CHANGE*
ugammaerr=0.13 ; upper error     *CHANGE*
lgammaerr=0.07 ; lower error     *CHANGE*
L2_10=4.352d40 ; 2-10 keV luminosity    *CHANGE*

h=0.6626d-26 ; Planck constant (CGS)
conv=1.602171d-9/h ; conversion factor keV -> Hz
nu0=2.*conv ; 2 keV -> Hz

; Creates a vector of frequencies
points=10
nui=nu0   ;nu0 ; 2 keV
nuf=10.*conv ; 10 keV
nux=findgen(points)/(points-1)*(nuf-nui) + nui

; Creates vectors containing the X-ray powerlaws, and plot them
alpha=gamma-1. ; best-fit
L0=nu0^(-alpha)*(1.-alpha)/((10.*conv)^(1.-alpha)-(2.*conv)^(1.-alpha))*L2_10
Lbest=L0*(nux/nu0)^(-alpha)
oplot, alog10(nux), alog10(nux*Lbest)
;plot, alog10(nux), alog10(nux*Lbest), yrange=[38,42]
alpha=gamma-lgammaerr-1. ; lower
L0=nu0^(-alpha)*(1.-alpha)/((10.*conv)^(1.-alpha)-(2.*conv)^(1.-alpha))*L2_10
Llow=L0*(nux/nu0)^(-alpha)
oplot, alog10(nux), alog10(nux*Llow)
alpha=gamma+ugammaerr-1. ; upper
L0=nu0^(-alpha)*(1.-alpha)/((10.*conv)^(1.-alpha)-(2.*conv)^(1.-alpha))*L2_10
Lup=L0*(nux/nu0)^(-alpha)
oplot, alog10(nux), alog10(nux*Lup)





device,/close
set_plot, 'X'

end
