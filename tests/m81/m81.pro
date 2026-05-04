; Compares my model with Quataert's model for M81

; Set output file
set_plot, 'ps'
device,/encapsulated,filename='m81models.eps'

; Uses astro routine plotsym to define filled circles
plotsym, 0, 0.4, /fill

; Defines aspect ratio, look for aspect.pro in Google
; http://www.dfanning.com/tips/ps_aspect.html
;plotPosition = ASPECT(1.)

; Plot
readcol,'/users/grads/nemmen/Work/doutorado/adaf_code/tests/templates/m81_adaf.csv',anu,anulum
readcol,'/users/grads/nemmen/Work/doutorado/adaf_code/tests/templates/m81_thindisk.csv',dnu,dnulum
plot, anu, anulum, xrange=[8.5,20],yrange=[36,43], xstyle=1, ystyle=1, linestyle=0
oplot, dnu, dnulum, linestyle=0

; Box with description of different lines
legend, ['Quataert+99 model','My model with wind'], lines =[0,5]

;plot, lognu, lognulum, psym=8, xrange=[-1.9,1.4],yrange=[-1.4,2.1], xstyle=1, ystyle=1, xtitle="log(P!djet!n/10!u43!n erg s!u-1!n)", ytitle="log(P!dBondi!n/10!u43!n erg s!u-1!n)", position=plotPosition

readcol,'spec_03', anu,anulum
readcol,'spec_03_ssd', dnu,dnulum
oplot, anu, anulum, linestyle=5
oplot, dnu, dnulum, linestyle=5

;readcol,'spec_md04',anu,anulum
;readcol,'spec_md04_ssd',dnu,dnulum
;oplot, anu, anulum, linestyle=3
;oplot, dnu, dnulum, linestyle=3

;readcol,'adafr04',anu,anulum
;readcol,'diskr04',dnu,dnulum
;oplot, anu, anulum, linestyle=5
;oplot, dnu, dnulum, linestyle=5

device,/close
set_plot, 'X'

end
