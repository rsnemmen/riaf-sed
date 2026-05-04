unset logscale
#set terminal x11
set terminal postscript color enhanced eps
unset arrow
set title "NGC 1097 test case (Nemmen et al. 2006)"
#set arrow from 11.9,40.13 to 12.4,40.13 nohead lt -1 lw 1.2
set output 'sed.eps'
plot [:21][36:44]'specgood.dat' with lines, \
'specvgood.dat' with lines, \
'specvgood01.dat' with lines, \
'specvgood02.dat' with lines, \
'../templates/n1097_spectrum.dat' title "Feng's ADAF", \
'../templates/n1097-thin.dat' title "Feng's SSD"

