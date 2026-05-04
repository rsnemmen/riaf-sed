#set terminal x11
set terminal postscript color enhanced eps
unset arrow
set title "NGC 1097 test case (Nemmen et al. 2006)"
set arrow from 11.9,40.13 to 12.4,40.13 nohead lt -1 lw 1.2
set output 'n1097.eps'
plot [][36:44]'tests/spectrum_sent.dat' with lines title "My ADAF", \
'tests/ssd_sent.dat' with lines title "My SSD (mdot=6E-4)", \
'tests/ssdstrong_sent.dat' with lines title "My SSD (mdot=6E-3)", \
'tests/n1097_spectrum.dat' title "Feng's ADAF", \
'tests/n1097-thin.dat' title "Feng's SSD"

