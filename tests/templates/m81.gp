set terminal x11
#set terminal postscript monochrome enhanced eps
unset arrow
set title "M81 test case (Quataert al. 1999)"
set output 'm81.eps'
#plot [:20][36:42]'../../perl/spectrum.dat' with lines title "My ADAF (STD)", \
#'../../perl/ssd.dat' with lines title "My SSD", \
#'../../perl/ses5.dat' with lines title "Without Comptonization", \
#'m81_adaf.csv' title "M81 ADAF", \
#'m81_thindisk.csv' title "M81 SSD"

plot [:20][36:42]'../../perl/spectrum.dat' with lines, \
'../../perl/ssd.dat' with lines, \
'../../perl/ses5.dat' with lines, \
'../../perl/specstd' with lines, \
'../../perl/ssdstd' with lines, \
'../../perl/comptstd' with lines

