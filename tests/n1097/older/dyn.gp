unset logscale
set logscale x
set terminal x11
set terminal postscript color enhanced eps
set output 'spec_angmom.eps'
unset arrow
x=12
plot 'dynvgood.dat' using 1:x, 'dynvgood01.dat' using 1:x, 'dynvgood02.dat' using 1:x
