#!/bin/sh
# Compiles the ADAF code.
#
# The -ffpe-trap option enables exception trapping and makes sure that NaN 
# occurrences result in core dumps, i.e. if a NaN happens the program will 
# stop. This fixes the problem of the code hanging with 100% CPU usage.
#
# Be careful with the optimization options you choose, as they may 
# potentially make things worse. Change the architecture according to your 
# own target machine. For Mac Pro and Ubuntu on Vaio I have been using 
# -march=i686 and -O with no problems. 

# Issues with Mac Pro:
# - I noticed that compiling spectrum_new with gfortran gives weird SEDs
#   with discontinuities. The problem vanishes with g77. (2/26/2008)
# - If I compile dynamics_new with gfortran, then I must compile ssd_new
#   with gfortran too, otherwise ssd_new does not understand the ouput of
#   dynamics_new. (2/28/2008)

#g77 -trapfpe dynamics_new.f -o dynamics_new
#g77 spectrum_new.f -o spectrum_new
gfortran -O -ffpe-trap=invalid,zero dynamics_new.f -o dynamics_new
gfortran -O spectrum_new.f -o spectrum_new
gfortran -O ssd_new.f -o ssd_new
gfortran -O ssd_alone.f -o ssd_alone

#g77 -O spectrum_new.f -o spectrum_new
#g77 -O ssd_new.f -o ssd_new
#g77 -O ssd_alone.f -o ssd_alone
