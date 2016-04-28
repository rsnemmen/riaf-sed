# The -ffpe-trap option enables exception trapping and makes sure that NaN 
# occurrences result in core dumps, i.e. if a NaN happens the program will 
# stop. This fixes the problem of the code hanging with 100% CPU usage.
#
# Old school way of compiling:
# gfortran -O -ffpe-trap=invalid,zero dynamics_new.f -o dynamics_new
# gfortran -O spectrum_new.f -o spectrum_new
# gfortran -O ssd_new.f -o ssd_new
# gfortran -O ssd_alone.f -o ssd_alone

FC = gfortran
FFLAGS = -O 
#OBJS = pram.o main.o bnd.o calcha.o calflx.o calsf.o caluu.o cdtcfl.o correction.o coord.o ct.o func.o grid.o hll.o mdgrmhd.o mpisub.o physfunc.o rec.o recov.o restar.o rkt.o

all: dynamics spectrum ssd ssd_alone

dynamics: dynamics_new.f
	$(FC) $(FFLAGS) -ffpe-trap=invalid,zero dynamics_new.f -o dynamics_new

spectrum: spectrum_new.f
	$(FC) $(FFLAGS) spectrum_new.f -o spectrum_new

ssd: ssd_new.f
	$(FC) $(FFLAGS) ssd_new.f -o ssd_new

ssd_alone: ssd_alone.f
	$(FC) $(FFLAGS) ssd_alone.f -o ssd_alone

clean:
	rm -f *.o *.mod *_new *_alone

