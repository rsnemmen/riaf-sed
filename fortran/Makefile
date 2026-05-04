# The -ffpe-trap option enables exception trapping and makes sure that NaN
# occurrences result in core dumps, i.e. if a NaN happens the program will
# stop. This fixes the problem of the code hanging with 100% CPU usage.

FC = gfortran
FFLAGS = -O
BINDIR = ../bin

all: $(BINDIR)/dynamics $(BINDIR)/spectrum $(BINDIR)/ssd_new $(BINDIR)/ssd_alone

$(BINDIR):
	mkdir -p $@

$(BINDIR)/dynamics: dynamics.f | $(BINDIR)
	$(FC) $(FFLAGS) -ffpe-trap=invalid,zero $< -o $@

$(BINDIR)/spectrum: spectrum.f | $(BINDIR)
	$(FC) $(FFLAGS) -fopenmp $< -o $@

$(BINDIR)/ssd_new: ssd_new.f | $(BINDIR)
	$(FC) $(FFLAGS) $< -o $@

$(BINDIR)/ssd_alone: ssd_alone.f | $(BINDIR)
	$(FC) $(FFLAGS) $< -o $@

clean:
	rm -f $(BINDIR)/dynamics $(BINDIR)/spectrum $(BINDIR)/ssd_new $(BINDIR)/ssd_alone
	rm -f *.o *.mod
