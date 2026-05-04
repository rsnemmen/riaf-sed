.PHONY: build clean clean-data smoke

build:
	$(MAKE) -C fortran

clean:
	$(MAKE) -C fortran clean

clean-data:
	find tests -maxdepth 1 -type f \( -name '*.dat' -o -name 'parameters.txt' \) \
		! -name 'largeR.dat' ! -name 'smallR.dat' -delete
	find perl -maxdepth 1 -type f -name '*.dat' -delete
	rm -rf perl/run0* perl/test

smoke:
	perl tests/smoke.pl
