.PHONY: build clean clean-data smoke

build:
	$(MAKE) -C fortran

clean:
	$(MAKE) -C fortran clean

clean-data:
	find examples -maxdepth 1 -type f -name '*.dat' \
		! -name 'largeR.dat' ! -name 'smallR.dat' -delete
	rm -f examples/parameters.txt
	rm -rf examples/out
	find perl -maxdepth 1 -type f -name '*.dat' -delete
	rm -rf perl/run0* perl/test

smoke:
	perl tests/smoke.pl
