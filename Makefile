.PHONY: build clean clean-data smoke

build:
	$(MAKE) -C fortran

clean:
	$(MAKE) -C fortran clean

clean-data:
	find . -type f -name '*.dat' \
		! -path './data/*' \
			! -path './examples/*' \
			! -path './tests/fixtures/*' \
			! -path './tests/testcases_for_code/*' \
		! -path './tests/templates/*' \
		! -path './tests/n1097/older/*' \
		! -path './fortran/hot.dat' \
		-delete
	find tests -maxdepth 1 -type f -name 'parameters.txt' -delete
	find . -type f -name '*.png' \
		! -path './docs/*' \
		! -path './tests/n1097/*' \
		-delete
	rm -rf perl/run0* perl/test

smoke:
	perl tests/smoke.pl
