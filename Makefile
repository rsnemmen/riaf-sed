.PHONY: build clean clean-data smoke

build:
	$(MAKE) -C fortran

clean:
	$(MAKE) -C fortran clean

clean-data:
	find tests -maxdepth 1 -type f \( -name '*.dat' -o -name 'parameters.txt' \) \
		! -name 'largeR.dat' ! -name 'smallR.dat' -delete
	find perl -maxdepth 1 -type f -name '*.dat' -delete
	find . -path ./data -prune -o -type f \( -name 'aomi*dat' -o -name 'romi*.dat' \) -delete
	find . \( -path ./docs -o -path ./tests/n1097 \) -prune -o -type f -name '*.png' -delete
	rm -rf perl/run0* perl/test

smoke:
	perl tests/smoke.pl
