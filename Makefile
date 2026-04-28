.PHONY: build clean smoke

build:
	$(MAKE) -C fortran

clean:
	$(MAKE) -C fortran clean

smoke:
	perl tests/smoke.pl
