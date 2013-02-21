V= 0.9.8-beta

all: efpmd tests

efpmd: libefp
	cd efpmd/src && $(MAKE)

tests: libefp
	cd tests && $(MAKE)

libefp:
	cd src && $(MAKE)

tags clean:
	cd src && $(MAKE) $@
	cd efpmd/src && $(MAKE) $@
	cd tests && $(MAKE) $@

check: tests
	@./tests/test

dist:
	git archive --format=tar.gz --prefix=libefp-$V/ -o libefp-$V.tar.gz HEAD

.PHONY: all efpmd tests libefp tags clean check dist
