include config.inc

all: efpmd

efpmd: libefp
	cd efpmd/libff && CC="$(CC)" CFLAGS="$(MYCFLAGS)" $(MAKE)
	cd efpmd/libopt && CC="$(CC)" CFLAGS="$(MYCFLAGS)" FC="$(FC)" FFLAGS="$(MYFFLAGS)" $(MAKE)
	cd efpmd/src && $(MAKE)

libefp:
	cd src && CC="$(CC)" CFLAGS="$(MYCFLAGS)" $(MAKE)

clean:
	cd src && $(MAKE) $@
	cd tests && $(MAKE) $@
	cd efpmd/libff && $(MAKE) $@
	cd efpmd/libopt && $(MAKE) $@
	cd efpmd/src && $(MAKE) $@

check checkomp checkmpi: efpmd
	cd tests && $(MAKE) $@

install: all
	install -d $(PREFIX)/bin
	install -d $(PREFIX)/include
	install -d $(PREFIX)/lib
	install -d $(FRAGLIB)/databases
	install -m 0755 efpmd/src/efpmd $(PREFIX)/bin
	install -m 0755 efpmd/tools/cubegen.pl $(PREFIX)/bin
	install -m 0755 efpmd/tools/trajectory.pl $(PREFIX)/bin
	install -m 0644 src/efp.h $(PREFIX)/include
	install -m 0644 src/libefp.a $(PREFIX)/lib
	install -m 0644 fraglib/*.efp $(FRAGLIB)
	install -m 0644 fraglib/makefp.inp $(FRAGLIB)
	install -m 0644 fraglib/databases/* $(FRAGLIB)/databases

dist:
	git archive --format=tar.gz --prefix=libefp/ -o libefp.tar.gz HEAD

.PHONY: all efpmd libefp clean check checkomp checkmpi install dist
