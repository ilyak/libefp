LIBEFP_A= libefp.a
LIBEFP_O= aidisp.o balance.o clapack.o disp.o efp.o elec.o \
	  electerms.o int.o log.o parse.o pol.o poldirect.o \
	  stream.o swf.o util.o xr.o

AR= ar rc
RANLIB= ranlib

$(LIBEFP_A): $(LIBEFP_O)
	$(AR) $@ $(LIBEFP_O)
	$(RANLIB) $@

clean:
	rm -f $(LIBEFP_A) $(LIBEFP_O)

.PHONY: clean
