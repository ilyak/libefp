TOP= ..

include $(TOP)/config.inc

CFLAGS= $(MYCFLAGS)

LIBEFP_A= libefp.a
ALL_O= ai_disp.o balance.o clapack.o disp.o efp.o elec.o elec_terms.o \
       int.o log.o parse.o pol.o pol_direct.o stream.o swf.o util.o xr.o

all: $(LIBEFP_A)

$(LIBEFP_A): $(ALL_O)
	$(AR) $@ $(ALL_O)
	$(RANLIB) $@

clean:
	$(RM) $(LIBEFP_A) $(ALL_O)

.PHONY: all clean
