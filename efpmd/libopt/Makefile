TOP= ../..

include $(TOP)/config.inc

CFLAGS= $(MYCFLAGS)
FFLAGS= $(MYFFLAGS)

LIBOPT_A= libopt.a
ALL_O= opt.o blas.o lbfgsb.o linpack.o timer.o

all: $(LIBOPT_A)

$(LIBOPT_A): $(ALL_O)
	$(AR) $@ $(ALL_O)
	$(RANLIB) $@

clean:
	$(RM) $(LIBOPT_A) $(ALL_O)

.PHONY: all clean
