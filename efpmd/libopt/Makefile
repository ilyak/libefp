LIBOPT_A= libopt.a
LIBOPT_O= opt.o blas.o lbfgsb.o linpack.o timer.o

AR= ar rc
RANLIB= ranlib

$(LIBOPT_A): $(LIBOPT_O)
	$(AR) $@ $(LIBOPT_O)
	$(RANLIB) $@

clean:
	rm -f $(LIBOPT_A) $(LIBOPT_O)

.PHONY: clean
