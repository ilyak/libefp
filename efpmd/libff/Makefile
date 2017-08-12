LIBFF_A= libff.a
LIBFF_O= ff.o

AR= ar rc
RANLIB= ranlib

$(LIBFF_A): $(LIBFF_O)
	$(AR) $@ $(LIBFF_O)
	$(RANLIB) $@

clean:
	rm -f $(LIBFF_A) $(LIBFF_O)

.PHONY: clean
