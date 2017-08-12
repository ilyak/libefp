include ../../config.inc

CFLAGS= -DFRAGLIB_PATH="\"$(FRAGLIB)\"" -I../../src -I../libopt -I../libff $(MYCFLAGS)
LDFLAGS= -L../../src -L../libopt -L../libff $(MYLDFLAGS)
LIBS= -lefp -lopt -lff $(MYLIBS) -lm

PROG= efpmd
ALL_O= cfg.o common.o efield.o energy.o grad.o gtest.o hess.o main.o \
       md.o msg.o opt.o parse.o rand.o sp.o

$(PROG): $(ALL_O)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $(ALL_O) $(LIBS)

clean:
	rm -f $(PROG) $(ALL_O)

.PHONY: clean
