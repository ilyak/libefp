TOP= ../..

include $(TOP)/config.inc

CFLAGS= -DFRAGLIB_PATH="\"$(FRAGLIB)\"" -I$(TOP)/src -I$(TOP)/efpmd/libopt -I$(TOP)/efpmd/libff $(MYCFLAGS)
LDFLAGS= -L$(TOP)/src -L$(TOP)/efpmd/libopt -L$(TOP)/efpmd/libff $(MYLDFLAGS)
LIBS= -lefp -lopt -lff $(MYLIBS) -lm

PROG= efpmd
ALL_O= cfg.o common.o efield.o energy.o grad.o gtest.o hess.o main.o \
       md.o msg.o opt.o parse.o rand.o sp.o

all: $(PROG)

$(PROG): $(ALL_O)
	$(CC) -o $@ $(CFLAGS) $(LDFLAGS) $(ALL_O) $(LIBS)

clean:
	$(RM) $(PROG) $(ALL_O)

.PHONY: all clean
