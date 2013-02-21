MYCFLAGS= -Wall -Wextra -Werror
MYLAPACK= -lalapack -lcblas -lf77blas -latlas

CC= gcc
FC= gfortran
CFLAGS= -std=c99 -g -I../../src $(MYCFLAGS) -DEFP_DATA_DIR="\"$(TOP)/fraglib\""
FFLAGS=
LDFLAGS= -L../../src
LIBS= $(MYLAPACK) -lefp -lm

RM= rm -f
CTAGS= ctags
TAGS= tags

PROG= efpmd
ALL_C= clapack.c common.c grad.c hess.c main.c \
       md.c opt.c optimizer.c parse.c rand.c sp.c
ALL_O= clapack.o common.o grad.o hess.o main.o \
       md.o opt.o optimizer.o parse.o rand.o sp.o \
       blas.o lbfgsb.o linpack.o timer.o

all: $(PROG)

$(PROG): $(ALL_O)
	$(FC) -o $@ $(LDFLAGS) $(ALL_O) $(LIBS)

clean:
	$(RM) $(PROG) $(ALL_O) $(TAGS)

$(TAGS): $(ALL_C)
	$(CTAGS) -f $(TAGS) $(ALL_C)

depend:
	@$(CC) $(CFLAGS) -MM $(ALL_C)

.PHONY: all clean depend

blas.o: lbfgsb/blas.f
	$(FC) -c $(FFLAGS) -o $@ $<

lbfgsb.o: lbfgsb/lbfgsb.f
	$(FC) -c $(FFLAGS) -o $@ $<

linpack.o: lbfgsb/linpack.f
	$(FC) -c $(FFLAGS) -o $@ $<

timer.o: lbfgsb/timer.f
	$(FC) -c $(FFLAGS) -o $@ $<

clapack.o: clapack.c clapack.h
common.o: common.c common.h cfg.h ../../src/efp.h phys.h
grad.o: grad.c common.h cfg.h ../../src/efp.h phys.h
hess.o: hess.c ../../src/math_util.h clapack.h common.h cfg.h \
 ../../src/efp.h phys.h
main.o: main.c common.h cfg.h ../../src/efp.h phys.h
md.o: md.c ../../src/math_util.h common.h cfg.h ../../src/efp.h phys.h \
 rand.h
opt.o: opt.c common.h cfg.h ../../src/efp.h phys.h optimizer.h
optimizer.o: optimizer.c optimizer.h
parse.o: parse.c common.h cfg.h ../../src/efp.h phys.h keywords.h
rand.o: rand.c rand.h
sp.o: sp.c common.h cfg.h ../../src/efp.h phys.h
