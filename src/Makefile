MYCFLAGS= -Wall -Wextra -Werror

CC= gcc
CFLAGS= -std=c99 -g $(MYCFLAGS)
LDFLAGS=
LIBS=

AR= ar rcu
RANLIB= ranlib
RM= rm -f
CTAGS= ctags
TAGS= tags

LIBEFP_A= libefp.a
ALL_C= disp.c efp.c elec.c elec_terms.c int.c parse.c pol.c swf.c util.c xr.c
ALL_O= disp.o efp.o elec.o elec_terms.o int.o parse.o pol.o swf.o util.o xr.o

all: $(LIBEFP_A)

$(LIBEFP_A): $(ALL_O)
	$(AR) $@ $(ALL_O)
	$(RANLIB) $@

clean:
	$(RM) $(LIBEFP_A) $(ALL_O) $(TAGS)

$(TAGS): $(ALL_C)
	$(CTAGS) -f $(TAGS) $(ALL_C)

depend:
	@$(CC) $(CFLAGS) -MM $(ALL_C)

.PHONY: all clean depend

disp.o: disp.c private.h efp.h int.h math_util.h swf.h terms.h
efp.o: efp.c elec.h math_util.h private.h efp.h int.h swf.h terms.h
elec.o: elec.c elec.h math_util.h private.h efp.h int.h swf.h terms.h
elec_terms.o: elec_terms.c elec.h math_util.h
int.o: int.c int.h math_util.h int_shift.h
parse.o: parse.c private.h efp.h int.h math_util.h swf.h terms.h
pol.o: pol.c elec.h math_util.h private.h efp.h int.h swf.h terms.h
swf.o: swf.c swf.h math_util.h
util.o: util.c private.h efp.h int.h math_util.h swf.h terms.h
xr.o: xr.c cblas.h private.h efp.h int.h math_util.h swf.h terms.h
