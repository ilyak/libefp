LIB = efp
WARNINGS_SET = HIGH
MKPROFILE = no
MKPICLIB = no

DEBUG_FLAGS ?= -O0 -g
VERSION != cat ${.CURDIR}/../VERSION

CFLAGS += ${DEBUG_FLAGS}
CFLAGS_efp.o += -DPACKAGE_VERSION="\"${VERSION}\""

SRCS = disp.c \
       efp.c \
       efp_util.c \
       elec.c \
       elec_grad.c \
       int.c \
       parse.c \
       pol.c \
       swf.c \
       xr.c

.include <lib.mk>
