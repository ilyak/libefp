LIB = efp
WARNINGS_SET = HIGH
MK_PROFILE = no
MK_PICLIB = no

LIBEFP_VERSION != cat ${.CURDIR}/../VERSION

CFLAGS += -std=c99
CFLAGS_efp.o += -DPACKAGE_VERSION="\"${LIBEFP_VERSION}\""

.if defined(WITH_DEBUG)
CFLAGS += -O0 -g
.endif

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
