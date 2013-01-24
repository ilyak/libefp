PROG = efpmd
WARNINGS_SET = HIGH

USER_LIBS ?= -lalapack -lcblas -lf77blas -latlas -lgfortran -lm
DEBUG_FLAGS ?= -O0 -g

VERSION != cat ${.CURDIR}/../../VERSION

CFLAGS += ${DEBUG_FLAGS}
CFLAGS_parse.o += -DEFP_DATA_DIR="\"${.CURDIR}/../../fraglib\""
CFLAGS_main.o += -DPACKAGE_VERSION="\"${VERSION}\""

LDADD += ${USER_LIBS}
LDADD += -L${.OBJDIR}/../../common -lcommon \
         -L${.OBJDIR}/../../src -lefp \
         -L${.OBJDIR}/../optimizer -loptimizer

SRCS = clapack.c \
       common.c \
       grad.c \
       hess.c \
       main.c \
       md.c \
       opt.c \
       parse.c \
       rand.c \
       sp.c

.include <prog.mk>
