PROG = efpmd
WARNINGS_SET = HIGH

LIBEFP_VERSION != cat ${.CURDIR}/../../VERSION

CFLAGS += -std=c99
CFLAGS_parse.o += -DEFP_DATA_DIR="\"${.CURDIR}/../../fraglib\""
CFLAGS_main.o += -DPACKAGE_VERSION="\"${LIBEFP_VERSION}\""

.if defined(WITH_DEBUG)
CFLAGS += -O0 -g
.endif

LDFLAGS += -L${.OBJDIR}/../../common \
           -L${.OBJDIR}/../../src \
           -L${.OBJDIR}/../optimizer

LDADD += -lcommon -lefp -loptimizer -lm
LDADD += -lalapack -lcblas -lf77blas -latlas -lgfortran

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
