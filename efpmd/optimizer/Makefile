LIB = optimizer
WARNINGS_SET = HIGH
MK_PROFILE = no
MK_PICLIB = no

CFLAGS += -std=c99

.if defined(WITH_DEBUG)
CFLAGS += -O0 -g
.endif

SRCS = optimizer.c \
       optimizer.h \
       blas.f \
       lbfgsb.f \
       linpack.f \
       timer.f

.include <lib.mk>
