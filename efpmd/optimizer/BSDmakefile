LIB = optimizer
WARNINGS_SET = HIGH
MKPROFILE = no
MKPICLIB = no

DEBUG_FLAGS ?= -O0 -g
CFLAGS += ${DEBUG_FLAGS}
FFLAGS += ${DEBUG_FLAGS}

SRCS = optimizer.c \
       optimizer.h \
       blas.f \
       lbfgsb.f \
       linpack.f \
       timer.f

.include <lib.mk>
