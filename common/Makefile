LIB = common
WARNINGS_SET = HIGH
MK_PROFILE = no
MK_PICLIB = no
CLEANFILES += tags

CFLAGS += -std=c99

.if defined(WITH_DEBUG)
CFLAGS += -O0 -g
.endif

SRCS = cblas.h \
       compat.c \
       compat.h \
       math_util.h \
       phys_const.h \
       util.c \
       util.h

.include <lib.mk>
