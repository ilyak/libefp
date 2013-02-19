PROG = test
WARNINGS_SET = HIGH
CLEANFILES += tags

CFLAGS += -std=c99
CFLAGS += -DABS_TOP_SRCDIR="\"${.CURDIR}/..\""

.if defined(WITH_DEBUG)
CFLAGS += -O0 -g
.endif

LDFLAGS += -L${.OBJDIR}/../common
LDFLAGS += -L${.OBJDIR}/../src

LDADD += -lefp -lcommon -lcheck -lm

.if defined(WITH_MKL)
LDADD += -lmkl_rt -lpthread
.elif defined(WITH_ATLAS)
LDADD += -lalapack -lcblas -lf77blas -latlas
.endif

SRCS = test.c \
       disp_1a.c \
       disp_1b.c \
       disp_1c.c \
       disp_1d.c \
       disp_2a.c \
       disp_2b.c \
       disp_3a.c \
       disp_3b.c \
       elec_1a.c \
       elec_1b.c \
       elec_1c.c \
       elec_2a.c \
       elec_2b.c \
       elec_3a.c \
       elec_3b.c \
       pol_1.c \
       pol_2.c \
       pol_3.c \
       qm_1.c \
       qm_2.c \
       total_1a.c \
       total_2a.c \
       total_3a.c \
       total_4a.c \
       total_4b.c \
       total_4c.c \
       total_5a.c \
       total_5b.c \
       total_5c.c \
       total_6a.c \
       total_6b.c \
       total_6c.c \
       xr_1a.c \
       xr_1b.c \
       xr_2.c \
       xr_3.c

run: ${PROG}
	./${PROG}

.include <prog.mk>
