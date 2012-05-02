EXTRA_DIST =
CLEANFILES =
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}
AM_MAKEFLAGS = --no-print-directory

AM_CPPFLAGS = \
	-include $(top_builddir)/config.h \
	-DSYSCONFDIR=\""$(sysconfdir)"\" \
	-DLIBEXECDIR=\""$(libexecdir)"\" \
	-I${top_srcdir}/src

AM_CFLAGS = ${my_CFLAGS} ${OPENMP_CFLAGS} \
	-fvisibility=hidden

AM_LDFLAGS = \
	-Wl,--gc-sections \
	-Wl,--as-needed

SED_PROCESS = \
	$(AM_V_GEN)$(MKDIR_P) $(dir $@) && $(SED) \
	-e 's,@VERSION\@,$(VERSION),g' \
	-e 's,@prefix\@,$(prefix),g' \
	-e 's,@exec_prefix\@,$(exec_prefix),g' \
	-e 's,@libdir\@,$(libdir),g' \
	-e 's,@includedir\@,$(includedir),g' \
	-e 's,@cblas_libs\@,$(cblas_libs),g' \
	< $< > $@ || rm $@

%.pc: %.pc.in Makefile
	$(SED_PROCESS)

LIBEFP_CURRENT = 0
LIBEFP_REVISION = 0
LIBEFP_AGE = 0

include_HEADERS = src/efp.h
lib_LTLIBRARIES = src/libefp.la

src_libefp_la_SOURCES = \
	src/ai_disp.c \
	src/ai_xr.c \
	src/cblas.h \
	src/chtr.c \
	src/disp.c \
	src/disp.h \
	src/efp.c \
	src/efp_private.h \
	src/elec.c \
	src/math_util.h \
	src/parse.c \
	src/pol.c \
	src/terms.h \
	src/xr.c

EXTRA_DIST += LICENSE
EXTRA_DIST += src/libefp.sym
EXTRA_DIST += \
	fraglib/makefp.inp \
	fraglib/h2o.efp \
	fraglib/nh3.efp

src_libefp_la_LDFLAGS = $(AM_LDFLAGS) \
	-version-info $(LIBEFP_CURRENT):$(LIBEFP_REVISION):$(LIBEFP_AGE) \
	-Wl,--version-script=$(top_srcdir)/src/libefp.sym
src_libefp_la_DEPENDENCIES = ${top_srcdir}/src/libefp.sym
src_libefp_la_LIBADD = -lm $(cblas_libs)

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = src/libefp.pc
EXTRA_DIST += src/libefp.pc.in
CLEANFILES += src/libefp.pc

TESTS = \
	tests/disp_1 \
	tests/elec_1 \
	tests/xr_1

check_PROGRAMS = $(TESTS)

TESTS_COMMON_FILES = tests/test_common.c tests/test_common.h
TESTS_LDADD = src/libefp.la
TESTS_CPPFLAGS = $(AM_CPPFLAGS) -DABS_TOP_SRCDIR=\"$(abs_top_srcdir)\"

tests_disp_1_SOURCES = tests/disp_1.c $(TESTS_COMMON_FILES)
tests_disp_1_LDADD = $(TESTS_LDADD)
tests_disp_1_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_elec_1_SOURCES = tests/elec_1.c $(TESTS_COMMON_FILES)
tests_elec_1_LDADD = $(TESTS_LDADD)
tests_elec_1_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_xr_1_SOURCES = tests/xr_1.c $(TESTS_COMMON_FILES)
tests_xr_1_LDADD = $(TESTS_LDADD)
tests_xr_1_CPPFLAGS = $(TESTS_CPPFLAGS)
