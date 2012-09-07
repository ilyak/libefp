VERSION_MAJOR = 0
VERSION_MINOR = 9

EXTRA_DIST =
CLEANFILES =
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}
AM_MAKEFLAGS = --no-print-directory
AM_LDFLAGS =

AM_CPPFLAGS = \
	-include $(top_builddir)/config.h \
	-DSYSCONFDIR=\""$(sysconfdir)"\" \
	-DLIBEXECDIR=\""$(libexecdir)"\" \
	-DLIBEFPVERSION=\""$(VERSION_MAJOR).$(VERSION_MINOR)"\" \
	-I${top_srcdir}/src

AM_CFLAGS = ${my_CFLAGS} ${OPENMP_CFLAGS} \
	-fvisibility=hidden

include_HEADERS = src/efp.h
lib_LTLIBRARIES = src/libefp.la

dist_pkgdata_DATA = \
	fraglib/makefp.inp \
	fraglib/acetone.efp \
	fraglib/c2h5oh.efp \
	fraglib/c6h6.efp \
	fraglib/ccl4.efp \
	fraglib/ch3oh.efp \
	fraglib/ch4.efp \
	fraglib/cl2.efp \
	fraglib/dcm.efp \
	fraglib/dmso.efp \
	fraglib/h2.efp \
	fraglib/h2o.efp \
	fraglib/nh3.efp

src_libefp_la_SOURCES = \
	common/compat.c \
	common/compat.h \
	common/math_util.h \
	src/ai_disp.c \
	src/ai_xr.c \
	src/cblas.h \
	src/chtr.c \
	src/disp.c \
	src/disp.h \
	src/efp.c \
	src/efp_private.h \
	src/elec.c \
	src/elec.h \
	src/elec_grad.c \
	src/int.c \
	src/int.h \
	src/int_shift.h \
	src/parse.c \
	src/parse.h \
	src/pol.c \
	src/terms.h \
	src/xr.c

EXTRA_DIST += AUTHORS LICENSE README.md

src_libefp_la_LDFLAGS = $(AM_LDFLAGS) \
	-version-info $(VERSION_MAJOR):$(VERSION_MINOR)

if BUILD_EFPMD
bin_PROGRAMS = efpmd/efpmd
noinst_PROGRAMS = efpmd/cubegen/cubegen
endif

efpmd_efpmd_SOURCES = \
	common/compat.c \
	common/compat.h \
	common/math_util.h \
	common/phys_const.h \
	efpmd/common.c \
	efpmd/common.h \
	efpmd/config.h \
	efpmd/grad.c \
	efpmd/hess.c \
	efpmd/main.c \
	efpmd/md.c \
	efpmd/opt.c \
	efpmd/parse.c \
	efpmd/sp.c \
	efpmd/optimizer/blas.f \
	efpmd/optimizer/lbfgsb.f \
	efpmd/optimizer/linpack.f \
	efpmd/optimizer/timer.f \
	efpmd/optimizer/optimizer.c \
	efpmd/optimizer/optimizer.h

efpmd_efpmd_LDADD = src/libefp.la
efpmd_efpmd_CFLAGS = $(AM_CFLAGS)
efpmd_efpmd_CPPFLAGS = $(AM_CPPFLAGS) \
	-DDATADIR=\""$(datarootdir)/libefp"\" \
	-DEFPMDVERSION=\""$(VERSION_MAJOR).$(VERSION_MINOR)"\"

efpmd_cubegen_cubegen_SOURCES = \
	efpmd/cubegen/cubegen.c

EXTRA_DIST += \
	efpmd/README.md \
	efpmd/optimizer/License.txt \
	efpmd/examples/c6h6-dimer_opt.in \
	efpmd/examples/h2o-ch3oh_grad.in \
	efpmd/examples/h2o-dimer_opt.in \
	efpmd/examples/h2o-nh3_grad.in \
	efpmd/examples/h2o-tetramer_opt.in \
	efpmd/examples/water-cube-125_grad.in \
	efpmd/examples/water-cube-125_sp.in \
	efpmd/examples/water-cube-3375_sp.in

TESTS = \
	tests/elec_1a \
	tests/elec_1b \
	tests/pol_1 \
	tests/disp_1a \
	tests/disp_1b \
	tests/xr_1 \
	tests/total_1a \
	tests/elec_2a \
	tests/elec_2b \
	tests/pol_2 \
	tests/disp_2a \
	tests/disp_2b \
	tests/xr_2 \
	tests/total_2a \
	tests/elec_3a \
	tests/elec_3b \
	tests/pol_3 \
	tests/disp_3a \
	tests/disp_3b \
	tests/xr_3 \
	tests/total_3a \
	tests/total_4a \
	tests/total_4b \
	tests/total_4c \
	tests/qm_5 \
	tests/qm_6

EXTRA_DIST += \
	tests/geometry_1.h \
	tests/geometry_2.h \
	tests/geometry_3.h \
	tests/geometry_4.h

check_PROGRAMS = $(TESTS)

TESTS_COMMON_FILES = tests/test_common.c tests/test_common.h
TESTS_LDADD = src/libefp.la
TESTS_CPPFLAGS = $(AM_CPPFLAGS) -DABS_TOP_SRCDIR=\"$(abs_top_srcdir)\"

tests_elec_1a_SOURCES = tests/elec_1a.c $(TESTS_COMMON_FILES)
tests_elec_1a_LDADD = $(TESTS_LDADD)
tests_elec_1a_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_elec_1b_SOURCES = tests/elec_1b.c $(TESTS_COMMON_FILES)
tests_elec_1b_LDADD = $(TESTS_LDADD)
tests_elec_1b_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_pol_1_SOURCES = tests/pol_1.c $(TESTS_COMMON_FILES)
tests_pol_1_LDADD = $(TESTS_LDADD)
tests_pol_1_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_disp_1a_SOURCES = tests/disp_1a.c $(TESTS_COMMON_FILES)
tests_disp_1a_LDADD = $(TESTS_LDADD)
tests_disp_1a_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_disp_1b_SOURCES = tests/disp_1b.c $(TESTS_COMMON_FILES)
tests_disp_1b_LDADD = $(TESTS_LDADD)
tests_disp_1b_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_xr_1_SOURCES = tests/xr_1.c $(TESTS_COMMON_FILES)
tests_xr_1_LDADD = $(TESTS_LDADD)
tests_xr_1_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_total_1a_SOURCES = tests/total_1a.c $(TESTS_COMMON_FILES)
tests_total_1a_LDADD = $(TESTS_LDADD)
tests_total_1a_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_elec_2a_SOURCES = tests/elec_2a.c $(TESTS_COMMON_FILES)
tests_elec_2a_LDADD = $(TESTS_LDADD)
tests_elec_2a_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_elec_2b_SOURCES = tests/elec_2b.c $(TESTS_COMMON_FILES)
tests_elec_2b_LDADD = $(TESTS_LDADD)
tests_elec_2b_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_pol_2_SOURCES = tests/pol_2.c $(TESTS_COMMON_FILES)
tests_pol_2_LDADD = $(TESTS_LDADD)
tests_pol_2_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_disp_2a_SOURCES = tests/disp_2a.c $(TESTS_COMMON_FILES)
tests_disp_2a_LDADD = $(TESTS_LDADD)
tests_disp_2a_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_disp_2b_SOURCES = tests/disp_2b.c $(TESTS_COMMON_FILES)
tests_disp_2b_LDADD = $(TESTS_LDADD)
tests_disp_2b_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_xr_2_SOURCES = tests/xr_2.c $(TESTS_COMMON_FILES)
tests_xr_2_LDADD = $(TESTS_LDADD)
tests_xr_2_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_total_2a_SOURCES = tests/total_2a.c $(TESTS_COMMON_FILES)
tests_total_2a_LDADD = $(TESTS_LDADD)
tests_total_2a_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_elec_3a_SOURCES = tests/elec_3a.c $(TESTS_COMMON_FILES)
tests_elec_3a_LDADD = $(TESTS_LDADD)
tests_elec_3a_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_elec_3b_SOURCES = tests/elec_3b.c $(TESTS_COMMON_FILES)
tests_elec_3b_LDADD = $(TESTS_LDADD)
tests_elec_3b_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_pol_3_SOURCES = tests/pol_3.c $(TESTS_COMMON_FILES)
tests_pol_3_LDADD = $(TESTS_LDADD)
tests_pol_3_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_disp_3a_SOURCES = tests/disp_3a.c $(TESTS_COMMON_FILES)
tests_disp_3a_LDADD = $(TESTS_LDADD)
tests_disp_3a_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_disp_3b_SOURCES = tests/disp_3b.c $(TESTS_COMMON_FILES)
tests_disp_3b_LDADD = $(TESTS_LDADD)
tests_disp_3b_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_xr_3_SOURCES = tests/xr_3.c $(TESTS_COMMON_FILES)
tests_xr_3_LDADD = $(TESTS_LDADD)
tests_xr_3_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_total_3a_SOURCES = tests/total_3a.c $(TESTS_COMMON_FILES)
tests_total_3a_LDADD = $(TESTS_LDADD)
tests_total_3a_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_total_4a_SOURCES = tests/total_4a.c $(TESTS_COMMON_FILES)
tests_total_4a_LDADD = $(TESTS_LDADD)
tests_total_4a_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_total_4b_SOURCES = tests/total_4b.c $(TESTS_COMMON_FILES)
tests_total_4b_LDADD = $(TESTS_LDADD)
tests_total_4b_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_total_4c_SOURCES = tests/total_4c.c $(TESTS_COMMON_FILES)
tests_total_4c_LDADD = $(TESTS_LDADD)
tests_total_4c_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_qm_5_SOURCES = tests/qm_5.c $(TESTS_COMMON_FILES)
tests_qm_5_LDADD = $(TESTS_LDADD)
tests_qm_5_CPPFLAGS = $(TESTS_CPPFLAGS)

tests_qm_6_SOURCES = tests/qm_6.c $(TESTS_COMMON_FILES)
tests_qm_6_LDADD = $(TESTS_LDADD)
tests_qm_6_CPPFLAGS = $(TESTS_CPPFLAGS)
