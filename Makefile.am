EXTRA_DIST =
CLEANFILES =
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}
AM_MAKEFLAGS = --no-print-directory
AM_LDFLAGS =

AM_CPPFLAGS = -include $(top_builddir)/config.h

AM_CFLAGS = ${my_CFLAGS} ${OPENMP_CFLAGS} \
	-fvisibility=hidden

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

noinst_LTLIBRARIES = common/libcommon.la

common_libcommon_la_SOURCES = \
	common/cblas.h \
	common/compat.c \
	common/compat.h \
	common/math_util.h \
	common/phys_const.h \
	common/util.c \
	common/util.h

lib_LTLIBRARIES = src/libefp.la
include_HEADERS = src/efp.h

src_libefp_la_SOURCES = \
	src/disp.c \
	src/disp.h \
	src/efp.c \
	src/efp_private.h \
	src/efp_util.c \
	src/elec.c \
	src/elec.h \
	src/elec_grad.c \
	src/int.c \
	src/int.h \
	src/int_shift.h \
	src/parse.c \
	src/pol.c \
	src/swf.c \
	src/swf.h \
	src/terms.h \
	src/xr.c

src_libefp_la_LIBADD = common/libcommon.la
src_libefp_la_LDFLAGS = $(AM_LDFLAGS)

EXTRA_DIST += AUTHORS LICENSE README.md VERSION

if BUILD_EFPMD
bin_PROGRAMS = efpmd/src/efpmd
endif

efpmd_src_efpmd_SOURCES = \
	efpmd/src/cfg.h \
	efpmd/src/clapack.c \
	efpmd/src/clapack.h \
	efpmd/src/common.c \
	efpmd/src/common.h \
	efpmd/src/grad.c \
	efpmd/src/hess.c \
	efpmd/src/keywords.h \
	efpmd/src/main.c \
	efpmd/src/md.c \
	efpmd/src/opt.c \
	efpmd/src/parse.c \
	efpmd/src/sp.c \
	efpmd/optimizer/blas.f \
	efpmd/optimizer/lbfgsb.f \
	efpmd/optimizer/linpack.f \
	efpmd/optimizer/timer.f \
	efpmd/optimizer/optimizer.c \
	efpmd/optimizer/optimizer.h

efpmd_src_efpmd_LDADD = common/libcommon.la src/libefp.la
efpmd_src_efpmd_CFLAGS = $(AM_CFLAGS)
efpmd_src_efpmd_CPPFLAGS = $(AM_CPPFLAGS) \
	-DEFP_DATA_DIR=\""$(datarootdir)/libefp"\"

EXTRA_DIST += \
	efpmd/README.md \
	efpmd/optimizer/License.txt \
	efpmd/scripts/cubegen.pl \
	efpmd/scripts/trajectory.pl \
	efpmd/examples/c6h6-dimer_opt.in \
	efpmd/examples/h2o-ch3oh_hess.in \
	efpmd/examples/h2o-dimer_opt.in \
	efpmd/examples/h2o-nh3_grad.in \
	efpmd/examples/h2o-tetramer_opt.in \
	efpmd/examples/water-cube-125_grad.in \
	efpmd/examples/water-cube-125_sp.in \
	efpmd/examples/water-cube-3375_sp.in

if HAVE_LIBCHECK
TESTS = tests/test
check_PROGRAMS = $(TESTS)
endif

tests_test_LDADD = common/libcommon.la src/libefp.la -lcheck
tests_test_CPPFLAGS = $(AM_CPPFLAGS) \
	-DABS_TOP_SRCDIR=\""$(abs_top_srcdir)"\"

tests_test_SOURCES = \
	tests/test.c \
	tests/test_common.h \
	tests/test_list.h \
	tests/geometry_1.h \
	tests/geometry_2.h \
	tests/geometry_3.h \
	tests/geometry_4.h \
	tests/geometry_5.h \
	tests/geometry_6.h \
	tests/disp_1a.c \
	tests/disp_1b.c \
	tests/disp_1c.c \
	tests/disp_1d.c \
	tests/disp_2a.c \
	tests/disp_2b.c \
	tests/disp_3a.c \
	tests/disp_3b.c \
	tests/elec_1a.c \
	tests/elec_1b.c \
	tests/elec_1c.c \
	tests/elec_2a.c \
	tests/elec_2b.c \
	tests/elec_3a.c \
	tests/elec_3b.c \
	tests/pol_1.c \
	tests/pol_2.c \
	tests/pol_3.c \
	tests/qm_1.c \
	tests/qm_2.c \
	tests/total_1a.c \
	tests/total_2a.c \
	tests/total_3a.c \
	tests/total_4a.c \
	tests/total_4b.c \
	tests/total_4c.c \
	tests/total_5a.c \
	tests/total_5b.c \
	tests/total_5c.c \
	tests/total_6a.c \
	tests/total_6b.c \
	tests/total_6c.c \
	tests/xr_1a.c \
	tests/xr_1b.c \
	tests/xr_2.c \
	tests/xr_3.c
