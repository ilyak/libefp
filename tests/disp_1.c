/*-
 * Copyright (c) 2012 Ilya Kaliman
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */

#include "test_common.h"

static const char *potential_files[] = {
	ABS_TOP_SRCDIR "/fraglib/h2o.efp",
	ABS_TOP_SRCDIR "/fraglib/nh3.efp",
	NULL
};

static const char *fragname[] = {
	"H2O_L",
	"NH3_L",
	 NULL
};

static const double xyzabc[] = {
	BOHR(0.0), BOHR(0.0), BOHR(0.0), 1.0, 2.0, 3.0,
	BOHR(5.0), BOHR(0.0), BOHR(0.0), 5.0, 2.0, 8.0
};

static const double ref_gradient[] = { /* from Q-Chem 4.0 */
	 6.2469370978534026e-05, -2.1601463511397569e-07,
	 4.0492586440871674e-07,  4.5003170183999674e-08,
	 2.5764689731580998e-06, -4.1445973600957137e-06,
	-6.2469370978534026e-05,  2.1601463511397569e-07,
	-4.0492586440871674e-07, -4.5003346926898037e-08,
	-6.4025064365653540e-06,  2.1035594363279749e-06
};

static const struct test_data test_data = {
	.potential_files = potential_files,
	.fragname = fragname,
	.xyzabc = xyzabc,
	.ref_energy = -0.000098903256, /* from Q-Chem 4.0 */
	.ref_gradient = ref_gradient,
	.opts = {
		.terms = EFP_TERM_DISP,
		.do_gradient = 1,
		.disp_damp = EFP_DISP_DAMP_TT
	}
};

DEFINE_TEST(test_data)
