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
	0.0, 0.0, 0.0, 1.0, 2.0, 3.0,
	6.0, 0.0, 0.0, 5.0, 2.0, 8.0
};

static const double ref_gradient[] = {
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

static const struct test_data test_data = {
	.potential_files = potential_files,
	.fragname = fragname,
	.xyzabc = xyzabc,
	.ref_energy = 0.0,
	.ref_gradient = ref_gradient,
	.opts = {
		.terms = EFP_TERM_ELEC,
		.do_gradient = 0
	}
};

DEFINE_TEST(test_data)
