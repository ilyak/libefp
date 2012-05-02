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

static const double ref_gradient[] = {
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

static enum efp_result
overlap_integrals_fn(struct efp_xr_block *block, double *s, double *sx,
		     void *user_data)
{
	static const double s_[] = {
0.0, 0.0, 0.0
	};

	static const double sx_[] = {
0.0, 0.0, 0.0
	};

	size_t size = block->basis_size_i * block->basis_size_j;

	if (size != ARRAY_SIZE(s_))
		return EFP_RESULT_INVALID_ARRAY_SIZE;

	memcpy(s, s_, size * sizeof(double));

	if (sx)
		memcpy(sx, sx_, 3 * size * sizeof(double));

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
kinetic_integrals_fn(struct efp_xr_block *block, double *t, double *tx,
		     void *user_data)
{
	static const double t_[] = {
0.0, 0.0, 0.0
	};

	static const double tx_[] = {
0.0, 0.0, 0.0
	};

	size_t size = block->basis_size_i * block->basis_size_j;

	if (size != ARRAY_SIZE(t_))
		return EFP_RESULT_INVALID_ARRAY_SIZE;

	memcpy(t, t_, size * sizeof(double));

	if (tx)
		memcpy(tx, tx_, 3 * size * sizeof(double));

	return EFP_RESULT_SUCCESS;
}

static const struct test_data test_data = {
	.potential_files = potential_files,
	.fragname = fragname,
	.xyzabc = xyzabc,
	.ref_energy = 0.0,
	.ref_gradient = ref_gradient,
	.opts = {
		.terms = EFP_TERM_XR,
		.do_gradient = 0
	},
	.callbacks = {
		.get_overlap_integrals = overlap_integrals_fn,
		.get_kinetic_integrals = kinetic_integrals_fn
	}
};

DEFINE_TEST(test_data)
