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
#include "geometry_2.h"

static enum efp_result
st_integrals_fn(struct efp_st_block *block, int compute_derivatives,
		struct efp_st_data *st, void *user_data)
{
	static const int basis_size = 345;

	if (compute_derivatives)
		return EFP_RESULT_NOT_IMPLEMENTED;

	if (block->basis_size_i != basis_size ||
	    block->basis_size_j != basis_size)
		return EFP_RESULT_INVALID_ARRAY_SIZE;

	FILE *fp;
	double *ptr;

	int size = block->basis_size_i * block->basis_size_j;

	fp = fopen(ABS_TOP_SRCDIR "/tests/data/sint_2", "r");
	if (!fp)
		return EFP_RESULT_FILE_NOT_FOUND;

	ptr = st->s;
	for (int i = 0; i < size; i++, ptr++)
		fscanf(fp, "%lf", ptr);

	fclose(fp);

	fp = fopen(ABS_TOP_SRCDIR "/tests/data/tint_2", "r");
	if (!fp)
		return EFP_RESULT_FILE_NOT_FOUND;

	ptr = st->t;
	for (int i = 0; i < size; i++, ptr++)
		fscanf(fp, "%lf", ptr);

	fclose(fp);
	return EFP_RESULT_SUCCESS;
}

static const double ref_gradient[] = { /* from Q-Chem 4.0 */
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

static const struct test_data test_data = {
	.potential_files = potential_files,
	.fragname = fragname,
	.xyzabc = xyzabc,
		/* elec + pol + disp + xr from Q-Chem 4.0 */
	.ref_energy =  0.001371996347 + -0.000190204501 +
		      -0.001468808757 +  0.000844260733,
	.ref_gradient = ref_gradient,
	.opts = {
		.terms = EFP_TERM_ELEC | EFP_TERM_POL |
			 EFP_TERM_DISP | EFP_TERM_XR,
		.disp_damp = EFP_DISP_DAMP_TT,
		.do_gradient = 0
	},
	.callbacks = {
		.get_st_integrals = st_integrals_fn
	}
};

DEFINE_TEST(test_data)
