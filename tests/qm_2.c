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

static const char files[] =
	ABS_TOP_SRCDIR "/fraglib/ch3oh.efp\n"
	ABS_TOP_SRCDIR "/fraglib/dmso.efp\n"
	ABS_TOP_SRCDIR "/fraglib/dcm.efp\n"
	ABS_TOP_SRCDIR "/fraglib/acetone.efp";

static const char names[] =
	"CH3OH_L\n"
	"DMSO_L\n"
	"DMSO_L\n"
	"ACETONE_L\n"
	"DCM_L\n"
	"ACETONE_L\n"
	"ACETONE_L";

static const double frag_xyzabc[] = {
	BOHR( 0.0), BOHR(-1.0), BOHR( 0.0),  0.0,  1.1,  2.0,
	BOHR(-5.0), BOHR(12.0), BOHR(-0.0),  3.0,  0.2,  5.0,
	BOHR(-0.0), BOHR(-3.0), BOHR( 5.0),  6.0,  2.3,  8.0,
	BOHR(-5.0), BOHR(-4.0), BOHR(-5.0),  9.0,  0.4,  1.0,
	BOHR(-9.0), BOHR(-5.0), BOHR( 1.0),  2.0,  1.5,  4.0,
	BOHR( 7.0), BOHR(-2.0), BOHR(11.0),  5.0,  0.6,  7.0,
	BOHR(-9.0), BOHR(-7.0), BOHR(-9.0),  8.0,  2.7,  0.0
};

static const double qm_znuc[] = {
	7.0, 10.0, 1.0, 2.0, 2.0
};

static const double qm_xyz[] = {
	BOHR( 4.0), BOHR( 5.0), BOHR( 5.0),
	BOHR( 8.0), BOHR( 5.0), BOHR( 6.0),
	BOHR( 5.0), BOHR( 8.0), BOHR( 5.0),
	BOHR( 8.0), BOHR( 9.0), BOHR( 8.0),
	BOHR( 5.0), BOHR( 8.0), BOHR( 8.0)
};

static enum efp_result get_electron_density_field(int n_pt,
		UNUSED const double *xyz, double *field,
		UNUSED void *user_data)
{
	/* no electrons */
	memset(field, 0, n_pt * 3 * sizeof(double));
	return EFP_RESULT_SUCCESS;
}

static const struct test_data test_data = {
	.files = files,
	.names = names,
	.geometry_xyzabc = frag_xyzabc,
	.n_ptc = ARRAY_SIZE(qm_znuc),
	.ptc_charges = qm_znuc,
	.ptc_xyz = qm_xyz,
	.ref_energy = -0.2314262632,
	.electron_density_field_fn = get_electron_density_field,
	.electron_density_field_user_data = NULL,
	.opts = {
		.terms = EFP_TERM_ELEC | EFP_TERM_POL |
			 EFP_TERM_DISP | EFP_TERM_XR |
			 EFP_TERM_AI_ELEC | EFP_TERM_AI_POL,
		.elec_damp = EFP_ELEC_DAMP_SCREEN,
		.disp_damp = EFP_DISP_DAMP_OVERLAP
	}
};

DEFINE_TEST(qm_2)
