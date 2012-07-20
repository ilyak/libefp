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

#include "efp_private.h"
#include "disp.h"

static double
get_damp_tt(double r)
{
	static const double a = 1.5; /* Tang-Toennies damping parameter */

	double ra = r * a;
	double ra2 = ra * ra;
	double ra3 = ra2 * ra;
	double ra4 = ra3 * ra;
	double ra5 = ra4 * ra;
	double ra6 = ra5 * ra;

	return 1.0 - exp(-ra) * (1.0 + ra + ra2 / 2.0 + ra3 / 6.0 +
			ra4 / 24.0 + ra5 / 120.0 + ra6 / 720.0);
}

static double
get_damp_tt_grad(double r)
{
	static const double a = 1.5; /* Tang-Toennies damping parameter */

	double ra = r * a;
	double ra2 = ra * ra;
	double ra6 = ra2 * ra2 * ra2;

	return a * exp(-ra) * ra6 / 720.0;
}

static inline double
get_damp_overlap(double s_ij)
{
	if (fabs(s_ij) < 1.0e-6)
		return 0.0;

	double ln_s = log(fabs(s_ij));

	return 1.0 - s_ij * s_ij * (1.0 - 2.0 * ln_s + 2.0 * ln_s * ln_s);
}

static double
point_point_disp(struct efp *efp, int fr_i_idx, int fr_j_idx,
		 int pt_i_idx, int pt_j_idx, int overlap_idx)
{
	struct frag *fr_i = efp->frags + fr_i_idx;
	struct frag *fr_j = efp->frags + fr_j_idx;

	const struct dynamic_polarizable_pt *pt_i =
				fr_i->dynamic_polarizable_pts + pt_i_idx;
	const struct dynamic_polarizable_pt *pt_j =
				fr_j->dynamic_polarizable_pts + pt_j_idx;

	double sum = 0.0;

	for (size_t k = 0; k < ARRAY_SIZE(disp_weights); k++)
		sum += disp_weights[k] * pt_i->trace[k] * pt_j->trace[k];

	double r = vec_dist(VEC(pt_i->x), VEC(pt_j->x));
	double r2 = r * r;
	double r6 = r2 * r2 * r2;

	double damp = 1.0;

	switch (efp->opts.disp_damp) {
	case EFP_DISP_DAMP_TT:
		damp = get_damp_tt(r);
		break;
	case EFP_DISP_DAMP_OVERLAP:
		damp = get_damp_overlap(fr_i->overlap_int[overlap_idx]);
		break;
	case EFP_DISP_DAMP_OFF:
		break;
	}

	double energy = -4.0 / 3.0 * sum * damp / r6;

	if (efp->do_gradient) {
		double gdamp = 0.0;

		switch (efp->opts.disp_damp) {
		case EFP_DISP_DAMP_TT:
			gdamp = get_damp_tt_grad(r);
			break;
		case EFP_DISP_DAMP_OVERLAP:
			/* XXX */
			break;
		case EFP_DISP_DAMP_OFF:
			break;
		}

		double g = 4.0 / 3.0 * sum * (gdamp / r - 6.0 * damp / r2) / r6;

		vec_t force = {
			g * (pt_j->x - pt_i->x),
			g * (pt_j->y - pt_i->y),
			g * (pt_j->z - pt_i->z)
		};

		add_force_torque(fr_i, fr_j, VEC(pt_i->x), VEC(pt_j->x),
				 &force);
	}
	return energy;
}

static double
frag_frag_disp(struct efp *efp, int frag_i, int frag_j, int overlap_idx)
{
	double sum = 0.0;

	int n_disp_i = efp->frags[frag_i].n_dynamic_polarizable_pts;
	int n_disp_j = efp->frags[frag_j].n_dynamic_polarizable_pts;

	for (int ii = 0, idx = overlap_idx; ii < n_disp_i; ii++)
		for (int jj = 0; jj < n_disp_j; jj++, idx++)
			sum += point_point_disp(efp, frag_i, frag_j, ii, jj, idx);

	return sum;
}

enum efp_result
efp_compute_disp(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_DISP))
		return EFP_RESULT_SUCCESS;

	double energy = 0.0;

	#pragma omp parallel for schedule(dynamic, 4) reduction(+:energy)
	for (int i = 0; i < efp->n_frag; i++) {
		for (int j = i + 1, idx = 0; j < efp->n_frag; j++) {
			energy += frag_frag_disp(efp, i, j, idx);
			idx += efp->frags[i].n_lmo * efp->frags[j].n_lmo;
		}
	}

	efp->energy.dispersion = energy;
	return EFP_RESULT_SUCCESS;
}

void
efp_update_disp(struct frag *frag)
{
	const mat_t *rotmat = &frag->rotmat;

	for (int i = 0; i < frag->n_dynamic_polarizable_pts; i++) {
		const struct dynamic_polarizable_pt *pt_in =
					frag->lib->dynamic_polarizable_pts + i;
		struct dynamic_polarizable_pt *pt_out =
					frag->dynamic_polarizable_pts + i;

		move_pt(VEC(frag->x), rotmat, VEC(frag->lib->x),
			VEC(pt_in->x), VEC(pt_out->x));
	}
}
