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
damp_tt(double r)
{
	static const double a = 1.5; /* Tang-Toennies damping parameter */

	double ra[7];
	powers(r * a, 7, ra);

	return 1.0 - exp(-ra[1]) * (1.0 + ra[1] + ra[2] / 2.0 + ra[3] / 6.0 +
			ra[4] / 24.0 + ra[5] / 120.0 + ra[6] / 720.0);
}

static double
damp_tt_grad(double r)
{
	static const double a = 1.5; /* Tang-Toennies damping parameter */

	double ra = r * a;
	double ra2 = ra * ra;
	double ra6 = ra2 * ra2 * ra2;

	return a * exp(-ra) * ra6 / 720.0;
}

static double
damp_overlap(struct efp *efp, int frag_i, int frag_j, int pt_i, int pt_j)
{
	int idx = disp_damp_overlap_idx(efp, frag_i, frag_j, pt_i, pt_j);
	return efp->disp_damp_overlap[idx];
}

static double
damp_overlap_grad(struct efp *efp, int frag_i, int frag_j, int pt_i, int pt_j)
{
	/* XXX */
	assert(0);
}

static double
compute_disp_pt(struct efp *efp, int i, int ii, int j, int jj, double *grad)
{
	const struct frag *fr_i = efp->frags + i;
	const struct frag *fr_j = efp->frags + j;

	const struct dynamic_polarizable_pt *pt_i =
				fr_i->dynamic_polarizable_pts + ii;
	const struct dynamic_polarizable_pt *pt_j =
				fr_j->dynamic_polarizable_pts + jj;

	double sum = 0.0;
	for (int k = 0; k < ARRAY_SIZE(disp_weights); k++)
		sum += disp_weights[k] * pt_i->trace[k] * pt_j->trace[k];

	double r = vec_dist(VEC(pt_i->x), VEC(pt_j->x));
	double r2 = r * r, r6 = r2 * r2 * r2;

	double damp;
	switch (efp->opts.disp_damp) {
	case EFP_DISP_DAMP_TT:
		damp = damp_tt(r);
		break;
	case EFP_DISP_DAMP_OVERLAP:
		damp = damp_overlap(efp, i, j, ii, jj);
		break;
	}

	double energy = 4.0 / 3.0 * sum * damp / r6;

	if (grad) {
		double dx = pt_i->x - pt_j->x;
		double dy = pt_i->y - pt_j->y;
		double dz = pt_i->z - pt_j->z;

		double gdamp;
		switch (efp->opts.disp_damp) {
		case EFP_DISP_DAMP_TT:
			gdamp = damp_tt_grad(r);
			break;
		case EFP_DISP_DAMP_OVERLAP:
			gdamp = damp_overlap_grad(efp, i, j, ii, jj);
			break;
		}

		double g = 4.0 / 3.0 * sum * (gdamp / r - 6.0 * damp / r2) / r6;
		double gx = g * dx;
		double gy = g * dy;
		double gz = g * dz;

		double *g_i = grad + 6 * i;
		double *g_j = grad + 6 * j;

		g_i[0] += gx;
		g_i[1] += gy;
		g_i[2] += gz;
		g_i[3] += gz * (pt_i->y - fr_i->y) - gy * (pt_i->z - fr_i->z);
		g_i[4] += gx * (pt_i->z - fr_i->z) - gz * (pt_i->x - fr_i->x);
		g_i[5] += gy * (pt_i->x - fr_i->x) - gx * (pt_i->y - fr_i->y);

		g_j[0] -= gx;
		g_j[1] -= gy;
		g_j[2] -= gz;
		g_j[3] -= gz * (pt_j->y - fr_j->y) - gy * (pt_j->z - fr_j->z);
		g_j[4] -= gx * (pt_j->z - fr_j->z) - gz * (pt_j->x - fr_j->x);
		g_j[5] -= gy * (pt_j->x - fr_j->x) - gx * (pt_j->y - fr_j->y);
	}
	return energy;
}

static double
compute_disp_frag(struct efp *efp, int i, int j, double *grad)
{
	double sum = 0.0;

	int n_disp_i = efp->frags[i].n_dynamic_polarizable_pts;
	int n_disp_j = efp->frags[j].n_dynamic_polarizable_pts;

	for (int ii = 0; ii < n_disp_i; ii++)
		for (int jj = 0; jj < n_disp_j; jj++)
			sum += compute_disp_pt(efp, i, ii, j, jj, grad);

	return sum;
}

enum efp_result
efp_compute_disp(struct efp *efp)
{
	double energy = 0.0;

	for (int i = 0; i < efp->n_frag; i++)
		for (int j = i + 1; j < efp->n_frag; j++)
			energy -= compute_disp_frag(efp, i, j, efp->grad);

	efp->energy[efp_get_term_index(EFP_TERM_DISP)] = energy;
	return EFP_RESULT_SUCCESS;
}

void
efp_update_disp(struct frag *frag, const struct mat *rotmat)
{
	for (int i = 0; i < frag->n_dynamic_polarizable_pts; i++)
		move_pt(VEC(frag->x), rotmat,
			VEC(frag->lib->dynamic_polarizable_pts[i].x),
			VEC(frag->dynamic_polarizable_pts[i].x));
}
