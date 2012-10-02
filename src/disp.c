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

#include <assert.h>

#include <util.h>

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

static double
disp_tt(struct frag *fr_i, struct frag *fr_j, int pt_i_idx,
		int pt_j_idx, double sum, int do_grad)
{
	const struct dynamic_polarizable_pt *pt_i =
				fr_i->dynamic_polarizable_pts + pt_i_idx;
	const struct dynamic_polarizable_pt *pt_j =
				fr_j->dynamic_polarizable_pts + pt_j_idx;

	double r = vec_dist(CVEC(pt_i->x), CVEC(pt_j->x));
	double r2 = r * r;
	double r6 = r2 * r2 * r2;

	double damp = get_damp_tt(r);
	double energy = -4.0 / 3.0 * sum * damp / r6;

	if (do_grad) {
		double gdamp = get_damp_tt_grad(r);
		double g = 4.0 / 3.0 * sum * (gdamp / r - 6.0 * damp / r2) / r6;

		vec_t force = {
			g * (pt_j->x - pt_i->x),
			g * (pt_j->y - pt_i->y),
			g * (pt_j->z - pt_i->z)
		};

		add_force_torque(fr_i, fr_j, CVEC(pt_i->x), CVEC(pt_j->x), &force);
	}
	return energy;
}

static double
disp_overlap(struct frag *fr_i, struct frag *fr_j, int pt_i_idx,
		int pt_j_idx, int overlap_idx, double sum, int do_grad)
{
	const struct dynamic_polarizable_pt *pt_i =
				fr_i->dynamic_polarizable_pts + pt_i_idx;
	const struct dynamic_polarizable_pt *pt_j =
				fr_j->dynamic_polarizable_pts + pt_j_idx;

	vec_t dr = vec_sub(CVEC(pt_j->x), CVEC(pt_i->x));

	double r = vec_len(&dr);
	double r2 = r * r;
	double r6 = r2 * r2 * r2;

	double s_ij = fr_i->overlap_int[overlap_idx];
	double ln_s = 0.0;
	double damp = 1.0;

	if (fabs(s_ij) > 1.0e-5) {
		ln_s = log(fabs(s_ij));
		damp = 1.0 - s_ij * s_ij * (1.0 - 2.0 * ln_s + 2.0 * ln_s * ln_s);
	}

	double energy = -4.0 / 3.0 * sum * damp / r6;

	if (do_grad) {
		six_t ds_ij = fr_i->overlap_int_deriv[overlap_idx];
		vec_t force, torque_i, torque_j;

		double t1 = -8.0 * sum / r6 / r2 * damp;
		double t2 = -16.0 / 3.0 * sum / r6 * ln_s * ln_s * s_ij;

		vec_t dr_i = vec_sub(CVEC(pt_i->x), VEC(fr_i->x));
		vec_t dr_j = vec_sub(CVEC(pt_j->x), VEC(fr_j->x));
		vec_t dr_com = vec_sub(VEC(fr_j->x), VEC(fr_i->x));

		force.x = t1 * dr.x - t2 * ds_ij.x;
		force.y = t1 * dr.y - t2 * ds_ij.y;
		force.z = t1 * dr.z - t2 * ds_ij.z;

		torque_i.x = t1 * (dr.z * dr_i.y - dr.y * dr_i.z) + t2 * ds_ij.a;
		torque_i.y = t1 * (dr.x * dr_i.z - dr.z * dr_i.x) + t2 * ds_ij.b;
		torque_i.z = t1 * (dr.y * dr_i.x - dr.x * dr_i.y) + t2 * ds_ij.c;

		torque_j.x = t1 * (dr.z * dr_j.y - dr.y * dr_j.z) +
			     t2 * (ds_ij.z * dr_com.y - ds_ij.y * dr_com.z) +
			     t2 * ds_ij.a;
		torque_j.y = t1 * (dr.x * dr_j.z - dr.z * dr_j.x) +
			     t2 * (ds_ij.x * dr_com.z - ds_ij.z * dr_com.x) +
			     t2 * ds_ij.b;
		torque_j.z = t1 * (dr.y * dr_j.x - dr.x * dr_j.y) +
			     t2 * (ds_ij.y * dr_com.x - ds_ij.x * dr_com.y) +
			     t2 * ds_ij.c;

		vec_atomic_add(&fr_i->force, &force);
		vec_atomic_add(&fr_i->torque, &torque_i);

		vec_atomic_sub(&fr_j->force, &force);
		vec_atomic_sub(&fr_j->torque, &torque_j);
	}
	return energy;
}

static double
disp_off(struct frag *fr_i, struct frag *fr_j, int pt_i_idx,
		int pt_j_idx, double sum, int do_grad)
{
	const struct dynamic_polarizable_pt *pt_i =
				fr_i->dynamic_polarizable_pts + pt_i_idx;
	const struct dynamic_polarizable_pt *pt_j =
				fr_j->dynamic_polarizable_pts + pt_j_idx;

	double r2 = vec_dist_2(CVEC(pt_i->x), CVEC(pt_j->x));
	double r6 = r2 * r2 * r2;

	double energy = -4.0 / 3.0 * sum / r6;

	if (do_grad) {
		double r8 = r6 * r2;
		double g = -8.0 * sum / r8;

		vec_t force = {
			g * (pt_j->x - pt_i->x),
			g * (pt_j->y - pt_i->y),
			g * (pt_j->z - pt_i->z)
		};

		add_force_torque(fr_i, fr_j, CVEC(pt_i->x), CVEC(pt_j->x), &force);
	}
	return energy;
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

	switch (efp->opts.disp_damp) {
		case EFP_DISP_DAMP_TT:
			return disp_tt(fr_i, fr_j, pt_i_idx, pt_j_idx,
				       sum, efp->do_gradient);
		case EFP_DISP_DAMP_OVERLAP:
			return disp_overlap(fr_i, fr_j, pt_i_idx, pt_j_idx,
					    overlap_idx, sum, efp->do_gradient);
		case EFP_DISP_DAMP_OFF:
			return disp_off(fr_i, fr_j, pt_i_idx, pt_j_idx,
					sum, efp->do_gradient);
	}
	assert(0);
}

static double
frag_frag_disp(struct efp *efp, int frag_i, int frag_j, int overlap_idx)
{
	double energy = 0.0;

	int n_disp_i = efp->frags[frag_i].n_dynamic_polarizable_pts;
	int n_disp_j = efp->frags[frag_j].n_dynamic_polarizable_pts;

	for (int ii = 0, idx = overlap_idx; ii < n_disp_i; ii++)
		for (int jj = 0; jj < n_disp_j; jj++, idx++)
			energy += point_point_disp(efp, frag_i, frag_j, ii, jj, idx);

	return energy;
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

		move_pt(CVEC(frag->x), rotmat, CVEC(pt_in->x), VEC(pt_out->x));
	}
}
