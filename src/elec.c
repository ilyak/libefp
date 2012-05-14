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

#include <string.h>

#include "efp_private.h"
#include "elec.h"

#define DR_A ((double *)&dr)[a]
#define DR_B ((double *)&dr)[b]
#define DR_C ((double *)&dr)[c]

#define D1_A   ((double *)&pt_i->dipole)[a]
#define D2_A   ((double *)&pt_j->dipole)[a]
#define D1_B   ((double *)&pt_i->dipole)[b]
#define D2_B   ((double *)&pt_j->dipole)[b]

#define Q1_AB  (pt_i->quadrupole[quad_idx(a, b)])
#define Q2_AB  (pt_j->quadrupole[quad_idx(a, b)])

#define O1_ABC (pt_i->octupole[oct_idx(a, b, c)])
#define O2_ABC (pt_j->octupole[oct_idx(a, b, c)])

#define SUM_A(x) (a = 0, sum_a = x, a = 1, sum_a += x, a = 2, sum_a += x)
#define SUM_B(x) (b = 0, sum_b = x, b = 1, sum_b += x, b = 2, sum_b += x)
#define SUM_C(x) (c = 0, sum_c = x, c = 1, sum_c += x, c = 2, sum_c += x)

#define SUM_AB(x)  (SUM_A(SUM_B(x)))
#define SUM_ABC(x) (SUM_A(SUM_B(SUM_C(x))))

static inline double
get_damp_screen(struct efp *efp, double r_ij, double pi, double pj)
{
	if (pj == HUGE_VAL) {   /* damping for nucleus */
		return 1.0 - exp(-pi * r_ij);
	}
	else if (fabs((pi - pj) * r_ij) < 1.0e-5) {
		return 1.0 - (1.0 + 0.5 * pi * r_ij) * exp(-pi * r_ij);
	}
	else {
		return 1.0 - exp(-pi * r_ij) * pj * pj / (pj * pj - pi * pi) -
			     exp(-pj * r_ij) * pi * pi / (pi * pi - pj * pj);
	}
}

static double
charge_mult_energy(struct efp *efp, double charge, const vec_t *pos,
		   int frag_idx, int pt_idx)
{
	struct frag *fr_i = efp->frags + frag_idx;
	struct multipole_pt *pt_i = fr_i->multipole_pts + pt_idx;

	double energy = 0.0;

	vec_t dr = vec_sub(VEC(pt_i->x), pos);
	double r = vec_len(&dr);

	double ri[8];
	powers(1.0 / r, 8, ri);

	int a, b, c;
	double sum_a, sum_b, sum_c;

	double damp = 1.0;
	if (efp->opts.elec_damp == EFP_ELEC_DAMP_SCREEN) {
		double scr_i = fr_i->screen_params[pt_idx];
		damp = get_damp_screen(efp, r, scr_i, HUGE_VAL);
	}

	energy += ri[1] * charge * pt_i->monopole * damp;
	energy -= ri[3] * charge * SUM_A(D1_A * DR_A);
	energy += ri[5] * charge * SUM_AB(Q1_AB * DR_A * DR_B);
	energy -= ri[7] * charge * SUM_ABC(O1_ABC * DR_A * DR_B * DR_C);

	return energy;
}

static double
compute_elec_pt(struct efp *efp, int i, int j, int ii, int jj)
{
	double energy = 0.0;

	struct frag *fr_i = efp->frags + i;
	struct frag *fr_j = efp->frags + j;

	struct multipole_pt *pt_i = fr_i->multipole_pts + ii;
	struct multipole_pt *pt_j = fr_j->multipole_pts + jj;

	vec_t dr = vec_sub(VEC(pt_j->x), VEC(pt_i->x));
	double r = vec_len(&dr);

	double ri[10];
	powers(1.0 / r, 10, ri);

	int a, b, c;
	double sum_a, sum_b, sum_c;
	double t1, t2, t3, t4;

	double q1 = pt_i->monopole;
	double q2 = pt_j->monopole;

	/* monopole - monopole */
	double damp = 1.0;
	if (efp->opts.elec_damp == EFP_ELEC_DAMP_SCREEN) {
		double scr_i = fr_i->screen_params[ii];
		double scr_j = fr_j->screen_params[jj];
		damp = get_damp_screen(efp, r, scr_i, scr_j);
	}
	energy += ri[1] * q1 * q2 * damp;

	/* monopole - dipole */
	energy += ri[3] * q2 * SUM_A(D1_A * DR_A);
	energy -= ri[3] * q1 * SUM_A(D2_A * DR_A);

	/* monopole - quadrupole */
	energy += ri[5] * q2 * SUM_AB(Q1_AB * DR_A * DR_B);
	energy += ri[5] * q1 * SUM_AB(Q2_AB * DR_A * DR_B);

	/* monopole - octupole */
	energy += ri[7] * q2 * SUM_ABC(O1_ABC * DR_A * DR_B * DR_C);
	energy -= ri[7] * q1 * SUM_ABC(O2_ABC * DR_A * DR_B * DR_C);

	/* dipole - dipole */
	energy += ri[3] * SUM_A(D1_A * D2_A);
	energy -= 3 * ri[5] * SUM_AB(D1_A * D2_B * DR_A * DR_B);

	/* dipole - quadrupole */
	energy += 2 * ri[5] * SUM_AB(Q1_AB * D2_A * DR_B);
	energy -= 2 * ri[5] * SUM_AB(Q2_AB * D1_A * DR_B);

	t1 = SUM_A(D2_A * DR_A);
	energy -= 5 * ri[7] * SUM_AB(Q1_AB * DR_A * DR_B) * t1;

	t1 = SUM_A(D1_A * DR_A);
	energy += 5 * ri[7] * SUM_AB(Q2_AB * DR_A * DR_B) * t1;

	/* quadrupole - quadrupole */
	t1 = 0.0;
	for (b = 0; b < 3; b++) {
		t2 = SUM_A(Q1_AB * DR_A);
		t3 = SUM_A(Q2_AB * DR_A);
		t1 += t2 * t3;
	}
	t2 = SUM_AB(Q1_AB * DR_A * DR_B);
	t3 = SUM_AB(Q2_AB * DR_A * DR_B);
	t4 = SUM_AB(Q1_AB * Q2_AB);

	energy -= 20 * ri[7] * t1 / 3.0;
	energy += 35 * ri[9] * t2 * t3 / 3.0;
	energy +=  2 * ri[5] * t4 / 3.0;

	/* gradient */
	if (efp->grad) {
	}
	return energy;
}

static double
compute_elec_frag(struct efp *efp, int fr_i_idx, int fr_j_idx)
{
	double energy = 0.0;

	struct frag *fr_i = efp->frags + fr_i_idx;
	struct frag *fr_j = efp->frags + fr_j_idx;

	/* nuclei - nuclei */
	for (int ii = 0; ii < fr_i->n_atoms; ii++) {
		for (int jj = 0; jj < fr_j->n_atoms; jj++) {
			struct efp_atom *at_i = fr_i->atoms + ii;
			struct efp_atom *at_j = fr_j->atoms + jj;

			double r = vec_dist(VEC(at_i->x), VEC(at_j->x));
			energy += at_i->znuc * at_j->znuc / r;
		}
	}

	/* nuclei - mult points */
	for (int ii = 0; ii < fr_i->n_atoms; ii++) {
		for (int jj = 0; jj < fr_j->n_multipole_pts; jj++) {
			struct efp_atom *at = fr_i->atoms + ii;

			energy += charge_mult_energy(efp, at->znuc, VEC(at->x),
						     fr_j_idx, jj);
		}
	}

	/* mult points - nuclei */
	for (int jj = 0; jj < fr_j->n_atoms; jj++) {
		for (int ii = 0; ii < fr_i->n_multipole_pts; ii++) {
			struct efp_atom *at = fr_j->atoms + jj;

			energy += charge_mult_energy(efp, at->znuc, VEC(at->x),
						     fr_i_idx, ii);
		}
	}

	/* mult points - mult points */
	for (int ii = 0; ii < fr_i->n_multipole_pts; ii++)
		for (int jj = 0; jj < fr_j->n_multipole_pts; jj++)
			energy += compute_elec_pt(efp, fr_i_idx, fr_j_idx,
					ii, jj);

	return energy;
}

enum efp_result
efp_compute_elec(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_ELEC))
		return EFP_RESULT_SUCCESS;

	if (efp->grad)
		return EFP_RESULT_NOT_IMPLEMENTED;

	double energy = 0.0;

	for (int i = 0; i < efp->n_frag; i++)
		for (int j = i + 1; j < efp->n_frag; j++)
			energy += compute_elec_frag(efp, i, j);

	efp->energy.electrostatic = energy;
	return EFP_RESULT_SUCCESS;
}

static void
rotate_quad(const mat_t *rotmat, const double *in, double *out)
{
	double full_in[9], full_out[9];

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			full_in[a * 3 + b] = in[quad_idx(a, b)];

	rotate_t2(rotmat, full_in, full_out);

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			out[quad_idx(a, b)] = full_out[a * 3 + b];
}

static void
rotate_oct(const mat_t *rotmat, const double *in, double *out)
{
	double full_in[27], full_out[27];

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			for (int c = 0; c < 3; c++) {
				int idx = 9 * a + 3 * b + c;
				full_in[idx] = in[oct_idx(a, b, c)];
			}

	rotate_t3(rotmat, full_in, full_out);

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			for (int c = 0; c < 3; c++) {
				int idx = 9 * a + 3 * b + c;
				out[oct_idx(a, b, c)] = full_out[idx];
			}
}

void
efp_update_elec(struct frag *frag, const mat_t *rotmat)
{
	for (int i = 0; i < frag->n_multipole_pts; i++) {
		/* move point position */
		move_pt(VEC(frag->x), rotmat,
			VEC(frag->lib->multipole_pts[i].x),
			VEC(frag->multipole_pts[i].x));

		/* rotate dipole */
		mat_vec(rotmat, &frag->lib->multipole_pts[i].dipole,
				&frag->multipole_pts[i].dipole);

		/* rotate quadrupole */
		rotate_quad(rotmat, frag->lib->multipole_pts[i].quadrupole,
				frag->multipole_pts[i].quadrupole);

		/* correction for Buckingham quadrupoles */
		double *quad = frag->multipole_pts[i].quadrupole;

		double qtr = quad[quad_idx(0, 0)] +
			     quad[quad_idx(1, 1)] +
			     quad[quad_idx(2, 2)];

		quad[0] = 1.5 * quad[0] - 0.5 * qtr;
		quad[1] = 1.5 * quad[1] - 0.5 * qtr;
		quad[2] = 1.5 * quad[2] - 0.5 * qtr;
		quad[3] = 1.5 * quad[3];
		quad[4] = 1.5 * quad[4];
		quad[5] = 1.5 * quad[5];

		/* rotate octupole */
		rotate_oct(rotmat, frag->lib->multipole_pts[i].octupole,
				frag->multipole_pts[i].octupole);

		/* correction for Buckingham octupoles */
		double *oct = frag->multipole_pts[i].octupole;

		double otrx = oct[oct_idx(0, 0, 0)] +
			      oct[oct_idx(0, 1, 1)] +
			      oct[oct_idx(0, 2, 2)];
		double otry = oct[oct_idx(0, 0, 1)] +
			      oct[oct_idx(1, 1, 1)] +
			      oct[oct_idx(1, 2, 2)];
		double otrz = oct[oct_idx(0, 0, 2)] +
			      oct[oct_idx(1, 1, 2)] +
			      oct[oct_idx(2, 2, 2)];

		oct[0] = 2.5 * oct[0] - 1.5 * otrx;
		oct[1] = 2.5 * oct[1] - 1.5 * otry;
		oct[2] = 2.5 * oct[2] - 1.5 * otrz;
		oct[3] = 2.5 * oct[3] - 0.5 * otry;
		oct[4] = 2.5 * oct[4] - 0.5 * otrz;
		oct[5] = 2.5 * oct[5] - 0.5 * otrx;
		oct[6] = 2.5 * oct[6] - 0.5 * otrz;
		oct[7] = 2.5 * oct[7] - 0.5 * otrx;
		oct[8] = 2.5 * oct[8] - 0.5 * otry;
		oct[9] = 2.5 * oct[9];
	}
}

static double
compute_ai_frag(struct efp *efp, int frag_idx)
{
	double energy = 0.0;
	struct frag *fr_i = efp->frags + frag_idx;

	for (int i = 0; i < fr_i->n_atoms; i++) {
		for (int j = 0; j < efp->qm_data.n_atoms; j++) {
			struct efp_atom *at_i = fr_i->atoms + i;
			struct efp_qm_atom *at_j = efp->qm_data.atoms + j;

			double r = vec_dist(VEC(at_i->x), VEC(at_j->x));
			energy += at_i->znuc * at_j->znuc / r;
		}
	}
	for (int i = 0; i < fr_i->n_multipole_pts; i++) {
		for (int j = 0; j < efp->qm_data.n_atoms; j++) {
			struct efp_qm_atom *at = efp->qm_data.atoms + j;

			energy += charge_mult_energy(efp, at->znuc, VEC(at->x),
						     frag_idx, i);
		}
	}
	return energy;
}

enum efp_result
efp_compute_ai_elec(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_AI_ELEC))
		return EFP_RESULT_SUCCESS;

	if (efp->grad)
		return EFP_RESULT_NOT_IMPLEMENTED;

	double energy = 0.0;

	for (int i = 0; i < efp->n_frag; i++)
		energy += compute_ai_frag(efp, i);

	efp->energy.ai_electrostatic = energy;
	return EFP_RESULT_SUCCESS;
}
