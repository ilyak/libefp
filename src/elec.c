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
#include "elec.h"

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

	vec_t dr = vec_sub(VEC(pt_i->x), pos);

	double r = vec_len(&dr);
	double r3 = r * r * r;
	double r5 = r3 * r * r;
	double r7 = r5 * r * r;

	double damp = 1.0;
	if (efp->opts.elec_damp == EFP_ELEC_DAMP_SCREEN) {
		double scr_i = fr_i->screen_params[pt_idx];
		damp = get_damp_screen(efp, r, scr_i, HUGE_VAL);
	}

	double energy = 0.0;

	energy += charge * pt_i->monopole * damp / r;
	energy -= charge * vec_dot(&pt_i->dipole, &dr) / r3;
	energy += charge * quadrupole_sum(pt_i->quadrupole, &dr) / r5;
	energy -= charge * octupole_sum(pt_i->octupole, &dr) / r7;

	return energy;
}

static void
mult_mult_grad(struct efp *efp, int fr_i_idx, int fr_j_idx,
	       int pt_i_idx, int pt_j_idx)
{
	/* monopole - monopole */

	/* monopole - dipole */

	/* dipole - monopole */

	/* monopole - quadrupole */

	/* quadrupole - monopole */

	/* monopole - octupole */

	/* octupole - monopole */

	/* dipole - dipole */

	/* dipole - quadrupole */

	/* quadrupole - dipole */

	/* quadrupole - quadrupole */
}

static double
mult_mult_energy(struct efp *efp, int fr_i_idx, int fr_j_idx,
		 int pt_i_idx, int pt_j_idx)
{
	double energy = 0.0;

	struct frag *fr_i = efp->frags + fr_i_idx;
	struct frag *fr_j = efp->frags + fr_j_idx;

	struct multipole_pt *pt_i = fr_i->multipole_pts + pt_i_idx;
	struct multipole_pt *pt_j = fr_j->multipole_pts + pt_j_idx;

	vec_t dr = vec_sub(VEC(pt_j->x), VEC(pt_i->x));

	double r = vec_len(&dr);
	double r3 = r * r * r;
	double r5 = r3 * r * r;
	double r7 = r5 * r * r;
	double r9 = r7 * r * r;

	double q1 = pt_i->monopole;
	double q2 = pt_j->monopole;
	const vec_t *d1 = &pt_i->dipole;
	const vec_t *d2 = &pt_j->dipole;
	const double *quad1 = pt_i->quadrupole;
	const double *quad2 = pt_j->quadrupole;
	const double *oct1 = pt_i->octupole;
	const double *oct2 = pt_j->octupole;

	/* monopole - monopole */
	double damp = 1.0;
	if (efp->opts.elec_damp == EFP_ELEC_DAMP_SCREEN) {
		double screen_i = fr_i->screen_params[pt_i_idx];
		double screen_j = fr_j->screen_params[pt_j_idx];
		damp = get_damp_screen(efp, r, screen_i, screen_j);
	}
	energy += q1 * q2 / r * damp;

	/* monopole - dipole */
	energy += q2 / r3 * vec_dot(d1, &dr);
	energy -= q1 / r3 * vec_dot(d2, &dr);

	/* monopole - quadrupole */
	energy += q2 / r5 * quadrupole_sum(quad1, &dr);
	energy += q1 / r5 * quadrupole_sum(quad2, &dr);

	/* monopole - octupole */
	energy += q2 / r7 * octupole_sum(oct1, &dr);
	energy -= q1 / r7 * octupole_sum(oct2, &dr);

	/* dipole - dipole */
	energy += vec_dot(d1, d2) / r3;
	energy -= 3.0 / r5 * vec_dot(d1, &dr) * vec_dot(d2, &dr);

	/* dipole - quadrupole */
	double q1d2dr = 0.0;
	double q2d1dr = 0.0;
	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			int idx = quad_idx(a, b);
			q1d2dr += quad1[idx] * vec_el(d2, a) * vec_el(&dr, b);
			q2d1dr += quad2[idx] * vec_el(d1, a) * vec_el(&dr, b);
		}
	}

	energy += 2.0 / r5 * q1d2dr - 2.0 / r5 * q2d1dr;
	energy -= 5.0 / r7 * quadrupole_sum(quad1, &dr) * vec_dot(d2, &dr);
	energy += 5.0 / r7 * quadrupole_sum(quad2, &dr) * vec_dot(d1, &dr);

	/* quadrupole - quadrupole */
	double q1q2 = 0.0;
	double q1q2dr = 0.0;
	for (int a = 0; a < 3; a++) {
		double t1 = 0.0;
		double t2 = 0.0;
		for (int b = 0; b < 3; b++) {
			int idx = quad_idx(a, b);
			t1 += quad1[idx] * vec_el(&dr, b);
			t2 += quad2[idx] * vec_el(&dr, b);
			q1q2 += quad1[idx] * quad2[idx];
		}
		q1q2dr += t1 * t2;
	}

	double q1dr = quadrupole_sum(quad1, &dr);
	double q2dr = quadrupole_sum(quad2, &dr);

	energy += (2.0 / r5 * q1q2 - 20.0 / r7 * q1q2dr +
			35.0 / r9 * q1dr * q2dr) / 3.0;

	return energy;
}

static double
frag_frag_elec(struct efp *efp, int fr_i_idx, int fr_j_idx)
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
		for (int jj = 0; jj < fr_j->n_multipole_pts; jj++) {
			energy += mult_mult_energy(efp, fr_i_idx, fr_j_idx,
						   ii, jj);

			if (efp->opts.do_gradient)
				mult_mult_grad(efp, fr_i_idx, fr_j_idx, ii, jj);
		}

	return energy;
}

enum efp_result
efp_compute_elec(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_ELEC))
		return EFP_RESULT_SUCCESS;

	if (efp->opts.do_gradient)
		return EFP_RESULT_NOT_IMPLEMENTED;

	double energy = 0.0;

	for (int i = 0; i < efp->n_frag; i++)
		for (int j = i + 1; j < efp->n_frag; j++)
			energy += frag_frag_elec(efp, i, j);

	efp->energy.electrostatic = energy;
	return EFP_RESULT_SUCCESS;
}

static void
rotate_quadrupole(const mat_t *rotmat, const double *in, double *out)
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
rotate_octupole(const mat_t *rotmat, const double *in, double *out)
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
		const struct multipole_pt *pt_in =
					frag->lib->multipole_pts + i;
		struct multipole_pt *pt_out =
					frag->multipole_pts + i;

		/* move point position */
		move_pt(VEC(frag->x), rotmat, VEC(pt_in->x), VEC(pt_out->x));

		/* rotate dipole */
		mat_vec(rotmat, &pt_in->dipole, &pt_out->dipole);

		/* rotate quadrupole */
		rotate_quadrupole(rotmat, pt_in->quadrupole,
				  pt_out->quadrupole);

		/* correction for Buckingham quadrupoles */
		double *quad = pt_out->quadrupole;

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
		rotate_octupole(rotmat, pt_in->octupole, pt_out->octupole);

		/* correction for Buckingham octupoles */
		double *oct = pt_out->octupole;

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
compute_ai_elec_frag(struct efp *efp, int frag_idx)
{
	double energy = 0.0;
	struct frag *fr_i = efp->frags + frag_idx;

	for (int i = 0; i < fr_i->n_atoms; i++) {
		for (int j = 0; j < efp->qm_data.n_atoms; j++) {
			struct efp_atom *at_i = fr_i->atoms + i;
			const double *xyz_j = efp->qm_data.atom_xyz + 3 * j;
			double znuc_j = efp->qm_data.atom_znuc[j];

			double r = vec_dist(VEC(at_i->x), (const vec_t *)xyz_j);
			energy += at_i->znuc * znuc_j / r;
		}
	}
	for (int i = 0; i < fr_i->n_multipole_pts; i++) {
		for (int j = 0; j < efp->qm_data.n_atoms; j++) {
			const double *xyz = efp->qm_data.atom_xyz + 3 * j;
			double znuc = efp->qm_data.atom_znuc[j];

			energy += charge_mult_energy(efp, znuc,
					(const vec_t *)xyz, frag_idx, i);
		}
	}
	return energy;
}

enum efp_result
efp_compute_ai_elec(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_AI_ELEC))
		return EFP_RESULT_SUCCESS;

	if (efp->opts.do_gradient)
		return EFP_RESULT_NOT_IMPLEMENTED;

	double energy = 0.0;

	for (int i = 0; i < efp->n_frag; i++)
		energy += compute_ai_elec_frag(efp, i);

	efp->energy.ai_electrostatic = energy;
	return EFP_RESULT_SUCCESS;
}
