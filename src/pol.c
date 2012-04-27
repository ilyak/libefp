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

static void
add_multipole_field(struct polarizable_pt *pt,
		    const struct multipole_pt *mult_pt)
{
	struct vec dr = {
		pt->x - mult_pt->x,
		pt->y - mult_pt->y,
		pt->z - mult_pt->z
	};

	double t1, t2;
	double r = vec_len(&dr);

	double ri[8];
	powers(1.0 / r, 8, ri);

	/* charge */
	pt->elec_field.x += mult_pt->monopole * dr.x * ri[3];
	pt->elec_field.y += mult_pt->monopole * dr.y * ri[3];
	pt->elec_field.z += mult_pt->monopole * dr.z * ri[3];

	/* dipole */
	t1 = mult_pt->dipole[0] * dr.x +
	     mult_pt->dipole[1] * dr.y +
	     mult_pt->dipole[2] * dr.z;

	pt->elec_field.x -= mult_pt->dipole[0] * ri[3] - 3 * ri[5] * t1 * dr.x;
	pt->elec_field.y -= mult_pt->dipole[1] * ri[3] - 3 * ri[5] * t1 * dr.y;
	pt->elec_field.z -= mult_pt->dipole[2] * ri[3] - 3 * ri[5] * t1 * dr.z;

	/* quadrupole */
	t1 = 0.0;
	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			t1 += mult_pt->quadrupole[t2_idx(a, b)] *
				ELV(&dr, a) * ELV(&dr, b);

	t2 = mult_pt->quadrupole[t2_idx(0, 0)] * dr.x +
	     mult_pt->quadrupole[t2_idx(1, 0)] * dr.y +
	     mult_pt->quadrupole[t2_idx(2, 0)] * dr.z;
	pt->elec_field.x -= 2 * ri[5] * t2 - 5 * ri[7] * t1 * dr.x;

	t2 = mult_pt->quadrupole[t2_idx(0, 1)] * dr.x +
	     mult_pt->quadrupole[t2_idx(1, 1)] * dr.y +
	     mult_pt->quadrupole[t2_idx(2, 1)] * dr.z;
	pt->elec_field.y -= 2 * ri[5] * t2 - 5 * ri[7] * t1 * dr.y;

	t2 = mult_pt->quadrupole[t2_idx(0, 2)] * dr.x +
	     mult_pt->quadrupole[t2_idx(1, 2)] * dr.y +
	     mult_pt->quadrupole[t2_idx(2, 2)] * dr.z;
	pt->elec_field.z -= 2 * ri[5] * t2 - 5 * ri[7] * t1 * dr.z;

	/* octupole-polarizability interactions are ignored */
}

static void
compute_elec_field_pt(struct efp *efp, int frag_idx, int pt_idx)
{
	struct frag *frag = efp->frags + frag_idx;
	struct polarizable_pt *pt = frag->polarizable_pts + pt_idx;

	vec_zero(&pt->elec_field);

	for (int i = 0; i < efp->n_frag; i++) {
		if (i == frag_idx)
			continue;

		/* field due to multipoles */
		for (int j = 0; j < efp->frags[i].n_multipole_pts; j++) {
			struct multipole_pt *mult_pt =
					efp->frags[i].multipole_pts + j;
			add_multipole_field(pt, mult_pt);
		}
	}

	if (efp->opts.terms & EFP_TERM_AI_POL) {
		/* field due to QM nuclei */
		for (int i = 0; i < efp->qm_data.n_atoms; i++) {
			struct efp_qm_atom *atom = efp->qm_data.atoms + i;

			struct vec dr = {
				pt->x - atom->x,
				pt->y - atom->y,
				pt->z - atom->z
			};

			double t = vec_len(&dr);
			t = t * t * t;
			t = 1.0 / t;

			pt->elec_field.x += atom->charge * dr.x * t;
			pt->elec_field.y += atom->charge * dr.y * t;
			pt->elec_field.z += atom->charge * dr.z * t;
		}
	}
}

static void
add_electron_density_field(struct efp *efp)
{
	int n_pt = 0;
	for (int i = 0; i < efp->n_frag; i++)
		n_pt += efp->frags[i].n_polarizable_pts;

	double xyz[3 * n_pt], field[3 * n_pt], *ptr;

	ptr = xyz;
	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;
			*ptr++ = pt->x;
			*ptr++ = pt->y;
			*ptr++ = pt->z;
		}
	}

	efp->callbacks.get_electron_density_field(n_pt, xyz, field,
			efp->callbacks.get_electron_density_field_user_data);

	ptr = field;
	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;
			pt->elec_field.x += *ptr++;
			pt->elec_field.y += *ptr++;
			pt->elec_field.z += *ptr++;
		}
	}
}

static void
compute_elec_field(struct efp *efp)
{
	for (int i = 0; i < efp->n_frag; i++)
		for (int j = 0; j < efp->frags[i].n_polarizable_pts; j++)
			compute_elec_field_pt(efp, i, j);

	if (efp->opts.terms & EFP_TERM_AI_POL)
		add_electron_density_field(efp);
}

static double
compute_pol_energy(struct efp *efp)
{
	double energy = 0.0;

	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;
			energy -= 0.5 * vec_dot(&pt->induced_dipole,
						&pt->elec_field);
		}
	}
	return energy;
}

static void
get_induced_dipole_field(struct efp *efp, int frag_idx,
			 struct polarizable_pt *pt,
			 struct vec *field, struct vec *field_conj)
{
	vec_zero(field);
	vec_zero(field_conj);

	for (int j = 0; j < efp->n_frag; j++) {
		if (j == frag_idx)
			continue;

		struct frag *fr_j = efp->frags + j;

		for (int jj = 0; jj < fr_j->n_polarizable_pts; jj++) {
			struct polarizable_pt *pt_j =
				fr_j->polarizable_pts + jj;

			struct vec dr = {
				pt->x - pt_j->x,
				pt->y - pt_j->y,
				pt->z - pt_j->z
			};

			double r, sum, ri[6];
			r = vec_len(&dr);
			powers(1.0 / r, 6, ri);

			sum = pt_j->induced_dipole.x * dr.x +
			      pt_j->induced_dipole.y * dr.y +
			      pt_j->induced_dipole.z * dr.z;

			field->x -= ri[3] * pt_j->induced_dipole.x -
						3.0 * ri[5] * sum * dr.x;
			field->y -= ri[3] * pt_j->induced_dipole.y -
						3.0 * ri[5] * sum * dr.y;
			field->z -= ri[3] * pt_j->induced_dipole.z -
						3.0 * ri[5] * sum * dr.z;

			sum = pt_j->induced_dipole_conj.x * dr.x +
			      pt_j->induced_dipole_conj.y * dr.y +
			      pt_j->induced_dipole_conj.z * dr.z;

			field_conj->x -= ri[3] * pt_j->induced_dipole_conj.x -
						3.0 * ri[5] * sum * dr.x;
			field_conj->y -= ri[3] * pt_j->induced_dipole_conj.y -
						3.0 * ri[5] * sum * dr.y;
			field_conj->z -= ri[3] * pt_j->induced_dipole_conj.z -
						3.0 * ri[5] * sum * dr.z;
		}
	}
}

static double
pol_scf_iter(struct efp *efp)
{
	/* compute new induced dipoles on polarizable points */
	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;

			/* electric field from other induced dipoles */
			struct vec field, field_conj;
			get_induced_dipole_field(efp, i, pt, &field,
					&field_conj);

			/* add field that doesn't change during scf */
			field.x += pt->elec_field.x;
			field.y += pt->elec_field.y;
			field.z += pt->elec_field.z;
			field_conj.x += pt->elec_field.x;
			field_conj.y += pt->elec_field.y;
			field_conj.z += pt->elec_field.z;

			mat_vec(&pt->tensor, &field, &pt->induced_dipole_new);
			mat_trans_vec(&pt->tensor, &field_conj,
					&pt->induced_dipole_conj_new);
		}
	}

	int n_pt = 0;
	double conv = 0.0;

	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		n_pt += frag->n_polarizable_pts;

		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;

			conv += vec_dist_2(&pt->induced_dipole_new,
					   &pt->induced_dipole);

			pt->induced_dipole = pt->induced_dipole_new;
			pt->induced_dipole_conj = pt->induced_dipole_conj_new;
		}
	}
	return conv / n_pt;
}

enum efp_result
efp_scf_update_pol(struct efp *efp, double *energy)
{
	static const double conv_epsilon = 1.0e-16;

	compute_elec_field(efp);

	/* compute induced dipoles self consistently */
	while (pol_scf_iter(efp) > conv_epsilon);

	*energy += compute_pol_energy(efp);
	return EFP_RESULT_SUCCESS;
}

static void
compute_grad(struct efp *efp)
{
	/* XXX */
}

enum efp_result
efp_compute_pol(struct efp *efp)
{
	/* induced dipoles are computed during SCF in efp_scf_update_pol */

	int idx = efp_get_term_index(EFP_TERM_POL);
	efp->energy[idx] = compute_pol_energy(efp);

	if (efp->grad)
		compute_grad(efp);

	return EFP_RESULT_SUCCESS;
}

enum efp_result
efp_compute_ai_pol(struct efp *efp)
{
	/* because polarization is computed self consistently we can't
	   separate it into EFP/EFP and EFP/AI parts. So everything goes
	   into EFP_TERM_POL and EFP_TERM_AI_POL is zero */
	efp->energy[efp_get_term_index(EFP_TERM_AI_POL)] = 0.0;
	return EFP_RESULT_SUCCESS;
}

static void
rotate_tensor(const struct mat *rotmat, const struct mat *in, struct mat *out)
{
	mat_zero(out);

	for (int a1 = 0; a1 < 3; a1++)
	for (int b1 = 0; b1 < 3; b1++)
		for (int a2 = 0; a2 < 3; a2++)
		for (int b2 = 0; b2 < 3; b2++)
			ELM(out, a2, b2) += ELM(in, a1, b1) *
				ELM(rotmat, a2, a1) * ELM(rotmat, b2, b1);
}

void
efp_update_pol(struct frag *frag, const struct mat *rotmat)
{
	for (int i = 0; i < frag->n_polarizable_pts; i++) {
		move_pt(VEC(frag->x), rotmat,
			VEC(frag->lib->polarizable_pts[i].x),
			VEC(frag->polarizable_pts[i].x));

		rotate_tensor(rotmat, &frag->lib->polarizable_pts[i].tensor,
				&frag->polarizable_pts[i].tensor);
	}
}
