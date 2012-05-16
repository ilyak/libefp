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

#define DR_A ((double *)&dr)[a]
#define DR_B ((double *)&dr)[b]

static void
add_multipole_field(struct polarizable_pt *pt,
		    const struct multipole_pt *mult_pt)
{
	vec_t dr = vec_sub(VEC(pt->x), VEC(mult_pt->x));

	double t1, t2;
	double r = vec_len(&dr);

	double ri[8];
	powers(1.0 / r, 8, ri);

	/* charge */
	pt->elec_field.x += mult_pt->monopole * dr.x * ri[3];
	pt->elec_field.y += mult_pt->monopole * dr.y * ri[3];
	pt->elec_field.z += mult_pt->monopole * dr.z * ri[3];

	/* dipole */
	t1 = mult_pt->dipole.x * dr.x +
	     mult_pt->dipole.y * dr.y +
	     mult_pt->dipole.z * dr.z;

	pt->elec_field.x -= mult_pt->dipole.x * ri[3] - 3 * ri[5] * t1 * dr.x;
	pt->elec_field.y -= mult_pt->dipole.y * ri[3] - 3 * ri[5] * t1 * dr.y;
	pt->elec_field.z -= mult_pt->dipole.z * ri[3] - 3 * ri[5] * t1 * dr.z;

	/* quadrupole */
	t1 = 0.0;
	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			t1 += mult_pt->quadrupole[quad_idx(a, b)] * DR_A * DR_B;

	t2 = mult_pt->quadrupole[quad_idx(0, 0)] * dr.x +
	     mult_pt->quadrupole[quad_idx(1, 0)] * dr.y +
	     mult_pt->quadrupole[quad_idx(2, 0)] * dr.z;
	pt->elec_field.x -= 2 * ri[5] * t2 - 5 * ri[7] * t1 * dr.x;

	t2 = mult_pt->quadrupole[quad_idx(0, 1)] * dr.x +
	     mult_pt->quadrupole[quad_idx(1, 1)] * dr.y +
	     mult_pt->quadrupole[quad_idx(2, 1)] * dr.z;
	pt->elec_field.y -= 2 * ri[5] * t2 - 5 * ri[7] * t1 * dr.y;

	t2 = mult_pt->quadrupole[quad_idx(0, 2)] * dr.x +
	     mult_pt->quadrupole[quad_idx(1, 2)] * dr.y +
	     mult_pt->quadrupole[quad_idx(2, 2)] * dr.z;
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

		/* field due to nuclei */
		for (int j = 0; j < efp->frags[i].n_atoms; j++) {
			struct efp_atom *at = efp->frags[i].atoms + j;

			vec_t dr = vec_sub(VEC(pt->x), VEC(at->x));
			double r = vec_len(&dr);
			double r3 = r * r * r;

			pt->elec_field.x += at->znuc * dr.x / r3;
			pt->elec_field.y += at->znuc * dr.y / r3;
			pt->elec_field.z += at->znuc * dr.z / r3;
		}

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
			const double *xyz = efp->qm_data.atom_xyz + 3 * i;
			double znuc = efp->qm_data.atom_znuc[i];

			vec_t dr = vec_sub(VEC(pt->x), (const vec_t *)xyz);
			double t = 1.0 / vec_len(&dr);
			t = t * t * t;

			pt->elec_field.x += znuc * dr.x * t;
			pt->elec_field.y += znuc * dr.y * t;
			pt->elec_field.z += znuc * dr.z * t;
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

static void
get_induced_dipole_field(struct efp *efp, int frag_idx,
			 struct polarizable_pt *pt,
			 vec_t *field, vec_t *field_conj)
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

			vec_t dr = vec_sub(VEC(pt->x), VEC(pt_j->x));
			double r = vec_len(&dr);

			double ri[6];
			powers(1.0 / r, 6, ri);

			double sum = vec_dot(&pt_j->induced_dipole, &dr);

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
			vec_t field, field_conj;
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

double
efp_compute_pol_energy(struct efp *efp)
{
	static const double scf_conv_epsilon = 1.0e-16;

	/* compute induced dipoles self consistently */
	while (pol_scf_iter(efp) > scf_conv_epsilon);

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

void
efp_pol_scf_init(struct efp *efp)
{
	compute_elec_field(efp);

	/* set initial approximation - all induced dipoles are zero */
	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;
		for (int j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;
			vec_zero(&pt->induced_dipole);
			vec_zero(&pt->induced_dipole_conj);
		}
	}
}

static void
charge_dipole_grad(const vec_t *p1, double q1, const vec_t *com1,
		   const vec_t *p2, const vec_t *d2, const vec_t *com2,
		   double *grad1, double *grad2)
{
	vec_t dr = vec_sub(p1, p2);

	double r = vec_len(&dr);
	double r3 = r * r * r;
	double r5 = r3 * r * r;

	double t1 = 3.0 * q1 / r5 * vec_dot(d2, &dr);
	double t2 = q1 / r3;

	double gx = t2 * d2->x - t1 * dr.x;
	double gy = t2 * d2->y - t1 * dr.y;
	double gz = t2 * d2->z - t1 * dr.z;

	grad1[0] -= gx;
	grad1[1] -= gy;
	grad1[2] -= gz;
	grad1[3] -= gz * (p2->y - com2->y) - gy * (p2->z - com2->z);
	grad1[4] -= gx * (p2->z - com2->z) - gz * (p2->x - com2->x);
	grad1[5] -= gy * (p2->x - com2->x) - gx * (p2->y - com2->y);

	grad2[0] += gx;
	grad2[1] += gy;
	grad2[2] += gz;
	grad2[3] += gz * (p1->y - com1->y) - gy * (p1->z - com1->z);
	grad2[4] += gx * (p1->z - com1->z) - gz * (p1->x - com1->x);
	grad2[5] += gy * (p1->x - com1->x) - gx * (p1->y - com1->y);
	grad2[3] += t2 * (d2->y * dr.z - d2->z * dr.y);
	grad2[4] += t2 * (d2->z * dr.x - d2->x * dr.z);
	grad2[5] += t2 * (d2->x * dr.y - d2->y * dr.x);
}

static void
dipole_dipole_grad(const vec_t *p1, const vec_t *d1, const vec_t *com1,
		   const vec_t *p2, const vec_t *d2, const vec_t *com2,
		   double *grad1, double *grad2)
{
	vec_t dr = vec_sub(p1, p2);

	double r = vec_len(&dr);
	double r3 = r * r * r;
	double r5 = r3 * r * r;
	double r7 = r5 * r * r;

	double d1dr = vec_dot(d1, &dr);
	double d2dr = vec_dot(d2, &dr);

	double t1 = 3.0 / r5 * vec_dot(d1, d2) - 15.0 / r7 * d1dr * d2dr;
	double t2 = 3.0 / r5;

	double gx = t1 * dr.x + t2 * (d2dr * d1->x + d1dr * d2->x);
	double gy = t1 * dr.x + t2 * (d2dr * d1->y + d1dr * d2->y);
	double gz = t1 * dr.x + t2 * (d2dr * d1->z + d1dr * d2->z);

	double dt1x = d1->y * (d2->z / r3 - t2 * dr.z * d2dr) -
		d1->z * (d2->y / r3 - t2 * dr.y * d2dr);
	double dt1y = d1->z * (d2->x / r3 - t2 * dr.x * d2dr) -
		d1->x * (d2->z / r3 - t2 * dr.z * d2dr);
	double dt1z = d1->x * (d2->y / r3 - t2 * dr.y * d2dr) -
		d1->y * (d2->x / r3 - t2 * dr.x * d2dr);

	double dt2x = d2->y * (d1->z / r3 - t2 * dr.z * d1dr) -
		d2->z * (d1->y / r3 - t2 * dr.y * d1dr);
	double dt2y = d2->z * (d1->x / r3 - t2 * dr.x * d1dr) -
		d2->x * (d1->z / r3 - t2 * dr.z * d1dr);
	double dt2z = d2->x * (d1->y / r3 - t2 * dr.y * d1dr) -
		d2->y * (d1->x / r3 - t2 * dr.x * d1dr);

	grad1[0] += gx;
	grad1[1] += gy;
	grad1[2] += gz;
	grad1[3] += gz * (p1->y - com1->y) - gy * (p1->z - com1->z) + dt1x;
	grad1[4] += gx * (p1->z - com1->z) - gz * (p1->x - com1->x) + dt1y;
	grad1[5] += gy * (p1->x - com1->x) - gx * (p1->y - com1->y) + dt1z;

	grad2[0] -= gx;
	grad2[1] -= gy;
	grad2[2] -= gz;
	grad2[3] -= gz * (p2->y - com2->y) - gy * (p2->z - com2->z) + dt2x;
	grad2[4] -= gx * (p2->z - com2->z) - gz * (p2->x - com2->x) + dt2y;
	grad2[5] -= gy * (p2->x - com2->x) - gx * (p2->y - com2->y) + dt2z;
}

static void
dipole_quadrupole_grad(const vec_t *p1, const vec_t *d1, const vec_t *com1,
		       const vec_t *p2, const double *quad2, const vec_t *com2,
		       double *grad1, double *grad2)
{
	vec_t dr = vec_sub(p1, p2);

	double r = vec_len(&dr);
	double r2 = r * r;
	double r5 = r2 * r2 * r;
	double r7 = r5 * r2;
	double r9 = r7 * r2;

	double q2s = 0.0;
	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			q2s += quad2[quad_idx(a, b)] * DR_A * DR_B;

	double q2sx = 0.0, q2sy = 0.0, q2sz = 0.0;
	for (int a = 0; a < 3; a++) {
		q2sx += quad2[quad_idx(0, a)] * DR_A;
		q2sy += quad2[quad_idx(1, a)] * DR_A;
		q2sz += quad2[quad_idx(2, a)] * DR_A;
	}

	double d1dr = vec_dot(d1, &dr);

	double t1 = d1->x * q2sx + d1->y * q2sy + d1->z * q2sz;
	double t2 = -10.0 / r7 * t1 + 35.0 / r9 * q2s * d1dr;

	double d1q2x = d1->x * quad2[quad_idx(0, 0)] +
		       d1->y * quad2[quad_idx(0, 1)] +
		       d1->z * quad2[quad_idx(0, 2)];
	double d1q2y = d1->x * quad2[quad_idx(1, 0)] +
		       d1->y * quad2[quad_idx(1, 1)] +
		       d1->z * quad2[quad_idx(1, 2)];
	double d1q2z = d1->x * quad2[quad_idx(2, 0)] +
		       d1->y * quad2[quad_idx(2, 1)] +
		       d1->z * quad2[quad_idx(2, 2)];

	double q2xdr = dr.x * quad2[quad_idx(0, 0)] +
		       dr.y * quad2[quad_idx(0, 1)] +
		       dr.z * quad2[quad_idx(0, 2)];
	double q2ydr = dr.x * quad2[quad_idx(1, 0)] +
		       dr.y * quad2[quad_idx(1, 1)] +
		       dr.z * quad2[quad_idx(1, 2)];
	double q2zdr = dr.x * quad2[quad_idx(2, 0)] +
		       dr.y * quad2[quad_idx(2, 1)] +
		       dr.z * quad2[quad_idx(2, 2)];

	double gx = t2 * dr.x + 2.0 / r5 * d1q2x - 5.0 / r7 * q2s * d1->x -
			10.0 / r7 * q2xdr * d1dr;
	double gy = t2 * dr.y + 2.0 / r5 * d1q2y - 5.0 / r7 * q2s * d1->y -
			10.0 / r7 * q2ydr * d1dr;
	double gz = t2 * dr.z + 2.0 / r5 * d1q2z - 5.0 / r7 * q2s * d1->z -
			10.0 / r7 * q2zdr * d1dr;

	double td1x = 2.0 / r5 * (d1->z * q2ydr - d1->y * q2zdr) +
			5.0 / r7 * q2s * (dr.z * d1->y - dr.y * d1->z);
	double td1y = 2.0 / r5 * (d1->x * q2zdr - d1->z * q2xdr) +
			5.0 / r7 * q2s * (dr.x * d1->z - dr.z * d1->x);
	double td1z = 2.0 / r5 * (d1->y * q2xdr - d1->x * q2ydr) +
			5.0 / r7 * q2s * (dr.y * d1->x - dr.x * d1->y);

	double tq2x = -10.0 / r7 * d1dr * (q2ydr * dr.z - q2zdr * dr.y) -
			2.0 / r5 * ((q2zdr * d1->y + dr.y * d1q2z) -
				    (q2ydr * d1->z + dr.z * d1q2y));
	double tq2y = -10.0 / r7 * d1dr * (q2zdr * dr.x - q2xdr * dr.z) -
			2.0 / r5 * ((q2xdr * d1->z + dr.z * d1q2x) -
				    (q2zdr * d1->x + dr.x * d1q2z));
	double tq2z = -10.0 / r7 * d1dr * (q2xdr * dr.y - q2ydr * dr.x) -
			2.0 / r5 * ((q2ydr * d1->x + dr.x * d1q2y) -
				    (q2xdr * d1->y + dr.y * d1q2x));

	grad1[0] += gx;
	grad1[1] += gy;
	grad1[2] += gz;
	grad1[3] += gz * (p1->y - com1->y) - gy * (p1->z - com1->z) - td1x;
	grad1[4] += gx * (p1->z - com1->z) - gz * (p1->x - com1->x) - td1y;
	grad1[5] += gy * (p1->x - com1->x) - gx * (p1->y - com1->y) - td1z;

	grad2[0] -= gx;
	grad2[1] -= gy;
	grad2[2] -= gz;
	grad2[3] -= gz * (p2->y - com2->y) - gy * (p2->z - com2->z) - tq2x;
	grad2[4] -= gx * (p2->z - com2->z) - gz * (p2->x - com2->x) - tq2y;
	grad2[5] -= gy * (p2->x - com2->x) - gx * (p2->y - com2->y) - tq2z;
}

static void
compute_grad_point(struct efp *efp, int frag_idx, int pt_idx)
{
	struct frag *fr_i = efp->frags + frag_idx;
	struct polarizable_pt *pt_i = fr_i->polarizable_pts + pt_idx;
	double *gr_i = efp->grad + 6 * frag_idx;

	vec_t dipole = {
		0.5 * (pt_i->induced_dipole.x + pt_i->induced_dipole_conj.x),
		0.5 * (pt_i->induced_dipole.y + pt_i->induced_dipole_conj.y),
		0.5 * (pt_i->induced_dipole.z + pt_i->induced_dipole_conj.z)
	};

	for (int j = 0; j < efp->n_frag; j++) {
		if (j == frag_idx)
			continue;

		struct frag *fr_j = efp->frags + j;
		double *gr_j = efp->grad + 6 * j;

		/* induced dipole - nuclei */
		for (int k = 0; k < fr_j->n_atoms; k++) {
			struct efp_atom *at_j = fr_j->atoms + k;
			charge_dipole_grad(VEC(at_j->x), at_j->znuc,
					   VEC(fr_j->x), VEC(pt_i->x), &dipole,
					   VEC(fr_i->x), gr_i, gr_j);
		}

		/* induced dipole - multipoles */
		for (int k = 0; k < fr_j->n_multipole_pts; k++) {
			struct multipole_pt *pt_j = fr_j->multipole_pts + k;

			/* induced dipole - charge */
			charge_dipole_grad(VEC(pt_j->x), pt_j->monopole,
					   VEC(fr_j->x), VEC(pt_i->x), &dipole,
					   VEC(fr_i->x), gr_i, gr_j);

			/* induced dipole - dipole */
			dipole_dipole_grad(VEC(pt_i->x), &dipole, VEC(fr_i->x),
					   VEC(pt_j->x), &pt_j->dipole,
					   VEC(fr_j->x), gr_i, gr_j);

			/* induced dipole - quadrupole */
			dipole_quadrupole_grad(VEC(pt_i->x), &dipole,
					       VEC(fr_i->x), VEC(pt_j->x),
					       pt_j->quadrupole, VEC(fr_j->x),
					       gr_i, gr_j);

			/* octupole-polarizability interactions are ignored */
		}

		/* induced dipole - induced dipoles */
		for (int k = 0; k < fr_j->n_polarizable_pts; k++) {
			struct polarizable_pt *pt_j =
					fr_j->polarizable_pts + k;

			dipole_dipole_grad(VEC(pt_i->x), &dipole, VEC(fr_i->x),
					   VEC(pt_j->x), &pt_j->induced_dipole,
					   VEC(fr_j->x), gr_i, gr_j);
		}
	}

	if (efp->opts.terms & EFP_TERM_AI_POL) {
		/* XXX */
	}
}

static void
compute_grad(struct efp *efp)
{
	for (int i = 0; i < efp->n_frag; i++)
		for (int j = 0; j < efp->frags[i].n_polarizable_pts; j++)
			compute_grad_point(efp, i, j);
}

enum efp_result
efp_compute_pol(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_POL))
		return EFP_RESULT_SUCCESS;

	efp->energy.polarization = efp_compute_pol_energy(efp);

	if (efp->grad)
		compute_grad(efp);

	return EFP_RESULT_SUCCESS;
}

void
efp_update_pol(struct frag *frag, const mat_t *rotmat)
{
	for (int i = 0; i < frag->n_polarizable_pts; i++) {
		move_pt(VEC(frag->x), rotmat,
			VEC(frag->lib->polarizable_pts[i].x),
			VEC(frag->polarizable_pts[i].x));

		const mat_t *in = &frag->lib->polarizable_pts[i].tensor;
		mat_t *out = &frag->polarizable_pts[i].tensor;

		rotate_t2(rotmat, (const double *)in, (double *)out);
	}
}
