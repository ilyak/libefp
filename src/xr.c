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

#include <stdlib.h>
#include <string.h>

#include "cblas.h"
#include "efp_private.h"
#include "disp.h"

static inline int
fock_idx(int i, int j)
{
	return i < j ?
		j * (j + 1) / 2 + i :
		i * (i + 1) / 2 + j;
}

static inline double
valence(double n)
{
	if (n > 2)
		n -= 2;
	if (n > 8)
		n -= 8;
	return n;
}

static double
charge_penetration_energy(double s_ij, double r_ij)
{
	if (fabs(s_ij) < 1.0e-6)
		return 0.0;

	double ln_s = log(fabs(s_ij));

	return -s_ij * s_ij / r_ij / sqrt(-2.0 * ln_s);
}

static void
charge_penetration_grad(struct frag *fr_i, struct frag *fr_j,
			int lmo_i_idx, int lmo_j_idx, double s_ij,
			const six_t *ds_ij)
{
	if (fabs(s_ij) < 1.0e-6)
		return;

	const vec_t *ct_i = fr_i->lmo_centroids + lmo_i_idx;
	const vec_t *ct_j = fr_j->lmo_centroids + lmo_j_idx;

	vec_t dr = vec_sub(ct_j, ct_i);

	double r_ij = vec_len(&dr);
	double ln_s = log(fabs(s_ij));

	double t1 = s_ij * s_ij / (r_ij * r_ij * r_ij) / sqrt(-2.0 * ln_s);
	double t2 = -s_ij / r_ij / sqrt(2.0) * (2.0 / sqrt(-ln_s) +
				0.5 / sqrt(-ln_s * ln_s * ln_s));

	vec_t force, torque_i, torque_j;

	force.x = t2 * ds_ij->x + t1 * dr.x;
	force.y = t2 * ds_ij->y + t1 * dr.y;
	force.z = t2 * ds_ij->z + t1 * dr.z;

	torque_i.x = t2 * ds_ij->a + t1 * (dr.y * (ct_i->z - fr_j->z) -
					   dr.z * (ct_i->y - fr_j->y));
	torque_i.y = t2 * ds_ij->b + t1 * (dr.z * (ct_i->x - fr_j->x) -
					   dr.x * (ct_i->z - fr_j->z));
	torque_i.z = t2 * ds_ij->c + t1 * (dr.x * (ct_i->y - fr_j->y) -
					   dr.y * (ct_i->x - fr_j->x));

	torque_j.x = force.y * (fr_j->z - fr_i->z) - force.z * (fr_j->y - fr_i->y);
	torque_j.y = force.z * (fr_j->x - fr_i->x) - force.x * (fr_j->z - fr_i->z);
	torque_j.z = force.x * (fr_j->y - fr_i->y) - force.y * (fr_j->x - fr_i->x);

	#pragma omp atomic
	fr_i->force.x += force.x;
	#pragma omp atomic
	fr_i->force.y += force.y;
	#pragma omp atomic
	fr_i->force.z += force.z;
	#pragma omp atomic
	fr_i->torque.x += torque_i.x;
	#pragma omp atomic
	fr_i->torque.y += torque_i.y;
	#pragma omp atomic
	fr_i->torque.z += torque_i.z;

	#pragma omp atomic
	fr_j->force.x -= force.x;
	#pragma omp atomic
	fr_j->force.y -= force.y;
	#pragma omp atomic
	fr_j->force.z -= force.z;
	#pragma omp atomic
	fr_j->torque.x -= torque_j.x;
	#pragma omp atomic
	fr_j->torque.y -= torque_j.y;
	#pragma omp atomic
	fr_j->torque.z -= torque_j.z;
}

static void
transform_integrals(int n_lmo_i, int n_lmo_j,
		    int wf_size_i, int wf_size_j,
		    const double *wf_i, const double *wf_j,
		    const double *s, double *lmo_s)
{
	double tmp[n_lmo_i * wf_size_j];

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		n_lmo_i, wf_size_j, wf_size_i, 1.0, wf_i, wf_size_i,
		s, wf_size_j, 0.0, tmp, wf_size_j);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
		n_lmo_i, n_lmo_j, wf_size_j, 1.0, tmp, wf_size_j,
		wf_j, wf_size_j, 0.0, lmo_s, n_lmo_j);
}

static void
transform_integral_derivatives(int n_lmo_i, int n_lmo_j,
			       int wf_size_i, int wf_size_j,
			       const double *wf_i, const double *wf_j,
			       const six_t *s, six_t *lmo_s)
{
	int wf_size = wf_size_i * wf_size_j;
	int lmo_size = n_lmo_i * n_lmo_j;

	double *ss = malloc(wf_size * 6 * sizeof(double));
	double *lmo_ss = malloc(lmo_size * 6 * sizeof(double));

	/* unpack */
	const six_t *s_ptr = s;

	for (int i = 0; i < wf_size; i++, s_ptr++) {
		ss[0 * wf_size + i] = s_ptr->x;
		ss[1 * wf_size + i] = s_ptr->y;
		ss[2 * wf_size + i] = s_ptr->z;
		ss[3 * wf_size + i] = s_ptr->a;
		ss[4 * wf_size + i] = s_ptr->b;
		ss[5 * wf_size + i] = s_ptr->c;
	}

	/* transform */
	for (int k = 0; k < 6; k++) {
		transform_integrals(n_lmo_i, n_lmo_j, wf_size_i, wf_size_j,
				wf_i, wf_j, ss + k * wf_size, lmo_ss + k * lmo_size);
	}

	/* pack */
	six_t *lmo_s_ptr = lmo_s;

	for (int i = 0; i < lmo_size; i++, lmo_s_ptr++) {
		lmo_s_ptr->x = lmo_ss[0 * lmo_size + i];
		lmo_s_ptr->y = lmo_ss[1 * lmo_size + i];
		lmo_s_ptr->z = lmo_ss[2 * lmo_size + i];
		lmo_s_ptr->a = lmo_ss[3 * lmo_size + i];
		lmo_s_ptr->b = lmo_ss[4 * lmo_size + i];
		lmo_s_ptr->c = lmo_ss[5 * lmo_size + i];
	}

	free(ss);
	free(lmo_ss);
}

static void
frag_frag_xr_grad(struct efp *efp, const double *lmo_s, const double *lmo_t,
		  int frag_i, int frag_j, int overlap_idx)
{
	struct frag *fr_i = efp->frags + frag_i;
	struct frag *fr_j = efp->frags + frag_j;

	six_t *ds = malloc(fr_i->xr_wf_size * fr_j->xr_wf_size * sizeof(six_t));
	six_t *dt = malloc(fr_i->xr_wf_size * fr_j->xr_wf_size * sizeof(six_t));

	six_t *lmo_ds = malloc(fr_i->n_lmo * fr_j->n_lmo * sizeof(six_t));
	six_t *lmo_dt = malloc(fr_i->n_lmo * fr_j->n_lmo * sizeof(six_t));

	/* compute derivatives of S and T integrals */
	efp_st_int_deriv(fr_i->n_xr_shells, fr_i->xr_shells,
			 fr_j->n_xr_shells, fr_j->xr_shells,
			 VEC(fr_i->x), fr_i->xr_wf_size, fr_j->xr_wf_size,
			 ds, dt);

	transform_integral_derivatives(fr_i->n_lmo, fr_j->n_lmo,
				       fr_i->xr_wf_size, fr_j->xr_wf_size,
				       fr_i->xr_wf, fr_j->xr_wf,
				       ds, lmo_ds);
	transform_integral_derivatives(fr_i->n_lmo, fr_j->n_lmo,
				       fr_i->xr_wf_size, fr_j->xr_wf_size,
				       fr_i->xr_wf, fr_j->xr_wf,
				       dt, lmo_dt);

	for (int i = 0, idx = overlap_idx; i < fr_i->n_lmo; i++) {
		for (int j = 0; j < fr_j->n_lmo; j++, idx++) {
			int ij = i * fr_j->n_lmo + j;

			if ((efp->opts.terms & EFP_TERM_DISP) &&
			    (efp->opts.disp_damp == EFP_DISP_DAMP_OVERLAP)) {
				fr_i->overlap_int_deriv[idx] = lmo_ds[ij];
			}

			if ((efp->opts.terms & EFP_TERM_ELEC) &&
			    (efp->opts.elec_damp == EFP_ELEC_DAMP_OVERLAP)) {
				charge_penetration_grad(fr_i, fr_j, i, j,
							lmo_s[ij], lmo_ds + ij);
			}
		}
	}

	free(ds), free(lmo_ds);
	free(dt), free(lmo_dt);
}

static void
frag_frag_xr(struct efp *efp, int frag_i, int frag_j, int overlap_idx,
	     double *exr_out, double *ecp_out)
{
	struct frag *fr_i = efp->frags + frag_i;
	struct frag *fr_j = efp->frags + frag_j;

	double *s = malloc(fr_i->xr_wf_size * fr_j->xr_wf_size * sizeof(double));
	double *t = malloc(fr_i->xr_wf_size * fr_j->xr_wf_size * sizeof(double));

	double lmo_s[fr_i->n_lmo * fr_j->n_lmo];
	double lmo_t[fr_i->n_lmo * fr_j->n_lmo];

	/* compute S and T integrals */
	efp_st_int(fr_i->n_xr_shells, fr_i->xr_shells,
		   fr_j->n_xr_shells, fr_j->xr_shells,
		   fr_j->xr_wf_size, s, t);

	transform_integrals(fr_i->n_lmo, fr_j->n_lmo,
			    fr_i->xr_wf_size, fr_j->xr_wf_size,
			    fr_i->xr_wf, fr_j->xr_wf,
			    s, lmo_s);
	transform_integrals(fr_i->n_lmo, fr_j->n_lmo,
			    fr_i->xr_wf_size, fr_j->xr_wf_size,
			    fr_i->xr_wf, fr_j->xr_wf,
			    t, lmo_t);

	double exr = 0.0;
	double ecp = 0.0;

	for (int i = 0, idx = overlap_idx; i < fr_i->n_lmo; i++) {
		for (int j = 0; j < fr_j->n_lmo; j++, idx++) {
			double s_ij = lmo_s[i * fr_j->n_lmo + j];
			double t_ij = lmo_t[i * fr_j->n_lmo + j];
			double r_ij = vec_dist(fr_i->lmo_centroids + i,
					       fr_j->lmo_centroids + j);

			if ((efp->opts.terms & EFP_TERM_DISP) &&
			    (efp->opts.disp_damp == EFP_DISP_DAMP_OVERLAP))
				fr_i->overlap_int[idx] = s_ij;

			if ((efp->opts.terms & EFP_TERM_ELEC) &&
			    (efp->opts.elec_damp == EFP_ELEC_DAMP_OVERLAP))
				ecp += charge_penetration_energy(s_ij, r_ij);

			/* xr - first part */
			if (fabs(s_ij) > 1.0e-6)
				exr += -2.0 * sqrt(-2.0 * log(fabs(s_ij)) / PI)
							* s_ij * s_ij / r_ij;

			/* xr - second part */
			for (int k = 0; k < fr_i->n_lmo; k++)
				exr -= s_ij * lmo_s[k * fr_j->n_lmo + j] *
					fr_i->xr_fock_mat[fock_idx(i, k)];
			for (int l = 0; l < fr_j->n_lmo; l++)
				exr -= s_ij * lmo_s[i * fr_j->n_lmo + l] *
					fr_j->xr_fock_mat[fock_idx(j, l)];
			exr += 2.0 * s_ij * t_ij;

			/* xr - third part */
			for (int jj = 0; jj < fr_j->n_atoms; jj++) {
				struct efp_atom *at = fr_j->atoms + jj;
				double r = vec_dist(fr_i->lmo_centroids + i,
							VEC(at->x));
				exr -= s_ij * s_ij * valence(at->znuc) / r;
			}
			for (int l = 0; l < fr_j->n_lmo; l++) {
				double r = vec_dist(fr_i->lmo_centroids + i,
						    fr_j->lmo_centroids + l);
				exr += 2.0 * s_ij * s_ij / r;
			}
			for (int ii = 0; ii < fr_i->n_atoms; ii++) {
				struct efp_atom *at = fr_i->atoms + ii;
				double r = vec_dist(fr_j->lmo_centroids + j,
							VEC(at->x));
				exr -= s_ij * s_ij * valence(at->znuc) / r;
			}
			for (int k = 0; k < fr_i->n_lmo; k++) {
				double r = vec_dist(fr_i->lmo_centroids + k,
						    fr_j->lmo_centroids + j);
				exr += 2.0 * s_ij * s_ij / r;
			}
			exr -= s_ij * s_ij / r_ij;
		}
	}

	*exr_out = 2.0 * exr;
	*ecp_out = ecp;

	if (efp->do_gradient)
		frag_frag_xr_grad(efp, lmo_s, lmo_t, frag_i, frag_j, overlap_idx);

	free(s), free(t);
}

enum efp_result
efp_compute_xr(struct efp *efp)
{
	if (!(efp->opts.terms & EFP_TERM_XR))
		return EFP_RESULT_SUCCESS;

	double exr = 0.0;
	double ecp = 0.0;

	#pragma omp parallel for schedule(dynamic, 4) reduction(+:exr,ecp)
	for (int i = 0; i < efp->n_frag; i++) {
		for (int j = i + 1, idx = 0; j < efp->n_frag; j++) {
			double exr_out, ecp_out;

			frag_frag_xr(efp, i, j, idx, &exr_out, &ecp_out);

			exr += exr_out;
			ecp += ecp_out;

			idx += efp->frags[i].n_lmo * efp->frags[j].n_lmo;
		}
	}

	efp->energy.exchange_repulsion = exr;
	efp->energy.charge_penetration = ecp;

	return EFP_RESULT_SUCCESS;
}

static inline int
func_d_idx(int a, int b)
{
	/* order in which GAMESS stores D functions */
	enum { xx = 0, yy, zz, xy, xz, yz };

	static const int idx[] = {
		xx, xy, xz, xy, yy, yz, xz, yz, zz
	};

	return idx[3 * a + b];
}

static inline int
func_f_idx(int a, int b, int c)
{
	/* order in which GAMESS stores F functions */
	enum { xxx = 0, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz };

	static const int idx[] = {
		xxx, xxy, xxz, xxy, xyy, xyz, xxz, xyz, xzz,
		xxy, xyy, xyz, xyy, yyy, yyz, xyz, yyz, yzz,
		xxz, xyz, xzz, xyz, yyz, yzz, xzz, yzz, zzz
	};

	return idx[9 * a + 3 * b + c];
}

static void
rotate_func_d(const mat_t *rotmat, const double *in, double *out)
{
	const double norm = sqrt(3.0) / 2.0;

	double full_in[9], full_out[9];

	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			full_in[3 * a + b] = in[func_d_idx(a, b)];
			if (a != b)
				full_in[3 * a + b] *= norm;
		}
	}

	rotate_t2(rotmat, full_in, full_out);

	for (int a = 0; a < 3; a++) {
		for (int b = 0; b < 3; b++) {
			if (a != b)
				full_out[3 * a + b] /= norm;
			out[func_d_idx(a, b)] = full_out[3 * a + b];
		}
	}
}

static void
rotate_func_f(const mat_t *rotmat, const double *in, double *out)
{
	const double norm1 = sqrt(5.0) / 3.0;
	const double norm2 = sqrt(3.0) / 2.0;

	double full_in[27], full_out[27];

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			for (int c = 0; c < 3; c++) {
				int full_idx = 9 * a + 3 * b + c;
				int in_idx = func_f_idx(a, b, c);

				full_in[full_idx] = in[in_idx];

				if (a != b || a != c)
					full_in[full_idx] *= norm1;

				if (a != b && a != c && b != c)
					full_in[full_idx] *= norm2;
			}

	rotate_t3(rotmat, full_in, full_out);

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			for (int c = 0; c < 3; c++) {
				int full_idx = 9 * a + 3 * b + c;
				int out_idx = func_f_idx(a, b, c);

				if (a != b || a != c)
					full_out[full_idx] /= norm1;

				if (a != b && a != c && b != c)
					full_out[full_idx] /= norm2;

				out[out_idx] = full_out[full_idx];
			}
}

void
efp_update_xr(struct frag *frag)
{
	const mat_t *rotmat = &frag->rotmat;

	/* update LMO centroids */
	for (int i = 0; i < frag->n_lmo; i++)
		move_pt(VEC(frag->x), rotmat, VEC(frag->lib->x),
			frag->lib->lmo_centroids + i, frag->lmo_centroids + i);

	/* update shells */
	for (int i = 0; i < frag->n_xr_shells; i++)
		move_pt(VEC(frag->x), rotmat, VEC(frag->lib->x),
			VEC(frag->lib->xr_shells[i].x),
			VEC(frag->xr_shells[i].x));

	/* rotate wavefunction */
	for (int lmo = 0; lmo < frag->n_lmo; lmo++) {
		const double *in = frag->lib->xr_wf + lmo * frag->xr_wf_size;
		double *out = frag->xr_wf + lmo * frag->xr_wf_size;

		for (int i = 0, func = 0; i < frag->n_xr_shells; i++) {
			switch (frag->xr_shells[i].type) {
			case 'S':
				func++;
				break;
			case 'L':
				func++;
				/* fall through */
			case 'P':
				mat_vec(rotmat, (const vec_t *)(in + func),
						(vec_t *)(out + func));
				func += 3;
				break;
			case 'D':
				rotate_func_d(rotmat, in + func, out + func);
				func += 6;
				break;
			case 'F':
				rotate_func_f(rotmat, in + func, out + func);
				func += 10;
				break;
			}
		}
	}
}
