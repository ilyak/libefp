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

static inline double
get_disp_damp_overlap(double s_ij)
{
	double ln_s = log(fabs(s_ij));
	return 1.0 - s_ij * s_ij * (1.0 - 2.0 * ln_s + 2.0 * ln_s * ln_s);
}

static void
calc_disp_damp_overlap(struct efp *efp, int frag_i, int frag_j,
		       int i, int j, double s_ij)
{
	if (!efp->disp_damp_overlap)
		return;

	int idx = disp_damp_overlap_idx(efp, frag_i, frag_j, i, j);
	efp->disp_damp_overlap[idx] =
		fabs(s_ij) > 1.0e-6 ? get_disp_damp_overlap(s_ij) : 0.0;
}

static inline double
get_charge_pen(double s_ij, double r_ij)
{
	double ln_s = log(fabs(s_ij));
	return -2.0 * s_ij * s_ij / r_ij / sqrt(-2.0 * ln_s);
}

static inline int
get_block_frag_count(struct efp *efp, int i)
{
	return efp->xr_block_frag_offset[i + 1] - efp->xr_block_frag_offset[i];
}

static double
compute_xr_frag(struct efp *efp, int frag_i, int frag_j, int stride,
		const double *s, const double *t)
{
	struct frag *fr_i = efp->frags + frag_i;
	struct frag *fr_j = efp->frags + frag_j;

	double lmo_s[fr_i->n_lmo][fr_j->n_lmo];
	double lmo_t[fr_i->n_lmo][fr_j->n_lmo];
	double tmp[fr_i->n_lmo * fr_j->xr_wf_size];

	/* lmo_s = wf_i * s * wf_j(t) */
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		fr_i->n_lmo, fr_j->xr_wf_size, fr_i->xr_wf_size, 1.0,
		fr_i->xr_wf, fr_i->xr_wf_size, s, stride, 0.0,
		tmp, fr_i->n_lmo);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
		fr_i->n_lmo, fr_j->n_lmo, fr_j->xr_wf_size, 1.0,
		tmp, fr_i->n_lmo, fr_j->xr_wf, fr_j->xr_wf_size, 0.0,
		&lmo_s[0][0], fr_i->n_lmo);

	/* lmo_t = wf_i * t * wf_j(t) */
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
		fr_i->n_lmo, fr_j->xr_wf_size, fr_i->xr_wf_size, 1.0,
		fr_i->xr_wf, fr_i->xr_wf_size, t, stride, 0.0,
		tmp, fr_i->n_lmo);
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
		fr_i->n_lmo, fr_j->n_lmo, fr_j->xr_wf_size, 1.0,
		tmp, fr_i->n_lmo, fr_j->xr_wf, fr_j->xr_wf_size, 0.0,
		&lmo_t[0][0], fr_i->n_lmo);

	double energy = 0.0;

	for (int i = 0; i < fr_i->n_lmo; i++) {
		for (int j = 0; j < fr_j->n_lmo; j++) {
			double s_ij = lmo_s[i][j];
			double t_ij = lmo_t[i][j];
			double r_ij = vec_dist(fr_i->lmo_centroids + i,
					       fr_j->lmo_centroids + j);

			calc_disp_damp_overlap(efp, frag_i, frag_j, i, j, s_ij);
			efp->charge_pen_energy += get_charge_pen(s_ij, r_ij);

			/* xr - first part */
			if (fabs(s_ij) > 1.0e-6)
				energy += -2.0 * sqrt(-2.0 * log(fabs(s_ij)) /
						PI) * s_ij * s_ij / r_ij;

			/* xr - second part */
			for (int k = 0; k < fr_i->n_lmo; k++)
				energy -= s_ij * lmo_s[k][j] *
					fr_i->xr_fock_mat[fock_idx(i, k)];
			for (int l = 0; l < fr_j->n_lmo; l++)
				energy -= s_ij * lmo_s[i][l] *
					fr_j->xr_fock_mat[fock_idx(j, l)];
			energy += 2.0 * s_ij * t_ij;

			/* xr - third part */
			for (int jj = 0; jj < fr_j->n_atoms; jj++) {
				struct efp_atom *at = fr_j->atoms + jj;
				double r = vec_dist(fr_i->lmo_centroids + i,
							VEC(at->x));
				energy -= s_ij * s_ij * valence(at->znuc) / r;
			}
			for (int l = 0; l < fr_j->n_lmo; l++) {
				double r = vec_dist(fr_i->lmo_centroids + i,
						    fr_j->lmo_centroids + l);
				energy += 2.0 * s_ij * s_ij / r;
			}
			for (int ii = 0; ii < fr_i->n_atoms; ii++) {
				struct efp_atom *at = fr_i->atoms + ii;
				double r = vec_dist(fr_j->lmo_centroids + j,
							VEC(at->x));
				energy -= s_ij * s_ij * valence(at->znuc) / r;
			}
			for (int k = 0; k < fr_i->n_lmo; k++) {
				double r = vec_dist(fr_i->lmo_centroids + k,
						    fr_j->lmo_centroids + j);
				energy += 2.0 * s_ij * s_ij / r;
			}
			energy -= s_ij * s_ij / r_ij;
		}
	}
	energy *= 2.0;
	return energy;
}

static enum efp_result
compute_xr_block(struct efp *efp, int block_i, int block_j, double *energy)
{
	enum efp_result res;

	struct frag *frags_i = efp->frags + efp->xr_block_frag_offset[block_i];
	struct frag *frags_j = efp->frags + efp->xr_block_frag_offset[block_j];
	int n_block_i = get_block_frag_count(efp, block_i);
	int n_block_j = get_block_frag_count(efp, block_j);

	struct efp_xr_block block;
	memset(&block, 0, sizeof(struct efp_xr_block));

	for (int i = 0; i < n_block_i; i++) {
		block.n_atoms_i += frags_i[i].n_atoms;
		block.basis_size_i += frags_i[i].xr_wf_size;
	}

	for (int j = 0; j < n_block_j; j++) {
		block.n_atoms_j += frags_j[j].n_atoms;
		block.basis_size_j += frags_j[j].xr_wf_size;
	}

	block.atoms_i = malloc(block.n_atoms_i * sizeof(struct efp_atom));
	block.atoms_j = malloc(block.n_atoms_j * sizeof(struct efp_atom));

	for (int i = 0, idx = 0; i < n_block_i; i++)
		for (int a = 0; a < frags_i[i].n_atoms; a++)
			block.atoms_i[idx++] = frags_i[i].atoms[a];

	for (int j = 0, idx = 0; j < n_block_j; j++)
		for (int a = 0; a < frags_j[j].n_atoms; a++)
			block.atoms_j[idx++] = frags_j[j].atoms[a];

	size_t size = block.basis_size_i * block.basis_size_j * sizeof(double);

	double *s = malloc(size), *t = malloc(size);
	double *sx = NULL, *tx = NULL;

	if (efp->grad) {
		sx = malloc(3 * size);
		tx = malloc(3 * size);
	}

	if (!s || !t || (efp->grad && (!sx || !tx))) {
		res = EFP_RESULT_NO_MEMORY;
		goto fail;
	}

	if ((res = efp->callbacks.get_overlap_integrals(&block, s, sx,
			efp->callbacks.get_overlap_integrals_user_data)))
		goto fail;

	if ((res = efp->callbacks.get_kinetic_integrals(&block, t, tx,
			efp->callbacks.get_kinetic_integrals_user_data)))
		goto fail;

	for (int i = 0; i < n_block_i; i++) {
		int j = 0, offset = i * block.basis_size_j;

		if (block_i == block_j)
			for (; j < i + 1; j++)
				offset += frags_j[j].xr_wf_size;

		for (; j < n_block_j; j++) {
			int frag_i = i + efp->xr_block_frag_offset[block_i];
			int frag_j = j + efp->xr_block_frag_offset[block_j];

			*energy += compute_xr_frag(efp, frag_i, frag_j,
						   block.basis_size_j,
						   s + offset, t + offset);

			offset += frags_j[j].xr_wf_size;
		}
	}

fail:
	free(block.atoms_i), free(block.atoms_j);
	free(s), free(sx), free(t), free(tx);
	return res;
}

enum efp_result
efp_compute_xr(struct efp *efp)
{
	if (efp->grad)
		return EFP_RESULT_NOT_IMPLEMENTED;

	enum efp_result res;
	double energy = 0.0;
	efp->charge_pen_energy = 0.0;

	/* Because of potentially huge number of fragments we can't just
	 * compute all overlap and kinetic energy integrals in one step.
	 * Instead we process fragments in blocks so that number of basis
	 * functions in each block is not greater than some reasonable number.
	 * Also see setup_xr function in efp.c for block setup details.
	 */
	for (int i = 0; i < efp->n_xr_blocks; i++)
		for (int j = i; j < efp->n_xr_blocks; j++)
			if ((res = compute_xr_block(efp, i, j, &energy)))
				return res;

	efp->energy[efp_get_term_index(EFP_TERM_XR)] = energy;
	return EFP_RESULT_SUCCESS;
}

static void
rotate_func_d(const struct mat *rotmat, const double *in, double *out)
{
	const double norm = sqrt(3.0) / 2.0;

	double normin[6];
	memcpy(normin, in, 6 * sizeof(double));

	normin[1] *= norm;
	normin[2] *= norm;
	normin[4] *= norm;

	rotate_t2(rotmat, normin, out);

	out[1] /= norm;
	out[2] /= norm;
	out[4] /= norm;
}

static void
rotate_func_f(const struct mat *rotmat, const double *in, double *out)
{
	const double norm1 = sqrt(5.0) / 3.0;
	const double norm2 = sqrt(3.0) / 2.0;

	double normin[10];
	memcpy(normin, in, 10 * sizeof(double));

	normin[1] *= norm1;
	normin[2] *= norm1;
	normin[3] *= norm1;
	normin[4] *= norm1 * norm2;
	normin[5] *= norm1;
	normin[7] *= norm1;
	normin[8] *= norm1;

	rotate_t3(rotmat, normin, out);

	out[1] /= norm1;
	out[2] /= norm1;
	out[3] /= norm1;
	out[4] /= norm1 * norm2;
	out[5] /= norm1;
	out[7] /= norm1;
	out[8] /= norm1;
}

void
efp_update_xr(struct frag *frag, const struct mat *rotmat)
{
	/* move LMO centroids */
	for (int i = 0; i < frag->n_lmo; i++)
		move_pt(VEC(frag->x), rotmat, frag->lib->lmo_centroids + i,
				frag->lmo_centroids + i);

	/* rotate wavefunction */
	for (int lmo = 0; lmo < frag->n_lmo; lmo++) {
		double *in = frag->lib->xr_wf + lmo * frag->xr_wf_size;
		double *out = frag->xr_wf + lmo * frag->xr_wf_size;

		for (int i = 0, func = 0; frag->shells[i]; i++) {
			switch (frag->shells[i]) {
			case 'S':
				func++;
				break;
			case 'L':
				func++;
				/* fall through */
			case 'P':
				rotate_t1(rotmat, in + func, out + func);
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
