/*-
 * Copyright (c) 2012-2017 Ilya Kaliman
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

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "balance.h"
#include "clapack.h"
#include "elec.h"
#include "private.h"
#include "stream.h"
#include "util.h"

static void
update_fragment(struct frag *frag)
{
	/* update atoms */
	for (size_t i = 0; i < frag->n_atoms; i++)
		efp_move_pt(CVEC(frag->x), &frag->rotmat,
			CVEC(frag->lib->atoms[i].x), VEC(frag->atoms[i].x));

	efp_update_elec(frag);
	efp_update_pol(frag);
	efp_update_disp(frag);
	efp_update_xr(frag);
}

static enum efp_result
set_coord_xyzabc(struct frag *frag, const double *coord)
{
	frag->x = coord[0];
	frag->y = coord[1];
	frag->z = coord[2];

	euler_to_matrix(coord[3], coord[4], coord[5], &frag->rotmat);
	update_fragment(frag);

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
set_coord_points(struct frag *frag, const double *coord)
{
	/* allow fragments with less than 3 atoms by using multipole points of
	 * ghost atoms; multipole points have the same coordinates as atoms */
	if (frag->n_multipole_pts < 3) {
		efp_log("fragment must contain at least three atoms");
		return EFP_RESULT_FATAL;
	}

	double ref[9] = {
		frag->lib->multipole_pts[0].x,
		frag->lib->multipole_pts[0].y,
		frag->lib->multipole_pts[0].z,
		frag->lib->multipole_pts[1].x,
		frag->lib->multipole_pts[1].y,
		frag->lib->multipole_pts[1].z,
		frag->lib->multipole_pts[2].x,
		frag->lib->multipole_pts[2].y,
		frag->lib->multipole_pts[2].z
	};

	vec_t p1;
	mat_t rot1, rot2;

	efp_points_to_matrix(coord, &rot1);
	efp_points_to_matrix(ref, &rot2);
	rot2 = mat_transpose(&rot2);
	frag->rotmat = mat_mat(&rot1, &rot2);
	p1 = mat_vec(&frag->rotmat, VEC(frag->lib->multipole_pts[0].x));

	/* center of mass */
	frag->x = coord[0] - p1.x;
	frag->y = coord[1] - p1.y;
	frag->z = coord[2] - p1.z;

	update_fragment(frag);

	return EFP_RESULT_SUCCESS;
}

static enum efp_result
set_coord_rotmat(struct frag *frag, const double *coord)
{
	if (!efp_check_rotation_matrix((const mat_t *)(coord + 3))) {
		efp_log("invalid rotation matrix specified");
		return EFP_RESULT_FATAL;
	}

	frag->x = coord[0];
	frag->y = coord[1];
	frag->z = coord[2];

	memcpy(&frag->rotmat, coord + 3, sizeof(frag->rotmat));
	update_fragment(frag);

	return EFP_RESULT_SUCCESS;
}

static void
free_frag(struct frag *frag)
{
	if (!frag)
		return;

	free(frag->atoms);
	free(frag->multipole_pts);
	free(frag->polarizable_pts);
	free(frag->dynamic_polarizable_pts);
    free(frag->dipquad_polarizable_pts);
	free(frag->lmo_centroids);
	free(frag->xr_fock_mat);
	free(frag->xr_wf);
	free(frag->xrfit);
	free(frag->screen_params);
	free(frag->ai_screen_params);

	for (size_t i = 0; i < 3; i++)
		free(frag->xr_wf_deriv[i]);

	for (size_t i = 0; i < frag->n_xr_atoms; i++) {
		for (size_t j = 0; j < frag->xr_atoms[i].n_shells; j++)
			free(frag->xr_atoms[i].shells[j].coef);
		free(frag->xr_atoms[i].shells);
	}

	free(frag->xr_atoms);

	/* don't do free(frag) here */
}

static enum efp_result
copy_frag(struct frag *dest, const struct frag *src)
{
	size_t size;

	memcpy(dest, src, sizeof(*dest));

	if (src->atoms) {
		size = src->n_atoms * sizeof(struct efp_atom);
		dest->atoms = (struct efp_atom *)malloc(size);
		if (!dest->atoms)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->atoms, src->atoms, size);
	}
	if (src->multipole_pts) {
		size = src->n_multipole_pts * sizeof(struct multipole_pt);
		dest->multipole_pts = (struct multipole_pt *)malloc(size);
		if (!dest->multipole_pts)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->multipole_pts, src->multipole_pts, size);
	}
	if (src->screen_params) {
		size = src->n_multipole_pts * sizeof(double);
		dest->screen_params = (double *)malloc(size);
		if (!dest->screen_params)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->screen_params, src->screen_params, size);
	}
	if (src->ai_screen_params) {
		size = src->n_multipole_pts * sizeof(double);
		dest->ai_screen_params = (double *)malloc(size);
		if (!dest->ai_screen_params)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->ai_screen_params, src->ai_screen_params, size);
	}
	if (src->polarizable_pts) {
		size = src->n_polarizable_pts * sizeof(struct polarizable_pt);
		dest->polarizable_pts = (struct polarizable_pt *)malloc(size);
		if (!dest->polarizable_pts)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->polarizable_pts, src->polarizable_pts, size);
	}
	if (src->dynamic_polarizable_pts) {
		size = src->n_dynamic_polarizable_pts *
				sizeof(struct dynamic_polarizable_pt);
		dest->dynamic_polarizable_pts =
				(struct dynamic_polarizable_pt *)malloc(size);
		if (!dest->dynamic_polarizable_pts)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->dynamic_polarizable_pts,
				src->dynamic_polarizable_pts, size);
	}
    if (src->dipquad_polarizable_pts) {
        size = src->n_dipquad_polarizable_pts *
               sizeof(struct dipquad_polarizable_pt);
        dest->dipquad_polarizable_pts =
                (struct dipquad_polarizable_pt *)malloc(size);
        if (!dest->dipquad_polarizable_pts)
            return EFP_RESULT_NO_MEMORY;
        memcpy(dest->dipquad_polarizable_pts,
               src->dipquad_polarizable_pts, size);
    }
	if (src->lmo_centroids) {
		size = src->n_lmo * sizeof(vec_t);
		dest->lmo_centroids = (vec_t *)malloc(size);
		if (!dest->lmo_centroids)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->lmo_centroids, src->lmo_centroids, size);
	}
	if (src->xr_atoms) {
		size = src->n_xr_atoms * sizeof(struct xr_atom);
		dest->xr_atoms = (struct xr_atom *)malloc(size);
		if (!dest->xr_atoms)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xr_atoms, src->xr_atoms, size);

		for (size_t j = 0; j < src->n_xr_atoms; j++) {
			const struct xr_atom *at_src = src->xr_atoms + j;
			struct xr_atom *at_dest = dest->xr_atoms + j;

			size = at_src->n_shells * sizeof(struct shell);
			at_dest->shells = (struct shell *)malloc(size);
			if (!at_dest->shells)
				return EFP_RESULT_NO_MEMORY;
			memcpy(at_dest->shells, at_src->shells, size);

			for (size_t i = 0; i < at_src->n_shells; i++) {
				size = (at_src->shells[i].type == 'L' ? 3 : 2) *
				    at_src->shells[i].n_funcs * sizeof(double);
				at_dest->shells[i].coef =
				    (double *)malloc(size);
				if (!at_dest->shells[i].coef)
					return EFP_RESULT_NO_MEMORY;
				memcpy(at_dest->shells[i].coef,
				    at_src->shells[i].coef, size);
			}
		}
	}
	if (src->xr_fock_mat) {
		size = src->n_lmo * (src->n_lmo + 1) / 2 * sizeof(double);
		dest->xr_fock_mat = (double *)malloc(size);
		if (!dest->xr_fock_mat)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xr_fock_mat, src->xr_fock_mat, size);
	}
	if (src->xr_wf) {
		size = src->n_lmo * src->xr_wf_size * sizeof(double);
		dest->xr_wf = (double *)malloc(size);
		if (!dest->xr_wf)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xr_wf, src->xr_wf, size);
	}
	if (src->xrfit) {
		size = src->n_lmo * 4 * sizeof(double);
		dest->xrfit = (double *)malloc(size);
		if (!dest->xrfit)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xrfit, src->xrfit, size);
	}
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
check_opts(const struct efp_opts *opts)
{
	if (opts->enable_pbc) {
		if ((opts->terms & EFP_TERM_AI_ELEC) ||
		    (opts->terms & EFP_TERM_AI_POL) ||
		    (opts->terms & EFP_TERM_AI_DISP) ||
		    (opts->terms & EFP_TERM_AI_XR) ||
		    (opts->terms & EFP_TERM_AI_CHTR)) {
			efp_log("periodic calculations are not supported "
			    "for QM/EFP");
			return EFP_RESULT_FATAL;
		}
		if (!opts->enable_cutoff) {
			efp_log("periodic calculations require interaction "
			    "cutoff to be enabled");
			return EFP_RESULT_FATAL;
		}
	}
	if (opts->enable_cutoff) {
		if (opts->swf_cutoff < 1.0) {
			efp_log("interaction cutoff is too small");
			return EFP_RESULT_FATAL;
		}
	}
    if (opts->enable_cutoff) {
        if (opts->swf_cutoff < opts->xr_cutoff) {
            efp_log("exchange-repulsion cutoff is smaller than interaction cutoff");
            return EFP_RESULT_FATAL;
        }
    }
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
check_frag_params(const struct efp_opts *opts, const struct frag *frag)
{
	if ((opts->terms & EFP_TERM_ELEC) || (opts->terms & EFP_TERM_AI_ELEC)) {
		if (!frag->multipole_pts) {
			efp_log("electrostatic parameters are missing");
			return EFP_RESULT_FATAL;
		}
		if (opts->elec_damp == EFP_ELEC_DAMP_SCREEN &&
		    frag->screen_params == NULL) {
			efp_log("screening parameters are missing");
			return EFP_RESULT_FATAL;
		}
	}
	if ((opts->terms & EFP_TERM_POL) || (opts->terms & EFP_TERM_AI_POL)) {
		if (!frag->polarizable_pts || !frag->multipole_pts) {
			efp_log("polarization parameters are missing");
			return EFP_RESULT_FATAL;
		}
	}
	if ((opts->terms & EFP_TERM_DISP) || (opts->terms & EFP_TERM_AI_DISP)) {
		if (frag->dynamic_polarizable_pts == NULL) {
			efp_log("dispersion parameters are missing");
			return EFP_RESULT_FATAL;
		}
		if (opts->disp_damp == EFP_DISP_DAMP_OVERLAP &&
		    frag->n_lmo != frag->n_dynamic_polarizable_pts) {
			efp_log("number of polarization points does not "
			    "match number of LMOs");
			return EFP_RESULT_FATAL;
		}
	}
	if ((opts->terms & EFP_TERM_XR) || (opts->terms & EFP_TERM_AI_XR)) {
		if (!frag->xr_atoms ||
		    !frag->xr_fock_mat ||
		    !frag->xr_wf ||
		    !frag->lmo_centroids) {
			efp_log("exchange repulsion parameters are missing");
			return EFP_RESULT_FATAL;
		}
	}
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
check_params(struct efp *efp)
{
	enum efp_result res;

	for (size_t i = 0; i < efp->n_frag; i++)
		if ((res = check_frag_params(&efp->opts, efp->frags + i)))
			return res;

	return EFP_RESULT_SUCCESS;
}

static int
do_elec(const struct efp_opts *opts)
{
	return opts->terms & EFP_TERM_ELEC;
}

static int
do_disp(const struct efp_opts *opts)
{
	return opts->terms & EFP_TERM_DISP;
}

static int
do_xr(const struct efp_opts *opts)
{
	int xr = (opts->terms & EFP_TERM_XR);
	int cp = (opts->terms & EFP_TERM_ELEC) &&
		 (opts->elec_damp == EFP_ELEC_DAMP_OVERLAP);
	int dd = (opts->terms & EFP_TERM_DISP) &&
		 (opts->disp_damp == EFP_DISP_DAMP_OVERLAP);
	return xr || cp || dd;
}

static void
compute_two_body_range(struct efp *efp, size_t frag_from, size_t frag_to,
    void *data)
{
	double e_elec = 0.0, e_disp = 0.0, e_xr = 0.0, e_cp = 0.0, e_elec_tmp = 0.0, e_disp_tmp = 0.0;

	(void)data;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) reduction(+:e_elec,e_disp,e_xr,e_cp)
#endif
	for (size_t i = frag_from; i < frag_to; i++) {
		size_t cnt = efp->n_frag % 2 ? (efp->n_frag - 1) / 2 :
		    i < efp->n_frag / 2 ? efp->n_frag / 2 :
		    efp->n_frag / 2 - 1;

		for (size_t j = i + 1; j < i + 1 + cnt; j++) {
			size_t fr_j = j % efp->n_frag;

			if (!efp_skip_frag_pair(efp, i, fr_j)) {
				double *s;
				six_t *ds;
				size_t n_lmo_ij = efp->frags[i].n_lmo *
				    efp->frags[fr_j].n_lmo;

				s = (double *)calloc(n_lmo_ij, sizeof(double));
				ds = (six_t *)calloc(n_lmo_ij, sizeof(six_t));

				if (do_xr(&efp->opts)) {
					double exr, ecp;

					efp_frag_frag_xr(efp, i, fr_j,
					    s, ds, &exr, &ecp);
					e_xr += exr;
					e_cp += ecp;

					/* */
					if (efp->opts.enable_pairwise) {
                        if (i == efp->opts.ligand) {
                            efp->pair_energies[fr_j].exchange_repulsion = exr;
                            efp->pair_energies[fr_j].charge_penetration = ecp;
                        }
                        if (fr_j == efp->opts.ligand) {
                            efp->pair_energies[i].exchange_repulsion = exr;
                            efp->pair_energies[i].charge_penetration = ecp;
                        }
                    }
				}
				if (do_elec(&efp->opts)) {
					e_elec_tmp = efp_frag_frag_elec(efp,
					    i, fr_j);
					e_elec += e_elec_tmp;
					/* */
					if (efp->opts.enable_pairwise) {
                        if (i == efp->opts.ligand) 
                            efp->pair_energies[fr_j].electrostatic = e_elec_tmp;                       
                        if (fr_j == efp->opts.ligand) 
                            efp->pair_energies[i].electrostatic = e_elec_tmp;
                    }
				}
				if (do_disp(&efp->opts)) {
					e_disp_tmp = efp_frag_frag_disp(efp,
					    i, fr_j, s, ds);
					e_disp += e_disp_tmp;
					/* */
					if (efp->opts.enable_pairwise) {
                        if (i == efp->opts.ligand) 
                            efp->pair_energies[fr_j].dispersion = e_disp_tmp;                       
                        if (fr_j == efp->opts.ligand) 
                            efp->pair_energies[i].dispersion = e_disp_tmp;
                    }
				}
				free(s);
				free(ds);
			}
		}
	}
	efp->energy.electrostatic += e_elec;
	efp->energy.dispersion += e_disp;
	efp->energy.exchange_repulsion += e_xr;
	efp->energy.charge_penetration += e_cp;
}

EFP_EXPORT enum efp_result
compute_two_body_crystal(struct efp *efp)
{
    double e_elec = 0.0, e_disp = 0.0, e_xr = 0.0, e_cp = 0.0, e_elec_tmp = 0.0, e_disp_tmp = 0.0;

// no parallelization
    int nsymm = efp->nsymm_frag;
    size_t *unique_frag = (size_t *)calloc(nsymm, sizeof(size_t));
    unique_symm_frag(efp, unique_frag);

    size_t *nsymm_frag = (size_t *)calloc(nsymm, sizeof(size_t));
    n_symm_frag(efp, nsymm_frag);

    for (size_t k = 0; k < nsymm; k++) {
        size_t i = unique_frag[k];
        struct frag *frag = efp->frags + i;

        // scaling factor that tells how many fragments like this are in the system
        size_t factor = nsymm_frag[k];

        for (size_t fr_j=0; fr_j<efp->n_frag; fr_j++){
            if ( fr_j != i && !efp_skip_frag_pair(efp, i, fr_j)) {

                struct frag *fragj = efp->frags + fr_j;
                double *s;
                six_t *ds;
                size_t n_lmo_ij = efp->frags[i].n_lmo *
                                  efp->frags[fr_j].n_lmo;

                s = (double *)calloc(n_lmo_ij, sizeof(double));
                ds = (six_t *)calloc(n_lmo_ij, sizeof(six_t));

                if (do_xr(&efp->opts)) {
                    double exr, ecp;

                    efp_frag_frag_xr(efp, i, fr_j,
                                     s, ds, &exr, &ecp);
                    e_xr += exr * factor;
                    e_cp += ecp * factor;

                    /* */
                    if (efp->opts.enable_pairwise) {
                        if (i == efp->opts.ligand) {
                            efp->pair_energies[fr_j].exchange_repulsion = exr;
                            efp->pair_energies[fr_j].charge_penetration = ecp;
                        }
                        if (fr_j == efp->opts.ligand) {
                            efp->pair_energies[i].exchange_repulsion = exr;
                            efp->pair_energies[i].charge_penetration = ecp;
                        }
                    }
                }
                if (do_elec(&efp->opts)) {
                    e_elec_tmp = efp_frag_frag_elec(efp,
                                                    i, fr_j);
                    e_elec += e_elec_tmp * factor;

                    /* */
                    if (efp->opts.enable_pairwise) {
                        if (i == efp->opts.ligand)
                            efp->pair_energies[fr_j].electrostatic = e_elec_tmp;
                        if (fr_j == efp->opts.ligand)
                            efp->pair_energies[i].electrostatic = e_elec_tmp;
                    }
                }
                if (do_disp(&efp->opts)) {
                    e_disp_tmp = efp_frag_frag_disp(efp,
                                                    i, fr_j, s, ds);
                    e_disp += e_disp_tmp * factor;
                    /* */
                    if (efp->opts.enable_pairwise) {
                        if (i == efp->opts.ligand)
                            efp->pair_energies[fr_j].dispersion = e_disp_tmp;
                        if (fr_j == efp->opts.ligand)
                            efp->pair_energies[i].dispersion = e_disp_tmp;
                    }
                }
                free(s);
                free(ds);

            }
        }
    }
    // really, we counted all pairwise interactions twice. Scaling back
    efp->energy.electrostatic += e_elec/2;
    efp->energy.dispersion += e_disp/2;
    efp->energy.exchange_repulsion += e_xr/2;
    efp->energy.charge_penetration += e_cp/2;

    return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_energy(struct efp *efp, struct efp_energy *energy)
{
	assert(efp);
	assert(energy);

	*energy = efp->energy;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_gradient(struct efp *efp, double *grad)
{
	assert(efp);
	assert(grad);

	if (!efp->do_gradient) {
		efp_log("gradient calculation was not requested");
		return EFP_RESULT_FATAL;
	}
	memcpy(grad, efp->grad, efp->n_frag * sizeof(six_t));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_atomic_gradient(struct efp *efp, double *grad)
{
	six_t *efpgrad = NULL; /* Calculated EFP gradient */
	vec_t *pgrad; /* Conversion of grad to vec_t type */
	size_t i, j, k, l;
	size_t nr; /* Number of atoms in the current fragment */
	size_t maxa; /* Maximum number of size of m, Ia, r arrays */
	vec_t *r = NULL; /* Radius-vector of each atom inside current fragment
			    with respect to CoM of that fragment */
	double mm, *m = NULL; /* Total Mass of fragments and individual atoms */
	double I, *Ia = NULL; /* Inertia along axis and contribution of each
				 individual atom */
	mat_t Id; /* Total inertia tensor of a fragment */
	vec_t v, g; /* Principal axis and Inertia along that axis */
	vec_t rbuf, rbuf2, tq, ri, rt;
	double dist, sina, ft, norm;
	enum efp_result res;

	assert(efp);
	assert(grad);

	if (!efp->do_gradient) {
		efp_log("gradient calculation was not requested");
		return EFP_RESULT_FATAL;
	}
	pgrad = (vec_t *)grad;

	/* Calculate maximum size of a fragment */
	maxa = 0;
	for (j = 0; j < efp->n_frag; j++) {
		if (efp->frags[j].n_atoms > maxa)
			maxa = efp->frags[j].n_atoms;
	}

	res = EFP_RESULT_NO_MEMORY;
	/* Create and initialize some arrays for work */
	if ((r = (vec_t *)malloc(maxa * sizeof(*r))) == NULL)
		goto error;
	if ((m = (double *)malloc(maxa * sizeof(*m))) == NULL)
		goto error;
	if ((Ia = (double *)malloc(maxa * sizeof(*Ia))) == NULL)
		goto error;

	/* Copy computed efp->grad */
	if ((efpgrad = (six_t *)malloc(efp->n_frag * sizeof(*efpgrad))) == NULL)
		goto error;
	memcpy(efpgrad, efp->grad, efp->n_frag * sizeof(*efpgrad));

	/* Main cycle (iterate fragments, distribute forces and torques) */
	for (k = 0, j = 0; j < efp->n_frag; j++) {
		nr = efp->frags[j].n_atoms;
		memset(r, 0, maxa * sizeof(*r));
		memset(m, 0, maxa * sizeof(*m));
		memset(Ia, 0, maxa * sizeof(*Ia));
		mm = 0.0;
		Id = mat_zero;
		v = vec_zero;
		g = vec_zero;

		for (i = 0; i < nr ; i++) {
			r[i].x = efp->frags[j].atoms[i].x - efp->frags[j].x;
			r[i].y = efp->frags[j].atoms[i].y - efp->frags[j].y;
			r[i].z = efp->frags[j].atoms[i].z - efp->frags[j].z;
			m[i] = efp->frags[j].atoms[i].mass;
			mm += m[i];

			/* Inertia tensor contribution calculations */
			Id.xx += m[i] * (r[i].y*r[i].y + r[i].z*r[i].z);
			Id.yy += m[i] * (r[i].x*r[i].x + r[i].z*r[i].z);
			Id.zz += m[i] * (r[i].x*r[i].x + r[i].y*r[i].y);
			Id.xy -= m[i] * r[i].x * r[i].y;
			Id.yx -= m[i] * r[i].x * r[i].y;
			Id.xz -= m[i] * r[i].x * r[i].z;
			Id.zx -= m[i] * r[i].x * r[i].z;
			Id.yz -= m[i] * r[i].y * r[i].z;
			Id.zy -= m[i] * r[i].y * r[i].z;
		}

		/* Try to diagonalize Id and get principal axis */
		if (efp_dsyev('V', 'U', 3, (double *)&Id, 3, (double *)&g)) {
			efp_log("inertia tensor diagonalization failed");
			res = EFP_RESULT_FATAL;
			goto error;
		}

		/* Add any additional forces from grad array to efpgrad array */
		for (i = 0; i < nr; i++) {
			efpgrad[j].x += pgrad[k+i].x;
			efpgrad[j].y += pgrad[k+i].y;
			efpgrad[j].z += pgrad[k+i].z;
			rbuf = vec_cross(&r[i], &pgrad[k+i]);
			efpgrad[j].a += rbuf.x;
			efpgrad[j].b += rbuf.y;
			efpgrad[j].c += rbuf.z;
			pgrad[k+i] = vec_zero;
		}

		/* Now we are ready to redistribute efpgrad over the atoms */

		/* Redistribute total translation grad[i]=m[i]/mm*efpgrad[j] */
		for (i = 0; i < nr; i++) {
			pgrad[k+i].x = efpgrad[j].x;
			pgrad[k+i].y = efpgrad[j].y;
			pgrad[k+i].z = efpgrad[j].z;
			vec_scale(&pgrad[k+i], m[i] / mm);
		}

		/* Redistribution of torque should be done over 3 principal
		   axes computed previously */
		for (l = 0; l < 3; l++) {
			v = ((vec_t *)&Id)[l];
			tq.x = efpgrad[j].a;
			tq.y = efpgrad[j].b;
			tq.z = efpgrad[j].c;

			/* Calculate contribution of each atom to moment of
			   inertia with respect to current axis */
			I = 0.0;
			for (i = 0; i < nr; i++) {
				rbuf = vec_cross(&v, &r[i]);
				dist = vec_len(&rbuf);
				Ia[i] = m[i] * dist * dist;
				I += Ia[i];
			}

			/* Project torque onto v axis */
			norm = vec_dot(&tq, &v);
			tq = v;
			vec_scale(&tq, norm);

			/* Now distribute torque using Ia[i]/I as a scale */
			for (i = 0; i < nr; i++) {
				if (eq(Ia[i], 0.0))
					continue;
				/* If atom is not on the current axis */
				rbuf = tq;
				vec_scale(&rbuf, Ia[i]/I);
				ft = vec_len(&rbuf);
				ri = r[i];
				vec_normalize(&ri);
				rt = tq;
				vec_normalize(&rt);
				rbuf2 = vec_cross(&rt, &ri);
				sina = vec_len(&rbuf2);
				vec_normalize(&rbuf2);
				vec_scale(&rbuf2, ft/sina/vec_len(&r[i]));
				/* Update grad with torque contribution of
				   atom i over axis v */
				pgrad[k+i] = vec_add(&pgrad[k+i], &rbuf2);
			}
		}
		k += nr;
	}
	res = EFP_RESULT_SUCCESS;
error:
	free(r);
	free(m);
	free(Ia);
	free(efpgrad);
	return res;
}

EFP_EXPORT enum efp_result
efp_set_point_charges(struct efp *efp, size_t n_ptc, const double *ptc,
    const double *xyz)
{
	assert(efp);
	efp->n_ptc = n_ptc;

	if (n_ptc == 0) {
		free(efp->ptc);
		free(efp->ptc_xyz);
		free(efp->ptc_grad);
		efp->ptc = NULL;
		efp->ptc_xyz = NULL;
		efp->ptc_grad = NULL;
		return EFP_RESULT_SUCCESS;
	}

	assert(ptc);
	assert(xyz);

	efp->ptc = (double *)realloc(efp->ptc, n_ptc * sizeof(double));
	efp->ptc_xyz = (vec_t *)realloc(efp->ptc_xyz, n_ptc * sizeof(vec_t));
	efp->ptc_grad = (vec_t *)realloc(efp->ptc_grad, n_ptc * sizeof(vec_t));

	memcpy(efp->ptc, ptc, n_ptc * sizeof(double));
	memcpy(efp->ptc_xyz, xyz, n_ptc * sizeof(vec_t));
	memset(efp->ptc_grad, 0, n_ptc * sizeof(vec_t));

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_point_charge_count(struct efp *efp, size_t *n_ptc)
{
	assert(efp);
	assert(n_ptc);

	*n_ptc = efp->n_ptc;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_point_charge_gradient(struct efp *efp, double *grad)
{
	assert(efp);
	assert(grad);

	if (!efp->do_gradient) {
		efp_log("gradient calculation was not requested");
		return EFP_RESULT_FATAL;
	}
	memcpy(grad, efp->ptc_grad, efp->n_ptc * sizeof(vec_t));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_point_charge_coordinates(struct efp *efp, double *xyz)
{
	assert(efp);
	assert(xyz);

	memcpy(xyz, efp->ptc_xyz, efp->n_ptc * sizeof(vec_t));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_set_point_charge_coordinates(struct efp *efp, const double *xyz)
{
	assert(efp);
	assert(xyz);

	memcpy(efp->ptc_xyz, xyz, efp->n_ptc * sizeof(vec_t));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_point_charge_values(struct efp *efp, double *ptc)
{
	assert(efp);
	assert(ptc);

	memcpy(ptc, efp->ptc, efp->n_ptc * sizeof(double));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_set_point_charge_values(struct efp *efp, const double *ptc)
{
	assert(efp);
	assert(ptc);

	memcpy(efp->ptc, ptc, efp->n_ptc * sizeof(double));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_set_coordinates(struct efp *efp, enum efp_coord_type coord_type,
    const double *coord)
{
	assert(efp);
	assert(coord);

	size_t stride;
	enum efp_result res;

	switch (coord_type) {
	case EFP_COORD_TYPE_XYZABC:
		stride = 6;
		break;
	case EFP_COORD_TYPE_POINTS:
		stride = 9;
		break;
	case EFP_COORD_TYPE_ROTMAT:
		stride = 12;
		break;
	}

	for (size_t i = 0; i < efp->n_frag; i++, coord += stride)
		if ((res = efp_set_frag_coordinates(efp, i, coord_type, coord)))
			return res;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_set_frag_coordinates(struct efp *efp, size_t frag_idx,
    enum efp_coord_type coord_type, const double *coord)
{
	struct frag *frag;

	assert(efp);
	assert(coord);
	assert(frag_idx < efp->n_frag);

	frag = efp->frags + frag_idx;

	switch (coord_type) {
	case EFP_COORD_TYPE_XYZABC:
		return set_coord_xyzabc(frag, coord);
	case EFP_COORD_TYPE_POINTS:
		return set_coord_points(frag, coord);
	case EFP_COORD_TYPE_ROTMAT:
		return set_coord_rotmat(frag, coord);
	}
	assert(0);
}

EFP_EXPORT enum efp_result
efp_get_coordinates(struct efp *efp, double *xyzabc)
{
	assert(efp);
	assert(xyzabc);

	for (size_t i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		double a, b, c;
		matrix_to_euler(&frag->rotmat, &a, &b, &c);

		*xyzabc++ = frag->x;
		*xyzabc++ = frag->y;
		*xyzabc++ = frag->z;
		*xyzabc++ = a;
		*xyzabc++ = b;
		*xyzabc++ = c;
	}
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_xyzabc(struct efp *efp, size_t frag_idx, double *xyzabc)
{
	struct frag *frag;
	double a, b, c;

	assert(efp);
	assert(frag_idx < efp->n_frag);
	assert(xyzabc);

	frag = efp->frags + frag_idx;
	matrix_to_euler(&frag->rotmat, &a, &b, &c);

	xyzabc[0] = frag->x;
	xyzabc[1] = frag->y;
	xyzabc[2] = frag->z;
	xyzabc[3] = a;
	xyzabc[4] = b;
	xyzabc[5] = c;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_set_periodic_box(struct efp *efp, double x, double y, double z, double alpha, double beta, double gamma)
{
	assert(efp);
	if (alpha < 0.01) {
        // assigning default angles of 90.0 degrees
        //printf("\n assigning angles to 90.0 \n");
        alpha = 90.0;
        beta = 90.0;
        gamma = 90.0;
    }

	efp->box.x = x;
	efp->box.y = y;
	efp->box.z = z;
	efp->box.a = alpha;
	efp->box.b = beta;
	efp->box.c = gamma;

	double max_cut = max_cutoff(efp->box);
	double cut = efp->opts.swf_cutoff;
	if (cut > max_cut) {
	    printf("\n Maximum allowed cutoff is %lf ", max_cut*0.52917721092);
        printf("\n Given cutoff is %lf \n", cut*0.52917721092);
        efp_log("periodic box dimensions must be at least twice "
                "the switching function cutoff");
        return EFP_RESULT_FATAL;
	}

    return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_periodic_box(struct efp *efp, double *xyzabc)
{
	assert(efp);
	assert(xyzabc);

    xyzabc[0] = efp->box.x;
	xyzabc[1] = efp->box.y;
	xyzabc[2] = efp->box.z;
    xyzabc[3] = efp->box.a;
    xyzabc[4] = efp->box.b;
    xyzabc[5] = efp->box.c;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_stress_tensor(struct efp *efp, double *stress)
{
	assert(efp);
	assert(stress);

	if (!efp->do_gradient) {
		efp_log("gradient calculation was not requested");
		return EFP_RESULT_FATAL;
	}

	*(mat_t *)stress = efp->stress;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_ai_screen(struct efp *efp, size_t frag_idx, double *screen)
{
	const struct frag *frag;
	size_t size;

	assert(efp);
	assert(screen);
	assert(frag_idx < efp->n_frag);

	frag = &efp->frags[frag_idx];

	if (frag->ai_screen_params == NULL) {
		efp_log("no screening parameters found for %s", frag->name);
		return EFP_RESULT_FATAL;
	}

	size = frag->n_multipole_pts * sizeof(double);
	memcpy(screen, frag->ai_screen_params, size);

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_prepare(struct efp *efp)
{
	assert(efp);

	efp->n_polarizable_pts = 0;

	for (size_t i = 0; i < efp->n_frag; i++) {
		efp->frags[i].polarizable_offset = efp->n_polarizable_pts;
		efp->n_polarizable_pts += efp->frags[i].n_polarizable_pts;
	}

	efp->n_fragment_field_pts = 0;
	if (efp->opts.enable_pairwise && efp->opts.ligand != -1) {
		size_t ligand_idx = efp->opts.ligand;
		for (size_t i = 0; i < efp->n_frag; i++) {
			efp->frags[i].fragment_field_offset = efp->n_fragment_field_pts;
			efp->n_fragment_field_pts += efp->frags[ligand_idx].n_polarizable_pts;
		}
	}
	efp->fragment_field = (vec_t *)calloc(efp->n_fragment_field_pts, sizeof(vec_t));

	efp->indip = (vec_t *)calloc(efp->n_polarizable_pts, sizeof(vec_t));
	efp->indipconj = (vec_t *)calloc(efp->n_polarizable_pts, sizeof(vec_t));
	efp->grad = (six_t *)calloc(efp->n_frag, sizeof(six_t));
	efp->skiplist = (char *)calloc(efp->n_frag * efp->n_frag, 1);
	efp->pair_energies = (struct efp_energy *)calloc(efp->n_frag, sizeof(struct efp_energy));
    efp->symmlist = (size_t *)calloc(efp->n_frag, sizeof(size_t));

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_set_orbital_energies(struct efp *efp, size_t n_core, size_t n_act,
    size_t n_vir, const double *oe)
{
	size_t size;

	assert(efp);
	assert(oe);

	efp->n_ai_core = n_core;
	efp->n_ai_act = n_act;
	efp->n_ai_vir = n_vir;

	size = (n_core + n_act + n_vir) * sizeof(double);

	efp->ai_orbital_energies = (double *)realloc(efp->ai_orbital_energies,
	    size);
	memcpy(efp->ai_orbital_energies, oe, size);

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_set_dipole_integrals(struct efp *efp, size_t n_core, size_t n_act,
    size_t n_vir, const double *dipint)
{
	size_t size;

	assert(efp);
	assert(dipint);

	efp->n_ai_core = n_core;
	efp->n_ai_act = n_act;
	efp->n_ai_vir = n_vir;

	size = 3 * (n_core + n_act + n_vir) * (n_core + n_act + n_vir);
	efp->ai_dipole_integrals = (double *)realloc(efp->ai_dipole_integrals,
	    size * sizeof(double));
	memcpy(efp->ai_dipole_integrals, dipint, size * sizeof(double));

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_wavefunction_dependent_energy(struct efp *efp, double *energy)
{
	assert(efp);
	assert(energy);

	if (!(efp->opts.terms & EFP_TERM_POL) &&
	    !(efp->opts.terms & EFP_TERM_AI_POL)) {
		*energy = 0.0;
		return EFP_RESULT_SUCCESS;
	}
	return efp_compute_pol_energy(efp, energy);
}

EFP_EXPORT enum efp_result
efp_compute(struct efp *efp, int do_gradient)
{
	enum efp_result res;

	assert(efp);

	if (efp->grad == NULL) {
		efp_log("call efp_prepare after all fragments are added");
		return EFP_RESULT_FATAL;
	}

	efp->do_gradient = do_gradient;

	if ((res = check_params(efp)))
		return res;

	memset(&efp->energy, 0, sizeof(efp->energy));
	memset(&efp->stress, 0, sizeof(efp->stress));
	memset(efp->grad, 0, efp->n_frag * sizeof(six_t));
	memset(efp->ptc_grad, 0, efp->n_ptc * sizeof(vec_t));
	memset(efp->pair_energies, 0, efp->n_frag * sizeof(efp->energy));

	if (efp->opts.symmetry == 0) { // standard case
        efp_balance_work(efp, compute_two_body_range, NULL);
	}
	else {  // high-symmetry crystals
	    if (res = compute_two_body_crystal(efp))
	        return res;
	}

	if (res = efp_compute_pol(efp))
	    return res;
	if (res = efp_compute_ai_elec(efp))
	    return res;
	if (res = efp_compute_ai_disp(efp))
	    return res;

#ifdef EFP_USE_MPI
	efp_allreduce(&efp->energy.electrostatic, 1);
	efp_allreduce(&efp->energy.dispersion, 1);
	efp_allreduce(&efp->energy.exchange_repulsion, 1);
	efp_allreduce(&efp->energy.charge_penetration, 1);

	if (efp->do_gradient) {
		efp_allreduce((double *)efp->grad, 6 * efp->n_frag);
		efp_allreduce((double *)efp->ptc_grad, 3 * efp->n_ptc);
		efp_allreduce((double *)&efp->stress, 9);
	}
#endif
	efp->energy.total = efp->energy.electrostatic +
			    efp->energy.charge_penetration +
			    efp->energy.electrostatic_point_charges +
			    efp->energy.polarization +
			    efp->energy.dispersion +
			    efp->energy.ai_dispersion +
			    efp->energy.exchange_repulsion;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_charge(struct efp *efp, size_t frag_idx, double *charge)
{
	struct frag *frag;
	double sum = 0.0;
	size_t i;

	assert(efp);
	assert(charge);
	assert(frag_idx < efp->n_frag);

	frag = efp->frags + frag_idx;

	for (i = 0; i < frag->n_atoms; i++)
		sum += frag->atoms[i].znuc;
	for (i = 0; i < frag->n_multipole_pts; i++)
		sum += frag->multipole_pts[i].monopole;

	*charge = sum;
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_multiplicity(struct efp *efp, size_t frag_idx, int *mult)
{
	assert(efp);
	assert(mult);
	assert(frag_idx < efp->n_frag);

	*mult = efp->frags[frag_idx].multiplicity;
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_multipole_count(struct efp *efp, size_t frag_idx, size_t *n_mult)
{
	assert(efp);
	assert(n_mult);
	assert(frag_idx < efp->n_frag);

	*n_mult = efp->frags[frag_idx].n_multipole_pts;
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_mult_rank(struct efp *efp, size_t frag_idx, size_t mult_idx, size_t *rank)
{
    assert(efp);
    assert(rank);
    assert(frag_idx < efp->n_frag);
    assert(mult_idx < efp->frags[frag_idx].n_multipole_pts);

    struct frag *frag = efp->frags + frag_idx;
    struct multipole_pt *pt = frag->multipole_pts + mult_idx;

    *rank = -1;

    for (size_t t = 0; t < 10; t++)
        if (pt->octupole[t] != 0) {
            *rank = 3;
            return EFP_RESULT_SUCCESS;
        }

    for (size_t t = 0; t < 6; t++)
        if (pt->quadrupole[t] != 0) {
            *rank = 2;
            return EFP_RESULT_SUCCESS;
        }

    if (pt->dipole.x != 0 || pt->dipole.y != 0 || pt->dipole.z != 0) {
        *rank = 1;
        return EFP_RESULT_SUCCESS;
    }

    if (pt->monopole != 0) {
        *rank = 0;
        return EFP_RESULT_SUCCESS;
    }

    return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_multipole_count(struct efp *efp, size_t *n_mult)
{
	size_t sum = 0;

	assert(efp);
	assert(n_mult);

	for (size_t i = 0; i < efp->n_frag; i++)
		sum += efp->frags[i].n_multipole_pts;

	*n_mult = sum;
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_multipole_coordinates(struct efp *efp, double *xyz)
{
	assert(efp);
	assert(xyz);

	for (size_t i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_multipole_pts; j++) {
			*xyz++ = frag->multipole_pts[j].x;
			*xyz++ = frag->multipole_pts[j].y;
			*xyz++ = frag->multipole_pts[j].z;
		}
	}
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_multipole_values(struct efp *efp, double *mult)
{
	assert(efp);
	assert(mult);

	for (size_t i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_multipole_pts; j++) {
			struct multipole_pt *pt = frag->multipole_pts + j;

			*mult++ = pt->monopole;

			*mult++ = pt->dipole.x;
			*mult++ = pt->dipole.y;
			*mult++ = pt->dipole.z;

			for (size_t t = 0; t < 6; t++)
				*mult++ = pt->quadrupole[t];
			for (size_t t = 0; t < 10; t++)
				*mult++ = pt->octupole[t];
		}
	}
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_induced_dipole_count(struct efp *efp, size_t frag_idx, size_t *n_dip)
{
    assert(efp);
    assert(n_dip);
    assert(frag_idx < efp->n_frag);

    *n_dip = efp->frags[frag_idx].n_polarizable_pts;
    return EFP_RESULT_SUCCESS;
}


EFP_EXPORT enum efp_result
efp_get_induced_dipole_count(struct efp *efp, size_t *n_dip)
{
	size_t sum = 0;

	assert(efp);
	assert(n_dip);

	for (size_t i = 0; i < efp->n_frag; i++)
		sum += efp->frags[i].n_polarizable_pts;

	*n_dip = sum;
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_induced_dipole_coordinates(struct efp *efp, double *xyz)
{
	assert(efp);
	assert(xyz);

	for (size_t i = 0; i < efp->n_frag; i++) {
		struct frag *frag = efp->frags + i;

		for (size_t j = 0; j < frag->n_polarizable_pts; j++) {
			struct polarizable_pt *pt = frag->polarizable_pts + j;

			*xyz++ = pt->x;
			*xyz++ = pt->y;
			*xyz++ = pt->z;
		}
	}
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_induced_dipole_values(struct efp *efp, double *dip)
{
	assert(efp);
	assert(dip);

	memcpy(dip, efp->indip, efp->n_polarizable_pts * sizeof(vec_t));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_induced_dipole_conj_values(struct efp *efp, double *dip)
{
	assert(efp);
	assert(dip);

	memcpy(dip, efp->indipconj, efp->n_polarizable_pts * sizeof(vec_t));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_lmo_count(struct efp *efp, size_t frag_idx, size_t *n_lmo)
{
	assert(efp != NULL);
	assert(frag_idx < efp->n_frag);
	assert(n_lmo != NULL);

	*n_lmo = efp->frags[frag_idx].n_lmo;
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_lmo_coordinates(struct efp *efp, size_t frag_idx, double *xyz)
{
	struct frag *frag;

	assert(efp != NULL);
	assert(frag_idx < efp->n_frag);
	assert(xyz != NULL);

	frag = efp->frags + frag_idx;

	if (frag->lmo_centroids == NULL) {
		efp_log("no LMO centroids for fragment %s", frag->name);
		return EFP_RESULT_FATAL;
	}
	memcpy(xyz, frag->lmo_centroids, frag->n_lmo * sizeof(vec_t));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_xrfit(struct efp *efp, size_t frag_idx, double *xrfit)
{
	struct frag *frag;

	assert(efp != NULL);
	assert(frag_idx < efp->n_frag);
	assert(xrfit != NULL);

	frag = efp->frags + frag_idx;

	if (frag->xrfit == NULL) {
		efp_log("no XRFIT parameters for fragment %s", frag->name);
		return EFP_RESULT_FATAL;
	}
	memcpy(xrfit, frag->xrfit, frag->n_lmo * 4 * sizeof(double));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT void
efp_shutdown(struct efp *efp)
{
	if (efp == NULL)
		return;
	for (size_t i = 0; i < efp->n_frag; i++)
		free_frag(efp->frags + i);
	for (size_t i = 0; i < efp->n_lib; i++) {
		free_frag(efp->lib[i]);
		free(efp->lib[i]);
	}
	free(efp->frags);
	free(efp->lib);
	free(efp->grad);
	free(efp->ptc);
	free(efp->ptc_xyz);
	free(efp->ptc_grad);
	free(efp->indip);
	free(efp->indipconj);
	free(efp->ai_orbital_energies);
	free(efp->ai_dipole_integrals);
	free(efp->skiplist);
	free(efp->fragment_field);
	free(efp->pair_energies);
    free(efp->symmlist);
	free(efp);
}

EFP_EXPORT enum efp_result
efp_set_opts(struct efp *efp, const struct efp_opts *opts)
{
	enum efp_result res;

	assert(efp);
	assert(opts);

	if ((res = check_opts(opts)))
		return res;

	efp->opts = *opts;
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_opts(struct efp *efp, struct efp_opts *opts)
{
	assert(efp);
	assert(opts);

	*opts = efp->opts;
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT void
efp_opts_default(struct efp_opts *opts)
{
	assert(opts);

	memset(opts, 0, sizeof(*opts));
	opts->terms = EFP_TERM_ELEC | EFP_TERM_POL | EFP_TERM_DISP |
	    EFP_TERM_XR | EFP_TERM_AI_ELEC | EFP_TERM_AI_POL;
}

EFP_EXPORT void
efp_set_error_log(void (*cb)(const char *))
{
	efp_set_log_cb(cb);
}

EFP_EXPORT enum efp_result
efp_add_fragment(struct efp *efp, const char *name)
{
	const struct frag *lib;

	assert(efp);
	assert(name);

	if (efp->skiplist) {
		efp_log("cannot add fragments after efp_prepare");
		return EFP_RESULT_FATAL;
	}
	if ((lib = efp_find_lib(efp, name)) == NULL) {
		efp_log("cannot find \"%s\" in any of .efp files", name);
		return EFP_RESULT_UNKNOWN_FRAGMENT;
	}

	efp->n_frag++;
	efp->frags = (struct frag *)realloc(efp->frags,
	    efp->n_frag * sizeof(struct frag));
	if (efp->frags == NULL)
		return EFP_RESULT_NO_MEMORY;

	enum efp_result res;
	struct frag *frag = efp->frags + efp->n_frag - 1;

	if ((res = copy_frag(frag, lib)))
		return res;

	for (size_t a = 0; a < 3; a++) {
		size_t size = frag->xr_wf_size * frag->n_lmo;

		frag->xr_wf_deriv[a] = (double *)calloc(size, sizeof(double));
		if (frag->xr_wf_deriv[a] == NULL)
			return EFP_RESULT_NO_MEMORY;
	}
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_skip_fragments(struct efp *efp, size_t i, size_t j, int value)
{
	assert(efp);
	assert(efp->skiplist); /* call efp_prepare first */
	assert(i < efp->n_frag);
	assert(j < efp->n_frag);

	efp->skiplist[i * efp->n_frag + j] = value ? 1 : 0;
	efp->skiplist[j * efp->n_frag + i] = value ? 1 : 0;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT struct efp *
efp_create(void)
{
	struct efp *efp = (struct efp *)calloc(1, sizeof(struct efp));

	if (efp == NULL)
		return NULL;

	efp_opts_default(&efp->opts);

	return efp;
}

EFP_EXPORT enum efp_result
efp_set_electron_density_field_fn(struct efp *efp,
    efp_electron_density_field_fn fn)
{
	assert(efp);

	efp->get_electron_density_field = fn;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_set_electron_density_field_user_data(struct efp *efp, void *user_data)
{
	assert(efp);

	efp->get_electron_density_field_user_data = user_data;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_count(struct efp *efp, size_t *n_frag)
{
	assert(efp);
	assert(n_frag);

	*n_frag = efp->n_frag;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_name(struct efp *efp, size_t frag_idx, size_t size,
    char *frag_name)
{
	assert(efp);
	assert(frag_name);
	assert(frag_idx < efp->n_frag);

	strncpy(frag_name, efp->frags[frag_idx].name, size);

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_mass(struct efp *efp, size_t frag_idx, double *mass_out)
{
	assert(efp);
	assert(mass_out);
	assert(frag_idx < efp->n_frag);

	const struct frag *frag = efp->frags + frag_idx;
	double mass = 0.0;

	for (size_t i = 0; i < frag->n_atoms; i++)
		mass += frag->atoms[i].mass;

	*mass_out = mass;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_inertia(struct efp *efp, size_t frag_idx, double *inertia_out)
{
	assert(efp);
	assert(inertia_out);
	assert(frag_idx < efp->n_frag);

	/* center of mass is in origin and axes are principal axes of inertia */

	const struct frag *frag = efp->frags[frag_idx].lib;
	vec_t inertia = vec_zero;

	for (size_t i = 0; i < frag->n_atoms; i++) {
		const struct efp_atom *atom = frag->atoms + i;

		inertia.x += atom->mass * (atom->y*atom->y + atom->z*atom->z);
		inertia.y += atom->mass * (atom->x*atom->x + atom->z*atom->z);
		inertia.z += atom->mass * (atom->x*atom->x + atom->y*atom->y);
	}

	inertia_out[0] = inertia.x;
	inertia_out[1] = inertia.y;
	inertia_out[2] = inertia.z;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_atom_count(struct efp *efp, size_t frag_idx, size_t *n_atoms)
{
	assert(efp);
	assert(n_atoms);
	assert(frag_idx < efp->n_frag);

	*n_atoms = efp->frags[frag_idx].n_atoms;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_frag_atoms(struct efp *efp, size_t frag_idx, size_t size,
    struct efp_atom *atoms)
{
	struct frag *frag;

	assert(efp);
	assert(atoms);
	assert(frag_idx < efp->n_frag);
	assert(size >= efp->frags[frag_idx].n_atoms);

	frag = efp->frags + frag_idx;
	memcpy(atoms, frag->atoms, frag->n_atoms * sizeof(struct efp_atom));

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT void
efp_torque_to_derivative(const double *euler, const double *torque,
    double *deriv)
{
	assert(euler);
	assert(torque);
	assert(deriv);

	double tx = torque[0];
	double ty = torque[1];
	double tz = torque[2];

	double sina = sin(euler[0]);
	double cosa = cos(euler[0]);
	double sinb = sin(euler[1]);
	double cosb = cos(euler[1]);

	deriv[0] = tz;
	deriv[1] = cosa * tx + sina * ty;
	deriv[2] = sinb * sina * tx - sinb * cosa * ty + cosb * tz;
}

EFP_EXPORT const char *
efp_banner(void)
{
	static const char banner[] =
		"LIBEFP ver. " LIBEFP_VERSION_STRING "\n"
		"Copyright (c) 2012-2017 Ilya Kaliman\n"
		"\n"
		"Journal References:\n"
		"  - Kaliman and Slipchenko, JCC 2013.\n"
		"    DOI: http://dx.doi.org/10.1002/jcc.23375\n"
		"  - Kaliman and Slipchenko, JCC 2015.\n"
		"    DOI: http://dx.doi.org/10.1002/jcc.23772\n"
		"\n"
		"Project web site: https://github.com/libefp2/libefp/\n";

	return banner;
}

EFP_EXPORT void
efp_print_banner(void)
{
	puts(efp_banner());
}

EFP_EXPORT const char *
efp_result_to_string(enum efp_result res)
{
	switch (res) {
	case EFP_RESULT_SUCCESS:
		return "Operation was successful.";
	case EFP_RESULT_FATAL:
		return "Fatal error has occurred.";
	case EFP_RESULT_NO_MEMORY:
		return "Insufficient memory.";
	case EFP_RESULT_FILE_NOT_FOUND:
		return "File not found.";
	case EFP_RESULT_SYNTAX_ERROR:
		return "Syntax error.";
	case EFP_RESULT_UNKNOWN_FRAGMENT:
		return "Unknown EFP fragment.";
	case EFP_RESULT_POL_NOT_CONVERGED:
		return "Polarization SCF procedure did not converge.";
	}
	assert(0);
}

EFP_EXPORT enum efp_result
efp_get_pairwise_energy(struct efp *efp, struct efp_energy *pair_energies){

        assert(efp);
        assert(pair_energies);

        memcpy(pair_energies, efp->pair_energies, efp->n_frag * sizeof(struct efp_energy));
        return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_set_pairwise_energy(struct efp *efp, struct efp_energy *pair_energies)
{
    assert(efp);
    assert(pair_energies);

    memcpy(efp->pair_energies, pair_energies, efp->n_frag * sizeof(struct efp_energy));
    return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_set_symmlist(struct efp *efp)
{
    assert(efp);
    assert(efp->symmlist);
    assert(efp->skiplist);

    if  (efp->opts.symmetry == 0)
    {
        for (size_t i = 0; i < efp->n_frag; i++) {
            efp->symmlist[i] = 0;
        }
    }

    else if (efp->opts.symm_frag == EFP_SYMM_FRAG_FRAG) {
        //printf("\n n_lib %d \n", efp->n_lib);

        // this needs to be changed for list settings of symmetry!!!
        efp->nsymm_frag = efp->n_lib;
        char name[32];
        char** unique_names=malloc(efp->n_lib * sizeof(name));
        for (int i = 0; i < efp->n_lib; i++){
                unique_names[i] = NULL;
        }

        int n_unique = 0;
        for (size_t i = 0; i < efp->n_frag; i++) {
            struct frag *frag = efp->frags + i;
            for (size_t m = 0; m < efp->n_lib; m++) {
                if (unique_names[m] != NULL && strcmp(frag->name, unique_names[m]) == 0) {
                    efp->symmlist[i] = m+1;
                    break;
                }
            }
            if (efp->symmlist[i] == 0) // did not find this fragment name in sofar unique names
            {
                unique_names[n_unique] = frag->name;
                n_unique++;
                efp->symmlist[i] = n_unique;
            }
            // printf("symm_list %d \n", efp->symmlist[i]);
        }

        // setup skiplist now...
        for (size_t i = 0; i < efp->n_frag; i++) {
            for (size_t j = 0; j < efp->n_frag; j++) {
                efp_skip_fragments(efp, i, j, 1);
            }
        }

        n_unique = 1;
        //memset(efp->skiplist,true,efp->n_frag*efp->n_frag);
        for (size_t i = 0; i < efp->n_frag; i++) {
            // this is the first occuracnce of the symmetry-unique fragment
            if (efp->symmlist[i] == n_unique) {
                for (size_t j = 0; j < efp->n_frag; j++) {
                    efp_skip_fragments(efp, i, j, 0);
                }
                n_unique ++;
            }
        }

        free(unique_names);
    }
    else {
        printf("\n DO NOT KNOW WHAT TO DO WITH THIS SYMMETRIC SYSTEM:  SYMM_FRAG IS UNKNOWN \n");
    }

    /*
    //printf("\n skiplist \n");
    for (int i = 0; i < efp->n_frag; i++){
        for (int j=0; j < efp->n_frag; j++){
            printf(" %s ", efp->skiplist[i*efp->n_frag + j] ? "true" : "false");
        }
    }
     */

    return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_symmlist(struct efp *efp, size_t frag_idx, size_t *symm){

    assert(efp);
    assert(efp->symmlist);

    *symm = efp->symmlist[frag_idx];
    return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_nsymm_frag(struct efp *efp, size_t *nsymm_frag){

    assert(efp);
    assert(efp->nsymm_frag);

    *nsymm_frag = efp->nsymm_frag;
    return EFP_RESULT_SUCCESS;
}

void
unique_symm_frag(struct efp *efp, size_t *unique_frag){
    //printf("\n Symmetry-unique fragments \n");
    int n = 0;
    int i = 0;
    do {
        if (efp->symmlist[i] > n) {
            unique_frag[n] = i;
            //printf(" %d ", unique_frag[n]);
            n++;
        }
        i++;
    } while (n < efp->nsymm_frag);
}

void
n_symm_frag(struct efp *efp, size_t *symm_frag) {

    for (size_t i = 0; i < efp->nsymm_frag; i++) {
        size_t counter = 0;
        for (size_t j = 0; j < efp->n_frag; j++) {
            if (efp->symmlist[i] == efp->symmlist[j])
                counter++;
        }
        symm_frag[i] = counter;
        // printf("\n symm_frag %d = %d", i, symm_frag[i]);
    }
}
