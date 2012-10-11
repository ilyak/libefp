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

#ifndef LIBEFP_EFP_PRIVATE_H
#define LIBEFP_EFP_PRIVATE_H

#include "../common/compat.h"
#include "../common/math_util.h"
#include "../common/util.h"

#include "efp.h"
#include "int.h"
#include "swf.h"
#include "terms.h"

#define EFP_EXPORT __attribute__((visibility("default")))
#define UNUSED __attribute__((unused))

#define ARRAY_SIZE(arr) (sizeof(arr)/sizeof(arr[0]))
#define EFP_INIT_MAGIC 0xEF2012AD

struct frag {
	/* fragment name */
	char *name;

	/* fragment center of mass */
	double x, y, z;

	/* rotation matrix representing orientation of a fragment */
	mat_t rotmat;

	/* pointer to the initial fragment state in library */
	const struct frag *lib;

	/* force on fragment center of mass */
	vec_t force;

	/* torque on fragment */
	vec_t torque;

	/* number of atoms in this fragment */
	int n_atoms;

	/* fragment atoms */
	struct efp_atom *atoms;

	/* distributed multipoles */
	struct multipole_pt {
		double x, y, z;
		double monopole;
		vec_t dipole;
		double quadrupole[6];
		double octupole[10];
	} *multipole_pts;

	/* number of distributed multipole points */
	int n_multipole_pts;

	/* electrostatic screening parameters */
	double *screen_params;

	/* ab initio electrostatic screening parameters */
	double *ai_screen_params;

	/* distributed polarizability points */
	struct polarizable_pt {
		double x, y, z;
		mat_t tensor;
		vec_t elec_field;
		vec_t induced_dipole;
		vec_t induced_dipole_new;
		vec_t induced_dipole_conj;
		vec_t induced_dipole_conj_new;
	} *polarizable_pts;

	/* number of distributed polarizability points */
	int n_polarizable_pts;

	/* dynamic polarizability points */
	struct dynamic_polarizable_pt {
		double x, y, z;
		double trace[12];
	} *dynamic_polarizable_pts;

	/* number of dynamic polarizability points */
	int n_dynamic_polarizable_pts;

	/* number of localized molecular orbitals */
	int n_lmo;

	/* localized molecular orbital centroids */
	vec_t *lmo_centroids;

	/* spin multiplicity */
	int multiplicity;

	/* number of exchange repulsion basis shells */
	int n_xr_shells;

	/* exchange repulsion basis shells */
	struct shell *xr_shells;

	/* upper triangle of fock matrix, size = n_lmo * (n_lmo + 1) / 2 */
	double *xr_fock_mat;

	/* exchange repulsion wavefunction size */
	int xr_wf_size;

	/* exchange repulsion wavefunction, size = n_lmo * xr_wf_size */
	double *xr_wf;

	/* rotational derivatives of MO coefficients */
	double *xr_wf_deriv[3];

	/* overlap integrals; used for overlap-based dispersion damping */
	double *overlap_int;

	/* derivatives of overlap integrals */
	six_t *overlap_int_deriv;
};

struct efp {
	/* number of fragments */
	int n_frag;

	/* array of fragments */
	struct frag *frags;

	/* number of fragments in the library */
	int n_lib;

	/* array with the library of fragment initial parameters */
	struct frag *lib;

	/* callbacks */
	struct efp_callbacks callbacks;

	/* user parameters for this EFP computation */
	struct efp_opts opts;

	/* gradient will also be computed if nonzero */
	int do_gradient;

	/* periodic simulation box size */
	vec_t box;

	/* number of ab initio atoms */
	int n_qm_atoms;

	struct qm_atom {
		double x, y, z;
		double znuc;
		vec_t grad;
	} *qm_atoms;

	/* EFP energy terms */
	struct efp_energy energy;

	/* initialization check */
	unsigned magic;
};

static inline struct swf
make_swf(struct efp *efp, const struct frag *fr_i, const struct frag *fr_j)
{
	struct swf swf = {
		.swf = 1.0,
		.dswf = { 0.0, 0.0, 0.0 },
		.cell = { 0.0, 0.0, 0.0 }
	};

	if (!efp->opts.enable_pbc)
		return swf;

	vec_t dr = vec_sub(CVEC(fr_j->x), CVEC(fr_i->x));

	swf.cell.x = efp->box.x * round(dr.x / efp->box.x);
	swf.cell.y = efp->box.y * round(dr.y / efp->box.y);
	swf.cell.z = efp->box.z * round(dr.z / efp->box.z);

	dr.x -= swf.cell.x;
	dr.y -= swf.cell.y;
	dr.z -= swf.cell.z;

	double r = vec_len(&dr);

	swf.swf = get_swf(r, efp->opts.swf_cutoff);
	double dswf = get_dswf(r, efp->opts.swf_cutoff);

	swf.dswf.x = -dswf * dr.x;
	swf.dswf.y = -dswf * dr.y;
	swf.dswf.z = -dswf * dr.z;

	return swf;
}

static inline void
add_force(struct frag *frag, const vec_t *pt, const vec_t *force, const vec_t *add)
{
	vec_t dr = vec_sub(CVEC(pt->x), CVEC(frag->x));
	vec_t torque = vec_cross(&dr, force);

	if (add) {
		torque.x += add->x;
		torque.y += add->y;
		torque.z += add->z;
	}

	vec_atomic_add(&frag->force, force);
	vec_atomic_add(&frag->torque, &torque);
}

static inline void
sub_force(struct frag *frag, const vec_t *pt, const vec_t *force, const vec_t *add)
{
	vec_t dr = vec_sub(CVEC(pt->x), CVEC(frag->x));
	vec_t torque = vec_cross(&dr, force);

	if (add) {
		torque.x += add->x;
		torque.y += add->y;
		torque.z += add->z;
	}

	vec_atomic_sub(&frag->force, force);
	vec_atomic_sub(&frag->torque, &torque);
}

#endif /* LIBEFP_EFP_PRIVATE_H */
