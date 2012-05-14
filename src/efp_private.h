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

#include "efp.h"
#include "math_util.h"
#include "terms.h"

#define EFP_EXPORT __attribute__((visibility("default")))
#define ARRAY_SIZE(arr) ((int)(sizeof(arr)/sizeof(arr[0])))
#define EFP_INIT_MAGIC 0xEF2012AD

struct frag {
	/* fragment name */
	char *name;

	/* fragment center of mass */
	double x, y, z;

	/* pointer to the initial fragment state in library */
	const struct frag *lib;

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

	/* name of the basis set used to generate fragment potential */
	char basis_name[32];

	/* shells of fragment basis (S,P,D,F,L); zero terminated */
	char *shells;

	/* upper triangle of fock matrix, size = n_lmo * (n_lmo + 1) / 2 */
	double *xr_fock_mat;

	/* exchange repulsion wavefunction size */
	int xr_wf_size;

	/* exchange repulsion wavefunction, size = n_lmo * xr_wf_size */
	double *xr_wf;
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

	/* number of blocks in xr integrals computation */
	int n_xr_blocks;

	/* offsets of fragments in their xr blocks */
	int *xr_block_frag_offset;

	/* contributions to dispersion damping from overlap integrals */
	double *disp_damp_overlap;

	/* fragment offsets in disp_damp_overlap array */
	int *disp_damp_overlap_offset;

	/* callbacks */
	struct efp_callbacks callbacks;

	/* user parameters for this EFP computation */
	struct efp_opts opts;

	/* information about ab initio region */
	struct efp_qm_data qm_data;

	/* EFP energy terms */
	struct efp_energy energy;

	/* total energy gradient */
	double *grad;

	/* gradient on QM atoms */
	double *qm_grad;

	/* initialization check */
	unsigned magic;
};

enum efp_result efp_read_potential(struct efp *efp, const char **files);
void efp_pol_scf_init(struct efp *efp);
double efp_compute_pol_energy(struct efp *efp);
void efp_update_elec(struct frag *frag, const mat_t *rotmat);
void efp_update_pol(struct frag *frag, const mat_t *rotmat);
void efp_update_disp(struct frag *frag, const mat_t *rotmat);
void efp_update_xr(struct frag *frag, const mat_t *rotmat);

#endif /* LIBEFP_EFP_PRIVATE_H */
