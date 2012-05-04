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

#include "efp_private.h"

typedef enum efp_result (*term_energy_fn)(struct efp *efp);

/* order must be the same as in efp_term enum */
static const term_energy_fn term_energy_funcs[] = {
	efp_compute_elec, efp_compute_pol, efp_compute_disp,
	efp_compute_xr, efp_compute_chtr,
	efp_compute_ai_elec, efp_compute_ai_pol, efp_compute_ai_disp,
	efp_compute_ai_xr, efp_compute_ai_chtr
};

static inline int
initialized(struct efp *efp)
{
	return efp && efp->magic == EFP_INIT_MAGIC;
}

EFP_EXPORT enum efp_result
efp_get_energy(struct efp *efp, double *energy)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!energy)
		return EFP_RESULT_INVALID_ARGUMENT;

	memcpy(energy, efp->energy, EFP_TERM_COUNT * sizeof(double));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_get_gradient(struct efp *efp, int n_grad, double *grad)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!grad)
		return EFP_RESULT_INVALID_ARGUMENT;

	if (!efp->grad)
		return EFP_RESULT_GRADIENT_NOT_REQUESTED;

	if (n_grad < 6 * efp->n_frag)
		return EFP_RESULT_INVALID_ARRAY_SIZE;

	memcpy(grad, efp->grad, 6 * efp->n_frag * sizeof(double));
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_update_qm_data(struct efp *efp, const struct efp_qm_data *qm_data)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!qm_data || !qm_data->atoms)
		return EFP_RESULT_INVALID_ARGUMENT;

	memcpy(&efp->qm_data, qm_data, sizeof(struct efp_qm_data));

	size_t size = qm_data->n_atoms * sizeof(struct efp_qm_atom);
	efp->qm_data.atoms = malloc(size);
	if (!efp->qm_data.atoms)
		return EFP_RESULT_NO_MEMORY;
	memcpy(efp->qm_data.atoms, qm_data->atoms, size);

	return EFP_RESULT_SUCCESS;
}

static void
update_fragment(struct efp *efp, int idx, double x, double y, double z,
			const struct mat *rotmat)
{
	struct frag *frag = efp->frags + idx;
	frag->x = x, frag->y = y, frag->z = z;

	/* update atoms */
	for (int i = 0; i < frag->n_atoms; i++)
		move_pt(VEC(frag->x), rotmat, VEC(frag->lib->atoms[i].x),
			VEC(frag->atoms[i].x));

	efp_update_elec(frag, rotmat);
	efp_update_pol(frag, rotmat);
	efp_update_disp(frag, rotmat);
	efp_update_xr(frag, rotmat);
}

static void
euler_to_matrix(double a, double b, double c, struct mat *out)
{
	double sina = sin(a), cosa = cos(a);
	double sinb = sin(b), cosb = cos(b);
	double sinc = sin(c), cosc = cos(c);
	out->xx =  cosa * cosc - sina * cosb * sinc;
	out->xy = -cosa * sinc - sina * cosb * cosc;
	out->xz =  sinb * sina;
	out->yx =  sina * cosc + cosa * cosb * sinc;
	out->yy = -sina * sinc + cosa * cosb * cosc;
	out->yz = -sinb * cosa;
	out->zx =  sinb * sinc;
	out->zy =  sinb * cosc;
	out->zz =  cosb;
}

EFP_EXPORT enum efp_result
efp_update_fragments(struct efp *efp, const double *xyzabc)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!xyzabc)
		return EFP_RESULT_INVALID_ARGUMENT;

	for (int i = 0; i < efp->n_frag; i++, xyzabc += 6) {
		double x = xyzabc[0], y = xyzabc[1], z = xyzabc[2];
		double a = xyzabc[3], b = xyzabc[4], c = xyzabc[5];

		struct mat rotmat;
		euler_to_matrix(a, b, c, &rotmat);

		update_fragment(efp, i, x, y, z, &rotmat);
	}
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_scf_update(struct efp *efp, double *energy)
{
	enum efp_result res;

	*energy = 0.0;
	if ((res = efp_scf_update_pol(efp, energy)))
		return res;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_compute(struct efp *efp)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (efp->grad)
		memset(efp->grad, 0, 6 * efp->n_frag * sizeof(double));

	enum efp_result res;

	for (int i = 0; i < ARRAY_SIZE(term_energy_funcs); i++)
		if (efp->opts.terms & (1 << i))
			if ((res = term_energy_funcs[i](efp)))
				return res;

	/* overlap-based charge penetration */
	if ((efp->opts.terms & EFP_TERM_ELEC) &&
	     efp->opts.elec_damp == EFP_ELEC_DAMP_OVERLAP) {
		int idx = efp_get_term_index(EFP_TERM_ELEC);
		efp->energy[idx] += efp->charge_pen_energy;
	}

	return EFP_RESULT_SUCCESS;
}

static void
free_frag(struct frag *frag)
{
	if (!frag)
		return;

	free(frag->name);
	free(frag->atoms);
	free(frag->multipole_pts);
	free(frag->polarizable_pts);
	free(frag->dynamic_polarizable_pts);
	free(frag->lmo_centroids);
	free(frag->xr_fock_mat);
	free(frag->xr_wf);
	free(frag->shells);
	free(frag->screen_params);
	free(frag->ai_screen_params);

	/* don't do free(frag) here */
}

EFP_EXPORT void
efp_shutdown(struct efp *efp)
{
	if (!efp)
		return;

	for (int i = 0; i < efp->n_frag; i++)
		free_frag(efp->frags + i);

	for (int i = 0; i < efp->n_lib; i++)
		free_frag(efp->lib + i);

	free(efp->frags);
	free(efp->lib);
	free(efp->xr_block_frag_offset);
	free(efp->disp_damp_overlap_offset);
	free(efp->disp_damp_overlap);
	free(efp->grad);
	free(efp->qm_grad);
	free(efp);
}

static enum efp_result
copy_frag(struct frag *dest, const struct frag *src)
{
	size_t size;
	memcpy(dest, src, sizeof(struct frag));

	if (src->name) {
		dest->name = strdup(src->name);
		if (!dest->name)
			return EFP_RESULT_NO_MEMORY;
	}
	if (src->shells) {
		dest->shells = strdup(src->shells);
		if (!dest->shells)
			return EFP_RESULT_NO_MEMORY;
	}
	if (src->atoms) {
		size = src->n_atoms * sizeof(struct efp_atom);
		dest->atoms = malloc(size);
		if (!dest->atoms)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->atoms, src->atoms, size);
	}
	if (src->multipole_pts) {
		size = src->n_multipole_pts * sizeof(struct multipole_pt);
		dest->multipole_pts = malloc(size);
		if (!dest->multipole_pts)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->multipole_pts, src->multipole_pts, size);
	}
	if (src->screen_params) {
		size = src->n_multipole_pts * sizeof(double);
		dest->screen_params = malloc(size);
		if (!dest->screen_params)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->screen_params, src->screen_params, size);
	}
	if (src->ai_screen_params) {
		size = src->n_multipole_pts * sizeof(double);
		dest->ai_screen_params = malloc(size);
		if (!dest->ai_screen_params)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->ai_screen_params, src->ai_screen_params, size);
	}
	if (src->polarizable_pts) {
		size = src->n_polarizable_pts * sizeof(struct polarizable_pt);
		dest->polarizable_pts = malloc(size);
		if (!dest->polarizable_pts)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->polarizable_pts, src->polarizable_pts, size);
	}
	if (src->dynamic_polarizable_pts) {
		size = src->n_dynamic_polarizable_pts *
				sizeof(struct dynamic_polarizable_pt);
		dest->dynamic_polarizable_pts = malloc(size);
		if (!dest->dynamic_polarizable_pts)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->dynamic_polarizable_pts,
				src->dynamic_polarizable_pts, size);
	}
	if (src->lmo_centroids) {
		size = src->n_lmo * sizeof(struct vec);
		dest->lmo_centroids = malloc(size);
		if (!dest->lmo_centroids)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->lmo_centroids, src->lmo_centroids, size);
	}
	if (src->xr_fock_mat) {
		size = src->n_lmo * (src->n_lmo + 1) / 2 * sizeof(double);
		dest->xr_fock_mat = malloc(size);
		if (!dest->xr_fock_mat)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xr_fock_mat, src->xr_fock_mat, size);
	}
	if (src->xr_wf) {
		size = src->n_lmo * src->xr_wf_size * sizeof(double);
		dest->xr_wf = malloc(size);
		if (!dest->xr_wf)
			return EFP_RESULT_NO_MEMORY;
		memcpy(dest->xr_wf, src->xr_wf, size);
	}
	return EFP_RESULT_SUCCESS;
}

static struct frag *
find_frag_in_library(struct efp *efp, const char *name)
{
	for (int i = 0; i < efp->n_lib; i++)
		if (!strcasecmp(efp->lib[i].name, name))
			return efp->lib + i;

	return NULL;
}

EFP_EXPORT void
efp_opts_default(struct efp_opts *opts)
{
	if (!opts)
		return;

	memset(opts, 0, sizeof(struct efp_opts));

	opts->terms = EFP_TERM_ELEC | EFP_TERM_POL | EFP_TERM_DISP |
		EFP_TERM_XR | EFP_TERM_AI_ELEC | EFP_TERM_AI_POL;
}

static enum efp_result
check_opts(const struct efp_opts *opts)
{
	enum efp_term terms = opts->terms;

	if (((terms & EFP_TERM_AI_ELEC) && !(terms & EFP_TERM_ELEC)) ||
	    ((terms & EFP_TERM_AI_POL) && !(terms & EFP_TERM_POL)) ||
	    ((terms & EFP_TERM_POL) && !(terms & EFP_TERM_ELEC)) ||
	    ((terms & EFP_TERM_AI_DISP) && !(terms & EFP_TERM_DISP)) ||
	    ((terms & EFP_TERM_AI_XR) && !(terms & EFP_TERM_XR)) ||
	    ((terms & EFP_TERM_AI_CHTR) && !(terms & EFP_TERM_CHTR)))
		return EFP_RESULT_INCONSISTENT_TERMS;

	if (terms & EFP_TERM_ELEC) {
		if (opts->elec_damp == EFP_ELEC_DAMP_OVERLAP &&
				!(terms & EFP_TERM_XR))
			return EFP_RESULT_OVERLAP_INTEGRALS_REQUIRED;
	}
	if (terms & EFP_TERM_DISP) {
		if (opts->disp_damp == EFP_DISP_DAMP_OVERLAP &&
				!(terms & EFP_TERM_XR))
			return EFP_RESULT_OVERLAP_INTEGRALS_REQUIRED;
	}
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
check_params(struct efp *efp)
{
	if (efp->opts.terms & EFP_TERM_ELEC) {
		for (int i = 0; i < efp->n_frag; i++)
			if (!efp->frags[i].multipole_pts)
				return EFP_RESULT_PARAMETERS_MISSING;
	}
	if (efp->opts.terms & EFP_TERM_POL) {
		for (int i = 0; i < efp->n_frag; i++)
			if (!efp->frags[i].polarizable_pts)
				return EFP_RESULT_PARAMETERS_MISSING;
	}
	if (efp->opts.terms & EFP_TERM_DISP) {
		for (int i = 0; i < efp->n_frag; i++)
			if (!efp->frags[i].dynamic_polarizable_pts)
				return EFP_RESULT_PARAMETERS_MISSING;
	}
	if (efp->opts.terms & EFP_TERM_XR) {
		if (!efp->callbacks.get_overlap_integrals ||
		    !efp->callbacks.get_kinetic_integrals)
			return EFP_RESULT_CALLBACK_NOT_SET;

		for (int i = 0; i < efp->n_frag; i++)
			if (!efp->frags[i].xr_fock_mat ||
			    !efp->frags[i].xr_wf ||
			    !efp->frags[i].lmo_centroids)
				return EFP_RESULT_PARAMETERS_MISSING;
	}
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
setup_xr(struct efp *efp)
{
	static const int max_basis_per_block = 2000;

	efp->n_xr_blocks = 0;
	for (int i = 0; i < efp->n_frag; efp->n_xr_blocks++)
		for (int n_basis = 0; i < efp->n_frag && n_basis +
			efp->frags[i].xr_wf_size <= max_basis_per_block; i++)
				n_basis += efp->frags[i].xr_wf_size;

	size_t size = (efp->n_xr_blocks + 1) * sizeof(int);
	efp->xr_block_frag_offset = malloc(size);
	if (!efp->xr_block_frag_offset)
		return EFP_RESULT_NO_MEMORY;

	efp->xr_block_frag_offset[0] = 0;
	for (int i = 0, block = 1; i < efp->n_frag; block++) {
		for (int n_basis = 0; i < efp->n_frag && n_basis +
			efp->frags[i].xr_wf_size <= max_basis_per_block; i++)
				n_basis += efp->frags[i].xr_wf_size;
		efp->xr_block_frag_offset[block] =
			efp->xr_block_frag_offset[block - 1] + i;
	}
	return EFP_RESULT_SUCCESS;
}

static enum efp_result
setup_disp(struct efp *efp)
{
	int n_disp = 0;
	for (int i = 0; i < efp->n_frag; i++)
		n_disp += efp->frags[i].n_dynamic_polarizable_pts;

	if (n_disp == 0)
		return EFP_RESULT_SUCCESS;

	efp->disp_damp_overlap_offset = malloc((efp->n_frag + 1) * sizeof(int));
	if (!efp->disp_damp_overlap_offset)
		return EFP_RESULT_NO_MEMORY;

	efp->disp_damp_overlap_offset[0] = 0;
	for (int i = 1; i <= efp->n_frag; i++)
		efp->disp_damp_overlap_offset[i] =
			efp->disp_damp_overlap_offset[i - 1] +
			efp->frags[i - 1].n_dynamic_polarizable_pts;

	efp->disp_damp_overlap = malloc(n_disp * n_disp * sizeof(double));
	if (!efp->disp_damp_overlap)
		return EFP_RESULT_NO_MEMORY;

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT enum efp_result
efp_init(struct efp **out, const struct efp_opts *opts,
	 const struct efp_callbacks *callbacks,
	 const char **potential_files, const char **fragname)
{
	if (!out || !opts || !callbacks || !potential_files || !fragname)
		return EFP_RESULT_INVALID_ARGUMENT;

	*out = calloc(1, sizeof(struct efp));
	if (!*out)
		return EFP_RESULT_NO_MEMORY;

	enum efp_result res;
	struct efp *efp = *out;

	if ((res = check_opts(opts)))
		return res;

	memcpy(&efp->opts, opts, sizeof(struct efp_opts));
	memcpy(&efp->callbacks, callbacks, sizeof(struct efp_callbacks));

	if ((res = efp_read_potential(efp, potential_files)))
		return res;

	while (fragname[efp->n_frag])
		efp->n_frag++;

	efp->frags = calloc(efp->n_frag, sizeof(struct frag));
	if (!efp->frags)
		return EFP_RESULT_NO_MEMORY;

	for (int i = 0; i < efp->n_frag; i++) {
		struct frag *frag = find_frag_in_library(efp, fragname[i]);
		if (!frag)
			return EFP_RESULT_UNKNOWN_FRAGMENT;

		if ((res = copy_frag(efp->frags + i, frag)))
			return res;
	}

	if (opts->do_gradient) {
		efp->grad = malloc(6 * efp->n_frag * sizeof(double));
		if (!efp->grad)
			return EFP_RESULT_NO_MEMORY;
	}

	if ((res = setup_xr(efp)))
		return res;

	if ((res = setup_disp(efp)))
		return res;

	if ((res = check_params(efp)))
		return res;

	efp->magic = EFP_INIT_MAGIC;
	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT int
efp_get_frag_count(struct efp *efp)
{
	return initialized(efp) ? efp->n_frag : -1;
}

EFP_EXPORT int
efp_get_frag_atom_count(struct efp *efp, int frag_idx)
{
	return initialized(efp) && frag_idx >= 0 && frag_idx < efp->n_frag ?
			efp->frags[frag_idx].n_atoms : -1;
}

EFP_EXPORT enum efp_result
efp_get_frag_atoms(struct efp *efp, int frag_idx, struct efp_atom *atoms)
{
	if (!initialized(efp))
		return EFP_RESULT_NOT_INITIALIZED;

	if (!atoms || frag_idx < 0 || frag_idx >= efp->n_frag)
		return EFP_RESULT_INVALID_ARGUMENT;

	struct frag *frag = efp->frags + frag_idx;
	memcpy(atoms, frag->atoms, frag->n_atoms * sizeof(struct efp_atom));

	return EFP_RESULT_SUCCESS;
}

EFP_EXPORT int
efp_get_term_index(enum efp_term term)
{
	switch (term) {
		case EFP_TERM_ELEC: return 0;
		case EFP_TERM_POL: return 1;
		case EFP_TERM_DISP: return 2;
		case EFP_TERM_XR: return 3;
		case EFP_TERM_CHTR: return 4;
		case EFP_TERM_AI_ELEC: return 5;
		case EFP_TERM_AI_POL: return 6;
		case EFP_TERM_AI_DISP: return 7;
		case EFP_TERM_AI_XR: return 8;
		case EFP_TERM_AI_CHTR: return 9;
	}
	return -1;
}

EFP_EXPORT const char *
efp_banner(void)
{
	static const char banner[] =
		"libefp\n"
		"The Effective Fragment Potential method implementation\n"
		"Copyright (c) 2012 Ilya Kaliman\n"
		"See LICENSE file for licensing terms\n"
		"Source code available at http://repo.or.cz/w/efp.git/\n";

	return banner;
}

EFP_EXPORT const char *
efp_result_to_string(enum efp_result res)
{
	switch (res) {
	case EFP_RESULT_SUCCESS:
return "no error";
	case EFP_RESULT_NO_MEMORY:
return "out of memory";
	case EFP_RESULT_NOT_IMPLEMENTED:
return "operation is not implemented";
	case EFP_RESULT_INVALID_ARGUMENT:
return "invalid argument to function";
	case EFP_RESULT_NOT_INITIALIZED:
return "efp was not properly initialized";
	case EFP_RESULT_FILE_NOT_FOUND:
return "file not found";
	case EFP_RESULT_SYNTAX_ERROR:
return "syntax error in potential data";
	case EFP_RESULT_BASIS_NOT_SPECIFIED:
return "no EFP basis specified";
	case EFP_RESULT_UNKNOWN_FRAGMENT:
return "unknown fragment type";
	case EFP_RESULT_DUPLICATE_PARAMETERS:
return "fragment parameters contain fragments with the same name";
	case EFP_RESULT_CALLBACK_NOT_SET:
return "required callback function is not set";
	case EFP_RESULT_CALLBACK_FAILED:
return "callback function failed";
	case EFP_RESULT_OVERLAP_INTEGRALS_REQUIRED:
return "overlap based damping requires exchange repulsion";
	case EFP_RESULT_GRADIENT_NOT_REQUESTED:
return "gradient computation was not requested";
	case EFP_RESULT_PARAMETERS_MISSING:
return "required EFP fragment parameters are missing";
	case EFP_RESULT_INVALID_ARRAY_SIZE:
return "invalid array size";
	case EFP_RESULT_UNSUPPORTED_SCREEN:
return "unsupported SCREEN group found in EFP data";
	case EFP_RESULT_INCONSISTENT_TERMS:
return "inconsistent energy terms selected";
	}
return "unknown result";
}
