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

#ifndef LIBEFP_EFP_H
#define LIBEFP_EFP_H

/* Public libefp interface. */

#ifdef __cplusplus
extern "C" {
#endif

/* Result of an operation. */
enum efp_result {
	EFP_RESULT_SUCCESS = 0,

	EFP_RESULT_NO_MEMORY,
	EFP_RESULT_NOT_IMPLEMENTED,
	EFP_RESULT_INVALID_ARGUMENT,
	EFP_RESULT_NOT_INITIALIZED,
	EFP_RESULT_FILE_NOT_FOUND,
	EFP_RESULT_SYNTAX_ERROR,
	EFP_RESULT_BASIS_NOT_SPECIFIED,
	EFP_RESULT_UNKNOWN_FRAGMENT,
	EFP_RESULT_DUPLICATE_PARAMETERS,
	EFP_RESULT_CALLBACK_NOT_SET,
	EFP_RESULT_CALLBACK_FAILED,
	EFP_RESULT_OVERLAP_INTEGRALS_REQUIRED,
	EFP_RESULT_GRADIENT_NOT_REQUESTED,
	EFP_RESULT_PARAMETERS_MISSING,
	EFP_RESULT_INVALID_ARRAY_SIZE,
	EFP_RESULT_UNSUPPORTED_SCREEN,
	EFP_RESULT_INCONSISTENT_TERMS
};

/* Total number of EFP terms. */
#define EFP_TERM_COUNT 10

/* Flags to specify EFP energy terms. */
enum efp_term {
	EFP_TERM_ELEC = 1 << 0,     /* EFP/EFP electrostatics */
	EFP_TERM_POL = 1 << 1,      /* EFP/EFP polarization */
	EFP_TERM_DISP = 1 << 2,     /* EFP/EFP dispersion */
	EFP_TERM_XR = 1 << 3,       /* EFP/EFP exchange repulsion */
	EFP_TERM_CHTR = 1 << 4,     /* EFP/EFP charge transfer */
	EFP_TERM_AI_ELEC = 1 << 5,  /* Ab-initio/EFP electrostatics */
	EFP_TERM_AI_POL = 1 << 6,   /* Ab-initio/EFP polarization */
	EFP_TERM_AI_DISP = 1 << 7,  /* Ab-initio/EFP dispersion */
	EFP_TERM_AI_XR = 1 << 8,    /* Ab-initio/EFP exchange repulsion */
	EFP_TERM_AI_CHTR = 1 << 9   /* Ab-initio/EFP charge transfer */
};

/* Fragment-fragment dispersion damping. */
enum efp_disp_damp {
	EFP_DISP_DAMP_OVERLAP = 0, /* overlap-based damping (default) */
	EFP_DISP_DAMP_TT           /* Tang-Toennies damping */
};

/* Fragment-fragment electrostatic damping. */
enum efp_elec_damp {
	EFP_ELEC_DAMP_SCREEN = 0,  /* controlled by SCREEN input (default) */
	EFP_ELEC_DAMP_OVERLAP      /* overlap-based damping */
};

/* Opaque EFP structure. */
struct efp;

/* Options controlling EFP computation. */
struct efp_opts {
	/* specifies which energy terms to compute */
	enum efp_term terms;

	/* set to nonzero to request gradient computation */
	int do_gradient;

	/* dispersion damping type */
	enum efp_disp_damp disp_damp;

	/* electrostatics damping type */
	enum efp_elec_damp elec_damp;
};

/* EFP atom info. */
struct efp_atom {
	char label[32];   /* atom label */
	char basis[32];   /* name of a basis on this atom */
	double x, y, z;   /* atom position */
	double mass;      /* atom mass */
	double znuc;      /* nuclear charge */
};

/* QM atom info. */
struct efp_qm_atom {
	double x, y, z;   /* atom position */
	double znuc;      /* effective nuclear charge */
};

/* Information about ab initio region. */
struct efp_qm_data {
	int n_atoms;                /* number of atoms in QM part */
	struct efp_qm_atom *atoms;  /* atom data */
};

/* XXX */
struct efp_xr_block {
	int n_atoms_i;             /* */
	struct efp_atom *atoms_i;  /* */
	int basis_size_i;          /* */
	int n_atoms_j;             /* */
	struct efp_atom *atoms_j;  /* */
	int basis_size_j;          /* */
};

/* Callback functions */

/* Compute overlap integrals over Gaussian basis functions.
 *
 * block (input): block information, see efp_xr_block.
 * s (output): array of size * size elements where overlap integrals should be
 *             written.
 * sx (output): if not NULL, an array of 3 * size * size elements where x,y,z
 *              derivatives of overlap integrals should be written.
 * user_data (input): user data which was specified during initialization.
 *
 * return: EFP_RESULT_CALLBACK_FAILED on error or EFP_RESULT_SUCCESS otherwise.
 */
typedef enum efp_result (*efp_overlap_integrals_fn)(struct efp_xr_block *block,
		double *s, double *sx, void *user_data);

/* Compute kinetic energy integrals over Gaussian basis functions.
 *
 * block (input): block information, see efp_xr_block.
 * t (output): array of size * size elements where kinetic energy integrals
 *             should be written.
 * tx (output): if not NULL, an array of 3 * size * size elements where x,y,z
 *              derivatives of kinetic energy integrals should be written.
 * user_data (input): user data which was specified during initialization.
 *
 * return: EFP_RESULT_CALLBACK_FAILED on error or EFP_RESULT_SUCCESS otherwise.
 */
typedef enum efp_result (*efp_kinetic_integrals_fn)(struct efp_xr_block *block,
		double *t, double *tx, void *user_data);

/* Compute field form electrons on specified points.
 *
 * n_pt (input): number of points in xyz array.
 * xyz (input): x,y,z coordinates of points where electric field is needed.
 * field (output): computed x,y,z components of electric field.
 * user_data (input): user data which was specified during initialization.
 *
 * return: EFP_RESULT_CALLBACK_FAILED on error or EFP_RESULT_SUCCESS otherwise.
 */
typedef enum efp_result (*efp_electron_density_field_fn)(int n_pt,
	const double *xyz, double *field, void *user_data);

/* Information about callback functions. */
struct efp_callbacks {
	/* callback function, see efp_overlap_integrals_fn */
	efp_overlap_integrals_fn get_overlap_integrals;

	/* will be passed as a last parameter to get_overlap_integrals */
	void *get_overlap_integrals_user_data;

	/* callback function, see efp_kinetic_integrals_fn */
	efp_kinetic_integrals_fn get_kinetic_integrals;

	/* will be passed as a last parameter to get_kinetic_integrals */
	void *get_kinetic_integrals_user_data;

	/* callback function, see efp_electron_density_field_fn */
	efp_electron_density_field_fn get_electron_density_field;

	/* will be passed as a last parameter to get_electron_density_field */
	void *get_electron_density_field_user_data;
};

/* Public methods */

/* Get a printable string with the information about this library. */
const char *efp_banner(void);

/* Set options to default values. */
void efp_opts_default(struct efp_opts *opts);

/* Initialize the EFP computation.
 *
 * out (output): initialized opaque efp structure.
 * opts (input): user-defined options controlling the computation.
 * callbacks (input): user-supplied callback functions (see efp_callbacks).
 * potential_files (input): NULL-terminated list with paths to the files with
 *                          Effective Fragment Potential data.
 * fragname (input): NULL-terminated list with names of the fragments
 *                   comprising the molecular system.
 *
 * return: EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_init(struct efp **out,
			 const struct efp_opts *opts,
			 const struct efp_callbacks *callbacks,
			 const char **potential_files,
			 const char **fragname);

/* Update information about ab initio region. */
enum efp_result efp_update_qm_data(struct efp *efp,
				   const struct efp_qm_data *qm_data);

/* Update positions and orientations of effective fragments.
 *
 * xyzabc (input): for each fragment specifies x,y,z components of fragment
 *                 center of mass position and three Euler rotation angles
 *                 representing orientation of a fragment.
 */
enum efp_result efp_update_fragments(struct efp *efp, const double *xyzabc);

/* Update wave function dependent terms. Must be called during QM SCF.
 *
 * energy (output): wave function dependent energy.
 */
enum efp_result efp_scf_update(struct efp *efp, double *energy);

/* Perform the EFP computation. This call can take long time to complete. */
enum efp_result efp_compute(struct efp *efp);

/* Get computed energies of corresponding EFP terms.
 *
 * energy (output): array of size EFP_TERM_COUNT where energy terms will be
 *                  stored. The sum of the elements of this array will give
 *                  total EFP/EFP and AI/EFP energy. Use efp_get_term_index to
 *                  get the term index in the energy array. NOTE that
 *                  polarization is computed self-consistently so we cannot
 *                  separate EFP/EFP and AI/EFP parts. If the AI/EFP part is
 *                  enabled it will be added to EFP/EFP part and corresponding
 *                  member in energy array will be zero.
 */
enum efp_result efp_get_energy(struct efp *efp, double *energy);

/* Get computed EFP energy gradient.
 *
 * n_grad (input): length of grad array. Must be greater than or equal to six
 *                 times number of fragments.
 *  grad (output): for each fragment contains x,y,z components of negative
 *                 force and torque.
 */
enum efp_result efp_get_gradient(struct efp *efp, int n_grad, double *grad);

/* Return number of fragments in this computation. */
int efp_get_frag_count(struct efp *efp);

/* Get number of atoms in the specified fragment. */
int efp_get_frag_atom_count(struct efp *efp, int frag_idx);

/* Get atoms comprising the specified fragment.
 *
 * frag_idx (input): index of a fragment between zero and total number of
 *                   fragments minus one.
 *   atoms (output): array of size returned by efp_get_frag_atom_count where
 *                   atom info will be stored.
 */
enum efp_result efp_get_frag_atoms(struct efp *efp, int frag_idx,
				   struct efp_atom *atoms);

/* Release all resources used by the efp. */
void efp_shutdown(struct efp *efp);

/* Convert efp_term flag to the index in array with energy terms. */
int efp_get_term_index(enum efp_term term);

/* Convert result to a human-readable message. */
const char *efp_result_to_string(enum efp_result res);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LIBEFP_EFP_H */
