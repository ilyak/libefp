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

/** \file efp.h
 * Public libefp interface.
 */

#ifdef __cplusplus
extern "C" {
#endif

/** Result of an operation. */
enum efp_result {
	/** Operation was successful. */
	EFP_RESULT_SUCCESS = 0,
	/** Insufficient memory. */
	EFP_RESULT_NO_MEMORY,
	/** Operation is not implemented. */
	EFP_RESULT_NOT_IMPLEMENTED,
	/** Argument to function was invalid. */
	EFP_RESULT_INVALID_ARGUMENT,
	/** EFP structure was not properly initialized. */
	EFP_RESULT_NOT_INITIALIZED,
	/** File not found on disk. */
	EFP_RESULT_FILE_NOT_FOUND,
	/** Syntax error in EFP potential data file. */
	EFP_RESULT_SYNTAX_ERROR,
	/** The basis name in EFP potential data file was not specified. */
	EFP_RESULT_BASIS_NOT_SPECIFIED,
	/** Unknown fragment type. */
	EFP_RESULT_UNKNOWN_FRAGMENT,
	/** Loaded EFP parameters contain duplicate fragments. */
	EFP_RESULT_DUPLICATE_PARAMETERS,
	/** Required callback was not set. */
	EFP_RESULT_CALLBACK_NOT_SET,
	/** Call to callback failed. */
	EFP_RESULT_CALLBACK_FAILED,
	/** Exchange-repulsion must be turned on.
	 * Required for overlap-based electrostatic and dispersion damping. */
	EFP_RESULT_OVERLAP_INTEGRALS_REQUIRED,
	/** Gradient computation was not requested. */
	EFP_RESULT_GRADIENT_NOT_REQUESTED,
	/** Certain EFP parameters are missing. */
	EFP_RESULT_PARAMETERS_MISSING,
	/** Wrong array length. */
	EFP_RESULT_INVALID_ARRAY_SIZE,
	/** Unsupported SCREEN group in EFP potential data file. */
	EFP_RESULT_UNSUPPORTED_SCREEN,
	/** Inconsistent selection of EFP terms. */
	EFP_RESULT_INCONSISTENT_TERMS
};

/** Total number of EFP terms. */
#define EFP_TERM_COUNT 10

/** Flags to specify EFP energy terms. */
enum efp_term {
	EFP_TERM_ELEC = 1 << 0,     /**< EFP/EFP electrostatics */
	EFP_TERM_POL = 1 << 1,      /**< EFP/EFP polarization */
	EFP_TERM_DISP = 1 << 2,     /**< EFP/EFP dispersion */
	EFP_TERM_XR = 1 << 3,       /**< EFP/EFP exchange repulsion */
	EFP_TERM_CHTR = 1 << 4,     /**< EFP/EFP charge transfer */
	EFP_TERM_AI_ELEC = 1 << 5,  /**< Ab initio/EFP electrostatics */
	EFP_TERM_AI_POL = 1 << 6,   /**< Ab initio/EFP polarization */
	EFP_TERM_AI_DISP = 1 << 7,  /**< Ab initio/EFP dispersion */
	EFP_TERM_AI_XR = 1 << 8,    /**< Ab initio/EFP exchange repulsion */
	EFP_TERM_AI_CHTR = 1 << 9   /**< Ab initio/EFP charge transfer */
};

/** Fragment-fragment dispersion damping type. */
enum efp_disp_damp {
	EFP_DISP_DAMP_OVERLAP = 0, /**< Overlap-based damping (default). */
	EFP_DISP_DAMP_TT           /**< Tang-Toennies damping. */
};

/** Fragment-fragment electrostatic damping type. */
enum efp_elec_damp {
	EFP_ELEC_DAMP_SCREEN = 0,  /**< SCREEN controlled damping (default). */
	EFP_ELEC_DAMP_OVERLAP      /**< Overlap-based damping. */
};

/** \struct efp
 * Main EFP opaque structure.
 */
struct efp;

/** Options controlling EFP computation. */
struct efp_opts {
	/** Specifies which energy terms to compute. */
	enum efp_term terms;

	/** Set to nonzero to request gradient computation. */
	int do_gradient;

	/** Dispersion damping type. */
	enum efp_disp_damp disp_damp;

	/** Electrostatic damping type. */
	enum efp_elec_damp elec_damp;
};

/** EFP atom info. */
struct efp_atom {
	char label[32];   /**< Atom label. */
	char basis[32];   /**< Name of a basis on this atom. */
	double x;         /**< X coordinate of atom position. */
	double y;         /**< Y coordinate of atom position. */
	double z;         /**< Z coordinate of atom position. */
	double mass;      /**< Atom mass. */
	double znuc;      /**< Nuclear charge. */
};

/** QM atom info. */
struct efp_qm_atom {
	double x;         /**< X coordinate of atom position. */
	double y;         /**< Y coordinate of atom position. */
	double z;         /**< Z coordinate of atom position. */
	double znuc;      /**< Effective nuclear charge. */
};

/** Information about ab initio region. */
struct efp_qm_data {
	int n_atoms;                /**< Number of atoms in QM part. */
	struct efp_qm_atom *atoms;  /**< Atom data. */
};

/**
 * Information about a block of atoms for overlap and kinetic energy integral
 * computation by efp_callbacks::get_st_integrals callback function.
 *
 * The block represents a rectangular matrix. The
 * efp_callbacks::get_st_integrals callback computes overlap and kinetic energy
 * integrals for each basis function from horizontal dimension of a block with
 * each basis function from vertical dimension of a block.
 */
struct efp_st_block {
	/** Number of atoms in block horizontal dimension. */
	int n_atoms_i;

	/** Atoms of horizontal block dimension. */
	struct efp_atom *atoms_i;

	/** Number of atoms in block vertical dimension. */
	int n_atoms_j;

	/** Atoms of vertical block dimension. */
	struct efp_atom *atoms_j;
};

/** Overlap and kinetic energy integrals data. */
struct efp_st_data {
	double *s;   /**< Overlap integrals matrix. */
	double *sx;  /**< X derivative of overlap integrals. */
	double *sy;  /**< Y derivative of overlap integrals. */
	double *sz;  /**< Z derivative of overlap integrals. */
	double *t;   /**< Kinetic energy integrals matrix. */
	double *tx;  /**< X derivative of kinetic energy integrals. */
	double *ty;  /**< Y derivative of kinetic energy integrals. */
	double *tz;  /**< Z derivative of kinetic energy integrals. */
	int size_i;  /**< Number of columns in integral matrices. */
	int size_j;  /**< Number of rows in integral matrices. */
};

/**
 * Compute overlap and kinetic energy integrals over Gaussian basis functions.
 *
 * \param[in] block Block information (see efp_st_block).
 * \param[in] compute_derivatives If nonzero also compute integral derivatives.
 * \param[out] st Structure where integrals should be stored (see efp_st_data).
 * \param[in] user_data User data which was specified during initialization.
 *
 * \return ::EFP_RESULT_CALLBACK_FAILED on error or ::EFP_RESULT_SUCCESS
 *         otherwise.
 */
typedef enum efp_result (*efp_st_integrals_fn)(const struct efp_st_block *block,
					       int compute_derivatives,
					       struct efp_st_data *st,
					       void *user_data);

/**
 * Compute electric field form electrons in the specified points.
 *
 * \param[in] n_pt Number of points in \a xyz array.
 * \param[in] xyz Coordinates of points where electric field should be computed.
 * \param[out] field Computed \a x \a y \a z components of electric field.
 * \param[in] user_data User data which was specified during initialization.
 *
 * \return ::EFP_RESULT_CALLBACK_FAILED on error or ::EFP_RESULT_SUCCESS
 *         otherwise.
 */
typedef enum efp_result (*efp_electron_density_field_fn)(int n_pt,
							 const double *xyz,
							 double *field,
							 void *user_data);

/** XXX */
typedef enum efp_result (*efp_ao_integrals_fn)(int n_charges,
					       const double *charges,
					       const double *xyz,
					       int n_basis,
					       double *ints,
					       void *user_data);

/** XXX */
typedef enum efp_result (*efp_field_integrals_fn)(double *z,
						  int order,
						  void *user_data);

/** Callback function information. */
struct efp_callbacks {
	/** Callback function, see ::efp_st_integrals_fn. */
	efp_st_integrals_fn get_st_integrals;

	/** Will be passed as a last parameter to
	 * efp_callbacks::get_st_integrals. */
	void *get_st_integrals_user_data;

	/** Callback function, see ::efp_electron_density_field_fn. */
	efp_electron_density_field_fn get_electron_density_field;

	/** Will be passed as a last parameter to
	 * efp_callbacks::get_electron_density_field. */
	void *get_electron_density_field_user_data;

	/** Callback function, see ::efp_ao_integrals_fn . */
	efp_ao_integrals_fn get_ao_integrals;

	/** Will be passed as a last parameter to
	 * efp_callbacks::get_ao_integrals. */
	void *get_ao_integrals_user_data;

	/** Callback function, see ::efp_field_integrals_fn . */
	efp_field_integrals_fn get_field_integrals;

	/** Will be passed as a last parameter to
	 * efp_callbacks::get_field_integrals. */
	void *get_field_integrals_user_data;
};

/**
 * Get a printable string with the information about the library.
 *
 * \return Banner string, zero-terminated.
 */
const char *efp_banner(void);

/**
 * Set options to default values.
 *
 * \param[out] opts Structure to store default option values.
 */
void efp_opts_default(struct efp_opts *opts);

/**
 * Initialize the EFP computation.
 *
 * \param[out] out Initialized efp structure.
 * \param[in] opts User defined options controlling the computation.
 * \param[in] callbacks User supplied callback functions (see efp_callbacks).
 * \param[in] potential_files NULL-terminated list with paths to the files with
 *                            Effective Fragment Potential data.
 * \param[in] fragname NULL-terminated list with names of the fragments
 *                     comprising the molecular system.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_init(struct efp **out,
			 const struct efp_opts *opts,
			 const struct efp_callbacks *callbacks,
			 const char **potential_files,
			 const char **fragname);

/**
 * Update information about ab initio region.
 *
 * \param[in] efp The efp structure.
 * \param[in] qm_data Information about ab initio region.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_update_qm_data(struct efp *efp,
				   const struct efp_qm_data *qm_data);

/**
 * Update positions and orientations of effective fragments.
 *
 * \param[in] efp The efp structure.
 * \param[in] xyzabc for each fragment specifies \a x \a y \a z components of
 *                   fragment center of mass position and three Euler rotation
 *                   angles representing orientation of a fragment.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_update_fragments(struct efp *efp, const double *xyzabc);

/**
 * Update positions and orientations of effective fragments.
 *
 * This is convenience function. It does the same as efp_update_fragments.
 * However to specify position and orientation of fragments it takes
 * coordinates of 3 points in space for each fragment. For each fragment point
 * 1 and first atom of fragment are made to coincide. The vector connecting
 * points 1 and 2 is aligned with the corresponding vector connecting fragment
 * atoms. The plane defined by points 1, 2, and 3 is made to coincide with the
 * corresponding fragment plane.
 *
 * \param[in] efp The efp structure.
 * \param[in] pts Array of 9 times the number of fragments numbers. For each
 *                fragment specifies \x \y \z coordinates of 3 points to
 *                determine the position and orientation of a corresponding
 *                fragment.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_update_fragments_2(struct efp *efp, const double *pts);

/**
 * Initialize SCF computation. Must be called before ab initio SCF cycle.
 *
 * \param[in] efp The efp structure.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_scf_init(struct efp *efp);

/**
 * Update wave function dependent terms. Must be called during ab initio SCF.
 *
 * \param[in] efp The efp structure.
 * \param[out] energy Wave function dependent EFP energy.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_scf_update(struct efp *efp, double *energy);

/**
 * Perform the EFP computation.
 *
 * This call can take long time to complete.
 *
 * \param[in] efp The efp structure.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_compute(struct efp *efp);

/**
 * Compute one-electron contributions to QM Hamiltonian from EFP fragments.
 *
 * \param[in] efp The efp structure.
 * \param[in] n_basis Number of basis functions.
 * \param[out] v One-electron contributions to QM Hamiltonian.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_qm_contribution(struct efp *efp, int n_basis, double *v);

/**
 * Get computed energies of corresponding EFP terms.
 *
 * \param[in] efp The efp structure.
 * \param[out] energy Array of size #EFP_TERM_COUNT where energy terms will be
 *                    stored. The sum of the elements of this array will give
 *                    total EFP/EFP and AI/EFP energy. Use efp_get_term_index
 *                    to get the term index in the energy array. NOTE that
 *                    polarization is computed self-consistently so we cannot
 *                    separate EFP/EFP and AI/EFP parts. If the AI/EFP part is
 *                    enabled it will be added to EFP/EFP part and
 *                    corresponding member in energy array will be zero.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_energy(struct efp *efp, double *energy);

/**
 * Get computed EFP energy gradient.
 *
 * \param[in] efp The efp structure.
 * \param[in] n_grad Length of the grad array. Must be greater than or equal to
 *                   six times number of fragments.
 * \param[out] grad For each fragment contains \a x \a y \a z components of
 *                  negative force and torque.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_gradient(struct efp *efp, int n_grad, double *grad);

/**
 * Return number of fragments in this computation.
 *
 * \param[in] efp The efp structure.
 *
 * \return Number of EFP fragments in this computation.
 */
int efp_get_frag_count(struct efp *efp);

/**
 * Get number of atoms in the specified fragment.
 *
 * \param[in] efp The efp structure.
 * \param[in] frag_idx Index of a fragment between zero and total number of
 *                     fragments minus one.
 *
 * \return Number of atoms in the specified fragment.
 */
int efp_get_frag_atom_count(struct efp *efp, int frag_idx);

/**
 * Get atoms comprising the specified fragment.
 *
 * \param[in] efp The efp structure.
 * \param[in] frag_idx Index of a fragment between zero and total number of
 *                     fragments minus one.
 * \param[out] atoms Array of size returned by efp_get_frag_atom_count where
 *                   atom info will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_atoms(struct efp *efp, int frag_idx,
				   struct efp_atom *atoms);

/** Release all resources used by this \a efp.
 *
 * \param[in] efp The efp structure to be released.
 */
void efp_shutdown(struct efp *efp);

/**
 * Convert #efp_term flag to the index in array of energy terms.
 *
 * \param term Input flag.
 *
 * \return Index of the flag in array of energy terms.
 */
int efp_get_term_index(enum efp_term term);

/**
 * Convert #efp_result to a human readable message.
 *
 * \param res Result value to be converted to string.
 *
 * \return Human readable string with description of the result,
 *         zero-terminated.
 */
const char *efp_result_to_string(enum efp_result res);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* LIBEFP_EFP_H */
