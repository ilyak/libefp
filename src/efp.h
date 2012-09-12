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
	/** Unexpected NULL argument to function was specified. */
	EFP_RESULT_ARGUMENT_NULL,
	/** EFP structure was not properly initialized. */
	EFP_RESULT_NOT_INITIALIZED,
	/** File not found on disk. */
	EFP_RESULT_FILE_NOT_FOUND,
	/** Syntax error in EFP parameters file. */
	EFP_RESULT_SYNTAX_ERROR,
	/** Unknown fragment type. */
	EFP_RESULT_UNKNOWN_FRAGMENT,
	/** EFP parameters contain duplicate fragments. */
	EFP_RESULT_DUPLICATE_PARAMETERS,
	/** Required callback function was not set. */
	EFP_RESULT_CALLBACK_NOT_SET,
	/** Call to callback function failed. */
	EFP_RESULT_CALLBACK_FAILED,
	/** Gradient computation was not requested. */
	EFP_RESULT_GRADIENT_NOT_REQUESTED,
	/** Polarization SCF did not converge. */
	EFP_RESULT_POL_NOT_CONVERGED,
	/** Certain EFP parameters are missing. */
	EFP_RESULT_PARAMETERS_MISSING,
	/** Euler angle beta is out of range [0,pi]. */
	EFP_RESULT_BAD_EULER_B,
	/** Incorrect enumeration value. */
	EFP_RESULT_INCORRECT_ENUM_VALUE,
	/** Index is out of range. */
	EFP_RESULT_INDEX_OUT_OF_RANGE,
	/** Wrong array length. */
	EFP_RESULT_INVALID_ARRAY_SIZE,
	/** Unsupported SCREEN group in EFP parameters file. */
	EFP_RESULT_UNSUPPORTED_SCREEN,
	/**
	 * Inconsistent selection of EFP terms.
	 *
	 * This means that AI/EFP terms were selected without selecting their
	 * EFP/EFP counterparts. Enabling polarization without electrostatics
	 * also produces this error. */
	EFP_RESULT_INCONSISTENT_TERMS
};

/** Flags to specify EFP energy terms. */
enum efp_term {
	/** EFP/EFP electrostatics. */
	EFP_TERM_ELEC = 1 << 0,
	/** EFP/EFP polarization. */
	EFP_TERM_POL = 1 << 1,
	/** EFP/EFP dispersion. */
	EFP_TERM_DISP = 1 << 2,
	/** EFP/EFP exchange repulsion. */
	EFP_TERM_XR = 1 << 3,
	/** EFP/EFP charge transfer, reserved for future use. */
	EFP_TERM_CHTR = 1 << 4,
	/** Ab initio/EFP electrostatics. */
	EFP_TERM_AI_ELEC = 1 << 5,
	/** Ab initio/EFP polarization. */
	EFP_TERM_AI_POL = 1 << 6,
	/** Ab initio/EFP dispersion, reserved for future use. */
	EFP_TERM_AI_DISP = 1 << 7,
	/** Ab initio/EFP exchange repulsion, reserved for future use. */
	EFP_TERM_AI_XR = 1 << 8,
	/** Ab initio/EFP charge transfer, reserved for future use. */
	EFP_TERM_AI_CHTR = 1 << 9
};

/** Fragment-fragment dispersion damping type. */
enum efp_disp_damp {
	EFP_DISP_DAMP_OVERLAP = 0, /**< Overlap-based damping (default). */
	EFP_DISP_DAMP_TT,          /**< Tang-Toennies damping. */
	EFP_DISP_DAMP_OFF          /**< No dispersion damping. */
};

/** Fragment-fragment electrostatic damping type. */
enum efp_elec_damp {
	EFP_ELEC_DAMP_SCREEN = 0,  /**< SCREEN-controlled damping (default). */
	EFP_ELEC_DAMP_OVERLAP,     /**< Overlap-based damping. */
	EFP_ELEC_DAMP_OFF          /**< No electrostatic damping. */
};

/** Fragment-fragment polarization damping type. */
enum efp_pol_damp {
	EFP_POL_DAMP_TT = 0,       /**< Tang-Toennies like damping (default). */
	EFP_POL_DAMP_OFF           /**< No polarization damping. */
};

/** Describes the way fragment coordinates are specified. */
enum efp_coord_type {
	/** Coordinates of center of mass of a fragment and Euler angles. */
	EFP_COORD_TYPE_XYZABC = 0,
	/** Coordinates of three points belonging to a fragment. */
	EFP_COORD_TYPE_POINTS
};

/** Specifies output format of the gradient. */
enum efp_grad_type {
	/** Store torque. */
	EFP_GRAD_TYPE_TORQUE = 0,
	/** Store derivatives of energy by Euler angles. */
	EFP_GRAD_TYPE_DERIVATIVE
};

/** \struct efp
 * Main EFP opaque structure.
 */
struct efp;

/** Options controlling EFP computation. */
struct efp_opts {
	/** Specifies which energy terms to compute. */
	unsigned terms;

	/** Dispersion damping type (see #efp_disp_damp). */
	enum efp_disp_damp disp_damp;

	/** Electrostatic damping type (see #efp_elec_damp). */
	enum efp_elec_damp elec_damp;

	/** Polarization damping type (see #efp_pol_damp). */
	enum efp_pol_damp pol_damp;
};

/** EFP energy terms. */
struct efp_energy {
	/**
	 * Electrostatic energy. */
	double electrostatic;
	/**
	 * Charge penetration energy from overlap-based electrostatic
	 * damping. Zero if overlap-based damping is turned off. */
	double charge_penetration;
	/**
	 * All polarization energy goes here. Polarization is computed
	 * self-consistently so it can't be separated into EFP/EFP and AI/EFP
	 * parts. */
	double polarization;
	/**
	 * Dispersion energy. */
	double dispersion;
	/**
	 * Exchange-repulsion energy. */
	double exchange_repulsion;
	/**
	 * Charge transfer energy. */
	double charge_transfer;
	/**
	 * AI/EFP electrostatic energy. */
	double ai_electrostatic;
	/**
	 * AI/EFP dispersion energy. */
	double ai_dispersion;
	/**
	 * AI/EFP exchange-repulsion energy. */
	double ai_exchange_repulsion;
	/**
	 * AI/EFP charge transfer energy. */
	double ai_charge_transfer;
	/**
	 * Sum of all the above EFP energy terms. */
	double total;
};

/** EFP atom info. */
struct efp_atom {
	char label[32];   /**< Atom label. */
	double x;         /**< X coordinate of atom position. */
	double y;         /**< Y coordinate of atom position. */
	double z;         /**< Z coordinate of atom position. */
	double mass;      /**< Atom mass. */
	double znuc;      /**< Nuclear charge. */
};

/**
 * Compute electric field form electrons in the specified points.
 *
 * \param[in] n_pt Number of points in \p xyz array.
 * \param[in] xyz Coordinates of points where electric field should be
 *                computed. The size of this array must be [3 * \p n_pt]
 *                elements.
 * \param[out] field Computed \a x \a y \a z components of electric field.
 *                   The size of this array must be at least [3 * \p n_pt]
 *                   elements.
 * \param[in] user_data User data which was specified during initialization.
 *
 * \return ::EFP_RESULT_CALLBACK_FAILED on error or ::EFP_RESULT_SUCCESS
 *         otherwise.
 */
typedef enum efp_result (*efp_electron_density_field_fn)(int n_pt,
							 const double *xyz,
							 double *field,
							 void *user_data);

/** Callback function information. */
struct efp_callbacks {
	/** Callback function, see ::efp_electron_density_field_fn. */
	efp_electron_density_field_fn get_electron_density_field;

	/** Will be passed as a last parameter to
	 * efp_callbacks::get_electron_density_field. */
	void *get_electron_density_field_user_data;
};

/**
 * Get a human readable banner string with information about the library.
 *
 * \return Banner string, zero-terminated.
 */
const char *efp_banner(void);

/**
 * Get default values of simulation options.
 *
 * \param[out] opts Structure to store the defaults.
 */
void efp_opts_default(struct efp_opts *opts);

/**
 * Initialize the EFP computation.
 *
 * \param[out] out Initialized efp structure.
 * \param[in] opts User defined options controlling the computation.
 * \param[in] callbacks User supplied callback functions (see efp_callbacks).
 * \param[in] potential_file_list Zero-terminated string with paths to the EFP
 *                                parameter files separated by the new-line
 *                                character. The last line must not include the
 *                                new-line character.
 * \param[in] frag_name_list Zero-terminated string with names of the fragments
 *                           separated by the new-line character. The last line
 *                           must not include the new-line character.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_init(struct efp **out,
			 const struct efp_opts *opts,
			 const struct efp_callbacks *callbacks,
			 const char *potential_file_list,
			 const char *frag_name_list);

/**
 * Update information about \a ab \a initio subsystem.
 *
 * \param[in] efp The efp structure.
 * \param[in] n_atoms Number of atoms in \a ab \a initio subsystem.
 * \param[in] znuc Array of \p n_atoms elements with atom nuclear charges.
 * \param[in] xyz Array of [3 * \p n_atoms] elements with \a x \a y \a z
 *                coordinates of atoms.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_set_qm_atoms(struct efp *efp,
				 int n_atoms,
				 const double *znuc,
				 const double *xyz);

/**
 * Get number of atoms in \a ab \a initio subsystem.
 *
 * \param[in] efp The efp structure.
 * \param[out] n_atoms Number of atoms in \a ab \a initio subsystem.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_qm_atom_count(struct efp *efp, int *n_atoms);

/**
 * Get information about \a ab \a initio subsystem.
 *
 * \param[in] efp The efp structure.
 * \param[in] n_atoms Expected number of atoms in \a ab \a initio subsystem.
 * \param[out] znuc Array of \p n_atoms elements where atom nuclear charges
 *                  will be stored.
 * \param[out] xyz Array of [3 * \p n_atoms] elements where \a x \a y \a z
 *                 coordinates of atom positions will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_qm_atoms(struct efp *efp,
				 int n_atoms,
				 double *znuc,
				 double *xyz);

/**
 * Update positions and orientations of effective fragments.
 *
 * \param[in] efp The efp structure.
 * \param[in] coord_type Specifies the type of coordinates in the \a coord
 *                       array (see #efp_coord_type).
 * \param[in] coord Array of fragment coordinates.
 *
 * If \p coord_type is \a EFP_COORD_TYPE_XYZABC then for each fragment the \p
 * coord array should contain \a x \a y \a z components of the center of mass
 * position and three Euler rotation angles representing orientation of a
 * fragment. The size of the \p coord array must be at least [6 * \a n]
 * elements, where \a n is the number of fragments.
 *
 * If \p coord_type is \a EFP_COORD_TYPE_POINTS then for each fragment the \p
 * coord array should contain the coordinates of 3 points in space. For each
 * fragment point 1 and first atom of fragment are made to coincide. The vector
 * connecting points 1 and 2 is aligned with the corresponding vector
 * connecting fragment atoms. The plane defined by points 1, 2, and 3 is made
 * to coincide with the corresponding fragment plane. The size of the \p coord
 * array must be at least [9 * \a n] elements, where \a n is the number of
 * fragments.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_set_coordinates(struct efp *efp,
				    enum efp_coord_type coord_type,
				    const double *coord);

/**
 * Get center of mass positions and Euler angles of the effective fragments.
 *
 * \param[in] efp The efp structure.
 * \param[in] n_frags Expected number of fragments.
 * \param[out] xyzabc Upon return [6 * \p n_frags] elements will be written to
 *                    this array. The coordinates of the center of mass and
 *                    Euler rotation angles for each fragment will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_coordinates(struct efp *efp,
				    int n_frags,
				    double *xyzabc);

/**
 * Update wave function dependent energy terms.
 *
 * This function must be called during \a ab \a initio SCF.
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
 * This call may take long time to complete.
 *
 * \param[in] efp The efp structure.
 * \param[in] do_gradient Also compute energy gradient if nonzero value is
 *                        specified.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_compute(struct efp *efp, int do_gradient);

/**
 * Get total number of EFP multipoles.
 *
 * \param[in] efp The efp structure.
 * \param[out] n_mult Array of 4 integers where the total number of charges,
 *                    dipoles, quadrupoles and octupoles will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_multipole_count(struct efp *efp, int *n_mult);

/**
 * Get all electrostatic multipoles from all fragments.
 *
 * \param[in] efp The efp structure.
 *
 * \param[out] xyz Array of four pointers to arrays where the corresponding
 *                 multipole point coordinates will be stored. The size of
 *                 each array must be at least [3 * \p n_mult] elements, where
 *                 \p n_mult is the corresponding value in the array returned
 *                 by the ::efp_get_multipole_count function.
 *
 * \param[out] z Array of four pointers to arrays where the charges, dipoles,
 *               quadrupoles, and octupoles will be stored.
 *
 *               The size of the first array \p z[0] (charges) must be at least
 *               [\p n_mult] elements, where \p n_mult is the corresponding
 *               first element of the array returned by the
 *               ::efp_get_multipole_count function.
 *
 *               The size of the second array \p z[1] (dipoles) must be at
 *               least [3 * \p n_mult] elements, where \p n_mult is the
 *               corresponding second element of the array returned by the
 *               ::efp_get_multipole_count function.
 *
 *               The size of the third array \p z[2] (quadrupoles) must be at
 *               least [6 * \p n_mult] elements, where \p n_mult is the
 *               corresponding third element of the array returned by the
 *               ::efp_get_multipole_count function.
 *
 *               Quadrupoles are stored in the following order:
 *                   \a xx, \a yy, \a zz, \a xy, \a xz, \a yz
 *
 *               The size of the fourth array \p z[3] (octupoles) must be at
 *               least [10 * \p n_mult] elements, where \p n_mult is the
 *               corresponding fourth element of the array returned by the
 *               ::efp_get_multipole_count function.
 *
 *               Octupoles are stored in the following order:
 *                   \a xxx, \a yyy, \a zzz, \a xxy, \a xxz,
 *                   \a xyy, \a yyz, \a xzz, \a yzz, \a xyz
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_multipoles(struct efp *efp, double **xyz, double **z);

/**
 * Get computed energy components.
 *
 * \param[in] efp The efp structure.
 * \param[out] energy Computed EFP energy components (see efp_energy).
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_energy(struct efp *efp, struct efp_energy *energy);

/**
 * Get computed EFP energy gradient.
 *
 * \param[in] efp The efp structure.
 * \param[in] grad_type Format of the returned gradient.
 * \param[in] n_frags Expected number of fragments.
 * \param[out] grad For each fragment \a x \a y \a z components of negative
 *                  force and torque will be written to this array. The size
 *                  of this array must be at least [3 * \p n_frags] elements.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_gradient(struct efp *efp,
				 enum efp_grad_type grad_type,
				 int n_frags,
				 double *grad);

/**
 * Get the gradient on atoms in \a ab \a initio subsystem due to EFP
 * interaction.
 *
 * This does not include additional correction due to changes in the
 * \a ab \a initio wave function affected by EFP subsystem.
 *
 * \param[in] efp The efp structure.
 * \param[in] n_atoms Expected number of atoms in the \a ab \a initio
 *                    subsystem.
 * \param[out] grad For each atom \a x \a y \a z components of energy
 *                  gradient are stored. The size of this array must be
 *                  at least [3 * \p n_atoms] elements.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_qm_gradient(struct efp *efp,
				    int n_atoms,
				    double *grad);

/**
 * Get the number of fragments in this computation.
 *
 * \param[in] efp The efp structure.
 * \param[out] n_frag Total number of fragments in this simulation.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_count(struct efp *efp, int *n_frag);

/**
 * Get the name of the specified effective fragment.
 *
 * \param[in] efp The efp structure.
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 *                     the total number of fragments minus one.
 * \param[in] size Size of a \p frag_name buffer.
 * \param[out] frag_name A buffer where name of the fragment will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_name(struct efp *efp,
				  int frag_idx,
				  int size,
				  char *frag_name);

/**
 * Get the number of atoms in the specified fragment.
 *
 * \param[in] efp The efp structure.
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 *                     the total number of fragments minus one.
 * \param[out] n_atoms Total number of atoms in the fragment.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_atom_count(struct efp *efp,
					int frag_idx,
					int *n_atoms);

/**
 * Get atoms comprising the specified fragment.
 *
 * \param[in] efp The efp structure.
 * \param[in] frag_idx Index of a fragment. Must be a value between zero and
 *                     the total number of fragments minus one.
 * \param[in] size Size of the \p atoms array. Must be greater than or equal to
 *                 the size returned by the ::efp_get_frag_atom_count function.
 * \param[out] atoms Array where atom information will be stored.
 *
 * \return ::EFP_RESULT_SUCCESS on success or error code otherwise.
 */
enum efp_result efp_get_frag_atoms(struct efp *efp,
				   int frag_idx,
				   int size,
				   struct efp_atom *atoms);

/**
 * Release all resources used by this \a efp.
 *
 * \param[in] efp The efp structure to be released.
 */
void efp_shutdown(struct efp *efp);

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
