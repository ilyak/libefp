This file contains a step-by-step instructions of how to integrate libefp into
an existing quantum chemistry package to get a working QM/EFP interface.

1. Initialization

	a. Print the libefp banner: efp_banner
	b. Create a EFP object: efp_create
	c. Setup your EFP options: efp_set_opts
	d. Load EFP parameters from ".efp" files: efp_add_potential
	e. Add your EFP fragments: efp_add_fragment
	f. Setup a callback function which computes electric field from
	   ab initio electrons: efp_set_electron_density_field_fn
	g. (optional)
	   Setup a callback function to get better error
	   messages: efp_set_error_log
	h. Call efp_prepare
	i. (optional, EFP-only)
	   Setup periodic box size. If you have enabled periodic
	   boundary conditions (options.enable_pbc = true), you should
	   call efp_set_periodic_box

2. Optimization or molecular dynamics step

	a. Set new fragment coordinates: efp_set_coordinates
	b. Update ab initio nuclei: efp_set_point_charges

2.1 Ab initio SCF cycle

	a. Start SCF cycle as usual
	b. Compute one-electron EFP contributions

		To compute one-electron contributions to quantum mechanical
		Hamiltonian you will need multipole moment integrals as well as
		electrostatics from EFP subsystem.

		To get the electrostatics use the following functions:

			atom charges from EFP fragments:
				efp_get_frag_atoms
			multipoles from EFP electrostatics:
				efp_get_multipoles

			polarization induced dipoles:
				efp_get_induced_dipole_values
			polarization conjugate induced dipoles:
				efp_get_induced_dipole_conj_values

		To compute QM/EFP interaction with induced dipoles use the
		average of induced dipoles and conjugate induced dipoles

	c. Compute wavefunction-dependent energy

		Because EFP polarization induced dipoles depend on the electric
		field from the wavefunction you have to converge them together.

		To compute wavefunction-dependent terms use:
			efp_get_wf_dependent_energy

	d. Add obtained wavefunction dependent energy to your SCF energy
	e. Go to step 2.1.a until SCF convergence criteria are met

2.2 Compute EFP energy and gradients

	a. Compute all EFP energy terms: efp_compute

2.3 Geometry update

	a. Get EFP energy: efp_get_energy

		Do not forget to subtract wavefunction dependent energy you have
		already added to your SCF energy from the total EFP energy

	b. Get gradient on EFP fragments: efp_get_gradient
	c. Get EFP gradient on ab initio nuclei: efp_get_point_charge_gradient
	d. Update system geometry

		Update all coordinates of your ab initio part and EFP fragments
		according to you optimization or molecular dynamics algorithm

	e. Go to the beginning of step 2

3. Cleanup resources

	a. Release all memory used by libefp: efp_shutdown
