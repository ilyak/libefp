! Fortran interface for libefp

module efp
use iso_c_binding
implicit none

type, bind(c) :: efp_opts
  integer(kind=c_int) terms
  integer(kind=c_int) disp_damp
  integer(kind=c_int) elec_damp
  integer(kind=c_int) pol_damp
  integer(kind=c_int) pol_driver
  integer(kind=c_int) enable_pbc
  integer(kind=c_int) enable_cutoff
  real(kind=c_double) swf_cutoff
end type efp_opts

type, bind(c) :: efp_energy
  real(kind=c_double) :: electrostatic
  real(kind=c_double) :: charge_penetration
  real(kind=c_double) :: electrostatic_point_charges
  real(kind=c_double) :: polarization
  real(kind=c_double) :: dispersion
  real(kind=c_double) :: ai_dispersion
  real(kind=c_double) :: exchange_repulsion
  real(kind=c_double) :: total
end type efp_energy

type, bind(c) :: efp_atom
  character(kind=c_char) :: label(32)
  real(kind=c_double) :: x, y, z
  real(kind=c_double) :: mass
  real(kind=c_double) :: znuc
end type efp_atom

interface

! const char *efp_banner(void);
function efp_banner() bind(c)
  use iso_c_binding, only: c_ptr
  type(c_ptr) :: efp_banner
end function

! void efp_print_banner(void);
subroutine efp_print_banner() bind(c)
end subroutine

! efp_t *efp_create(void);
function efp_create() bind(c)
  use iso_c_binding, only: c_ptr
  type(c_ptr) :: efp_create
end function

! void efp_opts_default(struct efp_opts *opts);
subroutine efp_opts_default(opts) bind(c)
  use iso_c_binding, only: c_ptr
  type(c_ptr), value :: opts
end subroutine

! efp_result_t efp_set_opts(struct efp *efp, const struct efp_opts *opts);
function efp_set_opts(efp, opts) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_set_opts
  type(c_ptr), value :: efp
  type(c_ptr), value :: opts
end function

! efp_result_t efp_get_opts(struct efp *efp, struct efp_opts *opts);
function efp_get_opts(efp, opts) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_opts
  type(c_ptr), value :: efp
  type(c_ptr), value :: opts
end function

! efp_result_t efp_add_potential(struct efp *efp, const char *path);
function efp_add_potential(efp, path) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_char
  integer(c_int) :: efp_add_potential
  type(c_ptr), value :: efp
  character(kind=c_char), dimension(*) :: path
end function

! efp_result_t efp_add_fragment(struct efp *efp, const char *name);
function efp_add_fragment(efp, name) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_char
  integer(c_int) :: efp_add_fragment
  type(c_ptr), value :: efp
  character(kind=c_char), dimension(*) :: name
end function

! efp_result_t efp_prepare(struct efp *efp);
function efp_prepare(efp) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_prepare
  type(c_ptr), value :: efp
end function

! efp_result_t efp_skip_fragments(struct efp *efp, size_t i, size_t j, int value);
function efp_skip_fragments(efp, i, j, value) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_skip_fragments
  type(c_ptr), value :: efp
  integer(c_size_t), value :: i
  integer(c_size_t), value :: j
  integer(c_int), value :: value
end function

! efp_result_t efp_set_electron_density_field_fn(struct efp *efp, efp_electron_density_field_fn fn);
function efp_set_electron_density_field_fn(efp, fn) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_funptr
  integer(c_int) :: efp_set_electron_density_field_fn
  type(c_ptr), value :: efp
  type(c_funptr), value :: fn
end function

! efp_result_t efp_set_electron_density_field_user_data(struct efp *efp, void *user_data);
function efp_set_electron_density_field_user_data(efp, user_data) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_set_electron_density_field_user_data
  type(c_ptr), value :: efp
  type(c_ptr), value :: user_data
end function

! efp_result_t efp_set_point_charges(struct efp *efp, size_t n_ptc, const double *ptc, const double *xyz);
function efp_set_point_charges(efp, n_ptc, ptc, xyz) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_set_point_charges
  type(c_ptr), value :: efp
  integer(c_size_t), value :: n_ptc
  type(c_ptr), value :: ptc
  type(c_ptr), value :: xyz
end function

! efp_result_t efp_get_point_charge_count(struct efp *efp, size_t *n_ptc);
function efp_get_point_charge_count(efp, n_ptc) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_point_charge_count
  type(c_ptr), value :: efp
  type(c_ptr), value :: n_ptc
end function

! efp_result_t efp_get_point_charge_values(struct efp *efp, double *ptc);
function efp_get_point_charge_values(efp, ptc) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_point_charge_values
  type(c_ptr), value :: efp
  type(c_ptr), value :: ptc
end function

! efp_result_t efp_set_point_charge_values(struct efp *efp, const double *ptc);
function efp_set_point_charge_values(efp, ptc) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_set_point_charge_values
  type(c_ptr), value :: efp
  type(c_ptr), value :: ptc
end function

! efp_result_t efp_get_point_charge_coordinates(struct efp *efp, double *xyz);
function efp_get_point_charge_coordinates(efp, xyz) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_point_charge_coordinates
  type(c_ptr), value :: efp
  type(c_ptr), value :: xyz
end function

! efp_result_t efp_set_point_charge_coordinates(struct efp *efp, const double *xyz);
function efp_set_point_charge_coordinates(efp, xyz) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_set_point_charge_coordinates
  type(c_ptr), value :: efp
  type(c_ptr), value :: xyz
end function

! efp_result_t efp_get_point_charge_gradient(struct efp *efp, double *grad);
function efp_get_point_charge_gradient(efp, grad) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_point_charge_gradient
  type(c_ptr), value :: efp
  type(c_ptr), value :: grad
end function

! efp_result_t efp_set_coordinates(struct efp *efp, efp_coord_type_t coord_type, const double *coord);
function efp_set_coordinates(efp, coord_type, coord) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_set_coordinates
  type(c_ptr), value :: efp
  integer(c_int), value :: coord_type
  type(c_ptr), value :: coord
end function

! efp_result_t efp_set_frag_coordinates(struct efp *efp, size_t frag_idx, efp_coord_type_t coord_type, const double *coord);
function efp_set_frag_coordinates(efp, frag_idx, coord_type, coord) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_set_frag_coordinates
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  integer(c_int), value :: coord_type
  type(c_ptr), value :: coord
end function

! efp_result_t efp_get_coordinates(struct efp *efp, double *xyzabc);
function efp_get_coordinates(efp, xyzabc) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_coordinates
  type(c_ptr), value :: efp
  type(c_ptr), value :: xyzabc
end function

! efp_result_t efp_get_frag_xyzabc(struct efp *efp, size_t frag_idx, double *xyzabc);
function efp_get_frag_xyzabc(efp, frag_idx, xyzabc) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_frag_xyzabc
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  type(c_ptr), value :: xyzabc
end function

! efp_result_t efp_set_periodic_box(struct efp *efp, double x, double y, double z);
function efp_set_periodic_box(efp, x, y, z) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_double
  integer(c_int) :: efp_set_periodic_box
  type(c_ptr), value :: efp
  real(c_double), value :: x
  real(c_double), value :: y
  real(c_double), value :: z
end function

! efp_result_t efp_get_stress_tensor(struct efp *efp, double *stress);
function efp_get_stress_tensor(efp, stress) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_stress_tensor
  type(c_ptr), value :: efp
  type(c_ptr), value :: stress
end function

! efp_result_t efp_get_ai_screen(struct efp *efp, size_t frag_idx, double *screen);
function efp_get_ai_screen(efp, frag_idx, screen) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_ai_screen
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  type(c_ptr), value :: screen
end function

! efp_result_t efp_set_orbital_energies(struct efp *efp, size_t n_core, size_t n_act, size_t n_vir, const double *oe);
function efp_set_orbital_energies(efp, n_core, n_act, n_vir, oe) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_set_orbital_energies
  type(c_ptr), value :: efp
  integer(c_size_t), value :: n_core
  integer(c_size_t), value :: n_act
  integer(c_size_t), value :: n_vir
  type(c_ptr), value :: oe
end function

! efp_result_t efp_set_dipole_integrals(struct efp *efp, size_t n_core, size_t n_act, size_t n_vir, const double *dipint);
function efp_set_dipole_integrals(efp, n_core, n_act, n_vir, dipint) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_set_dipole_integrals
  type(c_ptr), value :: efp
  integer(c_size_t), value :: n_core
  integer(c_size_t), value :: n_act
  integer(c_size_t), value :: n_vir
  type(c_ptr), value :: dipint
end function

! efp_result_t efp_get_wavefunction_dependent_energy(struct efp *efp, double *energy);
function efp_get_wavefunction_dependent_energy(efp, energy) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_wavefunction_dependent_energy
  type(c_ptr), value :: efp
  type(c_ptr), value :: energy
end function

! efp_result_t efp_compute(struct efp *efp, int do_gradient);
function efp_compute(efp, do_gradient) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_compute
  type(c_ptr), value :: efp
  integer(c_int), value :: do_gradient
end function

! efp_result_t efp_get_frag_charge(struct efp *efp, size_t frag_idx, double *charge);
function efp_get_frag_charge(efp, frag_idx, charge) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_frag_charge
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  type(c_ptr), value :: charge
end function

! efp_result_t efp_get_frag_multiplicity(struct efp *efp, size_t frag_idx, int *mult);
function efp_get_frag_multiplicity(efp, frag_idx, mult) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_frag_multiplicity
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  type(c_ptr), value :: mult
end function

! efp_result_t efp_get_frag_multipole_count(struct efp *efp, size_t frag_idx, size_t *n_mult);
function efp_get_frag_multipole_count(efp, frag_idx, n_mult) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_frag_multipole_count
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  type(c_ptr), value :: n_mult
end function

! efp_result_t efp_get_multipole_count(struct efp *efp, size_t *n_mult);
function efp_get_multipole_count(efp, n_mult) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_multipole_count
  type(c_ptr), value :: efp
  type(c_ptr), value :: n_mult
end function

! efp_result_t efp_get_multipole_coordinates(struct efp *efp, double *xyz);
function efp_get_multipole_coordinates(efp, xyz) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_multipole_coordinates
  type(c_ptr), value :: efp
  type(c_ptr), value :: xyz
end function

! efp_result_t efp_get_multipole_values(struct efp *efp, double *mult);
function efp_get_multipole_values(efp, mult) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_multipole_values
  type(c_ptr), value :: efp
  type(c_ptr), value :: mult
end function

! efp_result_t efp_get_induced_dipole_count(struct efp *efp, size_t *n_dip);
function efp_get_induced_dipole_count(efp, n_dip) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_induced_dipole_count
  type(c_ptr), value :: efp
  type(c_ptr), value :: n_dip
end function

! efp_result_t efp_get_induced_dipole_coordinates(struct efp *efp, double *xyz);
function efp_get_induced_dipole_coordinates(efp, xyz) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_induced_dipole_coordinates
  type(c_ptr), value :: efp
  type(c_ptr), value :: xyz
end function

! efp_result_t efp_get_induced_dipole_values(struct efp *efp, double *dip);
function efp_get_induced_dipole_values(efp, dip) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_induced_dipole_values
  type(c_ptr), value :: efp
  type(c_ptr), value :: dip
end function

! efp_result_t efp_get_induced_dipole_conj_values(struct efp *efp, double *dip);
function efp_get_induced_dipole_conj_values(efp, dip) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_induced_dipole_conj_values
  type(c_ptr), value :: efp
  type(c_ptr), value :: dip
end function

! efp_result_t efp_get_lmo_count(struct efp *efp, size_t frag_idx, size_t *n_lmo);
function efp_get_lmo_count(efp, frag_idx, n_lmo) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_lmo_count
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  type(c_ptr), value :: n_lmo
end function

! efp_result_t efp_get_lmo_coordinates(struct efp *efp, size_t frag_idx, double *xyz);
function efp_get_lmo_coordinates(efp, frag_idx, xyz) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_lmo_coordinates
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  type(c_ptr), value :: xyz
end function

! efp_result_t efp_get_xrfit(struct efp *efp, size_t frag_idx, double *xrfit);
function efp_get_xrfit(efp, frag_idx, xrfit) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_xrfit
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  type(c_ptr), value :: xrfit
end function

! efp_result_t efp_get_energy(struct efp *efp, struct efp_energy *energy);
function efp_get_energy(efp, energy) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_energy
  type(c_ptr), value :: efp
  type(c_ptr), value :: energy
end function

! efp_result_t efp_get_gradient(struct efp *efp, double *grad);
function efp_get_gradient(efp, grad) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_gradient
  type(c_ptr), value :: efp
  type(c_ptr), value :: grad
end function

! efp_result_t efp_get_frag_count(struct efp *efp, size_t *n_frag);
function efp_get_frag_count(efp, n_frag) bind(c)
  use iso_c_binding, only: c_int, c_ptr
  integer(c_int) :: efp_get_frag_count
  type(c_ptr), value :: efp
  type(c_ptr), value :: n_frag
end function

! efp_result_t efp_get_frag_name(struct efp *efp, size_t frag_idx, size_t size, char *frag_name);
function efp_get_frag_name(efp, frag_idx, size, frag_name) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t, c_char
  integer(c_int) :: efp_get_frag_name
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  integer(c_size_t), value :: size
  character(kind=c_char), dimension(*) :: frag_name
end function

! efp_result_t efp_get_frag_mass(struct efp *efp, size_t frag_idx, double *mass_out);
function efp_get_frag_mass(efp, frag_idx, mass_out) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_frag_mass
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  type(c_ptr), value :: mass_out
end function

! efp_result_t efp_get_frag_inertia(struct efp *efp, size_t frag_idx, double *inertia_out);
function efp_get_frag_inertia(efp, frag_idx, inertia_out) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_frag_inertia
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  type(c_ptr), value :: inertia_out
end function

! efp_result_t efp_get_frag_atom_count(struct efp *efp, size_t frag_idx, size_t *n_atoms);
function efp_get_frag_atom_count(efp, frag_idx, n_atoms) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_frag_atom_count
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  type(c_ptr), value :: n_atoms
end function

! efp_result_t efp_get_frag_atoms(struct efp *efp, size_t frag_idx, size_t size, struct efp_atom *atoms);
function efp_get_frag_atoms(efp, frag_idx, size, atoms) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_frag_atoms
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  integer(c_size_t), value :: size
  type(c_ptr), value :: atoms
end function

! efp_result_t efp_get_electric_field(struct efp *efp, size_t frag_idx, const double *xyz, double *field);
function efp_get_electric_field(efp, frag_idx, xyz, field) bind(c)
  use iso_c_binding, only: c_int, c_ptr, c_size_t
  integer(c_int) :: efp_get_electric_field
  type(c_ptr), value :: efp
  integer(c_size_t), value :: frag_idx
  type(c_ptr), value :: xyz
  type(c_ptr), value :: field
end function

! void efp_torque_to_derivative(const double *euler, const double *torque, double *deriv);
subroutine efp_torque_to_derivative(euler, torque, deriv) bind(c)
  use iso_c_binding, only: c_ptr
  type(c_ptr), value :: euler
  type(c_ptr), value :: torque
  type(c_ptr), value :: deriv
end subroutine

! void efp_shutdown(struct efp *efp);
subroutine efp_shutdown(efp) bind(c)
  use iso_c_binding, only: c_ptr
  type(c_ptr), value :: efp
end subroutine

end interface
end module efp
