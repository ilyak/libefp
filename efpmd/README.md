# EFPMD

EFPMD is a molecular simulation program based on LIBEFP. It supports EFP-only
single point energy and gradient calculations, semi-numerical Hessian, geometry
optimization, and molecular dynamics simulations in microcanonical and
canonical ensembles.

Simulations can be accelerated by running in parallel mode on multi-core CPUs.
To enable parallel computation set `OMP_NUM_THREADS` environmental variable to
the desired number of parallel threads. For example on 4-core machine do this
in shell before starting the program:

	export OMP_NUM_THREADS=4

This will enable parallel computation using 4 threads and should give almost 4x
speedup in all calculations.

## Input file format

Lines beginning with the `#` symbol are ignored during input parsing.

### Generic parameters

##### Type of the simulation

`run_type [sp|grad|hess|opt|md]`

`sp` - single point energy calculation.

`grad` - energy gradient calculation.

`hess` - semi-numerical Hessian calculation and normal mode analysis.

`opt` - geometry optimization.

`md` - molecular dynamics simulation.

Default value: `sp`

##### Format of fragment input

`coord [xyzabc|points|rotmat]`

`xyzabc` - Coordinates of the center of mass and Euler angles.

`points` - Coordinates of three atoms for each fragment.

`rotmat` - Coordinates of the center of mass and rotation matrix.

Default value: `xyzabc`

See fragment input specification for more details.

##### Units used for geometry specification

`units [angs|bohr]`

`angs` - Length unit is Angstrom.

`bohr` - Length unit is Bohr.

Default value: `angs`

##### Energy terms for EFP computation

`terms [elec [pol [disp [xr]]]]`

`elec` - Include electrostatics energy.

`pol` - Include polarization energy.

`disp` - Include dispersion energy.

`xr` - Include exchange repulsion energy.

Default value: `elec pol disp xr`

##### Electrostatic damping type

`elec_damp [screen|overlap|off]`

`screen` - Damping formula based on SCREEN group in the EFP potential.

`overlap` - Overlap-based damping formula. This damping correction is printed
as charge penetration energy.

`off` - No electrostatic damping.

Default value: `screen`

##### Dispersion damping type

`disp_damp [tt|overlap|off]`

`tt` - Damping based on the formula by Tang and Toennies.

`overlap` - Overlap-based dispersion damping.

`off` - No dispersion damping.

Default value: `tt`

##### Polarization damping type

`pol_damp [tt|off]`

`tt` - Tang and Toennies like damping formula.

`off` - No polarization damping.

Default value: `tt`

##### Enable/Disable the cutoff for fragment/fragment interactions

`enable_cutoff [true|false]`

Default value: `false`

##### Cutoff distance for fragment/fragment interactions

`swf_cutoff <value>`

Default value: `6.0`

Unit: Angstrom

##### Maximum number of steps to make

`max_steps <number>`

Default value: `100`

This specifies maximum number of steps for both geometry optimization and
molecular dynamics.

##### The path to the directory with fragment library

`fraglib_path <path>`

Default value: `"$(prefix)/share/libefp"` (data install directory)

The `<path>` parameter should not contain spaces or should be inside double
quotes otherwise.

##### The path to the directory with user-created fragments

`userlib_path <path>`

Default value: `"."` (current directory)

The `<path>` parameter should not contain spaces or should be in double quotes
otherwise.

### Periodic Boundary Conditions (PBC)

##### Enable/Disable PBC

`enable_pbc [true|false]`

Default value: `false`

Setting `enable_pbc` to `true` also sets `enable_cutoff` to `true`.

##### Periodic Box Size

`pbc_box <x> <y> <z>`

Default value: `15.0 15.0 15.0`

Unit: Angstrom

The smallest box dimension must be greater than `2 * swf_cutoff`.

### Geometry optimization related parameters

##### Optimization tolerance

`opt_tol <value>`

Default value: `1.0e-4`

Optimization will stop when maximum gradient component is less than `opt_tol`
and RMS gradient is less than one third of `opt_tol`.

### Molecular dynamics related parameters

##### Ensemble

`ensemble [nve|nvt]`

`nve` - Microcanonical ensemble.

`nvt` - Canonical ensemble with Nose-Hoover thermostat.

Default value: `nve`

##### Simulation temperature

`temperature <value>`

Unit: Kelvin

Default value: `300.0`

Target simulation temperature for NVT thermostat.

##### Time step

`time_step <value>`

Unit: Femtosecond

Default value: `1.0`

##### Print step

`print_step <value>`

Default value: `1`

Number of steps between outputs of the system state.

##### Thermostat parameter

`thermostat_tau <value>`

Unit: Femtosecond

Default value: `1.0e3`

Parameter of Nose-Hoover thermostat.

### Fragment input

One or more `fragment <name>` groups.

If `<name>` contains an `_l` suffix the parameters for this fragment will be
searched for in the `fraglib_path` directory. Otherwise the directory specified
by the `userlib_path` option will be used.

Each group is followed by the specification of fragment position and
orientation.

##### Format of input when `coord` is `xyzabc`

	fragment h2o
		0.0 0.0 0.0 0.0 0.0 0.0

The numbers are coordinates of the center of mass of a fragment and three Euler
rotation angels.

##### Format of input when `coord` is `points`

	fragment h2o
		0.0 0.0 0.0
		1.0 0.0 0.0
		0.0 1.0 0.0

The numbers are coordinates of three points belonging to a fragment.

##### Format of input when `coord` is `rotmat`

	fragment h2o
		0.0 0.0 0.0
		1.0 0.0 0.0
		0.0 1.0 0.0
		0.0 0.0 1.0

The numbers are coordinates of the center of mass of a fragment and a 3 x 3
rotation matrix.
