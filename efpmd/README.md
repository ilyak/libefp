# EFPMD

EFPMD is a molecular simulation program based on LIBEFP. It supports EFP-only
single point energy and gradient computations, geometry optimization, and
molecular dynamics simulations in microcanonical and canonical ensembles.

Simulations can be accelerated by running in parallel mode on multi-core CPUs.
To enable parallel computation set `OMP_NUM_THREADS` environmental variable to
the desired number of parallel threads. For example on 4-core machine in bash
shell do:

	export OMP_NUM_THREADS=4

before starting the program to enable parallel computation using 4 threads.
This should give almost 4x speedup in all computations.

## Input file format

Lines beginning with the `#` symbol are ignored during input parsing.

### Generic parameters

##### Type of the simulation

`run_type [sp|grad|cg|nve|nvt]`

`sp` - single point calculation.

`grad` - gradient calculation.

`cg` - geometry optimization using conjugate gradient method.

`nve` - molecular dynamics in microcanonical ensemble.

`nvt` - molecular dynamics in canonical ensemble.

Default value: `sp`

##### Format of fragment input

`coord [xyzabc|points]`

`xyzabc` - Coordinates of the center of mass and Euler angles.

`points` - Coordinates of three atoms for each fragment.

Default value: `xyzabc`

See fragment input specification for more details.

##### Geometry units

`units [angs|bohr]`

`angs` - Units are Angstroms.

`bohr` - Units are Bohr.

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

`overlap` - Overlap based damping which computes charge penetration energy.

`off` - No electrostatic damping.

Default value: `screen`

##### Dispersion damping type

`disp_damp [tt|overlap|off]`

`tt` - Damping based on the formula by Tang and Toennies.

`overlap` - Overlap-based dispersion damping.

`off` - No dispersion damping.

Default value: `tt`

##### Maximum number of steps to make

`max_steps <number>`

Default value: `100`

##### Number of steps between outputs

`print_step <number>`

Default value: `1`

##### The path to the directory with fragment library

`fraglib_path <path>`

Default value: `"$(prefix)/share/libefp"` (data install directory)

The `<path>` parameter should not contain spaces or should be in double quotes
otherwise.

##### The path to the directory with user-created fragments

`userlib_path <path>`

Default value: `"."` (current directory)

The `<path>` parameter should not contain spaces or should be in double quotes
otherwise.

### Geometry optimization related parameters

##### Line search step size

`ls_step_size <number>`

Default value: `20.0`

##### Optimization tolerance

`opt_tol <number>`

Default value: `1.0e-4`

Optimization will stop when maximum gradient component is less than `opt_tol`
and RMS gradient is less than one third of `opt_tol`.

### Molecular dynamics related parameters

##### Simulation temperature

`temperature <number>`

Units: `Kelvins`

Default value: `300.0`

##### Time step

`time_step <number>`

Units: `Femtoseconds`

Default value: `1.0`

### Fragment input

One or more `fragment <name>` groups.

Each `fragment <name>` line is followed by either a line of six numbers if
`coord` is `xyzabc` or three lines with three numbers on each line if `coord`
is `points`. If `<name>` contains an `_l` suffix the fragment parameters for
this fragment will be searched in the `fraglib_path` directory. Otherwise the
directory specified by the `userlib_path` option will be used.
