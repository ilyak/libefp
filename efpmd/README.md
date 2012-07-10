# efpmd

## Input file format

`run_type [sp|grad|cg|nve|nvt]`

Specifies the type of the simulation.

`sp` - single point calculation.

`grad` - gradient calculation.

`cg` - conjugate gradient geometry optimization.

`nve` - molecular dynamics in microcanonical ensemble.

`nvt` - molecular dynamics in canonical ensemble.


`coord [xyzabc|points]`

Specifies the format of fragment input.

`xyzabc` - Coordinates of the center of mass and Euler angles.

`points` - Coordinates of three atoms for each fragment.


`units [angs|bohr]`

Specifies which units to use.

`angs` - Units are Angstroms.

`bohr` - Units are Bohr.


`terms [elec [pol [disp [xr]]]]`

Specifies which energy terms to include in EFP computation.

`elec` - Include electrostatics energy.

`pol` - Include polarization energy.

`disp` - Include dispersion energy.

`xr` - Include exchange repulsion energy.


`elec_damp [screen|overlap|off]`

Electrostatic damping type.

`screen` - Damping formula based on SCREEN group in the EFP potential.

`overlap` - Overlap based damping which computes charge penetration energy.

`off` - No electrostatic damping.


`disp_damp [tt|overlap|off]`

Dispersion damping type.

`tt` - Damping based on the formula by Tang and Toennies.

`overlap` - Overlap-based dispersion damping.

`off` - No dispersion damping.


`fraglib_path <path>`

Specifies the path to the directory with fragment library.


`userlib_path <path>`

Specifies the path to the directory with user-created fragments.


`fragment <name>`

Specifies a fragment.
