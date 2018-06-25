#include "common.h"

void print_paired_energy_coord(struct state *state);
void sim_psp(struct state *state);

void print_paired_energy_coord(struct state *state){

        size_t n_frags;
        double coord[6 * n_frags];

        check_fail(efp_get_frag_count(state->efp, &n_frags));
        check_fail(efp_get_coordinates(state->efp, coord));
	
        for (size_t i=0; i <n_frags; i++){
                char name[64];
		char ligand[64];
		size_t n_atoms;
		size_t ligand_index = cfg_get_int(state->cfg, "ligand");
		
                check_fail(efp_get_frag_name(state->efp, i, sizeof(name),name));
		check_fail(efp_get_frag_name(state->efp, ligand_index, sizeof(ligand),ligand)); 
		check_fail(efp_get_frag_atom_count(state->efp, i, & n_atoms));

		struct efp_atom atoms[n_atoms];
                check_fail(efp_get_frag_atoms(state->efp, i, n_atoms, atoms));
		
		msg("fragment %s\n", name); 

		for (size_t a = 0; a < n_atoms; a++) {
                        double x = atoms[a].x * BOHR_RADIUS;
                        double y = atoms[a].y * BOHR_RADIUS;
                        double z = atoms[a].z * BOHR_RADIUS;

                        msg("   %-16s %12.6lf %12.6lf %12.6lf\n", atoms[a].label, x, y, z);
                }


                msg("   PAIRWISE-EFP ENERGY BETWEEN FRAGMENT %zu (%s) with LIGAND %zu (%s) \n", i+1, name, ligand_index, ligand);
		msg("   %16.8s %16.8s %16.8s %16.8s %16.8s %16.8s \n", "ELEC", "POL", "DISP", "EXCH REP", "CHR-PEN", "TOTAL"); 	
                print_six_t(state->energy_components + 6 * i);
                msg("\n\n");

        }
}
void sim_psp(struct state *state)
{
        msg("PAIRWISE SINGLE POINT ENERGY JOB\n\n\n");

        print_geometry(state->efp);
        compute_energy(state, false);
        print_paired_energy_coord(state);
        print_energy(state);

        msg("PAIRWISE SINGLE POINT ENERGY JOB COMPLETED SUCCESSFULLY\n");
}
