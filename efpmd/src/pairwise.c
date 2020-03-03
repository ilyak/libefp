/*-
 * Copyright (c) 2012-2015 Ilya Kaliman
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

#include "common.h"

void print_pair_energy(struct state *state);
void sim_pairwise(struct state *state);

void print_pair_energy(struct state *state){

        size_t n_frags;
        
        check_fail(efp_get_frag_count(state->efp, &n_frags));
        double coord[6 * n_frags];
        check_fail(efp_get_coordinates(state->efp, coord));

        struct efp_energy *energies;
        energies = xmalloc(n_frags * sizeof(struct efp_energy));
        check_fail(efp_get_pairwise_energy(state->efp, energies));
        
        char ligand[64];
        size_t lig_atoms;

        size_t ligand_index = cfg_get_int(state->cfg, "ligand");
        check_fail(efp_get_frag_name(state->efp, ligand_index, sizeof(ligand),ligand)); 
        check_fail(efp_get_frag_atom_count(state->efp, ligand_index, &lig_atoms));
        struct efp_atom latoms[lig_atoms];
        check_fail(efp_get_frag_atoms(state->efp, ligand_index, lig_atoms, latoms));
 
        char frag_name[64];
        size_t frag_atoms;               
        for (size_t i=0; i <n_frags; i++){
                check_fail(efp_get_frag_name(state->efp, i, sizeof(frag_name),frag_name));
                check_fail(efp_get_frag_atom_count(state->efp, i, &frag_atoms));

                struct efp_atom atoms[frag_atoms];
                check_fail(efp_get_frag_atoms(state->efp, i, frag_atoms, atoms));
                
                msg("   PAIRWISE ENERGY BETWEEN FRAGMENT %zu (%s) AND LIGAND %zu (%s) \n", i, frag_name, ligand_index, ligand);
                msg("fragment %s\n", name);
                for (size_t a = 0; a < n_atoms; a++) {
                        double x = atoms[a].x * BOHR_RADIUS;
                        double y = atoms[a].y * BOHR_RADIUS;
                        double z = atoms[a].z * BOHR_RADIUS;
                        msg("   %-16s %12.6lf %12.6lf %12.6lf\n", atoms[a].label, x, y, z);
                }
                msg("\n");

                msg("fragment %s\n", ligand);        
                for (size_t a = 0; a < lig_atoms; a++) {
                        double x = latoms[a].x * BOHR_RADIUS;
                        double y = latoms[a].y * BOHR_RADIUS;
                        double z = latoms[a].z * BOHR_RADIUS;
                        msg("   %-16s %12.6lf %12.6lf %12.6lf\n", latoms[a].label, x, y, z);
                }
                msg("\n\n");

                msg("    ENERGY COMPONENTS (ATOMIC UNITS)\n\n");
                msg("%30s %16.10lf\n", "ELECTROSTATIC ENERGY", energies.electrostatic);
                msg("%30s %16.10lf\n", "POLARIZATION ENERGY", energies.polarization);
                msg("%30s %16.10lf\n", "DISPERSION ENERGY", energies.dispersion);
                msg("%30s %16.10lf\n", "EXCHANGE REPULSION ENERGY",
                        energies.exchange_repulsion);
                msg("%30s %16.10lf\n", "CHARGE PENETRATION ENERGY",
                        energies.charge_penetration);
                msg("\n");

                msg("%30s %16.10lf\n", "TOTAL ENERGY", energies.total);
                msg("\n\n");

        }
}

void sim_pairwise(struct state *state)
{
        msg("PAIRWISE SINGLE POINT ENERGY JOB\n\n\n");

        print_geometry(state->efp);
        compute_energy(state, false);
        print_pair_energy(state);
        print_energy(state);

        msg("PAIRWISE SINGLE POINT ENERGY JOB COMPLETED SUCCESSFULLY\n");
}

