/*
Here are is my proposal for the the structure of the mc.c file:

void sim_mc(struct state *state);
      this is the main subroutine that calls the following functions:
void mc_create(struct state *state);
       allocates and assigns coordinates from the .efp and .inp input files for a libefp calculation.  
void initialize_positions(struct state *state);
        this is a necessary step to popular necessary data structures for MC. 
void update_step_mc(struct state *state);
       this is the function that your will be modifying for randomization of the center of mass (COM) position of fragments in a for loop. 

See below:
____________
*/

#include "common.h"
#include "rand.h"

#define MAX_ITER 10

struct body {
	vec_t pos; 
	vec_t mc_pos;
};

struct mc {
	size_t n_bodies;
	struct body *bodies;
	size_t n_freedom;
	vec_t box;
	int step; /* current md step */
	double potential_energy;
	double xr_energy; /* used in multistep md */
	double *xr_gradient; /* used in multistep md */
	double (*get_invariant)(const struct mc *);
	void (*update_step)(struct mc *);
	struct state *state;
	void *data; /* nvt/npt data */
};

void sim_mc(struct state *state);

/*
struct mc{
size_t n_bodies; 
struct body *bodies; 
size_t n_freedom; 
vec_t box; 
int step; 
double potential_energy; 
struct state *state; 
}

void sim_mc(struct state *state); 
void update_step_mc(struct mc *mc); 
*/


static struct mc *mc_create(struct state *state) {
	struct mc *mc = xcalloc(1, sizeof(struct mc));

	mc->state = state;
	mc->box = box_from_str(cfg_get_string(state->cfg, "periodic_box"));

	mc->update_step = update_step_mc;

	mc->n_bodies = state->sys->n_frags;
	mc->bodies = xcalloc(mc->n_bodies, sizeof(struct body));

	double coord[6 * mc->n_bodies];
	check_fail(efp_get_coordinates(state->efp, coord));

	for (size_t i = 0; i < mc->n_bodies; i++) {
		struct body *body = mc->bodies + i;

		body->pos.x = coord[6 * i + 0];
		body->pos.y = coord[6 * i + 1];
		body->pos.z = coord[6 * i + 2];

		double a = coord[6 * i + 3];
		double b = coord[6 * i + 4];
		double c = coord[6 * i + 5];

		euler_to_matrix(a, b, c, &body->rotmat);

	}

	return (mc);
}

static void initialize_positions(struct mc *mc) {
 
	for (size_t i = 0; i < mc->n_bodies; i++) {
		double crd[12]; 

		memcpy(crd, &mc->bodies[i].pos, 3 *sizeof(double)); 
		memcpy(crd + 3, &mc->bodies[i].rotmat, 9 * sizeof(double));

		check_fail(efp_set_frag_coordinates(mc->state->efp, i, EFP_COORD_TYPE_ROTMAT, crd));

	}

	compute_energy(mc->state, false); 

	mc->potential_energy = mc->state->energy; 

}

static void update_step_mc(struct mc *mc){
 

	for (size_t i = 0; i < mc->n_bodies; i++){

		struct body *body = mc->bodies + i; 

	//figure out how to scramble the positions for each fragment from body->pos.x; and store proposed new coordinates for each fragment at mc->bodies[i].proposed_pos. 

	}

	for (size_t i = 0; i < mc->n_bodies; i++){

		double crd[12]; 

		memcpy(crd, &mc->bodies[i].proposed_pos, 3 *sizeof(double));  
		memcpy(crd + 3, &mc->bodies[i].rotmat, 9 * sizeof(double)); //i'm not really sure if we need this..but we'll figure that out in the debugging process. 

		//this sets the coordinates in efp
		check_fail(efp_set_frag_coordinates(mc->state->efp, i, EFP_COORD_TYPE_ROTMAT, crd));

	}

	compute_energy(state, false); //compute_energy needs to be modified for the next part to work -> the comparison needs to be done in energy.c; 

	if accept_move == true {

		memcpy(&mc->bodies[i].pos, &mc->bodies[i].proposed_pos, 3 *sizeof(double));
	}

}

static void mc_shutdown(struct mc *mc)
{
	free(mc->bodies);
	free(mc->xr_gradient);
	free(mc->data);
	free(mc);
}


void sim_mc(struct state *state){
 
	msg ("MONTE CARLO DYNAMICS JOB\n\n\n"); 

	struct mc *mc = mc_create(state);
	initialize_positions(mc); 

	msg(" INITIAL CONFIGURATION\n\n"); 
	print_status(mc); 

	for (mc->step = 1; 
	mc->step <= cfg_get_int(state->cfg, "max_steps");
	mc->step++) {

		mc->update_step(mc);

		if (mc->step % cfg_get_int(state->cfg, "print_step") == 0){
			msg(" STATE AFTER %d STEPS \n\n", mc->step); 
			print_status(mc);
		}
	}

	mc_shutdown(mc);

	msg("MONTE CARLO JOB COMPLETED SUCCESSFULLY\n");

}

/*
__________

for energy.c there are only a few lines we need to add and call too, i'll work on this later; look for //begin where we need to modify 

#include "common.h"
*/

// current coordinates from efp struct are used 
/*
void compute_energy(struct state *state, bool do_grad)
{
struct efp_atom *atoms;
struct efp_energy efp_energy;
double xyz[3], xyzabc[6], *grad;
size_t ifrag, nfrag, iatom, natom;
int itotal;

// EFP part
check_fail(efp_compute(state->efp, do_grad));
check_fail(efp_get_energy(state->efp, &efp_energy));
check_fail(efp_get_frag_count(state->efp, &nfrag));

if (do_grad) {
check_fail(efp_get_gradient(state->efp, state->grad));
check_fail(efp_get_point_charge_gradient(state->efp,
    state->grad + 6 * nfrag));
}

//begin where we need to modify 
if (efp_energy.total - state->energy < 0){
return state->sys->accept_move = true; 
}else{

do the random thing for probability for sigma:

if true: 
state->sys->accept_move = true; 
else:
state->sys->accept_move = false; 
break; 
}
// end where we need to modify 

state->energy = efp_energy.total;

// constraints
for (ifrag = 0; ifrag < nfrag; ifrag++) {
const struct frag *frag = state->sys->frags + ifrag;

check_fail(efp_get_frag_xyzabc(state->efp, ifrag, xyzabc));

if (frag->constraint_enable) {
double dr2, drx, dry, drz;

drx = xyzabc[0] - frag->constraint_xyz.x;
dry = xyzabc[1] - frag->constraint_xyz.y;
drz = xyzabc[2] - frag->constraint_xyz.z;

dr2 = drx * drx + dry * dry + drz * drz;
state->energy += 0.5 * frag->constraint_k * dr2;

if (do_grad) {
grad = state->grad + 6 * ifrag;
grad[0] += frag->constraint_k * drx;
grad[1] += frag->constraint_k * dry;
grad[2] += frag->constraint_k * drz;
}
}
}

// MM force field part
if (state->ff == NULL)
return;

for (ifrag = 0, itotal = 0; ifrag < nfrag; ifrag++) {
check_fail(efp_get_frag_atom_count(state->efp, ifrag, &natom));
atoms = xmalloc(natom * sizeof(struct efp_atom));
check_fail(efp_get_frag_atoms(state->efp, ifrag, natom, atoms));

for (iatom = 0; iatom < natom; iatom++, itotal++)
ff_set_atom_xyz(state->ff, itotal, &atoms[iatom].x);

free(atoms);
}

ff_compute(state->ff, do_grad);

if (do_grad) {
for (ifrag = 0, itotal = 0, grad = state->grad; ifrag < nfrag; ifrag++, grad += 6) {
check_fail(efp_get_frag_xyzabc(state->efp, ifrag, xyzabc));
check_fail(efp_get_frag_atom_count(state->efp, ifrag, &natom));
atoms = xmalloc(natom * sizeof(struct efp_atom));
check_fail(efp_get_frag_atoms(state->efp, ifrag, natom, atoms));

for (iatom = 0; iatom < natom; iatom++, itotal++) {
ff_get_atom_gradient(state->ff, itotal, xyz);

grad[0] += xyz[0];
grad[1] += xyz[1];
grad[2] += xyz[2];

grad[3] += (atoms[iatom].y - xyzabc[1]) * xyz[2] -
   (atoms[iatom].z - xyzabc[2]) * xyz[1];
grad[4] += (atoms[iatom].z - xyzabc[2]) * xyz[0] -
   (atoms[iatom].x - xyzabc[0]) * xyz[2];
grad[5] += (atoms[iatom].x - xyzabc[0]) * xyz[1] -
   (atoms[iatom].y - xyzabc[1]) * xyz[0];
}

free(atoms);
}
}

state->energy += ff_get_energy(state->ff);
}

________

in efpmc/src/common.h

struct sys {
size_t n_frags;
struct frag *frags;
size_t n_charges;
struct charge *charges; bool/or int
accept_move; 
};
struct state {
struct efp *efp;
struct ff *ff;
struct cfg *cfg;
struct sys *sys;
double energy;
double *grad;
};
*/