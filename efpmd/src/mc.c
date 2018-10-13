#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "mc.h"
#include "rand.h"


#define DISPMAG_THRESHOLD 0.5
#define DISPMAG_MODIFIER 0.95
#define DISPMAG_MODIFY_STEPS 100

void sim_mc(struct state *state); 

struct mc_state {
	size_t n;
	size_t n_frags;
	double *x;
	double f;
	char task[60];
	mc_func_t func;
	void *data; // Why is this not a state*? It is cast to one later
	double *x_prop;
	double f_prop; 
	double n_charge; 
	int step;
	int n_accept; 
	int n_reject; 
	double dispmag;
	double dispmag_threshold; 
	double dispmag_modifier; 
	int dispmag_modify_steps; 
	double anglemag; 
	double temperature; 
};

struct mc_state *mc_create(size_t n)
{
	struct mc_state *state = calloc(1, sizeof(struct mc_state));
	assert(state);

	state->n = n;

	state->x = calloc(n, sizeof(double));
	assert(state->x);
	state->x_prop= calloc(n, sizeof(double));
	assert(state->x_prop);

	return state;
}

enum mc_result mc_init(struct mc_state *state, size_t n, const double *x)
{
//	fprintf(stdout, "\n----mc init----\n\n");
	assert(state);
	assert(n == state->n);
	assert(x);

	memcpy(state->x, x, n * sizeof(double));
	memset(state->task, '\0', sizeof(state->task));

//	state->dispmag = cfg_get_double(state->cfg, "dispmag_threshold"); 
	state->n_accept = 0; 
	state->n_reject = 0;
	state->step = 1;
	
	const char start[] = "START";
	strncpy(state->task, start, sizeof(start));

	if (strncmp(state->task, start, sizeof(start)) != 0) {
//		fprintf(stdout, "\n----failed strncmp %s----\n\n", state->task);
//return MC_RESULT_ERROR;
	}

//	fprintf(stdout, "Setting state->f\n");
	state->f = state->func(state->n, state->x, state->data);

	if (isnan(state->f)) {
	//	fprintf(stdout, "\n----failed isnan %lf----\n\n", state->f);
		return MC_RESULT_ERROR;
	}

	return MC_RESULT_SUCCESS;
}

void mc_set_func(struct mc_state *state, mc_func_t func)
{
	assert(state);
	state->func = func;
}

void mc_set_user_data(struct mc_state *state, void *data)
{
	assert(state);
	state->data = data;
}

double mc_get_fx(struct mc_state *state)
{
	assert(state);
	return state->f;
}

static void print_restart(struct efp *efp)
{
	size_t n_frags;
	check_fail(efp_get_frag_count(efp, &n_frags));

	double coord[6 * n_frags];
	check_fail(efp_get_coordinates(efp, coord));

	msg("    RESTART DATA\n\n");

	for (size_t i = 0; i < n_frags; i++) {
		char name[64];
		check_fail(efp_get_frag_name(efp, i, sizeof(name), name));

		coord[6 * i + 0] *= BOHR_RADIUS;
		coord[6 * i + 1] *= BOHR_RADIUS;
		coord[6 * i + 2] *= BOHR_RADIUS;

		print_fragment(name, coord + 6 * i, NULL);
	}

	size_t n_charges;
	check_fail(efp_get_point_charge_count(efp, &n_charges));

	if (n_charges > 0) {
		double q[n_charges];
		check_fail(efp_get_point_charge_values(efp, q));

		double xyz[3 * n_charges];
		check_fail(efp_get_point_charge_coordinates(efp, xyz));

		for (size_t i = 0; i < n_charges; i++) {
			double x = xyz[3 * i + 0] * BOHR_RADIUS;
			double y = xyz[3 * i + 1] * BOHR_RADIUS;
			double z = xyz[3 * i + 2] * BOHR_RADIUS;

			print_charge(q[i], x, y, z);
		}
	}

	msg("\n");
}

static void print_status(struct state *state, double e_diff)
{
        print_geometry(state->efp);
        print_restart(state->efp);
        print_energy(state);

        msg("%30s %16.10lf\n", "ENERGY CHANGE", e_diff);
        msg("\n\n");

        fflush(stdout);
}

static void mc_rand(struct mc_state *state){

	//#copy coord * gradient to n, x_prop, g_prop, state->data; 
	printf("n_frags in system: %li\n", state->n_frags);
	for (size_t i = 0; i < state->n_frags; i++){

		printf("dispmag %f", state->dispmag); 
		rand(); 
		printf("goes through mc_rand before\n");
		// Are all of these valid indices to state->x_prop?
		state->x_prop[6 * i + 0] = state->x[6 *i + 0] + (state->dispmag * rand_neg1_to_1()); 
		state->x_prop[6 * i + 1] = state->x[6 *i + 1] + (state->dispmag * rand_neg1_to_1()); 
		state->x_prop[6 * i + 2] = state->x[6 *i + 2] + (state->dispmag * rand_neg1_to_1()); 

		state->x_prop[6 * i + 3] = state->x[6 *i + 3] + (state->dispmag * rand_neg1_to_1());
	        state->x_prop[6 * i + 4] = state->x[6 *i + 4] + (state->dispmag * rand_neg1_to_1());
	        state->x_prop[6 * i + 5] = state->x[6 *i + 5] + (state->dispmag * rand_neg1_to_1());
		printf("com of frag %li x: %f y: %f z: %f \n", i, state->x_prop[6 * i + 0], state->x_prop[6 * i + 1], state->x_prop[6 * i + 2]); 
	}

	//#mess with translation

	printf("finishes mc_rand\n"); 
}

static bool check_acceptance_ratio(struct mc_state *state, double new_energy, double old_energy){
//	#Insert Jay's acceptance ratio script thing here;

	bool allow_move = false; 
//	fprintf(stdout, "Checking, new=%lf, old=%lf\n", new_energy, old_energy);

	if (new_energy < old_energy){
		allow_move = true; 
		state->f = new_energy; 
	}

	else {
	
		double temperature;
		temperature = state->temperature;  
		
		rand(); 
		double epsilon=	rand_normal(); 
//		double epsilon = rand_uniform_1();
		double delta_e = new_energy - old_energy; 
		double exp_value = -delta_e / (BOLTZMANN * temperature); 
		double bf = exp(exp_value);

//		fprintf(stdout, "Comparing, delta_e=%lf, exp_value=%lf, bf=%lf\n",
//			delta_e, exp_value, bf);

		if (epsilon < bf ) {
			allow_move = true;
			state->f = new_energy; 
		}
	}
	printf("check ratio\n"); 
	printf("new_energy %f\n", new_energy); 
	printf("old_energy %f\n", old_energy); 	
	return allow_move; 

}

static bool mc_step(struct mc_state *state){
	assert(state);

	bool allow_step = false;

	mc_rand(state);

	// state->task should be an integer, not a string (would be much easier 
	// to work with)
	
	// Nothing that I can see ever changes the value of state->task,
	// so both if statements will always fail.
	// Presumably there should be code somewhere that sets state-task
	// to FG or NEW_X.
	fprintf(stdout, "Setting state->f_prop\n");
	state->f_prop = state->func(state->n, state->x_prop, state->data);
	fprintf(stdout, "Set state->f_prop %lf, state->f is %lf\n",
		state->f_prop, state->f);

		// Maybe this code should always happen with no if around it?

	allow_step = check_acceptance_ratio(state, state->f_prop, state->f);

	return allow_step;
}

static double compute_efp(size_t n, const double *x, void *data)
{
	size_t n_frags, n_charge;
	//const double *coord = x;	
	// Why pass in void* and cast to state*? Far safer to have a state* parameter.
	struct state *state = (struct state *)data;
	struct efp_energy  energy; 
	// Why does this duplicate part of compute_energy,
	// i.e. efp_get_frag_count, then call compute_energy?
	check_fail(efp_get_frag_count(state->efp, &n_frags));

	/* Debugging */
//	fprintf(stdout, "compute_efp, n=%lu\n", n);
//	check_fail(efp_set_coordinates(struct efp *efp, enum efp_coord_type coord_type, const double *coord))
//	print_energy(state);
//	for (size_t i = 0; i < n_frags; i++, coord += 6) {
//		fprintf(stdout, "compute_efp, i=%lu, x[i]=%e\n", i, *coord);
//	}
	/* End debugging */

	check_fail(efp_get_point_charge_count(state->efp, &n_charge));

	assert(n == (6 * n_frags + 3 * n_charge));

	check_fail(efp_set_coordinates(state->efp, EFP_COORD_TYPE_XYZABC, x));
	check_fail(efp_set_point_charge_coordinates(state->efp, x + 6 * n_frags));

	compute_energy(state, false);
	check_fail(efp_get_energy(state->efp, &energy));
	/*
	efp_get_energy
	*/

	//fprintf(stdout, 'Energies 1 %f\n', energy.electrostatic);
	//fprintf(stdout, 'Energies 2 %lf\n', energy.polarization);
	//fprintf(stdout, 'Energies 3 %lf\n', energy.dispersion);
	//fprintf(stdout, 'Energies 4 %lf\n', energy.exchange_repulsion);
	state->energy = energy.electrostatic + energy.polarization + energy.dispersion + energy.exchange_repulsion;  
	
	//#we randomized gradients before so don't worry; 
	//#memset(gx, state->grad, (6 * n_frags + 3 * n_charge) * sizeof(double));
	return(state->energy);
}

void mc_shutdown(struct mc_state *state)
{
	if (state) {
		free(state->x);
		free(state->x_prop); 
		free(state);
	}
}

void sim_mc(struct state *state){
	
	msg ("MONTE CARLO JOB \n"); 

	size_t n_frags, n_charge, n_coord; 
	//double rms_grad, max_grad; 

	check_fail(efp_get_frag_count(state->efp, & n_frags));
	check_fail(efp_get_point_charge_count(state->efp, &n_charge));

	n_coord = 6 * n_frags + 3 * n_charge;// (xyzabc + charge coordinages = arrazy of size 9)

	struct mc_state *mc_state = mc_create(n_coord);

	mc_state->n_charge = n_charge; 
	mc_state->dispmag = cfg_get_double(state->cfg, "dispmag_modifier"); 

	if (!mc_state)
		error("unable to create an move!");

//	#sets compute_efp as opt_state->compute_efp
	mc_set_func(mc_state, compute_efp); 

//	#sets opt_state->data = state; 
	mc_set_user_data(mc_state, state);


//	#allocaates array of n_coord size and grad of n_coord]; 
	double coord[n_coord]; //, grad[n_coord];

//	#copies state->efp, into coord for each fragment; 
	check_fail(efp_get_coordinates(state->efp, coord));
	check_fail(efp_get_point_charge_coordinates(state->efp, coord + 6 * n_frags));

        mc_state->temperature = cfg_get_double(state->cfg, "temperature");
        mc_state->dispmag_threshold = cfg_get_double(state->cfg, "dispmag_threshold");
        mc_state->dispmag_modifier = cfg_get_double(state->cfg, "dispmag_modifier");
        mc_state->dispmag_modify_steps = cfg_get_int(state->cfg, "dispmag_modify_steps");

//	#does an efp_compute with coord
	if (mc_init(mc_state, n_coord, coord)) {
		error("unable to initialize an optimizer");
	}

	mc_state->n_frags = n_frags;

	double e_old = mc_get_fx(mc_state);
	printf("e_old %f\n", e_old); 
//	#sets up empty array; 
//	memset(grad,0, n_coord*sizeof(double)); 

	msg("    INITIAL STATE\n\n");
	print_status(state, 0.0);
	double e_new;
	bool allow_step;
	double acceptance_ratio;
	for (int step = 1; step <= cfg_get_int(state->cfg, "max_steps"); step++) {
		msg("   STATE AFTER %d STEPS\n\n", step);
		allow_step = mc_step(mc_state);
		fprintf(stdout, "Allow step: %d\n", allow_step);
		if (allow_step) {
			// Can only memcpy up to the size that was allocated
			memcpy(mc_state->x, mc_state->x_prop, mc_state->n * sizeof(double)); 
			mc_state->f = mc_state->f_prop;

			mc_state->n_accept++;

			e_new = mc_get_fx(mc_state);
			print_status(state, e_new-e_old);

		// if (state->step % DISPMAG_MODIFY_STEPS == 0) {
		// 	if (state->dispmag <= DISPMAG_THRESHOLD) {
		// 		fprintf(stdout, "dividing dispmag\n");
		// 		state->dispmag = state->dispmag / DISPMAG_MODIFIER;
		// 	} else {
		// 		fprintf(stdout, "multiplying dispmag\n");
		// 		state->dispmag = state->dispmag * DISPMAG_MODIFIER;	
		// 	}
		// 	msg("    NEW DISPMAG %e\n\n", state->dispmag);
		// }

		}else{
			mc_state->n_reject++; 
		}

		if ((step % mc_state->dispmag_modify_steps) == 0) {
			acceptance_ratio = (double)mc_state->n_accept/step;
			fprintf(stdout, "Step: %d, acceptance ratio: %lf\n", step, acceptance_ratio);
//			if ((acceptance_ratio < (mc_state->dispmag_threshold * 1.2))) {
			if (acceptance_ratio < 0.5){
				mc_state->dispmag = (mc_state->dispmag * mc_state->dispmag_modifier);
			printf("goes through increasing dispmag"); 
			} 
			if (acceptance_ratio > 0.5){
//			if ((acceptance_ratio > (mc_state->dispmag_threshold*0.8))) {
				mc_state->dispmag = (mc_state->dispmag / mc_state->dispmag_modifier);
			                        printf("goes through decreasing dispmag");
			}
			msg("    NEW DISPMAG %e\n\n", mc_state->dispmag);
		}

		// Step always errors:
		// step result is: 1, error is: 1

	}

	msg("	FINAL STATE\n\n");
	msg("	TOTAL STEP ACCEPTED: %d\n", mc_state->n_accept); 
	msg("	TOTAL STEPS REJECTED: %d\n", mc_state->n_reject);  
	print_status(state, e_new); 
	mc_shutdown(mc_state);

	msg("	MONTE CARLO JOB DONE\n"); 

}
