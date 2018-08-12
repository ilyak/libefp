#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "mc.h"
#include "rand.h"

void sim_mc(struct state *state); 

struct mc_state {
	size_t n;
	size_t n_frags;
	double *x;
	double f;
	double *g;
	char task[60];
	mc_func_t func;
	void *data;
	double *x_prop;
	double *g_prop; 
	double f_prop; 
	double n_charge; 
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
	state->g = calloc(n, sizeof(double));
	assert(state->g);
	state->g_prop = calloc(n, sizeof(double));
	assert(state->g_prop);

	return state;
}

enum mc_result mc_init(struct mc_state *state, size_t n, const double *x)
{
	fprintf(stdout, "\n----mc init----\n\n");
	assert(state);
	assert(n == state->n);
	assert(x);

	memcpy(state->x, x, n * sizeof(double));
	memset(state->task, '\0', sizeof(state->task));

	const char start[] = "START";
	strncpy(state->task, start, sizeof(start));

	if (strncmp(state->task, start, sizeof(start)) != 0) {
		fprintf(stdout, "\n----failed strncmp %s----\n\n", state->task);
		return MC_RESULT_ERROR;
	}

	state->f = state->func(state->n, state->x, state->data);

	if (isnan(state->f)) {
		fprintf(stdout, "\n----failed isnan %lf----\n\n", state->f);
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

	for (size_t i = 0; i < state->n_frags; i++){

		state->x_prop[6 * i + 0] = state->x[6 *i + 0] + rand_neg1_to_1(); 
		state->x_prop[6 * i + 1] = state->x[6 *i + 1] + rand_neg1_to_1(); 
		state->x_prop[6 * i + 2] = state->x[6 *i + 2] + rand_neg1_to_1(); 

	}

	//#mess with translation

	for (size_t i = 0; i < state->n_frags; i++){

		state->x_prop[6 * i + 3] = state->x[6 *i + 3] + rand_neg1_to_1(); 
		state->x_prop[6 * i + 4] = state->x[6 *i + 4] + rand_neg1_to_1(); 
		state->x_prop[6 * i + 5] = state->x[6 *i + 5] + rand_neg1_to_1(); 

	}
}

static void check_acceptance_ratio(struct mc_state *state, double new_energy, double old_energy){
//	#Insert Jay's acceptance ratio script thing here;

	bool allow_move = false; 

	if(new_energy < old_energy){
		allow_move = true; 
	}

	else{
	
		double temperature = 298.0; 
		double epsilon = rand_uniform_1();
		double delta_e = new_energy - old_energy; 
		double exp_value = -delta_e / (BOLTZMANN * temperature); 
		double bf = exp(exp_value);

		if(epsilon < bf ){
			allow_move = true;
		}
	}

	if (allow_move){
		memcpy(state->x, state->x_prop, (6 * state->n + 3 * state->n_charge) * sizeof(double)); 
	}; 

}

enum mc_result mc_step(struct mc_state *state){
	assert(state);

next:
	mc_rand(state);
	
	if (strncmp(state->task, "FG", strlen("FG")) == 0) {
		state->f_prop = state->func(state->n, state->x_prop, state->data);

		check_acceptance_ratio(state, state->f, state->f_prop); 

		if (isnan(state->f))
			return MC_RESULT_ERROR;

		goto next;
	}

	if (strncmp(state->task, "NEW_X", strlen("NEW_X")) == 0)
		return MC_RESULT_SUCCESS;

	return MC_RESULT_ERROR;
}

static double compute_efp(size_t n, const double *x, void *data)
{
	size_t n_frags, n_charge;
	struct state *state = (struct state *)data;

	check_fail(efp_get_frag_count(state->efp, &n_frags));
	check_fail(efp_get_point_charge_count(state->efp, &n_charge));

	assert(n == (6 * n_frags + 3 * n_charge));

	check_fail(efp_set_coordinates(state->efp, EFP_COORD_TYPE_XYZABC, x));
	check_fail(efp_set_point_charge_coordinates(state->efp, x + 6 * n_frags));

	compute_energy(state, false);
	
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
	
	msg ("YO MONTE CARLO JOB BITCHES"); 

	size_t n_frags, n_charge, n_coord; 
	double rms_grad, max_grad; 

	check_fail(efp_get_frag_count(state->efp, & n_frags));
	check_fail(efp_get_point_charge_count(state->efp, &n_charge));

	n_coord = 6 * n_frags + 3 * n_charge;// (xyzabc + charge coordinages = arrazy of size 9)

	struct mc_state *mc_state = mc_create(n_coord);

	mc_state->n_charge = n_charge; 

	if (!mc_state)
		error("unable to create an move!");

//	#sets compute_efp as opt_state->compute_efp
	mc_set_func(mc_state, compute_efp); 

//	#sets opt_state->data = state; 
	mc_set_user_data(mc_state, state);


//	#allocaates array of n_coord size and grad of n_coord]; 
	double coord[n_coord], grad[n_coord];

//	#copies state->efp, into coord for each fragment; 
	check_fail(efp_get_coordinates(state->efp, coord));
	check_fail(efp_get_point_charge_coordinates(state->efp, coord + 6 * n_frags));


//	#does an efp_compute with coord
	if (mc_init(mc_state, n_coord, coord)) {
		error("unable to initialize an optimizer");
	}

	double e_old = mc_get_fx(mc_state);

//	#sets up empty array; 
//	memset(grad,0, n_coord*sizeof(double)); 

	msg("    INITIAL STATE\n\n");
	print_status(state, 0.0);
	double e_new; 
	for (int step = 1; step <= cfg_get_int(state->cfg, "max_steps"); step++) {
		if (mc_step(mc_state)){
			e_new = mc_get_fx(mc_state);
			msg("	STATE AFTER %d STEPS\n\n", step);
			print_status(state, e_new-e_old);
		}
	}
	msg("	FINAL STATE\n\n");
	print_status(state, e_new); 

	mc_shutdown(mc_state);

	msg("	MONTE CARLO JOB DONE BITCHES\n"); 

}
