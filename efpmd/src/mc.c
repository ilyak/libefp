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
#include "rand.h"

#define MAX_ITER 10

struct body {
	mat_t rotmat;
	vec_t pos;
	vec_t vel;
	vec_t vel_old;
	vec_t angmom;
	vec_t angmom_old;
	vec_t force;
	vec_t torque;
	vec_t inertia;
	vec_t inertia_inv;
	double mass;
};

struct mc {
	size_t n_bodies;
	struct body *bodies;
	size_t n_freedom;
	vec_t box;
	int step; /* current mc step */
	double potential_energy;
	double mc_energy; /* energy calculated using monte carlo method */
	double xr_energy; /* used in multistep mc */
	double *xr_gradient; /* used in multistep mc */
	double (*get_invariant)(const struct mc *);
	void (*update_step)(struct mc *);
	struct state *state;
	void *data; /* nvt/npt data */
};

void sim_mc(struct state *state);

static vec_t wrap(const struct mc *mc, const vec_t *pos)
{
	if (!cfg_get_bool(mc->state->cfg, "enable_pbc"))
		return *pos;

	vec_t sub = {
		mc->box.x * floor(pos->x / mc->box.x),
		mc->box.y * floor(pos->y / mc->box.y),
		mc->box.z * floor(pos->z / mc->box.z)
	};

	return vec_sub(pos, &sub);
}

static double get_kinetic_energy(const struct mc *mc)
{
	double ke = 0.0;

	for (size_t i = 0; i < mc->n_bodies; i++) {
		struct body *body = mc->bodies + i;

		ke += body->mass * body->vel.x * body->vel.x;
		ke += body->mass * body->vel.y * body->vel.y;
		ke += body->mass * body->vel.z * body->vel.z;

		ke += body->angmom.x * body->angmom.x * body->inertia_inv.x;
		ke += body->angmom.y * body->angmom.y * body->inertia_inv.y;
		ke += body->angmom.z * body->angmom.z * body->inertia_inv.z;
	}

	return 0.5 * ke;
}

static double get_temperature(const struct mc *mc)
{
	double ke = get_kinetic_energy(mc);

	return 2.0 * ke / BOLTZMANN / mc->n_freedom;
}

static double get_volume(const struct mc *mc)
{
	return mc->box.x * mc->box.y * mc->box.z;
}

static double get_pressure(const struct mc *mc)
{
	double volume = get_volume(mc);
	vec_t pressure = vec_zero;

	for (size_t i = 0; i < mc->n_bodies; i++) {
		const struct body *body = mc->bodies + i;

		pressure.x += body->mass * body->vel.x * body->vel.x;
		pressure.y += body->mass * body->vel.y * body->vel.y;
		pressure.z += body->mass * body->vel.z * body->vel.z;
	}

	mat_t stress;
	check_fail(efp_get_stress_tensor(mc->state->efp, (double *)&stress));

	pressure.x = (pressure.x + stress.xx) / volume;
	pressure.y = (pressure.y + stress.yy) / volume;
	pressure.z = (pressure.z + stress.zz) / volume;

	return (pressure.x + pressure.y + pressure.z) / 3.0;
}

static double get_invariant_nve(const struct mc *mc)
{
	return mc->potential_energy + get_kinetic_energy(mc);
}

static vec_t get_system_com(const struct mc *mc)
{
	double mass = 0.0;
	vec_t com = vec_zero;

	for (size_t i = 0; i < mc->n_bodies; i++) {
		struct body *body = mc->bodies + i;
		vec_t pos = wrap(mc, &body->pos);

		com.x += body->mass * pos.x;
		com.y += body->mass * pos.y;
		com.z += body->mass * pos.z;

		mass += body->mass;
	}

	com.x /= mass;
	com.y /= mass;
	com.z /= mass;

	return com;
}

static vec_t get_system_com_velocity(const struct mc *mc)
{
	double mass = 0.0;
	vec_t cv = vec_zero;

	for (size_t i = 0; i < mc->n_bodies; i++) {
		struct body *body = mc->bodies + i;

		cv.x += body->vel.x * body->mass;
		cv.y += body->vel.y * body->mass;
		cv.z += body->vel.z * body->mass;

		mass += body->mass;
	}

	cv.x /= mass;
	cv.y /= mass;
	cv.z /= mass;

	return cv;
}

static vec_t get_system_angular_momentum(const struct mc *mc)
{
	vec_t cp = get_system_com(mc);
	vec_t cv = get_system_com_velocity(mc);

	vec_t am = vec_zero;

	for (size_t i = 0; i < mc->n_bodies; i++) {
		struct body *body = mc->bodies + i;

		vec_t pos = wrap(mc, &body->pos);
		vec_t dr = vec_sub(&pos, &cp);
		vec_t dv = vec_sub(&body->vel, &cv);

		am.x += (dr.y * dv.z - dr.z * dv.y) * body->mass;
		am.y += (dr.z * dv.x - dr.x * dv.z) * body->mass;
		am.z += (dr.x * dv.y - dr.y * dv.x) * body->mass;
	}

	return am;
}

static mat_t get_system_inertia_tensor(const struct mc *mc)
{
	mat_t inertia = mat_zero;
	vec_t com = get_system_com(mc);

	for (size_t i = 0; i < mc->n_bodies; i++) {
		struct body *body = mc->bodies + i;

		vec_t pos = wrap(mc, &body->pos);
		vec_t dr = vec_sub(&pos, &com);

		inertia.xx += body->mass * (dr.y * dr.y + dr.z * dr.z);
		inertia.yy += body->mass * (dr.x * dr.x + dr.z * dr.z);
		inertia.zz += body->mass * (dr.x * dr.x + dr.y * dr.y);
		inertia.xy -= body->mass * dr.x * dr.y;
		inertia.xz -= body->mass * dr.x * dr.z;
		inertia.yz -= body->mass * dr.y * dr.z;
	}

	inertia.yx = inertia.xy;
	inertia.zx = inertia.xz;
	inertia.zy = inertia.yz;

	return inertia;
}

static void remove_system_drift(struct mc *mc)
{
	vec_t cp = get_system_com(mc);
	vec_t cv = get_system_com_velocity(mc);
	vec_t am = get_system_angular_momentum(mc);

	mat_t inertia = get_system_inertia_tensor(mc);
	mat_t inertia_inv = mat_zero;

	if (inertia.xx < EPSILON ||
	    inertia.yy < EPSILON ||
	    inertia.zz < EPSILON) {
		inertia_inv.xx = inertia.xx < EPSILON ? 0.0 : 1.0 / inertia.xx;
		inertia_inv.yy = inertia.yy < EPSILON ? 0.0 : 1.0 / inertia.yy;
		inertia_inv.zz = inertia.zz < EPSILON ? 0.0 : 1.0 / inertia.zz;
	}
	else {
		inertia_inv = mat_inv(&inertia);
	}

	vec_t av = mat_vec(&inertia_inv, &am);

	for (size_t i = 0; i < mc->n_bodies; i++) {
		struct body *body = mc->bodies + i;
		vec_t pos = wrap(mc, &body->pos);

		vec_t cross = {
			av.y * (pos.z - cp.z) - av.z * (pos.y - cp.y),
			av.z * (pos.x - cp.x) - av.x * (pos.z - cp.z),
			av.x * (pos.y - cp.y) - av.y * (pos.x - cp.x)
		};

		body->vel.x -= cv.x + cross.x;
		body->vel.y -= cv.y + cross.y;
		body->vel.z -= cv.z + cross.z;
	}

	vec_t cv2 = get_system_com_velocity(mc);
	vec_t am2 = get_system_angular_momentum(mc);

	assert(vec_len(&cv2) < EPSILON && vec_len(&am2) < EPSILON);
}

static void compute_forces(struct mc *mc)
{
	for (size_t i = 0; i < mc->n_bodies; i++) {
		double crd[12];

		memcpy(crd, &mc->bodies[i].pos, 3 * sizeof(double));
		memcpy(crd + 3, &mc->bodies[i].rotmat, 9 * sizeof(double));

		check_fail(efp_set_frag_coordinates(mc->state->efp, i,
		    EFP_COORD_TYPE_ROTMAT, crd));
	}

	if (cfg_get_bool(mc->state->cfg, "enable_multistep")) {
		struct efp_opts opts, opts_save;
		int multistep_steps;

		multistep_steps = cfg_get_int(mc->state->cfg,
		    "multistep_steps");
		check_fail(efp_get_opts(mc->state->efp, &opts));
		if (mc->step % multistep_steps == 0) {
			opts_save = opts;
			opts.terms = EFP_TERM_XR; /* xr only */
			check_fail(efp_set_opts(mc->state->efp, &opts));
			compute_energy(mc->state, true);
			mc->xr_energy = mc->state->energy;
			memcpy(mc->xr_gradient, mc->state->grad,
			    6 * mc->n_bodies * sizeof(double));
			opts = opts_save;
		}
		opts.terms &= ~EFP_TERM_XR; /* turn off xr */
		check_fail(efp_set_opts(mc->state->efp, &opts));
		compute_energy(mc->state, true);
		mc->state->energy += mc->xr_energy;
		for (size_t i = 0; i < 6 * mc->n_bodies; i++)
			mc->state->grad[i] += mc->xr_gradient[i];
	} else
		compute_energy(mc->state, true);

	mc->potential_energy = mc->state->energy;

	for (size_t i = 0; i < mc->n_bodies; i++) {
		struct body *body = mc->bodies + i;

		body->force.x = -mc->state->grad[6 * i + 0];
		body->force.y = -mc->state->grad[6 * i + 1];
		body->force.z = -mc->state->grad[6 * i + 2];

		body->torque.x = -mc->state->grad[6 * i + 3];
		body->torque.y = -mc->state->grad[6 * i + 4];
		body->torque.z = -mc->state->grad[6 * i + 5];

		/* convert torque to body frame */
		body->torque = mat_trans_vec(&body->rotmat, &body->torque);
	}
}

static void set_body_mass_and_inertia(struct efp *efp, size_t idx,
    struct body *body)
{
	double mass, inertia[3];

	check_fail(efp_get_frag_mass(efp, idx, &mass));
	check_fail(efp_get_frag_inertia(efp, idx, inertia));

	body->mass = AMU_TO_AU * mass;

	body->inertia.x = AMU_TO_AU * inertia[0];
	body->inertia.y = AMU_TO_AU * inertia[1];
	body->inertia.z = AMU_TO_AU * inertia[2];

	body->inertia_inv.x = body->inertia.x < EPSILON ? 0.0 :
	    1.0 / body->inertia.x;
	body->inertia_inv.y = body->inertia.y < EPSILON ? 0.0 :
	    1.0 / body->inertia.y;
	body->inertia_inv.z = body->inertia.z < EPSILON ? 0.0 :
	    1.0 / body->inertia.z;
}

static void rotate_step(size_t a1, size_t a2, double angle, vec_t *angmom,
    mat_t *rotmat)
{
	mat_t rot = { 1.0, 0.0, 0.0,
		      0.0, 1.0, 0.0,
		      0.0, 0.0, 1.0 };

	double cosa = cos(angle);
	double sina = sin(angle);

	mat_set(&rot, a1, a1,  cosa);
	mat_set(&rot, a2, a2,  cosa);
	mat_set(&rot, a1, a2,  sina);
	mat_set(&rot, a2, a1, -sina);

	*angmom = mat_vec(&rot, angmom);

	/* transpose */
	mat_set(&rot, a1, a2, -sina);
	mat_set(&rot, a2, a1,  sina);

	*rotmat = mat_mat(rotmat, &rot);
}

/*
 * Rotation algorithm reference:
 *
 * Andreas Dullweber, Benedict Leimkuhler, and Robert McLachlan
 *
 * Symplectic splitting methods for rigid body molecular dynamics
 *
 * J. Chem. Phys. 107, 5840 (1997)
 */
static void rotate_body(struct body *body, double dt)
{
	double angle;

	/* rotate about x axis */
	angle = 0.5 * dt * body->angmom.x * body->inertia_inv.x;
	rotate_step(1, 2, angle, &body->angmom, &body->rotmat);

	/* rotate about y axis */
	angle = 0.5 * dt * body->angmom.y * body->inertia_inv.y;
	rotate_step(2, 0, angle, &body->angmom, &body->rotmat);

	/* rotate about z axis */
	angle = dt * body->angmom.z * body->inertia_inv.z;
	rotate_step(0, 1, angle, &body->angmom, &body->rotmat);

	/* rotate about y axis */
	angle = 0.5 * dt * body->angmom.y * body->inertia_inv.y;
	rotate_step(2, 0, angle, &body->angmom, &body->rotmat);

	/* rotate about x axis */
	angle = 0.5 * dt * body->angmom.x * body->inertia_inv.x;
	rotate_step(1, 2, angle, &body->angmom, &body->rotmat);
}

static void update_step_mc(struct mc *mc)
{
	/* TODO compute mc energy change at each step */
	double dt = cfg_get_double(mc->state->cfg, "time_step");

	for (size_t i = 0; i < mc->n_bodies; i++) {
		struct body *body = mc->bodies + i;

		body->vel.x += 0.5 * body->force.x * dt / body->mass;
		body->vel.y += 0.5 * body->force.y * dt / body->mass;
		body->vel.z += 0.5 * body->force.z * dt / body->mass;

		body->angmom.x += 0.5 * body->torque.x * dt;
		body->angmom.y += 0.5 * body->torque.y * dt;
		body->angmom.z += 0.5 * body->torque.z * dt;

		body->pos.x += body->vel.x * dt;
		body->pos.y += body->vel.y * dt;
		body->pos.z += body->vel.z * dt;

		rotate_body(body, dt);
	}

	compute_forces(mc);

	for (size_t i = 0; i < mc->n_bodies; i++) {
		struct body *body = mc->bodies + i;

		body->vel.x += 0.5 * body->force.x * dt / body->mass;
		body->vel.y += 0.5 * body->force.y * dt / body->mass;
		body->vel.z += 0.5 * body->force.z * dt / body->mass;

		body->angmom.x += 0.5 * body->torque.x * dt;
		body->angmom.y += 0.5 * body->torque.y * dt;
		body->angmom.z += 0.5 * body->torque.z * dt;
	}
}

/*
Thus we
define dAB, the minimum distance between particles A
and B, as the shortest distance between A and any of
the particles B, of which there is one in each of the
squares which comprise the complete substance.

If we know the positions of the N particles in the
square, we can easily calculate, for example, the potential
energy of the system, 
 	0.5 * sum i = 1 to N sum j = 1 to N (i != j) V(d sub ij)
(Here V is the potential between molecules, and dij is
the minimum distance between particles i and j as
defined above.) 

This we do as follows: We place the N particles in any
configuration, for example, in a regular lattice. Then
we move each of the particles in succession according
to the following prescription: 
	x becomes x + a * epsilon1
	y becomes y + a * epsilon2
where a is the maximum allowed displacement, which
for the sake of this argument is arbitrary, and epsilon1 and epsilon2
are random numbers between (-1) and 1. Then, after
we move a particle, it is equally likely to be anywhere
within a square of side 2a centered about its original
position. (In accord with the periodicity assumption,
if the indicated move would put the particle outside the
square, this only means that it re-enters the square from
the opposite side.) 
*/

static void print_info(const struct mc *mc)
{
	double e_kin = get_kinetic_energy(mc);
	double invariant = mc->get_invariant(mc);
	double temperature = get_temperature(mc);

	print_energy(mc->state);

	msg("%30s %16.10lf\n\n", "MC ENERGY", mc->mc_energy);

	msg("%30s %16.10lf\n", "KINETIC ENERGY", e_kin);
	msg("%30s %16.10lf\n", "INVARIANT", invariant);
	msg("%30s %16.10lf\n", "TEMPERATURE (K)", temperature);

	if (cfg_get_enum(mc->state->cfg, "ensemble") == ENSEMBLE_TYPE_NPT) {
		double pressure = get_pressure(mc) / BAR_TO_AU;

		msg("%30s %16.10lf bar\n", "PRESSURE", pressure);
	}

	if (cfg_get_bool(mc->state->cfg, "enable_pbc")) {
		double x = mc->box.x * BOHR_RADIUS;
		double y = mc->box.y * BOHR_RADIUS;
		double z = mc->box.z * BOHR_RADIUS;

		msg("%30s %9.3lf %9.3lf %9.3lf A^3\n",
		    "PERIODIC BOX SIZE", x, y, z);
	}

	msg("\n\n");
}

static void print_restart(const struct mc *mc)
{
	msg("    RESTART DATA\n\n");

	for (size_t i = 0; i < mc->n_bodies; i++) {
		struct body *body = mc->bodies + i;

		char name[64];
		check_fail(efp_get_frag_name(mc->state->efp, i,
		    sizeof(name), name));

		double xyzabc[6] = { body->pos.x * BOHR_RADIUS,
				     body->pos.y * BOHR_RADIUS,
				     body->pos.z * BOHR_RADIUS };

		matrix_to_euler(&body->rotmat,
		    xyzabc + 3, xyzabc + 4, xyzabc + 5);

		double vel[6] = { body->vel.x,
				  body->vel.y,
				  body->vel.z,
				  body->angmom.x * body->inertia_inv.x,
				  body->angmom.y * body->inertia_inv.y,
				  body->angmom.z * body->inertia_inv.z };

		print_fragment(name, xyzabc, vel);
	}

	msg("\n");
}

static struct mc *mc_create(struct state *state)
{
	struct mc *mc = xcalloc(1, sizeof(struct mc));

	mc->state = state;
	mc->box = box_from_str(cfg_get_string(state->cfg, "periodic_box"));

	mc->get_invariant = get_invariant_nve;
	mc->update_step = update_step_mc;

	mc->n_bodies = state->sys->n_frags;
	mc->bodies = xcalloc(mc->n_bodies, sizeof(struct body));
	mc->xr_gradient = xcalloc(6 * mc->n_bodies, sizeof(double));

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

		body->vel.x = mc->state->sys->frags[i].vel[0];
		body->vel.y = mc->state->sys->frags[i].vel[1];
		body->vel.z = mc->state->sys->frags[i].vel[2];

		set_body_mass_and_inertia(state->efp, i, body);

		body->angmom.x = mc->state->sys->frags[i].vel[3] *
		    body->inertia.x;
		body->angmom.y = mc->state->sys->frags[i].vel[4] *
		    body->inertia.y;
		body->angmom.z = mc->state->sys->frags[i].vel[5] *
		    body->inertia.z;

		mc->n_freedom += 3;

		if (body->inertia.x > EPSILON)
			mc->n_freedom++;
		if (body->inertia.y > EPSILON)
			mc->n_freedom++;
		if (body->inertia.z > EPSILON)
			mc->n_freedom++;
	}

	return (mc);
}

static void velocitize(struct mc *mc)
{
	rand_init();

	double temperature = cfg_get_double(mc->state->cfg, "temperature");
	double ke = temperature * BOLTZMANN *
	    mc->n_freedom / (2.0 * 6.0 * mc->n_bodies);

	for (size_t i = 0; i < mc->n_bodies; i++) {
		struct body *body = mc->bodies + i;

		double vel = sqrt(2.0 * ke / body->mass);

		body->vel.x = vel * rand_normal();
		body->vel.y = vel * rand_normal();
		body->vel.z = vel * rand_normal();

		body->angmom.x = sqrt(2.0 * ke * body->inertia.x) *
		    rand_normal();
		body->angmom.y = sqrt(2.0 * ke * body->inertia.y) *
		    rand_normal();
		body->angmom.z = sqrt(2.0 * ke * body->inertia.z) *
		    rand_normal();
	}
}

static void print_status(const struct mc *mc)
{
	print_geometry(mc->state->efp);
	print_restart(mc);
	print_info(mc);

	fflush(stdout);
}

static void mc_shutdown(struct mc *mc)
{
	free(mc->bodies);
	free(mc->xr_gradient);
	free(mc->data);
	free(mc);
}

void sim_mc(struct state *state)
{
	msg("MONTE CARLO JOB\n\n\n");

	struct mc *mc = mc_create(state);

	if (cfg_get_bool(state->cfg, "velocitize"))
		velocitize(mc);

	remove_system_drift(mc);
	compute_forces(mc);

	msg("    INITIAL STATE\n\n");
	print_status(mc);

	for (mc->step = 1;
	    mc->step <= cfg_get_int(state->cfg, "max_steps");
	    mc->step++) {
		mc->update_step(mc);

		if (mc->step % cfg_get_int(state->cfg, "print_step") == 0) {
			msg("    STATE AFTER %d STEPS\n\n", mc->step);
			print_status(mc);
		}
	}

	mc_shutdown(mc);

	msg("MONTE CARLO JOB COMPLETED SUCCESSFULLY\n");
}
