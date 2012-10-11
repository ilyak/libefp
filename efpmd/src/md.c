/*-
 * Copyright (c) 2012 Ilya Kaliman
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

#include "../../common/util.h"

#include "common.h"

#define HOOVER_MAX_ITER 20

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

struct nvt_data {
	double chi;
	double chi_dt;
};

struct md {
	int n_bodies;
	struct body *bodies;
	int n_freedom;
	double potential_energy;
	double (*get_invariant)(const struct md *);
	void (*update_step)(struct md *);
	struct efp *efp;
	const struct config *config;
	void *data;
};

void sim_md(struct efp *, const struct config *);

static double get_kinetic_energy(const struct md *md)
{
	double ke = 0.0;

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		ke += body->mass * body->vel.x * body->vel.x;
		ke += body->mass * body->vel.y * body->vel.y;
		ke += body->mass * body->vel.z * body->vel.z;

		ke += body->angmom.x * body->angmom.x * body->inertia_inv.x;
		ke += body->angmom.y * body->angmom.y * body->inertia_inv.y;
		ke += body->angmom.z * body->angmom.z * body->inertia_inv.z;
	}

	return 0.5 * ke;
}

static double get_temperature(const struct md *md)
{
	double ke = get_kinetic_energy(md);

	return 2.0 * ke / BOLTZMANN / md->n_freedom;
}

static double get_invariant_nve(const struct md *md)
{
	return md->potential_energy + get_kinetic_energy(md);
}

static double get_invariant_nvt(const struct md *md)
{
	struct nvt_data *data = (struct nvt_data *)md->data;

	double tau = md->config->thermostat_tau;
	double target = md->config->target_temperature;

	double virt = BOLTZMANN * target * md->n_freedom * (data->chi_dt +
			0.5 * data->chi * data->chi * tau * tau);

	return md->potential_energy + get_kinetic_energy(md) + virt;
}

static vec_t get_system_com(const struct md *md)
{
	double mass = 0.0;
	vec_t com = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		com.x += body->mass * body->pos.x;
		com.y += body->mass * body->pos.y;
		com.z += body->mass * body->pos.z;

		mass += body->mass;
	}

	com.x /= mass;
	com.y /= mass;
	com.z /= mass;

	return com;
}

static vec_t get_system_com_velocity(const struct md *md)
{
	double mass = 0.0;
	vec_t cv = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

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

static vec_t get_system_angular_momentum(const struct md *md)
{
	vec_t cp = get_system_com(md);
	vec_t cv = get_system_com_velocity(md);

	vec_t am = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		vec_t dr = vec_sub(&body->pos, &cp);
		vec_t dv = vec_sub(&body->vel, &cv);

		am.x += (dr.y * dv.z - dr.z * dv.y) * body->mass;
		am.y += (dr.z * dv.x - dr.x * dv.z) * body->mass;
		am.z += (dr.x * dv.y - dr.y * dv.x) * body->mass;
	}

	return am;
}

static mat_t get_system_inertia_tensor(const struct md *md)
{
	mat_t inertia = { 0.0, 0.0, 0.0,
			  0.0, 0.0, 0.0,
			  0.0, 0.0, 0.0 };

	vec_t com = get_system_com(md);

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		vec_t dr = { body->pos.x - com.x,
			     body->pos.y - com.y,
			     body->pos.z - com.z };

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

static void remove_system_drift(struct md *md)
{
	vec_t cp = get_system_com(md);
	vec_t cv = get_system_com_velocity(md);
	vec_t am = get_system_angular_momentum(md);

	mat_t inertia = get_system_inertia_tensor(md);
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

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		vec_t cross = {
			av.y * (body->pos.z - cp.z) - av.z * (body->pos.y - cp.y),
			av.z * (body->pos.x - cp.x) - av.x * (body->pos.z - cp.z),
			av.x * (body->pos.y - cp.y) - av.y * (body->pos.x - cp.x)
		};

		body->vel.x -= cv.x + cross.x;
		body->vel.y -= cv.y + cross.y;
		body->vel.z -= cv.z + cross.z;
	}

	vec_t cv2 = get_system_com_velocity(md);
	vec_t am2 = get_system_angular_momentum(md);

	assert(vec_len(&cv2) < EPSILON && vec_len(&am2) < EPSILON);
}

static void make_coord_array(const struct md *md, double *coord)
{
	for (int i = 0; i < md->n_bodies; i++, coord += 12) {
		struct body *body = md->bodies + i;

		coord[0] = body->pos.x;
		coord[1] = body->pos.y;
		coord[2] = body->pos.z;

		memcpy(coord + 3, &body->rotmat, sizeof(mat_t));
	}
}

static void compute_forces(struct md *md)
{
	struct efp_energy energy;
	double coord[12 * md->n_bodies];
	double grad[6 * md->n_bodies];

	make_coord_array(md, coord);

	check_fail(efp_set_coordinates(md->efp, EFP_COORD_TYPE_ROTMAT, coord));
	check_fail(efp_compute(md->efp, 1));
	check_fail(efp_get_energy(md->efp, &energy));
	check_fail(efp_get_gradient(md->efp, md->n_bodies, grad));

	md->potential_energy = energy.total;

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->force.x = -grad[6 * i + 0];
		body->force.y = -grad[6 * i + 1];
		body->force.z = -grad[6 * i + 2];

		body->torque.x = -grad[6 * i + 3];
		body->torque.y = -grad[6 * i + 4];
		body->torque.z = -grad[6 * i + 5];

		/* convert torque to body frame */
		body->torque = mat_trans_vec(&body->rotmat, &body->torque);
	}
}

static void set_body_mass_and_inertia(struct efp *efp, int idx, struct body *body)
{
	double mass, inertia[3];

	check_fail(efp_get_frag_mass(efp, idx, &mass));
	check_fail(efp_get_frag_inertia(efp, idx, inertia));

	body->mass = AMU_TO_AU * mass;

	body->inertia.x = AMU_TO_AU * inertia[0];
	body->inertia.y = AMU_TO_AU * inertia[1];
	body->inertia.z = AMU_TO_AU * inertia[2];

	body->inertia_inv.x = body->inertia.x < EPSILON ? 0.0 : 1.0 / body->inertia.x;
	body->inertia_inv.y = body->inertia.y < EPSILON ? 0.0 : 1.0 / body->inertia.y;
	body->inertia_inv.z = body->inertia.z < EPSILON ? 0.0 : 1.0 / body->inertia.z;
}

static void rotate_step(int a1, int a2, double angle, vec_t *angmom, mat_t *rotmat)
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

static void update_step_nve(struct md *md)
{
	double dt = md->config->time_step;

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

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

	compute_forces(md);

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x += 0.5 * body->force.x * dt / body->mass;
		body->vel.y += 0.5 * body->force.y * dt / body->mass;
		body->vel.z += 0.5 * body->force.z * dt / body->mass;

		body->angmom.x += 0.5 * body->torque.x * dt;
		body->angmom.y += 0.5 * body->torque.y * dt;
		body->angmom.z += 0.5 * body->torque.z * dt;
	}
}

/*
 * NVT with Nose-Hoover thermostat:
 *
 * William G. Hoover
 *
 * Canonical dynamics: Equilibrium phase-space distributions
 *
 * Phys. Rev. A 31, 1695 (1985)
 */
static void update_step_nvt(struct md *md)
{
	struct nvt_data *data = (struct nvt_data *)md->data;
	double dt = md->config->time_step;

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x += 0.5 * dt * (body->force.x / body->mass -
					body->vel.x * data->chi);
		body->vel.y += 0.5 * dt * (body->force.y / body->mass -
					body->vel.y * data->chi);
		body->vel.z += 0.5 * dt * (body->force.z / body->mass -
					body->vel.z * data->chi);

		body->angmom.x += 0.5 * dt * (body->torque.x -
					body->angmom.x * data->chi);
		body->angmom.y += 0.5 * dt * (body->torque.y -
					body->angmom.y * data->chi);
		body->angmom.z += 0.5 * dt * (body->torque.z -
					body->angmom.z * data->chi);

		body->pos.x += body->vel.x * dt;
		body->pos.y += body->vel.y * dt;
		body->pos.z += body->vel.z * dt;

		rotate_body(body, dt);
	}

	double tau = md->config->thermostat_tau;
	double target = md->config->target_temperature;

	data->chi += 0.5 * dt * (get_temperature(md) / target - 1.0) / tau / tau;
	data->chi_dt += 0.5 * dt * data->chi;

	compute_forces(md);

	double chi_init = data->chi;
	vec_t angmom_init[md->n_bodies], vel_init[md->n_bodies];

	for (int i = 0; i < md->n_bodies; i++) {
		angmom_init[i] = md->bodies[i].angmom;
		vel_init[i] = md->bodies[i].vel;
	}

	for (int iter = 0; iter < HOOVER_MAX_ITER; iter++) {
		double chi_prev = data->chi;
		double ratio = get_temperature(md) / target;

		data->chi = chi_init + 0.5 * dt * (ratio - 1.0) / tau / tau;

		for (int i = 0; i < md->n_bodies; i++) {
			struct body *body = md->bodies + i;

			body->vel.x = vel_init[i].x + 0.5 * dt *
				(body->force.x / body->mass - vel_init[i].x * data->chi);
			body->vel.y = vel_init[i].y + 0.5 * dt *
				(body->force.y / body->mass - vel_init[i].y * data->chi);
			body->vel.z = vel_init[i].z + 0.5 * dt *
				(body->force.z / body->mass - vel_init[i].z * data->chi);

			body->angmom.x = angmom_init[i].x + 0.5 * dt *
				(body->torque.x - angmom_init[i].x * data->chi);
			body->angmom.y = angmom_init[i].y + 0.5 * dt *
				(body->torque.y - angmom_init[i].y * data->chi);
			body->angmom.z = angmom_init[i].z + 0.5 * dt *
				(body->torque.z - angmom_init[i].z * data->chi);
		}

		if (fabs(data->chi - chi_prev) < EPSILON)
			break;

		if (iter == HOOVER_MAX_ITER - 1)
			error("Nose-Hoover did not converge");
	}

	data->chi_dt += 0.5 * dt * data->chi;
}

static void print_info(const struct md *md)
{
	printf("             POTENTIAL ENERGY %16.10lf\n", md->potential_energy);
	printf("                    INVARIANT %16.10lf\n", md->get_invariant(md));
	printf("           SYSTEM TEMPERATURE %16.10lf\n", get_temperature(md));
	printf("\n\n");
}

static void print_restart(const struct md *md)
{
	printf("    RESTART DATA (ATOMIC UNITS)\n\n");

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		char name[64];
		check_fail(efp_get_frag_name(md->efp, i, sizeof(name), name));

		double xyzabc[6] = { body->pos.x,
				     body->pos.y,
				     body->pos.z };

		matrix_to_euler(&body->rotmat, xyzabc + 3, xyzabc + 4, xyzabc + 5);

		double vel[6] = { body->vel.x,
				  body->vel.y,
				  body->vel.z,
				  body->angmom.x * body->inertia_inv.x,
				  body->angmom.y * body->inertia_inv.y,
				  body->angmom.z * body->inertia_inv.z };

		print_fragment(name, xyzabc, vel);
	}

	printf("\n");
}

static struct md *md_create(struct efp *efp, const struct config *config)
{
	struct md *md = xcalloc(1, sizeof(struct md));

	md->efp = efp;
	md->config = config;

	switch (config->ensemble_type) {
		case ENSEMBLE_TYPE_NVE:
			md->get_invariant = get_invariant_nve;
			md->update_step = update_step_nve;
			break;
		case ENSEMBLE_TYPE_NVT:
			md->get_invariant = get_invariant_nvt;
			md->update_step = update_step_nvt;
			md->data = xcalloc(1, sizeof(struct nvt_data));
			break;
		default:
			assert(0);
	}

	md->n_bodies = config->n_frags;
	md->bodies = xcalloc(md->n_bodies, sizeof(struct body));

	double coord[6 * md->n_bodies];
	check_fail(efp_get_coordinates(efp, md->n_bodies, coord));

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->pos.x = coord[6 * i + 0];
		body->pos.y = coord[6 * i + 1];
		body->pos.z = coord[6 * i + 2];

		double a = coord[6 * i + 3];
		double b = coord[6 * i + 4];
		double c = coord[6 * i + 5];

		euler_to_matrix(a, b, c, &body->rotmat);

		body->vel.x = md->config->frags[i].vel[0];
		body->vel.y = md->config->frags[i].vel[1];
		body->vel.z = md->config->frags[i].vel[2];

		set_body_mass_and_inertia(efp, i, body);

		body->angmom.x = md->config->frags[i].vel[3] * body->inertia.x;
		body->angmom.y = md->config->frags[i].vel[4] * body->inertia.y;
		body->angmom.z = md->config->frags[i].vel[5] * body->inertia.z;

		md->n_freedom += 3;

		if (body->inertia.x > EPSILON)
			md->n_freedom++;
		if (body->inertia.y > EPSILON)
			md->n_freedom++;
		if (body->inertia.z > EPSILON)
			md->n_freedom++;
	}

	return md;
}

static void print_status(const struct md *md)
{
	print_geometry(md->efp);
	print_restart(md);
	print_info(md);

	fflush(stdout);
}

static void md_shutdown(struct md *md)
{
	free(md->bodies);
	free(md->data);
	free(md);
}

void sim_md(struct efp *efp, const struct config *config)
{
	printf("MOLECULAR DYNAMICS JOB\n\n\n");

	struct md *md = md_create(efp, config);

	remove_system_drift(md);
	compute_forces(md);

	printf("    INITIAL STATE\n\n");
	print_status(md);

	for (int i = 1; i <= config->max_steps; i++) {
		md->update_step(md);

		if (i % config->print_step == 0) {
			printf("    STATE AFTER %d STEPS\n\n", i);
			print_status(md);
		}
	}

	md_shutdown(md);

	printf("MOLECULAR DYNAMICS JOB COMPLETED SUCCESSFULLY\n");
}
