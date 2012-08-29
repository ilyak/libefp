#include "../common/math_util.h"

#include "common.h"
#include "sim.h"

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
	double mass;
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
};

static double get_kinetic_energy(const struct md *md)
{
	double ke = 0.0;

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		ke += body->mass * body->vel.x * body->vel.x;
		ke += body->mass * body->vel.y * body->vel.y;
		ke += body->mass * body->vel.z * body->vel.z;

		if (body->inertia.x > EPSILON) /* non-linear body */
			ke += body->angmom.x * body->angmom.x / body->inertia.x;

		ke += body->angmom.y * body->angmom.y / body->inertia.y;
		ke += body->angmom.z * body->angmom.z / body->inertia.z;
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

static double get_invariant_nvt(UNUSED const struct md *md)
{
	assert(0);
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
	mat_t inertia_inv;

	if (inertia.xx < EPSILON ||
	    inertia.yy < EPSILON ||
	    inertia.zz < EPSILON) {
		mat_zero(&inertia_inv);
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

static void compute_forces(struct md *md)
{
	enum efp_result res;
	struct efp_energy energy;
	double coord[6 * md->n_bodies];
	double grad[6 * md->n_bodies];

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		double a, b, c;
		matrix_to_euler(&body->rotmat, &a, &b, &c);

		coord[6 * i + 0] = body->pos.x;
		coord[6 * i + 1] = body->pos.y;
		coord[6 * i + 2] = body->pos.z;

		coord[6 * i + 3] = a;
		coord[6 * i + 4] = b;
		coord[6 * i + 5] = c;
	}

	if ((res = efp_set_coordinates(md->efp, EFP_COORD_TYPE_XYZABC, coord)))
		lib_error(res);

	if ((res = efp_compute(md->efp, 1)))
		lib_error(res);

	if ((res = efp_get_energy(md->efp, &energy)))
		lib_error(res);

	md->potential_energy = energy.total;

	if ((res = efp_get_gradient(md->efp, md->n_bodies, grad)))
		lib_error(res);

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->force.x = -grad[6 * i + 0];
		body->force.y = -grad[6 * i + 1];
		body->force.z = -grad[6 * i + 2];

		body->torque.x = -grad[6 * i + 3];
		body->torque.y = -grad[6 * i + 4];
		body->torque.z = -grad[6 * i + 5];

		mat_t rotmat = mat_transpose(&body->rotmat);
		body->torque = mat_vec(&rotmat, &body->torque);
	}
}

static void set_body_mass_and_inertia(struct body *body, int n_atoms,
		const struct efp_atom *atoms)
{
	mat_t inertia = { 0.0, 0.0, 0.0,
			  0.0, 0.0, 0.0,
			  0.0, 0.0, 0.0 };
	body->mass = 0.0;

	for (int i = 0; i < n_atoms; i++) {
		const struct efp_atom *atom = atoms + i;

		vec_t dr = { atom->x - body->pos.x,
			     atom->y - body->pos.y,
			     atom->z - body->pos.z };

		mat_t rotmat = mat_transpose(&body->rotmat);
		dr = mat_vec(&rotmat, &dr);

		inertia.xx += atom->mass * (dr.y * dr.y + dr.z * dr.z);
		inertia.yy += atom->mass * (dr.x * dr.x + dr.z * dr.z);
		inertia.zz += atom->mass * (dr.x * dr.x + dr.y * dr.y);
		inertia.xy -= atom->mass * dr.x * dr.y;
		inertia.xz -= atom->mass * dr.x * dr.z;
		inertia.yz -= atom->mass * dr.y * dr.z;

		body->mass += AMU_TO_AU * atom->mass;
	}

	if (fabs(inertia.xy) > EPSILON ||
	    fabs(inertia.xz) > EPSILON ||
	    fabs(inertia.yz) > EPSILON ||
	    inertia.xx > inertia.yy ||
	    inertia.yy > inertia.zz ||
	    inertia.yy < EPSILON)
		error("INCORRECT FRAGMENT FRAME");

	body->inertia.x = AMU_TO_AU * inertia.xx;
	body->inertia.y = AMU_TO_AU * inertia.yy;
	body->inertia.z = AMU_TO_AU * inertia.zz;
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
	*rotmat = mat_mat(&rot, rotmat);
}

static void rotate_body(struct body *body, double dt)
{
	double angle;
	mat_t rotmat = mat_transpose(&body->rotmat);

	/* rotate about x axis */
	if (body->inertia.x > EPSILON) { /* non-linear body */
		angle = 0.5 * dt * body->angmom.x / body->inertia.x;
		rotate_step(1, 2, angle, &body->angmom, &rotmat);
	}

	/* rotate about y axis */
	angle = 0.5 * dt * body->angmom.y / body->inertia.y;
	rotate_step(2, 0, angle, &body->angmom, &rotmat);

	/* rotate about z axis */
	angle = dt * body->angmom.z / body->inertia.z;
	rotate_step(0, 1, angle, &body->angmom, &rotmat);

	/* rotate about y axis */
	angle = 0.5 * dt * body->angmom.y / body->inertia.y;
	rotate_step(2, 0, angle, &body->angmom, &rotmat);

	/* rotate about x axis */
	if (body->inertia.x > EPSILON) { /* non-linear body */
		angle = 0.5 * dt * body->angmom.x / body->inertia.x;
		rotate_step(1, 2, angle, &body->angmom, &rotmat);
	}

	body->rotmat = mat_transpose(&rotmat);
}

static void update_step_nve(struct md *md)
{
	double dt = FS_TO_AU * md->config->time_step;

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

static void update_step_nvt(UNUSED struct md *md)
{
	assert(0);
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
	char name[64];
	enum efp_result res;

	printf("    RESTART DATA (ATOMIC UNITS)\n\n");

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		if ((res = efp_get_frag_name(md->efp, i, sizeof(name), name)))
			lib_error(res);

		double xyzabc[6] = { body->pos.x,
				     body->pos.y,
				     body->pos.z };

		matrix_to_euler(&body->rotmat, xyzabc + 3, xyzabc + 4, xyzabc + 5);

		double vel[6] = { body->vel.x,
				  body->vel.y,
				  body->vel.z,
				  body->inertia.x > EPSILON ?
					body->angmom.x / body->inertia.x : 0.0,
				  body->angmom.y / body->inertia.y,
				  body->angmom.z / body->inertia.z };

		print_fragment(name, xyzabc, vel);
	}

	printf("\n");
}

static struct md *md_create(struct efp *efp, const struct config *config)
{
	struct md *md = calloc(1, sizeof(struct md));

	if (!md)
		error("NO MEMORY");

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
			break;
		default:
			assert(0);
	}

	enum efp_result res;

	if ((res = efp_get_frag_count(efp, &md->n_bodies)))
		lib_error(res);

	md->bodies = calloc(md->n_bodies, sizeof(struct body));

	if (!md->bodies)
		error("NO MEMORY");

	double coord[6 * md->n_bodies];

	if ((res = efp_get_coordinates(efp, md->n_bodies, coord)))
		lib_error(res);

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

		int n_atoms;
		if ((res = efp_get_frag_atom_count(efp, i, &n_atoms)))
			lib_error(res);

		struct efp_atom atoms[n_atoms];
		if ((res = efp_get_frag_atoms(efp, i, n_atoms, atoms)))
			lib_error(res);

		set_body_mass_and_inertia(body, n_atoms, atoms);

		body->angmom.x = md->config->frags[i].vel[3] * body->inertia.x;
		body->angmom.y = md->config->frags[i].vel[4] * body->inertia.y;
		body->angmom.z = md->config->frags[i].vel[5] * body->inertia.z;

		md->n_freedom += body->inertia.x > EPSILON ? 6 : 5;
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
	free(md);
}

void sim_md(struct efp *efp, const struct config *config)
{
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
}
