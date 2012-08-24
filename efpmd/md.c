#include <assert.h>

#include "../common/math_util.h"

#include "common.h"
#include "md.h"

struct body {
	mat_t rotmat;
	vec_t pos;
	vec_t vel;
	vec_t vel_old;
	vec_t angvel;
	vec_t angvel_old;
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
	double invariant;
	struct efp *efp;
	const struct config *config;
};

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

static vec_t get_system_angular_velocity(const struct md *md)
{
	vec_t cp = get_system_com(md);
	vec_t cv = get_system_com_velocity(md);

	vec_t av = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		vec_t dr = vec_sub(&body->pos, &cp);
		vec_t dv = vec_sub(&body->vel, &cv);

		av.x += dr.y * dv.z - dr.z * dv.y;
		av.y += dr.z * dv.x - dr.x * dv.z;
		av.z += dr.x * dv.y - dr.y * dv.x;
	}

	return av;
}

static void remove_system_drift(struct md *md)
{
	vec_t cp = get_system_com(md);
	vec_t cv = get_system_com_velocity(md);
	vec_t av = get_system_angular_velocity(md);

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
}

static double get_atoms_mass(int n_atoms, const struct efp_atom *atoms)
{
	double mass = 0.0;

	for (int i = 0; i < n_atoms; i++)
		mass += atoms[i].mass;

	return mass;
}

static vec_t get_atoms_com(int n_atoms, const struct efp_atom *atoms)
{
	double mass = 0.0;
	vec_t com = { 0.0, 0.0, 0.0 };

	for (int i = 0; i < n_atoms; i++) {
		const struct efp_atom *atom = atoms + i;

		com.x += atom->mass * atom->x;
		com.y += atom->mass * atom->y;
		com.z += atom->mass * atom->z;

		mass += atom->mass;
	}

	com.x /= mass;
	com.y /= mass;
	com.z /= mass;

	return com;
}

static mat_t get_atoms_inertia_tensor(int n_atoms, const struct efp_atom *atoms)
{
	mat_t inertia = { 0.0, 0.0, 0.0,
			  0.0, 0.0, 0.0,
			  0.0, 0.0, 0.0 };

	vec_t com = get_atoms_com(n_atoms, atoms);

	for (int i = 0; i < n_atoms; i++) {
		const struct efp_atom *atom = atoms + i;

		vec_t dr = { atom->x - com.x,
			     atom->y - com.y,
			     atom->z - com.z };

		inertia.xx += atom->mass * (dr.y * dr.y + dr.z * dr.z);
		inertia.yy += atom->mass * (dr.x * dr.x + dr.z * dr.z);
		inertia.zz += atom->mass * (dr.x * dr.x + dr.y * dr.y);
		inertia.xy -= atom->mass * dr.x * dr.y;
		inertia.xz -= atom->mass * dr.x * dr.z;
		inertia.yz -= atom->mass * dr.y * dr.z;
	}

	inertia.yx = inertia.xy;
	inertia.zx = inertia.xz;
	inertia.zy = inertia.yz;

	return inertia;
}

/* XXX */
#include <mkl_lapacke.h>

static void mat_eigen(const mat_t *mat, mat_t *evec, vec_t *eval)
{
	*evec = *mat;

	LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', 3, (double *)evec, 3, (double *)eval);
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
	}
}

struct md *md_create(struct efp *efp, const struct config *config)
{
	struct md *md = calloc(1, sizeof(struct md));

	if (!md)
		memory_error();

	md->efp = efp;
	md->config = config;

	enum efp_result res;

	if ((res = efp_get_frag_count(efp, &md->n_bodies)))
		lib_error(res);

	md->bodies = calloc(md->n_bodies, sizeof(struct body));

	if (!md->bodies)
		memory_error();

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x = md->config->frags[i].vel[0];
		body->vel.y = md->config->frags[i].vel[1];
		body->vel.z = md->config->frags[i].vel[2];

		body->angvel.x = md->config->frags[i].vel[3];
		body->angvel.y = md->config->frags[i].vel[4];
		body->angvel.z = md->config->frags[i].vel[5];

		int n_atoms;
		if ((res = efp_get_frag_atom_count(efp, i, &n_atoms)))
			lib_error(res);

		struct efp_atom atoms[n_atoms];
		if ((res = efp_get_frag_atoms(efp, i, n_atoms, atoms)))
			lib_error(res);

		body->pos = get_atoms_com(n_atoms, atoms);
		body->mass = get_atoms_mass(n_atoms, atoms);

		mat_t inertia = get_atoms_inertia_tensor(n_atoms, atoms);
		mat_eigen(&inertia, &body->rotmat, &body->inertia);

		double det = mat_det(&body->rotmat);
		if (det < 0.0)
			mat_negate(&body->rotmat);

		body->inertia.x *= AMU_TO_AU;
		body->inertia.y *= AMU_TO_AU;
		body->inertia.z *= AMU_TO_AU;

		if (fabs(body->inertia.x) < EPSILON)
			md->n_freedom += 5;
		else
			md->n_freedom += 6;
	}

	remove_system_drift(md);

	vec_t cv = get_system_com_velocity(md);
	vec_t av = get_system_angular_velocity(md);

	double cv2 = vec_len(&cv);
	double av2 = vec_len(&av);

	assert(cv2 < EPSILON && av2 < EPSILON);

	compute_forces(md);

	return md;
}

static void rotate_step(int a1, int a2, double angle, vec_t *angvel, mat_t *rotmat)
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

	vec_t av = *angvel;
	mat_vec(&rot, &av, angvel);

	mat_t rm = *rotmat;
	mat_mat(&rot, &rm, rotmat);
}

static void rotate_body(struct body *body, double dt)
{
	double angle;

	/* rotate about z axis */
	angle = 0.5 * dt * body->angvel.z;
	rotate_step(0, 1, angle, &body->angvel, &body->rotmat);

	/* rotate about y axis */
	angle = 0.5 * dt * body->angvel.y;
	rotate_step(2, 0, angle, &body->angvel, &body->rotmat);

	/* rotate about x axis */
	angle = dt * body->angvel.x;
	rotate_step(1, 2, angle, &body->angvel, &body->rotmat);

	/* rotate about y axis */
	angle = 0.5 * dt * body->angvel.y;
	rotate_step(2, 0, angle, &body->angvel, &body->rotmat);

	/* rotate about z axis */
	angle = 0.5 * dt * body->angvel.z;
	rotate_step(0, 1, angle, &body->angvel, &body->rotmat);
}

static double get_kinetic_energy(const struct md *md)
{
	double ke = 0.0;

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		ke += body->mass * body->vel.x * body->vel.x;
		ke += body->mass * body->vel.y * body->vel.y;
		ke += body->mass * body->vel.z * body->vel.z;

		ke += body->inertia.x * body->angvel.x * body->angvel.x;
		ke += body->inertia.y * body->angvel.y * body->angvel.y;
		ke += body->inertia.z * body->angvel.z * body->angvel.z;
	}

	return 0.5 * ke;
}

static double get_temperature(const struct md *md)
{
	double ke = get_kinetic_energy(md);
	return 2.0 * ke / BOLTZMANN / md->n_freedom;
}

void md_update_step_nve(struct md *md)
{
	double dt = FS_TO_AU * md->config->time_step;

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		body->vel.x += 0.5 * body->force.x * dt / body->mass;
		body->vel.y += 0.5 * body->force.y * dt / body->mass;
		body->vel.z += 0.5 * body->force.z * dt / body->mass;

		body->angvel.x += 0.5 * body->torque.x * dt / body->mass;
		body->angvel.y += 0.5 * body->torque.y * dt / body->mass;
		body->angvel.z += 0.5 * body->torque.z * dt / body->mass;

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

		body->angvel.x += 0.5 * body->torque.x * dt / body->mass;
		body->angvel.y += 0.5 * body->torque.y * dt / body->mass;
		body->angvel.z += 0.5 * body->torque.z * dt / body->mass;
	}

	md->invariant = md->potential_energy + get_kinetic_energy(md);
}

void md_update_step_nvt(UNUSED struct md *md)
{
}

void md_print_info(struct md *md)
{
	printf("    POTENTIAL ENERGY: %20.10lf\n", md->potential_energy);
	printf("           INVARIANT: %20.10lf\n", md->invariant);
	printf("  SYSTEM TEMPERATURE: %20.10lf\n", get_temperature(md));
	printf("\n\n");
}

static void print_fragment(const char *name, const struct body *body)
{
	double a, b, c;
	matrix_to_euler(&body->rotmat, &a, &b, &c);

	printf("fragment %s\n", name);
	printf(" %12.6lf %12.6lf %12.6lf", body->pos.x, body->pos.y, body->pos.z);
	printf(" %12.6lf %12.6lf %12.6lf", a, b, c);
	printf("\nvelocity\n");
	printf(" %12.6lf %12.6lf %12.6lf", body->vel.x, body->vel.y, body->vel.z);
	printf(" %12.6lf %12.6lf %12.6lf", body->angvel.x, body->angvel.y, body->angvel.z);
	printf("\n\n");
}

void md_print_geometry(struct md *md)
{
	char name[64];
	enum efp_result res;

	print_geometry(md->efp);
	printf("    RESTART DATA\n\n");

	for (int i = 0; i < md->n_bodies; i++) {
		struct body *body = md->bodies + i;

		if ((res = efp_get_frag_name(md->efp, i, sizeof(name), name)))
			lib_error(res);

		print_fragment(name, body);
	}
}

void md_shutdown(struct md *md)
{
	free(md->bodies);
	free(md);
}
