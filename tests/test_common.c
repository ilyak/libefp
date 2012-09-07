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

#include "test_common.h"

#define NUM_GRAD_DELTA 0.001

static void message(const char *msg)
{
	fprintf(stderr, "%s\n", msg);
}

static void error(const char *title, enum efp_result res)
{
	fprintf(stderr, "%s:\n", title);
	fprintf(stderr, "    %s\n", efp_result_to_string(res));
}

static int eq(double a, double b)
{
	static const double eps = 5.0e-6;
	return fabs(a - b) < eps;
}

static enum efp_result test_qm_numerical_grad(struct efp *efp,
			const double *grad, int *fail)
{
	enum efp_result res;

	int n_qm_atoms;
	if ((res = efp_get_qm_atom_count(efp, &n_qm_atoms)))
		return res;

	double znuc[n_qm_atoms], xyz[3 * n_qm_atoms];
	if ((res = efp_get_qm_atoms(efp, n_qm_atoms, znuc, xyz)))
		return res;

	*fail = 0;

	for (int i = 0; i < n_qm_atoms; i++) {
		double grad_num[3];

		for (int j = 0; j < 3; j++) {
			double coord = xyz[3 * i + j];

			struct efp_energy e1;
			xyz[3 * i + j] = coord - NUM_GRAD_DELTA;
			if ((res = efp_set_qm_atoms(efp, n_qm_atoms, znuc, xyz)))
				return res;
			if ((res = efp_compute(efp, 0)))
				return res;
			if ((res = efp_get_energy(efp, &e1)))
				return res;

			struct efp_energy e2;
			xyz[3 * i + j] = coord + NUM_GRAD_DELTA;
			if ((res = efp_set_qm_atoms(efp, n_qm_atoms, znuc, xyz)))
				return res;
			if ((res = efp_compute(efp, 0)))
				return res;
			if ((res = efp_get_energy(efp, &e2)))
				return res;

			xyz[3 * i + j] = coord;
			grad_num[j] = (e2.total - e1.total) / (2.0 * NUM_GRAD_DELTA);
		}

		for (int j = 0; j < 3; j++)
			if (!eq(grad_num[j], grad[3 * i + j]))
				*fail = 1;
	}

	if ((res = efp_set_qm_atoms(efp, n_qm_atoms, znuc, xyz)))
		return res;

	return EFP_RESULT_SUCCESS;
}

static enum efp_result test_frag_numerical_grad(struct efp *efp,
			const double *grad, int *fail)
{
	enum efp_result res;

	int n_frag;
	if ((res = efp_get_frag_count(efp, &n_frag)))
		return res;

	double xyzabc[6 * n_frag];
	if ((res = efp_get_coordinates(efp, n_frag, xyzabc)))
		return res;

	*fail = 0;

	for (int i = 0; i < n_frag; i++) {
		double grad_num[6], grad_euler[6];

		for (int j = 0; j < 6; j++) {
			double coord = xyzabc[6 * i + j];

			struct efp_energy e1;
			xyzabc[6 * i + j] = coord - NUM_GRAD_DELTA;
			if ((res = efp_set_coordinates(efp, EFP_COORD_TYPE_XYZABC, xyzabc)))
				return res;
			if ((res = efp_compute(efp, 0)))
				return res;
			if ((res = efp_get_energy(efp, &e1)))
				return res;

			struct efp_energy e2;
			xyzabc[6 * i + j] = coord + NUM_GRAD_DELTA;
			if ((res = efp_set_coordinates(efp, EFP_COORD_TYPE_XYZABC, xyzabc)))
				return res;
			if ((res = efp_compute(efp, 0)))
				return res;
			if ((res = efp_get_energy(efp, &e2)))
				return res;

			xyzabc[6 * i + j] = coord;
			grad_num[j] = (e2.total - e1.total) / (2.0 * NUM_GRAD_DELTA);
		}

		/* convert torque to energy derivatives by Euler angles */

		double tx = grad[6 * i + 3];
		double ty = grad[6 * i + 4];
		double tz = grad[6 * i + 5];

		double a = xyzabc[6 * i + 3];
		double b = xyzabc[6 * i + 4];

		double sina = sin(a);
		double cosa = cos(a);
		double sinb = sin(b);
		double cosb = cos(b);

		grad_euler[0] = grad[6 * i + 0];
		grad_euler[1] = grad[6 * i + 1];
		grad_euler[2] = grad[6 * i + 2];
		grad_euler[3] = tz;
		grad_euler[4] = cosa * tx + sina * ty;
		grad_euler[5] = sinb * sina * tx - sinb * cosa * ty + cosb * tz;

		for (int j = 0; j < 6; j++)
			if (!eq(grad_num[j], grad_euler[j]))
				*fail = 1;
	}

	if ((res = efp_set_coordinates(efp, EFP_COORD_TYPE_XYZABC, xyzabc)))
		return res;

	return EFP_RESULT_SUCCESS;
}

static enum efp_result print_atoms(struct efp *efp)
{
	enum efp_result res;

	int n_qm_atoms;
	if ((res = efp_get_qm_atom_count(efp, &n_qm_atoms)))
		return res;

	double qm_znuc[n_qm_atoms];
	double qm_xyz[3 * n_qm_atoms];

	if ((res = efp_get_qm_atoms(efp, n_qm_atoms, qm_znuc, qm_xyz)))
		return res;

	for (int i = 0; i < n_qm_atoms; i++) {
		double x = qm_xyz[3 * i + 0];
		double y = qm_xyz[3 * i + 1];
		double z = qm_xyz[3 * i + 2];

		printf("QQ %12.8lf %12.8lf %12.8lf\n",
			ANGSTROM(x), ANGSTROM(y), ANGSTROM(z));
	}

	int n_frag;
	if ((res = efp_get_frag_count(efp, &n_frag)))
		return res;

	for (int i = 0; i < n_frag; i++) {
		int n_atoms;
		if ((res = efp_get_frag_atom_count(efp, i, &n_atoms)))
			return res;

		struct efp_atom atoms[n_atoms];
		if ((res = efp_get_frag_atoms(efp, i, n_atoms, atoms)))
			return res;

		for (int a = 0; a < n_atoms; a++) {
			struct efp_atom *atom = atoms + a;
			double x = atom->x, y = atom->y, z = atom->z;

			printf("%s %12.8lf %12.8lf %12.8lf\n", atom->label,
				ANGSTROM(x), ANGSTROM(y), ANGSTROM(z));
		}
	}
	return EFP_RESULT_SUCCESS;
}

int run_test(const struct test_data *test_data)
{
	enum efp_result res;
	struct efp *efp;
	int status = EXIT_SUCCESS;

	if ((res = efp_init(&efp, &test_data->opts, &test_data->callbacks,
			test_data->potential_files, test_data->fragname))) {
		error("efp_init", res);
		goto fail;
	}

	const double *geometry;
	enum efp_coord_type coord_type;

	if (test_data->geometry_xyzabc) {
		geometry = test_data->geometry_xyzabc;
		coord_type = EFP_COORD_TYPE_XYZABC;
	}
	else {
		geometry = test_data->geometry_points;
		coord_type = EFP_COORD_TYPE_POINTS;
	}

	if ((res = efp_set_coordinates(efp, coord_type, geometry))) {
		error("efp_update_fragments", res);
		goto fail;
	}

	int do_qm = test_data->qm_znuc && test_data->qm_xyz;

	if (do_qm) {
		if ((res = efp_set_qm_atoms(efp, test_data->n_qm_atoms,
				test_data->qm_znuc, test_data->qm_xyz))) {
			error("efp_set_qm_atoms", res);
			goto fail;
		}
	}

	if ((res = print_atoms(efp))) {
		error("print_atoms", res);
		goto fail;
	}

	/* Begin imaginary ab initio SCF */
	double scf_energy;
	if ((res = efp_scf_update(efp, &scf_energy))) {
		error("efp_scf_update", res);
		goto fail;
	}
	/* End imaginary ab initio SCF */

	if ((res = efp_compute(efp, test_data->test_gradient))) {
		error("efp_compute", res);
		goto fail;
	}

	struct efp_energy energy;

	if ((res = efp_get_energy(efp, &energy))) {
		error("efp_get_energy", res);
		goto fail;
	}

	if (!eq(energy.total, test_data->ref_energy)) {
		message("wrong energy");
		status = EXIT_FAILURE;
	}

	if (test_data->test_gradient) {
		int failed;

		int n_frag;
		if ((res = efp_get_frag_count(efp, &n_frag))) {
			error("efp_get_frag_count", res);
			goto fail;
		}

		double frag_grad[6 * n_frag];
		if ((res = efp_get_gradient(efp, n_frag, frag_grad))) {
			error("efp_get_gradient", res);
			goto fail;
		}

		if (do_qm) {
			int n_qm_atoms;
			if ((res = efp_get_qm_atom_count(efp, &n_qm_atoms))) {
				error("efp_get_qm_atom_count", res);
				goto fail;
			}

			double qm_grad[3 * n_qm_atoms];
			if ((res = efp_get_qm_gradient(efp, n_qm_atoms, qm_grad))) {
				error("efp_get_qm_gradient", res);
				goto fail;
			}

			if ((res = test_qm_numerical_grad(efp, qm_grad, &failed))) {
				error("test_qm_numerical_grad", res);
				goto fail;
			}
			if (failed) {
				message("wrong qm gradient");
				status = EXIT_FAILURE;
			}
		}

		if ((res = test_frag_numerical_grad(efp, frag_grad, &failed))) {
			error("test_frag_numerical_grad", res);
			goto fail;
		}
		if (failed) {
			message("wrong fragment gradient");
			status = EXIT_FAILURE;
		}
	}

fail:
	efp_shutdown(efp);
	return res ? EXIT_FAILURE : status;
}
