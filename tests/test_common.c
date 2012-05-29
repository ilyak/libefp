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

static void
message(const char *msg)
{
	fprintf(stderr, "%s\n", msg);
}

static void
error(const char *title, enum efp_result res)
{
	fprintf(stderr, "%s:\n", title);
	fprintf(stderr, "    %s\n", efp_result_to_string(res));
}

static inline int
eq(double a, double b)
{
	static const double eps = 1.0e-6;
	return fabs(a - b) < eps;
}

static int
test_numerical_gradient(struct efp *efp,
			const double *xyzabc,
			const double *grad)
{
	static const double grad_delta = 0.0001;

	int n_frag = efp_get_frag_count(efp);
	int n_grad = 6 * n_frag;

	double xyzabc_new[n_grad];
	memcpy(xyzabc_new, xyzabc, n_grad * sizeof(double));

	int result = 0;

	for (int i = 0; i < n_frag; i++) {
		double grad_num[6], grad_euler[6];

		for (int j = 0; j < 6; j++) {
			struct efp_energy e1;
			xyzabc_new[6 * i + j] = xyzabc[6 * i + j] - grad_delta;
			efp_set_coordinates(efp, xyzabc_new);
			efp_compute(efp, 0);
			efp_get_energy(efp, &e1);

			struct efp_energy e2;
			xyzabc_new[6 * i + j] = xyzabc[6 * i + j] + grad_delta;
			efp_set_coordinates(efp, xyzabc_new);
			efp_compute(efp, 0);
			efp_get_energy(efp, &e2);

			xyzabc_new[6 * i + j] = xyzabc[6 * i + j];
			grad_num[j] = (e2.total - e1.total) / (2 * grad_delta);
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
				result = 1;
	}
	return result;
}

static inline enum efp_result
print_atoms(struct efp *efp)
{
	enum efp_result res;
	int n_frags = efp_get_frag_count(efp);

	for (int i = 0; i < n_frags; i++) {
		int n_atoms = efp_get_frag_atom_count(efp, i);

		struct efp_atom atoms[n_atoms];
		if ((res = efp_get_frag_atoms(efp, i, atoms)))
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

enum efp_result
st_integrals_from_file(const struct efp_st_block *block,
		       int compute_derivatives,
		       struct efp_st_data *st, void *user_data,
		       int expected_size_i, int expected_size_j,
		       const char *s_path, const char *t_path)
{
	if (compute_derivatives)
		return EFP_RESULT_NOT_IMPLEMENTED;

	if (st->size_i != expected_size_i ||
	    st->size_j != expected_size_j)
		return EFP_RESULT_INVALID_ARRAY_SIZE;

	FILE *fp;
	double *ptr;

	int size = st->size_i * st->size_j;

	fp = fopen(s_path, "r");
	if (!fp)
		return EFP_RESULT_FILE_NOT_FOUND;

	ptr = st->s;
	for (int i = 0; i < size; i++, ptr++)
		fscanf(fp, "%lf", ptr);

	fclose(fp);

	fp = fopen(t_path, "r");
	if (!fp)
		return EFP_RESULT_FILE_NOT_FOUND;

	ptr = st->t;
	for (int i = 0; i < size; i++, ptr++)
		fscanf(fp, "%lf", ptr);

	fclose(fp);
	return EFP_RESULT_SUCCESS;
}

int
run_test(const struct test_data *test_data)
{
	enum efp_result res;
	int status = EXIT_SUCCESS;
	struct efp *efp;

	if ((res = efp_init(&efp, &test_data->opts, &test_data->callbacks,
			test_data->potential_files, test_data->fragname))) {
		error("efp_init", res);
		goto fail;
	}

	if (test_data->geometry_xyzabc) {
		if ((res = efp_set_coordinates(efp,
					test_data->geometry_xyzabc))) {
			error("efp_update_fragments", res);
			goto fail;
		}
	}
	else {
		if ((res = efp_set_coordinates_2(efp,
					test_data->geometry_points))) {
			error("efp_update_fragments_2", res);
			goto fail;
		}
	}

	if ((res = print_atoms(efp))) {
		error("print_atoms", res);
		goto fail;
	}

	if ((res = efp_scf_init(efp))) {
		error("efp_scf_init", res);
		goto fail;
	}

	/* Begin imaginary SCF */
	double scf_energy;
	if ((res = efp_scf_update(efp, &scf_energy))) {
		error("efp_scf_update", res);
		goto fail;
	}
	/* End imaginary SCF */

	int do_gradient = test_data->ref_gradient ||
			  test_data->test_numerical_gradient;

	if ((res = efp_compute(efp, do_gradient))) {
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

	if (do_gradient) {
		int n_frag = efp_get_frag_count(efp);
		int n_grad = 6 * n_frag;
		double grad[n_grad];

		if ((res = efp_get_gradient(efp, n_grad, grad))) {
			error("efp_get_gradient", res);
			goto fail;
		}

		if (test_data->ref_gradient) {
			int wrong_grad = 0;

			for (int i = 0; i < n_grad; i++)
				if (!eq(grad[i], test_data->ref_gradient[i]))
					wrong_grad = 1;

			if (wrong_grad) {
				message("wrong gradient");
				status = EXIT_FAILURE;
			}
		}

		if (test_data->test_numerical_gradient) {
			if(test_numerical_gradient(efp,
					test_data->geometry_xyzabc, grad)) {
				message("wrong numerical gradient");
				status = EXIT_FAILURE;
			}
		}
	}

fail:
	efp_shutdown(efp);
	return res ? EXIT_FAILURE : status;
}
