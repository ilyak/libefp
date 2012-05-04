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
	static const double eps = 1.0e-5;
	return fabs(a - b) < eps;
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
			printf("%s %10.6lf %10.6lf %10.6lf\n", atom->label,
					atom->x, atom->y, atom->z);
		}
	}
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

	if ((res = efp_update_fragments(efp, test_data->xyzabc))) {
		error("efp_update_fragments", res);
		goto fail;
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

	if ((res = efp_compute(efp))) {
		error("efp_compute", res);
		goto fail;
	}

	double energy = 0.0;
	double energy_terms[EFP_TERM_COUNT];

	if ((res = efp_get_energy(efp, energy_terms))) {
		error("efp_get_energy", res);
		goto fail;
	}

	for (int i = 0; i < EFP_TERM_COUNT; i++)
		energy += energy_terms[i];

	if (!eq(energy, test_data->ref_energy)) {
		message("wrong energy");
		status = EXIT_FAILURE;
	}

	if (test_data->opts.do_gradient) {
		int n_grad = 6 * efp_get_frag_count(efp);
		double grad[n_grad];

		if ((res = efp_get_gradient(efp, n_grad, grad))) {
			error("efp_get_gradient", res);
			goto fail;
		}

		int wrong_grad = 0;

		for (int i = 0; i < n_grad; i++)
			if (!eq(grad[i], test_data->ref_gradient[i]))
				wrong_grad = 1;

		if (wrong_grad) {
			message("wrong gradient");
			status = EXIT_FAILURE;
		}
	}

fail:
	efp_shutdown(efp);
	return res ? EXIT_FAILURE : status;
}
