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

#include <clapack.h>
#include <util.h>

#include "common.h"

void sim_hess(struct efp *, const struct config *);

static void get_inertia_fact(const double *inertia, const double *rotmat,
					double *inertia_fact)
{
	double fact[3];

	for (int i = 0; i < 3; i++)
		fact[i] = inertia[i] < EPSILON ? 0.0 : 1.0 / sqrt(inertia[i]);

	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			for (int k = 0; k < 3; k++)
				inertia_fact[3 * i + j] = fact[k] *
					rotmat[3 * i + k] * rotmat[3 * j + k];
}

static void get_weight_fact(struct efp *efp, double *mass_fact, mat_t *inertia_fact)
{
	int n_frags;
	check_fail(efp_get_frag_count(efp, &n_frags));

	double xyzabc[6 * n_frags];
	check_fail(efp_get_coordinates(efp, n_frags, xyzabc));

	for (int i = 0; i < n_frags; i++) {
		double mass;
		check_fail(efp_get_frag_mass(efp, i, &mass));

		double inertia[3];
		check_fail(efp_get_frag_inertia(efp, i, inertia));

		double a = xyzabc[6 * i + 3];
		double b = xyzabc[6 * i + 4];
		double c = xyzabc[6 * i + 5];

		mat_t rotmat;
		euler_to_matrix(a, b, c, &rotmat);

		mass_fact[i] = 1.0 / sqrt(mass);
		get_inertia_fact(inertia, (double *)&rotmat, (double *)&inertia_fact[i]);
	}
}

static void w_tr_tr(double fact1, double fact2, int stride, double *hess)
{
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			hess[stride * i + j] *= fact1 * fact2;
}

static void w_tr_rot(double fact1, const mat_t *fact2, int stride, double *hess)
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			double w = 0.0;

			for (int k = 0; k < 3; k++)
				w += mat_get(fact2, i, k) * hess[stride * i + k];

			hess[stride * i + j] = w * fact1;
		}
	}
}

static void w_rot_rot(const mat_t *fact1, const mat_t *fact2, int stride, double *hess)
{
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			double w1 = 0.0;

			for (int ii = 0; ii < 3; ii++) {
				double w2 = 0.0;

				for (int jj = 0; jj < 3; jj++)
					w2 += mat_get(fact2, j, jj) * hess[stride * i + jj];

				w1 += w2 * mat_get(fact1, i, ii);
			}

			hess[i * stride + j] = w1;
		}
	}
}

static void hess_weight(int n_frags, const double *mass_fact,
			const mat_t *inertia_fact, double *hess)
{
	int stride = 6 * n_frags;

	for (int i = 0; i < n_frags; i++) {
		for (int j = 0; j < n_frags; j++) {
			double *block = hess + 6 * stride * i + 6 * j;

			w_tr_tr(mass_fact[i], mass_fact[j], stride,
					block);

			w_tr_rot(mass_fact[i], inertia_fact + j, stride,
					block + 3);

			w_tr_rot(mass_fact[j], inertia_fact + i, stride,
					block + 3 * stride);

			w_rot_rot(inertia_fact + i, inertia_fact + j, stride,
					block + 3 * stride + 3);
		}
	}
}

static double to_rcm(double eval)
{
	double freq = copysign(sqrt(fabs(eval)), eval);

	return AMU_TO_AU * FINE_CONST / 2.0 / PI / BOHR_RADIUS * 1.0e8 * freq;
}

static void print_freq(int n_coord, const double *eval, const double *evec)
{
	for (int i = 0; i < n_coord; i++) {
		printf("    MODE %4d    FREQUENCY %10.3lf cm-1\n\n", i + 1, to_rcm(eval[i]));
		print_vector(n_coord, evec + i * n_coord);
	}
}

void sim_hess(struct efp *efp, const struct config *config)
{
	printf("HESSIAN JOB\n\n\n");

	int n_frags = config->n_frags;
	int n_coord = 6 * n_frags;

	double *xyzabc = xmalloc(n_coord * sizeof(double));
	double *grad = xmalloc(n_coord * sizeof(double));
	double *grad_0 = xmalloc(n_coord * sizeof(double));
	double *hess = xmalloc(n_coord * n_coord * sizeof(double));

	print_geometry(efp);
	check_fail(efp_compute(efp, 1));
	print_energy(efp);
	print_gradient(efp);

	check_fail(efp_get_coordinates(efp, n_frags, xyzabc));
	check_fail(efp_get_gradient(efp, n_frags, grad_0));
	torque_to_deriv(n_frags, xyzabc, grad_0);

	for (int i = 0; i < n_coord; i++) {
		printf("COMPUTING DISPLACEMENT %5d OF %d\n", i + 1, n_coord);
		fflush(stdout);

		double save = xyzabc[i];
		xyzabc[i] = save + config->hess_delta;

		check_fail(efp_set_coordinates(efp, EFP_COORD_TYPE_XYZABC, xyzabc));
		check_fail(efp_compute(efp, 1));

		check_fail(efp_get_gradient(efp, n_frags, grad));
		torque_to_deriv(n_frags, xyzabc, grad);

		for (int j = 0; j < n_coord; j++)
			hess[i * n_coord + j] = (grad[j] - grad_0[j]) / config->hess_delta;

		xyzabc[i] = save;
	}

	/* restore original coordinates */
	check_fail(efp_set_coordinates(efp, EFP_COORD_TYPE_XYZABC, xyzabc));

	/* reduce errors by computing the average of H(i,j) and H(j,i) */
	for (int i = 0; i < n_coord; i++) {
		for (int j = i + 1; j < n_coord; j++) {
			double sum = hess[i * n_coord + j] + hess[j * n_coord + i];

			hess[i * n_coord + j] = sum / 2.0;
			hess[j * n_coord + i] = sum / 2.0;
		}
	}

	printf("\n\n    HESSIAN MATRIX\n\n");
	print_matrix(n_coord, n_coord, hess);

	printf("    VIBRATIONAL NORMAL MODE ANALYSIS\n\n");

	double *eigen = xmalloc(n_coord * sizeof(double));
	double *work = xmalloc(4 * n_coord * sizeof(double));

	double mass_fact[n_frags];
	mat_t inertia_fact[n_frags];

	get_weight_fact(efp, mass_fact, inertia_fact);
	hess_weight(n_frags, mass_fact, inertia_fact, hess);

	if (u_dsyev('V', 'U', n_coord, hess, n_coord, eigen, work, 4 * n_coord))
		error("UNABLE TO DIAGONALIZE MASS-WEIGHTED HESSIAN MATRIX");

	print_freq(n_coord, eigen, hess);

	free(eigen);
	free(work);
	free(xyzabc);
	free(grad);
	free(grad_0);
	free(hess);

	printf("HESSIAN JOB COMPLETED SUCCESSFULLY\n");
}
