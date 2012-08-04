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

#include "common.h"
#include "opt.h"
#include "sim.h"

static double energy_func(int n, const double *x, double *dfx, void *data)
{
	enum efp_result res;
	struct efp *efp = (struct efp *)data;
	struct efp_energy energy;

	if ((res = efp_set_coordinates(efp, EFP_COORD_TYPE_XYZABC, x)))
		lib_error(res);

	if ((res = efp_compute(efp, 1)))
		lib_error(res);

	if ((res = efp_get_energy(efp, &energy)))
		lib_error(res);

	if ((res = efp_get_gradient(efp, n, dfx)))
		lib_error(res);

	return energy.total;
}

static int check_conv(double rms_grad, double max_grad, double opt_tol)
{
	return fabs(max_grad) < opt_tol && fabs(rms_grad) < opt_tol / 3.0;
}

static void print_status(struct efp *efp, double e_diff, double rms_grad,
			 double max_grad)
{
	print_geometry(efp);
	print_energy(efp);

	printf("ENERGY DIFFERENCE FROM PREVIOUS STEP %20.10lf\n", e_diff);
	printf("        ROOT MEAN SQUARE OF GRADIENT %20.10lf\n", rms_grad);
	printf("          MAXIMUM GRADIENT COMPONENT %20.10lf\n", max_grad);
}

void sim_cg(struct efp *efp, const struct config *config)
{
	int n_frag;
	enum efp_result res;

	if ((res = efp_get_frag_count(efp, &n_frag)))
		lib_error(res);

	int n_grad = 6 * n_frag;
	double coord[n_grad], grad[n_grad];
	double e_old = 0.0;

	for (int step = 1; step <= config->max_steps; step++) {
		double e_new = opt_cg_step(energy_func, n_grad, coord, grad, efp);
		double e_diff = e_new - e_old;
		double rms_grad = 0.0;
		double max_grad = 0.0;

		for (int i = 0; i < n_grad; i++) {
			rms_grad += grad[i] * grad[i];

			if (fabs(grad[i]) > max_grad)
				max_grad = grad[i];
		}

		rms_grad = sqrt(rms_grad / n_grad);

		if (check_conv(rms_grad, max_grad, config->opt_tol)) {
			printf("\nOPTIMIZATION CONVERGED\n\n");
			printf("FINAL ENERGY AND GEOMETRY\n\n");
			print_status(efp, e_diff, rms_grad, max_grad);
			return;
		}

		if (step % config->print_step == 0) {
			printf("ENERGY AND GEOMETRY AFTER %d STEPS\n\n", step);
			print_status(efp, e_diff, rms_grad, max_grad);
		}

		e_old = e_new;
	}
}
