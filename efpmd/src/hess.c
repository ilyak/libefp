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

#include <util.h>

#include "common.h"

void sim_hess(struct efp *, const struct config *);

void sim_hess(struct efp *efp, const struct config *config)
{
	enum efp_result res;
	int n_frags = config->n_frags;
	int n_coord = 6 * n_frags;

	double *xyzabc = xmalloc(n_coord * sizeof(double));
	double *grad = xmalloc(n_coord * sizeof(double));
	double *grad_0 = xmalloc(n_coord * sizeof(double));
	double *hess = xmalloc(n_coord * n_coord * sizeof(double));

	if ((res = efp_compute(efp, 1)))
		lib_error(res);

	if ((res = efp_get_coordinates(efp, n_frags, xyzabc)))
		lib_error(res);

	if ((res = efp_get_gradient(efp, n_frags, grad_0)))
		lib_error(res);

	torque_to_deriv(n_frags, xyzabc, grad_0);

	print_geometry(efp);
	print_energy(efp);
	print_gradient(efp);

	for (int i = 0; i < n_coord; i++) {
		printf("COMPUTING DISPLACEMENT %5d OF %d\n", i + 1, n_coord);
		fflush(stdout);

		double save = xyzabc[i];
		xyzabc[i] = save + config->hess_delta;

		if ((res = efp_set_coordinates(efp, EFP_COORD_TYPE_XYZABC, xyzabc)))
			lib_error(res);

		if ((res = efp_compute(efp, 1)))
			lib_error(res);

		if ((res = efp_get_gradient(efp, n_frags, grad)))
			lib_error(res);

		torque_to_deriv(n_frags, xyzabc, grad);

		for (int j = 0; j < n_coord; j++)
			hess[i * n_coord + j] = (grad[j] - grad_0[j]) / config->hess_delta;

		xyzabc[i] = save;
	}

	printf("\n\n    HESSIAN MATRIX\n\n");
	print_matrix(n_coord, n_coord, hess);

	free(xyzabc);
	free(grad);
	free(grad_0);
	free(hess);
}
