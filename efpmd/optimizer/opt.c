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

#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "opt.h"

#define copy(n, a, b) memcpy(a, b, n * sizeof(double))

struct opt_state {
	size_t n_dim;
	opt_fn fn;
	double ls_tol;
	double ls_fn_tol;
	double ls_step_size;
	size_t ls_max_iter;
	double *x_cur;
	double fx_cur;
	double *gx_cur;
	double *gx_prev;
	double *dir;
	void *user_data;
};

struct opt_state *opt_create(size_t n_dim)
{
	struct opt_state *state = calloc(1, sizeof(struct opt_state));

	if (state) {
		state->ls_tol = 1.0e-8;
		state->ls_fn_tol = 1.0e-8;
		state->ls_step_size = 1.0;
		state->ls_max_iter = 200;
		state->n_dim = n_dim;
		state->x_cur = malloc(state->n_dim * sizeof(double));
		state->gx_cur = malloc(state->n_dim * sizeof(double));
		state->gx_prev = malloc(state->n_dim * sizeof(double));
		state->dir = malloc(state->n_dim * sizeof(double));
	}

	return state;
}

void opt_shutdown(struct opt_state *state)
{
	if (state) {
		free(state->x_cur);
		free(state->gx_cur);
		free(state->gx_prev);
		free(state->dir);
		free(state);
	}
}

void opt_set_fn(struct opt_state *state, opt_fn fn)
{
	if (state)
		state->fn = fn;
}

void opt_set_user_data(struct opt_state *state, void *data)
{
	if (state)
		state->user_data = data;
}

void opt_set_ls_max_iter(struct opt_state *state, size_t max_iter)
{
	if (state)
		state->ls_max_iter = max_iter;
}

void opt_set_ls_tol(struct opt_state *state, double tol)
{
	if (state)
		state->ls_tol = tol;
}

void opt_set_ls_step_size(struct opt_state *state, double step_size)
{
	if (state)
		state->ls_step_size = step_size;
}

void opt_set_ls_fn_tol(struct opt_state *state, double tol)
{
	if (state)
		state->ls_fn_tol = tol;
}

enum opt_result opt_init(struct opt_state *state, const double *x)
{
	enum opt_result res;

	if (!state || !state->fn)
		return OPT_RESULT_ERROR;

	for (size_t i = 0; i < state->n_dim; i++)
		state->x_cur[i] = x[i];

	if ((res = state->fn(state->n_dim, state->x_cur, &state->fx_cur,
					state->gx_cur, state->user_data)))
		return res;

	for (size_t i = 0; i < state->n_dim; i++) {
		state->dir[i] = -state->gx_cur[i];
		state->gx_prev[i] = state->gx_cur[i];
	}

	return OPT_RESULT_SUCCESS;
}

static double dot(size_t n, const double *a, const double *b)
{
	double sum = 0.0;

	while (n--)
		sum += (*a++) * (*b++);

	return sum;
}

static enum opt_result do_line_search(struct opt_state *state)
{
	size_t iter = 0;
	enum opt_result res;

	double a, fa, xa[state->n_dim], ga[state->n_dim], sla;
	double b, fb, xb[state->n_dim], gb[state->n_dim], slb;
	double c, fc, xc[state->n_dim], gc[state->n_dim], slc;

	double init_sl;

	a = 0.0;
	fa = state->fx_cur;
	copy(state->n_dim, xa, state->x_cur);
	copy(state->n_dim, ga, state->gx_cur);
	sla = dot(state->n_dim, state->gx_cur, state->dir);

	if (sla > 0.0) {
		for (size_t i = 0; i < state->n_dim; i++)
			state->dir[i] = -state->gx_cur[i];

		sla = dot(state->n_dim, state->gx_cur, state->dir);
	}

	init_sl = sla;
	c = state->ls_step_size;

	for (size_t i = 0; i < state->n_dim; i++)
		xc[i] = state->x_cur[i] + state->dir[i] * c;

	if ((res = state->fn(state->n_dim, xc, &fc, gc, state->user_data)))
		return res;

	if (fc < fa) {
		copy(state->n_dim, state->x_cur, xc);
		copy(state->n_dim, state->gx_cur, gc);
		state->fx_cur = fc;
		return OPT_RESULT_SUCCESS;
	}

	slc = dot(state->n_dim, gc, state->dir);

	do {
		if (fabs(c - a) < state->ls_tol)
			return OPT_RESULT_ERROR;

		double eta = 3.0 * (fa - fc) / (c - a) + sla + slc;
		double phi = sqrt(eta * eta - sla * slc);

		b = a + (c - a) * (1.0 - (slc + phi - eta) / (slc - sla + 2.0 * phi));

		for (size_t i = 0; i < state->n_dim; i++)
			xb[i] = state->x_cur[i] + state->dir[i] * b;

		if ((res = state->fn(state->n_dim, xb, &fb, gb, state->user_data)))
			return res;

		slb = dot(state->n_dim, gb, state->dir);

		if (fb < state->fx_cur + init_sl * state->ls_fn_tol * b) {
			copy(state->n_dim, state->x_cur, xb);
			copy(state->n_dim, state->gx_cur, gb);
			state->fx_cur = fb;
			return OPT_RESULT_SUCCESS;
		}

		if (slb < 0.0 && fb < fa) {
			fa = fb;
			a = b;
			sla = slb;
			copy(state->n_dim, xa, xb);
			copy(state->n_dim, ga, gb);
		}
		else {
			fc = fb;
			c = b;
			slc = slb;
			copy(state->n_dim, xc, xb);
			copy(state->n_dim, gc, gb);
		}
	} while ((fb > fa || fb > fc) && (++iter < state->ls_max_iter));

	if (iter >= state->ls_max_iter)
		return OPT_RESULT_ERROR;

	if (fa < fc) {
		copy(state->n_dim, state->x_cur, xa);
		copy(state->n_dim, state->gx_cur, ga);
		state->fx_cur = fa;
	}
	else {
		copy(state->n_dim, state->x_cur, xc);
		copy(state->n_dim, state->gx_cur, gc);
		state->fx_cur = fc;
	}

	return OPT_RESULT_SUCCESS;
}

enum opt_result opt_step(struct opt_state *state)
{
	if (!state)
		return OPT_RESULT_ERROR;

	double gx_delta[state->n_dim];

	/* Polak-Ribiere conjugate gradient formula */

	for (size_t i = 0; i < state->n_dim; i++)
		gx_delta[i] = state->gx_cur[i] - state->gx_prev[i];

	double beta = dot(state->n_dim, gx_delta, state->gx_cur) /
			dot(state->n_dim, state->gx_prev, state->gx_prev);

	for (size_t i = 0; i < state->n_dim; i++) {
		state->dir[i] = -state->gx_cur[i] + beta * state->dir[i];
		state->gx_prev[i] = state->gx_cur[i];
	}

	return do_line_search(state);
}

double opt_get_fx(struct opt_state *state)
{
	if (state)
		return state->fx_cur;

	return NAN;
}

enum opt_result opt_get_x(struct opt_state *state, size_t size, double *out)
{
	if (!state || !out || state->n_dim > size)
		return OPT_RESULT_ERROR;

	copy(state->n_dim, out, state->x_cur);
	return OPT_RESULT_SUCCESS;
}

enum opt_result opt_get_gx(struct opt_state *state, size_t size, double *out)
{
	if (!state || !out || state->n_dim > size)
		return OPT_RESULT_ERROR;

	copy(state->n_dim, out, state->gx_cur);
	return OPT_RESULT_SUCCESS;
}
