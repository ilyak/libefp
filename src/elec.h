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

#ifndef LIBEFP_ELEC_H
#define LIBEFP_ELEC_H

static inline int
quad_idx(int a, int b)
{
	/* order in which GAMESS stores quadrupoles */
	enum { xx = 0, yy, zz, xy, xz, yz };

	static const int idx[] = {
		xx, xy, xz, xy, yy, yz, xz, yz, zz
	};

	return idx[a * 3 + b];
}

static inline int
oct_idx(int a, int b, int c)
{
	/* order in which GAMESS stores octupoles */
	enum { xxx = 0, yyy, zzz, xxy, xxz, xyy, yyz, xzz, yzz, xyz };

	static const int idx[] = {
		xxx, xxy, xxz, xxy, xyy, xyz, xxz, xyz, xzz,
		xxy, xyy, xyz, xyy, yyy, yyz, xyz, yyz, yzz,
		xxz, xyz, xzz, xyz, yyz, yzz, xzz, yzz, zzz
	};

	return idx[a * 9 + b * 3 + c];
}

static inline double
quadrupole_sum(const double *quad, const vec_t *dr)
{
	const double *pdr = (const double *)dr;
	double sum = 0.0;

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			sum += quad[quad_idx(a, b)] * pdr[a] * pdr[b];

	return sum;
}

static inline double
octupole_sum(const double *oct, const vec_t *dr)
{
	const double *pdr = (const double *)dr;
	double sum = 0.0;

	for (int a = 0; a < 3; a++)
		for (int b = 0; b < 3; b++)
			for (int c = 0; c < 3; c++)
				sum += oct[oct_idx(a, b, c)] * pdr[a] *
						pdr[b] * pdr[c];

	return sum;
}

void efp_charge_dipole_grad(double q1, const vec_t *d2, const vec_t *dr,
			    vec_t *force, vec_t *add1, vec_t *add2);

void efp_charge_quadrupole_grad(double q1, const double *quad2, const vec_t *dr,
				vec_t *force, vec_t *add1, vec_t *add2);

void efp_charge_octupole_grad(double q1, const double *oct2, const vec_t *dr,
			      vec_t *force, vec_t *add1, vec_t *add2);

void efp_dipole_dipole_grad(const vec_t *d1, const vec_t *d2, const vec_t *dr,
			    vec_t *force, vec_t *add1, vec_t *add2);

void efp_dipole_quadrupole_grad(const vec_t *d1, const double *quad2,
				const vec_t *dr, vec_t *force, vec_t *add1,
				vec_t *add2);

void efp_quadrupole_quadrupole_grad(const double *quad1, const double *quad2,
				    const vec_t *dr, vec_t *force, vec_t *add1,
				    vec_t *add2);

#endif /* LIBEFP_ELEC_H */
