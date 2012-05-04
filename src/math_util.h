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

#ifndef LIBEFP_MATH_UTIL_H
#define LIBEFP_MATH_UTIL_H

#include <math.h>
#include <string.h>

#include "cblas.h"

#define PI 3.14159265358979323846

struct vec {
	double x, y, z;
};

struct mat {
	double xx, xy, xz, yx, yy, yz, zx, zy, zz;
};

#define VEC(x) ((struct vec *)(&(x)))
#define ELV(v, i) (((double *)v)[i])
#define ELM(m, i, j) (((double *)m)[3 * i + j])

static inline void
vec_zero(struct vec *vec)
{
	memset(vec, 0, sizeof(struct vec));
}

//static inline void
//vec_add(const struct vec *a, const struct vec *b, struct vec *out)
//{
//	out->x = a->x + b->x, out->y = a->y + b->y, out->z = a->z + b->z;
//}

//static inline void
//vec_sub(const struct vec *a, const struct vec *b, struct vec *out)
//{
//	out->x = a->x - b->x, out->y = a->y - b->y, out->z = a->z - b->z;
//}

static inline double
vec_dot(const struct vec *a, const struct vec *b)
{
	return a->x * b->x + a->y * b->y + a->z * b->z;
}

static inline double
vec_len_2(const struct vec *a)
{
	return vec_dot(a, a);
}

static inline double
vec_len(const struct vec *a)
{
	return sqrt(vec_len_2(a));
}

static inline double
vec_dist_2(const struct vec *a, const struct vec *b)
{
	struct vec dr = {
		a->x - b->x, a->y - b->y, a->z - b->z
	};
	return vec_len_2(&dr);
}

static inline double
vec_dist(const struct vec *a, const struct vec *b)
{
	return sqrt(vec_dist_2(a, b));
}

static inline int
eq(double a, double b)
{
	static const double eps = 1.0e-8;
	return fabs(a - b) < eps;
}

static inline void
mat_zero(struct mat *mat)
{
	memset(mat, 0, sizeof(struct mat));
}

static inline void
mat_vec(const struct mat *mat, const struct vec *vec, struct vec *out)
{
	out->x = mat->xx * vec->x + mat->xy * vec->y + mat->xz * vec->z;
	out->y = mat->yx * vec->x + mat->yy * vec->y + mat->yz * vec->z;
	out->z = mat->zx * vec->x + mat->zy * vec->y + mat->zz * vec->z;
}

static inline void
mat_trans_vec(const struct mat *mat, const struct vec *vec, struct vec *out)
{
	out->x = mat->xx * vec->x + mat->yx * vec->y + mat->zx * vec->z;
	out->y = mat->xy * vec->x + mat->yy * vec->y + mat->zy * vec->z;
	out->z = mat->xz * vec->x + mat->yz * vec->y + mat->zz * vec->z;
}

static inline void
move_pt(const struct vec *pos, const struct mat *rotmat,
	const struct vec *init, struct vec *out)
{
	mat_vec(rotmat, init, out);
	out->x += pos->x, out->y += pos->y, out->z += pos->z;
}

static inline void
powers(double x, int n, double *p)
{
	p[0] = 1.0;

	while (--n > 0) {
		p[1] = p[0] * x;
		p++;
	}
}

static inline void
rotate_t1(const struct mat *rotmat, const double *in, double *out)
{
	mat_vec(rotmat, (const struct vec *)in, (struct vec *)out);
}

static inline void
rotate_t2(const struct mat *rotmat, const double *in, double *out)
{
	memset(out, 0, 9 * sizeof(double));

	for (int a1 = 0; a1 < 3; a1++)
	for (int b1 = 0; b1 < 3; b1++)
		for (int a2 = 0; a2 < 3; a2++)
		for (int b2 = 0; b2 < 3; b2++)
			out[a2 * 3 + b2] += in[a1 * 3 + b1] *
				ELM(rotmat, a2, a1) * ELM(rotmat, b2, b1);
}

static inline void
rotate_t3(const struct mat *rotmat, const double *in, double *out)
{
	memset(out, 0, 27 * sizeof(double));

	for (int a1 = 0; a1 < 3; a1++)
	for (int b1 = 0; b1 < 3; b1++)
	for (int c1 = 0; c1 < 3; c1++)
		for (int a2 = 0; a2 < 3; a2++)
		for (int b2 = 0; b2 < 3; b2++)
		for (int c2 = 0; c2 < 3; c2++)
			out[a2 * 9 + b2 * 3 + c2] += in[a1 * 9 + b1 * 3 + c1] *
				ELM(rotmat, a2, a1) *
				ELM(rotmat, b2, b1) *
				ELM(rotmat, c2, c1);
}

#endif /* LIBEFP_MATH_UTIL_H */
