/*-
 * Copyright (c) 2012-2017 Ilya Kaliman
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

#include <ctype.h>

#include "private.h"
#include "util.h"

int
efp_skip_frag_pair(const struct efp *efp, size_t fr_i_idx, size_t fr_j_idx)
{
	size_t idx = fr_i_idx * efp->n_frag + fr_j_idx;

	if (efp->skiplist[idx])
		return 1;
	if (!efp->opts.enable_cutoff)
		return 0;

	const struct frag *fr_i = efp->frags + fr_i_idx;
	const struct frag *fr_j = efp->frags + fr_j_idx;

	double cutoff2 = efp->opts.swf_cutoff * efp->opts.swf_cutoff;
	vec_t dr = vec_sub(CVEC(fr_j->x), CVEC(fr_i->x));

	if (efp->opts.enable_pbc) {
	    // orthogonal box
	    if (efp->box.a == 90.0 && efp->box.b == 90.0 && efp->box.c == 90.0) {
            vec_t cell = {efp->box.x * round(dr.x / efp->box.x),
                          efp->box.y * round(dr.y / efp->box.y),
                          efp->box.z * round(dr.z / efp->box.z)};
            dr = vec_sub(&dr, &cell);
        }
	    else {
	        cart_to_frac(efp->box, &dr);
	        vec_t cell = {round(dr.x), round(dr.y), round(dr.z)};
	        dr = vec_sub(&dr, &cell);
	        frac_to_cart(efp->box, &dr);
	    }
	}
	return vec_len_2(&dr) > cutoff2;
}

struct swf
efp_make_swf(const struct efp *efp, const struct frag *fr_i,
    const struct frag *fr_j, int short_range_cutoff)
{
	struct swf swf;

	memset(&swf, 0, sizeof(swf));
	swf.swf = 1.0;
	swf.dr = vec_sub(CVEC(fr_j->x), CVEC(fr_i->x));

	if (!efp->opts.enable_cutoff)
		return swf;
	if (efp->opts.enable_pbc) {
        if (efp->box.a == 90.0 && efp->box.b == 90.0 && efp->box.c == 90.0) {
            swf.cell.x = efp->box.x * round(swf.dr.x / efp->box.x);
            swf.cell.y = efp->box.y * round(swf.dr.y / efp->box.y);
            swf.cell.z = efp->box.z * round(swf.dr.z / efp->box.z);
            swf.dr.x -= swf.cell.x;
            swf.dr.y -= swf.cell.y;
            swf.dr.z -= swf.cell.z;
        }
        else {
            cart_to_frac(efp->box, &swf.dr);
            swf.cell.x = round(swf.dr.x);
            swf.cell.y = round(swf.dr.y);
            swf.cell.z = round(swf.dr.z);
            swf.dr = vec_sub(&swf.dr, &swf.cell);
            frac_to_cart(efp->box, &swf.cell);
            frac_to_cart(efp->box, &swf.dr);
        }
	}

	double r = vec_len(&swf.dr);

	// differentiating between cutoff for exrep and other terms
	double dswf = 0.0;
	if (short_range_cutoff == 1) {
        swf.swf = efp_get_swf(r, efp->opts.xr_cutoff);
        dswf = efp_get_dswf(r, efp->opts.xr_cutoff);
	}
	else {
        swf.swf = efp_get_swf(r, efp->opts.swf_cutoff);
        dswf = efp_get_dswf(r, efp->opts.swf_cutoff);
    }
	swf.dswf.x = -dswf * swf.dr.x;
	swf.dswf.y = -dswf * swf.dr.y;
	swf.dswf.z = -dswf * swf.dr.z;

	return swf;
}

int
efp_check_rotation_matrix(const mat_t *rotmat)
{
	vec_t ax = { rotmat->xx, rotmat->yx, rotmat->zx };
	vec_t ay = { rotmat->xy, rotmat->yy, rotmat->zy };
	vec_t az = { rotmat->xz, rotmat->yz, rotmat->zz };

	if (!eq(vec_len(&ax), 1.0) ||
	    !eq(vec_len(&ay), 1.0) ||
	    !eq(vec_len(&az), 1.0))
		return 0;
	if (!eq(vec_dot(&ax, &ay), 0.0))
		return 0;

	vec_t cross = vec_cross(&ax, &ay);

	if (!eq(cross.x, az.x) ||
	    !eq(cross.y, az.y) ||
	    !eq(cross.z, az.z))
		return 0;

	return 1;
}

void
efp_points_to_matrix(const double *pts, mat_t *rotmat)
{
	vec_t p1 = { pts[0], pts[1], pts[2] };
	vec_t p2 = { pts[3], pts[4], pts[5] };
	vec_t p3 = { pts[6], pts[7], pts[8] };

	vec_t r12 = vec_sub(&p2, &p1);
	vec_t r13 = vec_sub(&p3, &p1);

	vec_normalize(&r12);
	vec_normalize(&r13);

	double dot = vec_dot(&r12, &r13);

	r13.x -= dot * r12.x;
	r13.y -= dot * r12.y;
	r13.z -= dot * r12.z;

	vec_t cross = vec_cross(&r12, &r13);

	vec_normalize(&r13);
	vec_normalize(&cross);

	rotmat->xx = r12.x;
	rotmat->yx = r12.y;
	rotmat->zx = r12.z;
	rotmat->xy = r13.x;
	rotmat->yy = r13.y;
	rotmat->zy = r13.z;
	rotmat->xz = cross.x;
	rotmat->yz = cross.y;
	rotmat->zz = cross.z;
}

const struct frag *
efp_find_lib(struct efp *efp, const char *name)
{
	for (size_t i = 0; i < efp->n_lib; i++)
		if (efp_strcasecmp(efp->lib[i]->name, name) == 0)
			return efp->lib[i];
	return NULL;
}

void
efp_add_stress(const vec_t *dr, const vec_t *force, mat_t *stress)
{
#ifdef _OPENMP
#pragma omp critical
#endif
	{
		stress->xx += dr->x * force->x;
		stress->xy += dr->x * force->y;
		stress->xz += dr->x * force->z;
		stress->yx += dr->y * force->x;
		stress->yy += dr->y * force->y;
		stress->yz += dr->y * force->z;
		stress->zx += dr->z * force->x;
		stress->zy += dr->z * force->y;
		stress->zz += dr->z * force->z;
	}
}

void
efp_add_force(six_t *grad, const vec_t *com, const vec_t *pt,
    const vec_t *force, const vec_t *add)
{
	vec_t dr = vec_sub(CVEC(pt->x), com);
	vec_t torque = vec_cross(&dr, force);

	if (add) {
		torque.x += add->x;
		torque.y += add->y;
		torque.z += add->z;
	}
	six_atomic_add_xyz(grad, force);
	six_atomic_add_abc(grad, &torque);
}

void
efp_sub_force(six_t *grad, const vec_t *com, const vec_t *pt,
    const vec_t *force, const vec_t *add)
{
	vec_t dr = vec_sub(CVEC(pt->x), com);
	vec_t torque = vec_cross(&dr, force);

	if (add) {
		torque.x += add->x;
		torque.y += add->y;
		torque.z += add->z;
	}
	six_atomic_sub_xyz(grad, force);
	six_atomic_sub_abc(grad, &torque);
}

void
efp_move_pt(const vec_t *com, const mat_t *rotmat,
    const vec_t *pos_int, vec_t *out)
{
	*out = mat_vec(rotmat, pos_int);
	out->x += com->x;
	out->y += com->y;
	out->z += com->z;
}

void
efp_rotate_t2(const mat_t *rotmat, const double *in, double *out)
{
	for (size_t i = 0; i < 3 * 3; i++)
		out[i] = 0.0;

	for (size_t a1 = 0; a1 < 3; a1++)
	for (size_t b1 = 0; b1 < 3; b1++)
		for (size_t a2 = 0; a2 < 3; a2++)
		for (size_t b2 = 0; b2 < 3; b2++)
			out[a2 * 3 + b2] += in[a1 * 3 + b1] *
			    mat_get(rotmat, a2, a1) *
			    mat_get(rotmat, b2, b1);
}

void
efp_rotate_t3(const mat_t *rotmat, const double *in, double *out)
{
	double a, b, c;

	for (size_t i = 0; i < 3 * 3 * 3; i++)
		out[i] = 0.0;

	for (size_t a1 = 0; a1 < 3; a1++)
	for (size_t b1 = 0; b1 < 3; b1++)
	for (size_t c1 = 0; c1 < 3; c1++)
		for (size_t a2 = 0; a2 < 3; a2++)
		for (size_t b2 = 0; b2 < 3; b2++)
		for (size_t c2 = 0; c2 < 3; c2++) {
			a = mat_get(rotmat, a2, a1);
			b = mat_get(rotmat, b2, b1);
			c = mat_get(rotmat, c2, c1);
			out[a2 * 9 + b2 * 3 + c2] +=
			    in[a1 * 9 + b1 * 3 + c1] * a * b * c;
		}
}

mat_t
rotmat_2frags(const mat_t *rotmat1, const mat_t *rotmat2)
{
    // rotmat2*rotmat1^T
    mat_t rotmat1_t = mat_transpose(rotmat1);
    return mat_mat(rotmat2,rotmat1_t);
}

int
efp_strcasecmp(const char *s1, const char *s2)
{
	while (tolower(*s1) == tolower(*s2++))
		if (*s1++ == '\0')
			return 0;
	return tolower(*s1) - tolower(*--s2);
}

int
efp_strncasecmp(const char *s1, const char *s2, size_t n)
{
	if (n != 0) {
		do {
			if (tolower(*s1) != tolower(*s2++))
				return tolower(*s1) - tolower(*--s2);
			if (*s1++ == '\0')
				break;
		} while (--n != 0);
	}
	return 0;
}

void
find_plane(const vec_t pt1, const vec_t pt2, const vec_t pt3, vec_t *normal, double d) {
    // determines plane by three points: ax + by + cz + d = 0, where a,b,c are given by a normal vector
    vec_t vec21 = vec_sub(&pt2,&pt1);
    vec_t vec31 = vec_sub(&pt3,&pt1);
    *normal = vec_cross(&vec21,&vec31);
    d = - vec_dot(normal,&pt1);
}

double
max_cutoff(const six_t box) {
    vec_t point000 = {0.0, 0.0, 0.0};
    vec_t point100 = {1.0, 0.0, 0.0};
    vec_t point010 = {0.0, 1.0, 0.0};
    vec_t point001 = {0.0, 0.0, 1.0};
    vec_t point = {0.5,0.5,0.5};
    frac_to_cart(box, &point100);
    frac_to_cart(box, &point010);
    frac_to_cart(box, &point001);
    frac_to_cart(box, &point);
    // plane xy: points 000, 100, 010
    vec_t normal;
    double d;
    double dist;
    double min_dist;
    // plane xy: points 000, 100, 010
    find_plane(point100, point000, point010, &normal, d);
    dist = (vec_dot(&normal, &point) + d) / vec_len(&normal);
    dist = sqrt(dist*dist);
    min_dist = dist;
    // plane xz: points 000, 100, 001
    find_plane(point100, point000, point001, &normal, d);
    dist = (vec_dot(&normal, &point) + d) / vec_len(&normal);
    dist = sqrt(dist*dist);
    if (min_dist > dist) min_dist = dist;
    // plane yz: points 000, 010, 001
    find_plane(point010, point000, point001, &normal, d);
    dist = (vec_dot(&normal, &point) + d) / vec_len(&normal);
    dist = sqrt(dist*dist);
    if (min_dist > dist) min_dist = dist;
    return min_dist;
}