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

#include "efp_private.h"

int skip_frag_pair(struct efp *efp, int fr_i_idx, int fr_j_idx)
{
	if (!efp->opts.enable_cutoff)
		return 0;

	const struct frag *fr_i = efp->frags + fr_i_idx;
	const struct frag *fr_j = efp->frags + fr_j_idx;

	double cutoff2 = efp->opts.swf_cutoff * efp->opts.swf_cutoff;
	vec_t dr = vec_sub(CVEC(fr_j->x), CVEC(fr_i->x));

	if (efp->opts.enable_pbc) {
		vec_t cell = { efp->box.x * round(dr.x / efp->box.x),
			       efp->box.y * round(dr.y / efp->box.y),
			       efp->box.z * round(dr.z / efp->box.z) };

		dr = vec_sub(&dr, &cell);
	}

	return vec_len_2(&dr) > cutoff2;
}

struct swf make_swf(struct efp *efp, const struct frag *fr_i, const struct frag *fr_j)
{
	struct swf swf = {
		.swf = 1.0,
		.dswf = vec_zero,
		.cell = vec_zero
	};

	if (!efp->opts.enable_cutoff)
		return swf;

	vec_t dr = vec_sub(CVEC(fr_j->x), CVEC(fr_i->x));

	if (efp->opts.enable_pbc) {
		swf.cell.x = efp->box.x * round(dr.x / efp->box.x);
		swf.cell.y = efp->box.y * round(dr.y / efp->box.y);
		swf.cell.z = efp->box.z * round(dr.z / efp->box.z);

		dr.x -= swf.cell.x;
		dr.y -= swf.cell.y;
		dr.z -= swf.cell.z;
	}

	double r = vec_len(&dr);

	swf.swf = get_swf(r, efp->opts.swf_cutoff);
	double dswf = get_dswf(r, efp->opts.swf_cutoff);

	swf.dswf.x = -dswf * dr.x;
	swf.dswf.y = -dswf * dr.y;
	swf.dswf.z = -dswf * dr.z;

	return swf;
}

void add_stress(const vec_t *dr, const vec_t *force, mat_t *stress)
{
	#pragma omp critical
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

void add_force(struct frag *frag, const vec_t *pt, const vec_t *force, const vec_t *add)
{
	vec_t dr = vec_sub(CVEC(pt->x), CVEC(frag->x));
	vec_t torque = vec_cross(&dr, force);

	if (add) {
		torque.x += add->x;
		torque.y += add->y;
		torque.z += add->z;
	}

	vec_atomic_add(&frag->force, force);
	vec_atomic_add(&frag->torque, &torque);
}

void sub_force(struct frag *frag, const vec_t *pt, const vec_t *force, const vec_t *add)
{
	vec_t dr = vec_sub(CVEC(pt->x), CVEC(frag->x));
	vec_t torque = vec_cross(&dr, force);

	if (add) {
		torque.x += add->x;
		torque.y += add->y;
		torque.z += add->z;
	}

	vec_atomic_sub(&frag->force, force);
	vec_atomic_sub(&frag->torque, &torque);
}
