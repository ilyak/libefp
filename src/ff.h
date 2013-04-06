/*-
 * Copyright (c) 2012-2013 Ilya Kaliman
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

#ifndef LIBEFP_FF_H
#define LIBEFP_FF_H

#include "math_util.h"

enum ff_res {
	FF_OK = 0,
	FF_FILE_NOT_FOUND,
	FF_BAD_FORMAT,
	FF_STRING_TOO_LONG,
	FF_NO_PARAMETERS
};

struct ff;

struct ff *efp_ff_create(void);
enum ff_res efp_ff_parse(struct ff *, const char *);
enum ff_res efp_ff_add_atom(struct ff *, const char *);
enum ff_res efp_ff_add_bond(struct ff *, size_t, size_t);
enum ff_res efp_ff_add_angle(struct ff *, size_t, size_t, size_t);
enum ff_res efp_ff_add_torsion(struct ff *, size_t, size_t, size_t, size_t);
enum ff_res efp_ff_auto_angles(struct ff *);
enum ff_res efp_ff_auto_torsions(struct ff *);
size_t efp_ff_get_atom_count(const struct ff *);
void efp_ff_set_atom_pos(struct ff *, size_t, vec_t);
vec_t efp_ff_get_atom_pos(const struct ff *, size_t);
double efp_ff_compute(const struct ff *, vec_t *);
void efp_ff_free(struct ff *);

#endif /* LIBEFP_FF_H */
