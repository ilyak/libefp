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

#include <assert.h>
#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "stream.h"
#include "ff.h"

#define DEG2RAD (PI / 180.0)
#define KCALMOL2AU (1.0 / 627.50947)
#define ANG2BOHR (1.0 / 0.52917721092)

struct ff {
	size_t n_atoms;

	struct atom {
		char type[16];
		vec_t pos;
	} *atoms;

	struct bond_fp {
		char type1[16];
		char type2[16];
		double force;
		double ideal;
		struct bond_fp *next;
	} *bond_fps;

	struct bond {
		size_t idx1;
		size_t idx2;
		struct bond_fp *fp;
		struct bond *next;
	} *bonds;

	struct angle_fp {
		char type1[16];
		char type2[16];
		char type3[16];
		double force;
		double ideal;
		struct angle_fp *next;
	} *angle_fps;

	struct angle {
		size_t idx1;
		size_t idx2;
		size_t idx3;
		struct angle_fp *fp;
		struct angle *next;
	} *angles;

	struct torsion_fp {
		char type1[16];
		char type2[16];
		char type3[16];
		char type4[16];
		double vi[6];
		double ci[6];
		double si[6];
		struct torsion_fp *next;
	} *torsion_fps;

	struct torsion {
		size_t idx1;
		size_t idx2;
		size_t idx3;
		size_t idx4;
		struct torsion_fp *fp;
		struct torsion *next;
	} *torsions;
};

static const char *str_skip(const char *str, size_t cnt)
{
	if (cnt == 0 || *str == '\0')
		return str;

	while (*++str) {
		if (isspace(*str) && !isspace(str[-1]))
			if (--cnt == 0)
				break;
	}

	while (*str && isspace(*str))
		str++;

	return str;
}

static struct bond_fp *find_bond_fp(struct bond_fp *fps,
		const char *type1, const char *type2)
{
	while (fps) {
		if ((strcmp(type1, fps->type1) == 0 &&
		     strcmp(type2, fps->type2) == 0) ||
		    (strcmp(type2, fps->type1) == 0 &&
		     strcmp(type1, fps->type2) == 0))
			return fps;

		fps = fps->next;
	}

	return NULL;
}

static struct angle_fp *find_angle_fp(struct angle_fp *fps,
		const char *type1, const char *type2, const char *type3)
{
	while (fps) {
		if ((strcmp(type1, fps->type1) == 0 &&
		     strcmp(type2, fps->type2) == 0 &&
		     strcmp(type3, fps->type3) == 0) ||
		    (strcmp(type3, fps->type1) == 0 &&
		     strcmp(type2, fps->type2) == 0 &&
		     strcmp(type1, fps->type3) == 0))
			return fps;

		fps = fps->next;
	}

	return NULL;
}

static struct torsion_fp *find_torsion_fp(struct torsion_fp *fps,
		const char *type1, const char *type2,
		const char *type3, const char *type4)
{
	bool is_a, is_b, is_c, is_d;

	while (fps) {
		is_a = strcmp(fps->type1, "X") == 0 || strcmp(type1, fps->type1) == 0;
		is_b = strcmp(fps->type2, "X") == 0 || strcmp(type2, fps->type2) == 0;
		is_c = strcmp(fps->type3, "X") == 0 || strcmp(type3, fps->type3) == 0;
		is_d = strcmp(fps->type4, "X") == 0 || strcmp(type4, fps->type4) == 0;

		if (is_a && is_b && is_c && is_d)
			return fps;

		is_a = strcmp(fps->type1, "X") == 0 || strcmp(type4, fps->type1) == 0;
		is_b = strcmp(fps->type2, "X") == 0 || strcmp(type3, fps->type2) == 0;
		is_c = strcmp(fps->type3, "X") == 0 || strcmp(type2, fps->type3) == 0;
		is_d = strcmp(fps->type4, "X") == 0 || strcmp(type1, fps->type4) == 0;

		if (is_a && is_b && is_c && is_d)
			return fps;

		fps = fps->next;
	}

	return NULL;
}

static double bond_energy(const struct bond *bond, const struct atom *atoms,
		vec_t *grad)
{
	const struct atom *at1 = atoms + bond->idx1;
	const struct atom *at2 = atoms + bond->idx2;
	const struct bond_fp *fp = bond->fp;

	vec_t dr = vec_sub(&at1->pos, &at2->pos);

	double r = vec_len(&dr);
	double dt = r - fp->ideal;

	double energy = fp->force * dt * dt;

	if (grad) {
		double dedt = 2.0 * fp->force * dt;

		double dedx1 = dedt * dr.x / r;
		double dedy1 = dedt * dr.y / r;
		double dedz1 = dedt * dr.z / r;

		grad[bond->idx1].x += dedx1;
		grad[bond->idx1].y += dedy1;
		grad[bond->idx1].z += dedz1;
		grad[bond->idx2].x -= dedx1;
		grad[bond->idx2].y -= dedy1;
		grad[bond->idx2].z -= dedz1;
	}

	return energy;
}

static double angle_energy(const struct angle *angle, const struct atom *atoms,
		vec_t *grad)
{
	const struct atom *at1 = atoms + angle->idx1;
	const struct atom *at2 = atoms + angle->idx2;
	const struct atom *at3 = atoms + angle->idx3;
	const struct angle_fp *fp = angle->fp;

	vec_t dr12 = vec_sub(&at1->pos, &at2->pos);
	vec_t dr32 = vec_sub(&at3->pos, &at2->pos);

	double r12 = vec_len_2(&dr12);
	double r32 = vec_len_2(&dr32);
	double cosa = vec_dot(&dr12, &dr32) / sqrt(r12 * r32);
	double dt = acos(cosa) - fp->ideal;

	double energy = fp->force * dt * dt;

	if (grad) {
		vec_t p = vec_cross(&dr32, &dr12);

		double rp = vec_len(&p);
		double dedt = 2.0 * fp->force * dt;
		double ta = -dedt / (r12 * rp);
		double tc = dedt / (r32 * rp);

		double dedx1 = ta * (dr12.y * p.z - dr12.z * p.y);
		double dedy1 = ta * (dr12.z * p.x - dr12.x * p.z);
		double dedz1 = ta * (dr12.x * p.y - dr12.y * p.x);
		double dedx3 = tc * (dr32.y * p.z - dr32.z * p.y);
		double dedy3 = tc * (dr32.z * p.x - dr32.x * p.z);
		double dedz3 = tc * (dr32.x * p.y - dr32.y * p.x);

		grad[angle->idx1].x += dedx1;
		grad[angle->idx1].y += dedy1;
		grad[angle->idx1].z += dedz1;
		grad[angle->idx2].x -= dedx1 + dedx3;
		grad[angle->idx2].y -= dedy1 + dedy3;
		grad[angle->idx2].z -= dedz1 + dedz3;
		grad[angle->idx3].x += dedx3;
		grad[angle->idx3].y += dedy3;
		grad[angle->idx3].z += dedz3;
	}

	return energy;
}

static double torsion_energy(const struct torsion *torsion, const struct atom *atoms,
		vec_t *grad)
{
	size_t idx1 = torsion->idx1;
	size_t idx2 = torsion->idx2;
	size_t idx3 = torsion->idx3;
	size_t idx4 = torsion->idx4;
	const struct atom *at1 = atoms + idx1;
	const struct atom *at2 = atoms + idx2;
	const struct atom *at3 = atoms + idx3;
	const struct atom *at4 = atoms + idx4;
	const struct torsion_fp *fp = torsion->fp;

	vec_t dr21 = vec_sub(&at2->pos, &at1->pos);
	vec_t dr32 = vec_sub(&at3->pos, &at2->pos);
	vec_t dr43 = vec_sub(&at4->pos, &at3->pos);
	vec_t t = vec_cross(&dr21, &dr32);
	vec_t u = vec_cross(&dr32, &dr43);
	vec_t tu = vec_cross(&t, &u);

	double rt = vec_len(&t);
	double ru = vec_len(&u);
	double r32 = vec_len(&dr32);

	double cosine[6], sine[6];
	double energy = 0.0, dedphi = 0.0;

	cosine[0] = vec_dot(&t, &u) / (rt * ru);
	sine[0] = vec_dot(&tu, &dr32) / (r32 * rt * ru);

	for (size_t i = 1; i < 6; i++) {
		cosine[i] = cosine[0] * cosine[i - 1] - sine[0] * sine[i - 1];
		sine[i] = cosine[0] * sine[i - 1] + sine[0] * cosine[i - 1];
	}

	for (size_t i = 0; i < 6; i++) {
		energy += fp->vi[i] * (1.0 + (cosine[i] * fp->ci[i] +
				sine[i] * fp->si[i]));
		dedphi += fp->vi[i] * (1.0 + i) * (cosine[i] * fp->si[i] -
				sine[i] * fp->ci[i]);
	}

	if (grad) {
		vec_t dr31 = vec_sub(&at3->pos, &at1->pos);
		vec_t dr42 = vec_sub(&at4->pos, &at2->pos);
		vec_t dedt = vec_cross(&t, &dr32);
		vec_t dedu = vec_cross(&u, &dr32);

		vec_scale(&dedt,  dedphi / (r32 * rt * rt));
		vec_scale(&dedu, -dedphi / (r32 * ru * ru));

		grad[idx1].x += dr32.z * dedt.y - dr32.y * dedt.z;
		grad[idx1].y += dr32.x * dedt.z - dr32.z * dedt.x;
		grad[idx1].z += dr32.y * dedt.x - dr32.x * dedt.y;
		grad[idx2].x += dr31.y * dedt.z - dr31.z * dedt.y;
		grad[idx2].y += dr31.z * dedt.x - dr31.x * dedt.z;
		grad[idx2].z += dr31.x * dedt.y - dr31.y * dedt.x;
		grad[idx2].x += dr43.z * dedu.y - dr43.y * dedu.z;
		grad[idx2].y += dr43.x * dedu.z - dr43.z * dedu.x;
		grad[idx2].z += dr43.y * dedu.x - dr43.x * dedu.y;
		grad[idx3].x += dr21.z * dedt.y - dr21.y * dedt.z;
		grad[idx3].y += dr21.x * dedt.z - dr21.z * dedt.x;
		grad[idx3].z += dr21.y * dedt.x - dr21.x * dedt.y;
		grad[idx3].x += dr42.y * dedu.z - dr42.z * dedu.y;
		grad[idx3].y += dr42.z * dedu.x - dr42.x * dedu.z;
		grad[idx3].z += dr42.x * dedu.y - dr42.y * dedu.x;
		grad[idx4].x += dr32.z * dedu.y - dr32.y * dedu.z;
		grad[idx4].y += dr32.x * dedu.z - dr32.z * dedu.x;
		grad[idx4].z += dr32.y * dedu.x - dr32.x * dedu.y;
	}

	return energy;
}

static enum ff_res parse_bond_fp(struct ff *ff, const char *str)
{
	struct bond_fp *fp = malloc(sizeof(struct bond_fp));

	int nsc = sscanf(str, "bond %10s %10s %lf %lf", fp->type1, fp->type2,
				&fp->force, &fp->ideal);

	if (nsc < 4) {
		free(fp);
		return FF_BAD_FORMAT;
	}

	fp->next = ff->bond_fps;
	ff->bond_fps = fp;
	return FF_OK;
}

static enum ff_res parse_angle_fp(struct ff *ff, const char *str)
{
	struct angle_fp *fp = malloc(sizeof(struct angle_fp));

	int nsc = sscanf(str, "angle %10s %10s %10s %lf %lf", fp->type1, fp->type2,
				fp->type3, &fp->force, &fp->ideal);

	if (nsc < 5) {
		free(fp);
		return FF_BAD_FORMAT;
	}

	fp->ideal = DEG2RAD * fp->ideal;
	fp->next = ff->angle_fps;
	ff->angle_fps = fp;
	return FF_OK;
}

static enum ff_res parse_torsion_fp(struct ff *ff, const char *str)
{
	struct torsion_fp *fp = malloc(sizeof(struct torsion_fp));

	int nsc = sscanf(str, "torsion %10s %10s %10s %10s", fp->type1, fp->type2,
				fp->type3, fp->type4);

	if (nsc < 4) {
		free(fp);
		return FF_BAD_FORMAT;
	}

	for (size_t i = 0; i < 6; i++) {
		fp->vi[i] = 0.0;
		fp->ci[i] = 1.0;
		fp->si[i] = 0.0;
	}

	const char *ptr = str_skip(str, 5);

	while (*ptr) {
		int mult;
		double vi, phi;

		nsc = sscanf(ptr, "%lf %lf %d", &vi, &phi, &mult);

		if (nsc < 3 || mult < 1 || mult > 6) {
			free(fp);
			return FF_BAD_FORMAT;
		}

		phi = DEG2RAD * phi;
		mult = mult - 1;

		fp->vi[mult] = vi;
		fp->ci[mult] = cos(phi);
		fp->si[mult] = sin(phi);

		ptr = str_skip(ptr, 3);
	}

	fp->next = ff->torsion_fps;
	ff->torsion_fps = fp;
	return FF_OK;
}

static bool check_angle(const struct bond *b1, const struct bond *b2,
		size_t *idx1, size_t *idx2, size_t *idx3)
{
	if (b1->idx1 == b2->idx1) {
		*idx1 = b1->idx2;
		*idx2 = b1->idx1;
		*idx3 = b2->idx2;
		return true;
	}

	if (b1->idx1 == b2->idx2) {
		*idx1 = b1->idx2;
		*idx2 = b1->idx1;
		*idx3 = b2->idx1;
		return true;
	}

	if (b1->idx2 == b2->idx1) {
		*idx1 = b1->idx1;
		*idx2 = b1->idx2;
		*idx3 = b2->idx2;
		return true;
	}

	if (b1->idx2 == b2->idx2) {
		*idx1 = b1->idx1;
		*idx2 = b1->idx2;
		*idx3 = b2->idx1;
		return true;
	}

	return false;
}

static bool check_torsion(const struct angle *a1, const struct angle *a2,
		size_t *idx1, size_t *idx2, size_t *idx3, size_t *idx4)
{
	if (a1->idx2 == a2->idx1 && a1->idx3 == a2->idx2) {
		*idx1 = a1->idx1;
		*idx2 = a1->idx2;
		*idx3 = a1->idx3;
		*idx4 = a2->idx3;
		return true;
	}

	if (a1->idx2 == a2->idx1 && a1->idx1 == a2->idx2) {
		*idx1 = a1->idx3;
		*idx2 = a1->idx2;
		*idx3 = a1->idx1;
		*idx4 = a2->idx3;
		return true;
	}

	if (a1->idx2 == a2->idx3 && a1->idx3 == a2->idx2) {
		*idx1 = a1->idx1;
		*idx2 = a1->idx2;
		*idx3 = a1->idx3;
		*idx4 = a2->idx1;
		return true;
	}

	if (a1->idx2 == a2->idx3 && a1->idx1 == a2->idx2) {
		*idx1 = a1->idx3;
		*idx2 = a1->idx2;
		*idx3 = a1->idx1;
		*idx4 = a2->idx1;
		return true;
	}

	return false;
}

struct ff *efp_ff_create(void)
{
	return calloc(1, sizeof(struct ff));
}

#define REVERSE_LIST(list, type)                  \
	do {                                      \
		type *prev = NULL;                \
                                                  \
		while (list) {                    \
			type *next = list->next;  \
			list->next = prev;        \
			prev = list;              \
			list = next;              \
		}                                 \
                                                  \
		list = prev;                      \
	} while (0)

enum ff_res efp_ff_parse(struct ff *ff, const char *path)
{
	enum ff_res res = FF_OK;
	struct stream *stream;

	assert(ff);
	assert(path);
	assert(ff->bond_fps == NULL);
	assert(ff->angle_fps == NULL);
	assert(ff->torsion_fps == NULL);

	if ((stream = efp_stream_open(path)) == NULL)
		return FF_FILE_NOT_FOUND;

	while (!efp_stream_eof(stream)) {
		efp_stream_next_line(stream);
		efp_stream_skip_space(stream);

		if (efp_stream_eol(stream))
			continue;

		if (strncmp(efp_stream_get_ptr(stream), "bond",
				strlen("bond")) == 0) {
			if ((res = parse_bond_fp(ff, efp_stream_get_ptr(stream))))
				break;
		}
		else if (strncmp(efp_stream_get_ptr(stream), "angle",
				strlen("angle")) == 0) {
			if ((res = parse_angle_fp(ff, efp_stream_get_ptr(stream))))
				break;
		}
		else if (strncmp(efp_stream_get_ptr(stream), "torsion",
				strlen("torsion")) == 0) {
			if ((res = parse_torsion_fp(ff, efp_stream_get_ptr(stream))))
				break;
		}
		else {
			res = FF_BAD_FORMAT;
			break;
		}
	}

	REVERSE_LIST(ff->bond_fps, struct bond_fp);
	REVERSE_LIST(ff->angle_fps, struct angle_fp);
	REVERSE_LIST(ff->torsion_fps, struct torsion_fp);

	efp_stream_close(stream);
	return res;
}

enum ff_res efp_ff_add_atom(struct ff *ff, const char *type)
{
	struct atom *atom;

	assert(ff);
	assert(type);

	if (strlen(type) >= sizeof(atom->type))
		return FF_STRING_TOO_LONG;

	ff->n_atoms++;
	ff->atoms = realloc(ff->atoms, ff->n_atoms * sizeof(struct atom));
	atom = ff->atoms + ff->n_atoms - 1;
	memset(atom, 0, sizeof(struct atom));
	strcpy(atom->type, type);

	return FF_OK;
}

enum ff_res efp_ff_add_bond(struct ff *ff, size_t idx1, size_t idx2)
{
	assert(ff);
	assert(idx1 < ff->n_atoms);
	assert(idx2 < ff->n_atoms);
	assert(idx1 != idx2);

	struct bond_fp *fp = find_bond_fp(ff->bond_fps,
			ff->atoms[idx1].type, ff->atoms[idx2].type);

	if (!fp)
		return FF_NO_PARAMETERS;

	struct bond *bond = malloc(sizeof(struct bond));

	bond->fp = fp;
	bond->idx1 = idx1;
	bond->idx2 = idx2;
	bond->next = ff->bonds;
	ff->bonds = bond;

	return FF_OK;
}

enum ff_res efp_ff_add_angle(struct ff *ff, size_t idx1, size_t idx2, size_t idx3)
{
	assert(ff);
	assert(idx1 < ff->n_atoms);
	assert(idx2 < ff->n_atoms);
	assert(idx3 < ff->n_atoms);
	assert(idx1 != idx2 && idx1 != idx3 && idx2 != idx3);

	struct angle_fp *fp = find_angle_fp(ff->angle_fps,
			ff->atoms[idx1].type, ff->atoms[idx2].type,
			ff->atoms[idx3].type);

	if (!fp)
		return FF_NO_PARAMETERS;

	struct angle *angle = malloc(sizeof(struct angle));

	angle->fp = fp;
	angle->idx1 = idx1;
	angle->idx2 = idx2;
	angle->idx3 = idx3;
	angle->next = ff->angles;
	ff->angles = angle;

	return FF_OK;
}

enum ff_res efp_ff_add_torsion(struct ff *ff, size_t idx1, size_t idx2,
		size_t idx3, size_t idx4)
{
	assert(ff);
	assert(idx1 < ff->n_atoms);
	assert(idx2 < ff->n_atoms);
	assert(idx3 < ff->n_atoms);
	assert(idx4 < ff->n_atoms);
	assert(idx1 != idx2 && idx1 != idx3 && idx1 != idx4);
	assert(idx2 != idx3 && idx2 != idx4 && idx3 != idx4);

	struct torsion_fp *fp = find_torsion_fp(ff->torsion_fps,
			ff->atoms[idx1].type, ff->atoms[idx2].type,
			ff->atoms[idx3].type, ff->atoms[idx4].type);

	if (!fp)
		return FF_NO_PARAMETERS;

	struct torsion *torsion = malloc(sizeof(struct torsion));

	torsion->fp = fp;
	torsion->idx1 = idx1;
	torsion->idx2 = idx2;
	torsion->idx3 = idx3;
	torsion->idx4 = idx4;
	torsion->next = ff->torsions;
	ff->torsions = torsion;

	return FF_OK;
}

enum ff_res efp_ff_auto_angles(struct ff *ff)
{
	enum ff_res res;
	size_t idx1, idx2, idx3;

	assert(ff);

	for (const struct bond *b1 = ff->bonds; b1; b1 = b1->next) {
		for (const struct bond *b2 = b1->next; b2; b2 = b2->next) {
			if (check_angle(b1, b2, &idx1, &idx2, &idx3))
				if ((res = efp_ff_add_angle(ff, idx1, idx2, idx3)))
					return res;
		}
	}

	return FF_OK;
}

enum ff_res efp_ff_auto_torsions(struct ff *ff)
{
	enum ff_res res;
	size_t idx1, idx2, idx3, idx4;

	assert(ff);

	for (const struct angle *a1 = ff->angles; a1; a1 = a1->next) {
		for (const struct angle *a2 = a1->next; a2; a2 = a2->next) {
			if (check_torsion(a1, a2, &idx1, &idx2, &idx3, &idx4))
				if ((res = efp_ff_add_torsion(ff, idx1, idx2, idx3, idx4)))
					return res;
		}
	}

	return FF_OK;
}

size_t efp_ff_get_atom_count(const struct ff *ff)
{
	assert(ff);

	return ff->n_atoms;
}

void efp_ff_set_atom_pos(struct ff *ff, size_t idx, vec_t pos)
{
	assert(ff);
	assert(idx < ff->n_atoms);

	pos.x /= ANG2BOHR;
	pos.y /= ANG2BOHR;
	pos.z /= ANG2BOHR;

	ff->atoms[idx].pos = pos;
}

vec_t efp_ff_get_atom_pos(const struct ff *ff, size_t idx)
{
	assert(ff);
	assert(idx < ff->n_atoms);

	vec_t pos = {
		ff->atoms[idx].pos.x * ANG2BOHR,
		ff->atoms[idx].pos.y * ANG2BOHR,
		ff->atoms[idx].pos.z * ANG2BOHR
	};

	return pos;
}

double efp_ff_compute(const struct ff *ff, vec_t *grad)
{
	double energy = 0.0;

	assert(ff);

	if (grad)
		memset(grad, 0, ff->n_atoms * sizeof(vec_t));

	for (struct bond *bond = ff->bonds; bond; bond = bond->next)
		energy += bond_energy(bond, ff->atoms, grad);

	for (struct angle *angle = ff->angles; angle; angle = angle->next)
		energy += angle_energy(angle, ff->atoms, grad);

	for (struct torsion *torsion = ff->torsions; torsion; torsion = torsion->next)
		energy += torsion_energy(torsion, ff->atoms, grad);

	if (grad)
		for (size_t i = 0; i < ff->n_atoms; i++)
			vec_scale(grad + i, KCALMOL2AU / ANG2BOHR);

	return energy * KCALMOL2AU;
}

#define FREE_LIST(list, type)                                        \
	while (list) {                                               \
		type *ptr = list;                                    \
		list = ptr->next;                                    \
		free(ptr);                                           \
	}

void efp_ff_free(struct ff *ff)
{
	if (ff) {
		FREE_LIST(ff->bonds, struct bond);
		FREE_LIST(ff->angles, struct angle);
		FREE_LIST(ff->torsions, struct torsion);
		FREE_LIST(ff->bond_fps, struct bond_fp);
		FREE_LIST(ff->angle_fps, struct angle_fp);
		FREE_LIST(ff->torsion_fps, struct torsion_fp);
		free(ff->atoms);
		free(ff);
	}
}
