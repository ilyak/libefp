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

#include "int.h"

/* Overlap and kinetic energy integral computation routines.
 * Used by exchange-repulsion code. */

static void
set_coef(double *con, char type, const double *coef)
{
	memset(con, 0, 20 * sizeof(double));

	switch (type) {
		case 'S':
			con[0] = *coef;
			return;
		case 'L':
			con[0] = *coef++;
			/* fall through */
		case 'P':
			for (int i = 1; i < 4; i++)
				con[i] = *coef;
			return;
		case 'D':
			for (int i = 4; i < 10; i++)
				con[i] = *coef;
			return;
		case 'F':
			for (int i = 10; i < 20; i++)
				con[i] = *coef;
			return;
	}
}

static void
make_int(int ni, int nj, double tt, double x, double y, double z,
	 double xi, double yi, double zi, double xj, double yj, double zj,
	 double *outx, double *outy, double *outz)
{
	static const int imin[] = { 0, 1, 3,  6, 10, 15, 21, 28, 36, 45 };
	static const int imax[] = { 1, 3, 6, 10, 15, 21, 28, 36, 45, 55 };

	static const double h[] = {
		 0.0000000000000000, -0.7071067811865475,  0.7071067811865475,
		-1.2247448713915889,  0.0000000000000000,  1.2247448713915889,
		-1.6506801238857844, -0.5246476232752903,  0.5246476232752903,
		 1.6506801238857844, -2.0201828704560856, -0.9585724646138185,
		 0.0000000000000000,  0.9585724646138185,  2.0201828704560856,
		-2.3506049736744923, -1.3358490740136970, -0.4360774119276165,
		 0.4360774119276165,  1.3358490740136970,  2.3506049736744923,
		-2.6519613568352334, -1.6735516287674714, -0.8162878828589647,
		 0.0000000000000000,  0.8162878828589647,  1.6735516287674714,
		 2.6519613568352334
	};

	static const double w[] = {
		 1.7724538509055161,  0.8862269254527580,  0.8862269254527580,
		 0.2954089751509193,  1.1816359006036774,  0.2954089751509193,
		 0.0813128354472451,  0.8049140900055128,  0.8049140900055128,
		 0.0813128354472451,  0.0199532420590459,  0.3936193231522411,
		 0.9453087204829419,  0.3936193231522411,  0.0199532420590459,
		 0.0045300099055088,  0.1570673203228566,  0.7246295952243925,
		 0.7246295952243925,  0.1570673203228566,  0.0045300099055088,
		 0.0009717812450995,  0.0545155828191270,  0.4256072526101277,
		 0.8102646175568073,  0.4256072526101277,  0.0545155828191270,
		 0.0009717812450995
	};

	int npts = (ni + nj) / 2;

	double xint = 0.0, yint = 0.0, zint = 0.0;

	for (int i = imin[npts]; i < imax[npts]; i++) {
		double px = w[i];
		double py = w[i];
		double pz = w[i];

		double tmp = h[i] * tt;

		if(ni > 0) {
			double ax = tmp + x - xi;
			double ay = tmp + y - yi;
			double az = tmp + z - zi;

			for (int j = 0; j < ni; j++) {
				px *= ax;
				py *= ay;
				pz *= az;
			}
		}

		if(nj > 0) {
			double bx = tmp + x - xj;
			double by = tmp + y - yj;
			double bz = tmp + z - zj;

			for (int j = 0; j < nj; j++) {
				px *= bx;
				py *= by;
				pz *= bz;
			}
		}

		xint += px;
		yint += py;
		zint += pz;
	}

	*outx = xint;
	*outy = yint;
	*outz = zint;
}

static int
shell_type(char type)
{
	switch (type) {
		case 'S': return 0;
		case 'L': return 1;
		case 'P': return 2;
		case 'D': return 3;
		case 'F': return 4;
	}
	return -1;
}

void
efp_st_int(struct efp *efp, int frag_i, int frag_j, double *s, double *t)
{
	static const int shift_ix[] = {
		0, 5, 0,  0, 10, 0, 0, 5, 5, 0,
		5, 0, 0, 10, 10, 5, 0, 5, 0, 5
	};
	static const int shift_iy[] = {
		0,  0, 5, 0, 0, 10,  0, 5, 0, 5,
		0, 15, 0, 5, 0, 10, 10, 0, 5, 5
	};
	static const int shift_iz[] = {
		0, 0,  0, 5, 0, 0, 10,  0,  5, 5,
		0, 0, 15, 0, 5, 0,  5, 10, 10, 5
	};
	static const int shift_jx[] = {
		0, 1, 0, 0, 2, 0, 0, 1, 1, 0,
		3, 0, 0, 2, 2, 1, 0, 1, 0, 1
	};
	static const int shift_jy[] = {
		0, 0, 1, 0, 0, 2, 0, 1, 0, 1,
		0, 3, 0, 1, 0, 2, 2, 0, 1, 1
	};
	static const int shift_jz[] = {
		0, 0, 0, 1, 0, 0, 2, 0, 1, 1,
		0, 0, 3, 0, 1, 0, 1, 2, 2, 1
	};

	static const double norm[] = {
		1.0000000000000000, 1.0000000000000000, 1.0000000000000000,
		1.0000000000000000, 1.0000000000000000, 1.0000000000000000,
		1.0000000000000000, 1.7320508075688801, 1.7320508075688801,
		1.7320508075688801, 1.0000000000000000, 1.0000000000000000,
		1.0000000000000000, 2.2360679774997898, 2.2360679774997898,
		2.2360679774997898, 2.2360679774997898, 2.2360679774997898,
		2.2360679774997898, 3.8729833462074232
	};

	static const double tolerance = 46.051701859881; /* 20 ln10 */

	static const int primitive_start[] = {
		0, 0, 1, 4, 10
	};
	static const int primitive_end[] = {
		1, 4, 4, 10, 20
	};
	static const int shell_l[] = {
		1, 2, 2, 3, 4
	};

	struct frag *fr_i = efp->frags + frag_i;
	struct frag *fr_j = efp->frags + frag_j;

	int shift_x[100];
	int shift_y[100];
	int shift_z[100];

	double xin[90];
	double yin[90];
	double zin[90];

	double ft[100], dij[100];
	double sblk[100], tblk[100];

	/* shell i */
	for (int ii = 0, loc_i = 0; ii < fr_i->n_xr_shells; ii++) {
		struct shell *sh_i = fr_i->xr_shells + ii;
		const vec_t *ri = VEC(fr_i->atoms[sh_i->atom_index].x);

		int type_i = shell_type(sh_i->type);
		int start_i = primitive_start[type_i];
		int end_i = primitive_end[type_i];
		int count_i = end_i - start_i;

		/* shell j */
		for (int jj = 0, loc_j = 0; jj < fr_j->n_xr_shells; jj++) {
			struct shell *sh_j = fr_j->xr_shells + jj;
			const vec_t *rj = VEC(fr_j->atoms[sh_j->atom_index].x);

			int type_j = shell_type(sh_j->type);
			int start_j = primitive_start[type_j];
			int end_j = primitive_end[type_j];
			int count_j = end_j - start_j;

			int count = 0;

			for (int i = start_i; i < end_i; i++) {
				for (int j = start_j; j < end_j; j++, count++) {
					shift_x[count] = shift_ix[i] + shift_jx[j];
					shift_y[count] = shift_iy[i] + shift_jy[j];
					shift_z[count] = shift_iz[i] + shift_jz[j];

					if (j < 1)
						ft[count] = 3.0;
					else if (j >= 1 && j < 4)
						ft[count] = 5.0;
					else if (j >= 4 && j < 10)
						ft[count] = 7.0;
					else
						ft[count] = 9.0;
				}
			}

			memset(sblk, 0, count * sizeof(double));
			memset(tblk, 0, count * sizeof(double));

			const double *coef_i = sh_i->coef;

			/* primitive i */
			for (int ig = 0; ig < sh_i->n_funcs; ig++) {
				double ai = *coef_i++;

				double con_i[20];
				set_coef(con_i, sh_i->type, coef_i);

				coef_i++;
				if (sh_i->type == 'L')
					coef_i++;

				const double *coef_j = sh_j->coef;

				/* primitive j */
				for (int jg = 0; jg < sh_j->n_funcs; jg++) {
					double aj = *coef_j++;

					double aa = 1.0 / (ai + aj);
					double dum = aj * ai * vec_dist_2(ri, rj) * aa;

					if (dum > tolerance) {
						coef_j++;
						if (sh_j->type == 'L')
							coef_j++;

						continue;
					}

					double con_j[20];
					set_coef(con_j, sh_j->type, coef_j);

					coef_j++;
					if (sh_j->type == 'L')
						coef_j++;

					double ax = (ai * ri->x + aj * rj->x) * aa;
					double ay = (ai * ri->y + aj * rj->y) * aa;
					double az = (ai * ri->z + aj * rj->z) * aa;

					double fac = exp(-dum);

					for (int i = start_i, idx = 0; i < end_i; i++)
						for (int j = start_j; j < end_j; j++, idx++)
							dij[idx] = fac * con_i[i] * norm[i] * con_j[j] * norm[j];

					double taa = sqrt(aa);
					double t1 = -2.0 * aj * aj * taa;
					double t2 = -0.5 * taa;

					for (int i = 0, idx = 0; i < shell_l[type_i]; i++, idx += 5) {
						for (int j = 0; j < shell_l[type_j]; j++) {
							double xint, yint, zint;

							make_int(i, j, taa, ax, ay, az,
								 ri->x, ri->y, ri->z,
								 rj->x, rj->y, rj->z,
								 &xint, &yint, &zint);
							xin[idx + j] = xint * taa;
							yin[idx + j] = yint * taa;
							zin[idx + j] = zint * taa;

							make_int(i, j + 2, taa, ax, ay, az,
								 ri->x, ri->y, ri->z,
								 rj->x, rj->y, rj->z,
								 &xint, &yint, &zint);
							xin[idx + j + 30] = xint * t1;
							yin[idx + j + 30] = yint * t1;
							zin[idx + j + 30] = zint * t1;

							if (j >= 2) {
								make_int(i, j - 2, taa, ax, ay, az,
									 ri->x, ri->y, ri->z,
									 rj->x, rj->y, rj->z,
									 &xint, &yint, &zint);
							}
							else {
								xint = 0.0;
								yint = 0.0;
								zint = 0.0;
							}
							double t3 = j * (j - 1) * t2;
							xin[idx + j + 60] = xint * t3;
							yin[idx + j + 60] = yint * t3;
							zin[idx + j + 60] = zint * t3;
						}
					}
					for (int i = 0; i < count; i++) {
						int nx = shift_x[i];
						int ny = shift_y[i];
						int nz = shift_z[i];
						double tmp1 = xin[nx] * yin[ny] * zin[nz];
						double tmp2 = (xin[nx + 30] + xin[nx + 60]) * yin[ny] * zin[nz] +
							      (yin[ny + 30] + yin[ny + 60]) * xin[nx] * zin[nz] +
							      (zin[nz + 30] + zin[nz + 60]) * xin[nx] * yin[ny];
						sblk[i] = sblk[i] + dij[i] * tmp1;
						tblk[i] = tblk[i] + dij[i] * (tmp1 * aj * ft[i] + tmp2);
					}
				}
			}

			/* store integrals */
			for (int i = 0, idx = 0; i < count_i; i++) {
				int idx2 = (loc_i + i) * fr_j->xr_wf_size + loc_j;

				for (int j = 0; j < count_j; j++, idx++, idx2++) {
					s[idx2] = sblk[idx];
					t[idx2] = tblk[idx];
				}
			}
			loc_j += count_j;
		}
		loc_i += count_i;
	}
}

void
efp_st_int_deriv(struct efp *efp, int frag_i, int frag_j,
		 double *sx, double *tx)
{
}
