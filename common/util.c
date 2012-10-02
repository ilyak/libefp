#include "util.h"
#include "clapack.h"

int check_rotation_matrix(const mat_t *rotmat)
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

int matrix_eigen(const mat_t *mat, vec_t *eval, mat_t *evec)
{
	double work[32];

	*evec = *mat;

	return u_dsyev('V', 'U', 3, (double *)evec, 3, (double *)eval, work, 32);
}

void euler_to_matrix(double a, double b, double c, mat_t *out)
{
	double sina = sin(a), cosa = cos(a);
	double sinb = sin(b), cosb = cos(b);
	double sinc = sin(c), cosc = cos(c);

	out->xx =  cosa * cosc - sina * cosb * sinc;
	out->xy = -cosa * sinc - sina * cosb * cosc;
	out->xz =  sinb * sina;
	out->yx =  sina * cosc + cosa * cosb * sinc;
	out->yy = -sina * sinc + cosa * cosb * cosc;
	out->yz = -sinb * cosa;
	out->zx =  sinb * sinc;
	out->zy =  sinb * cosc;
	out->zz =  cosb;
}

void matrix_to_euler(const mat_t *rotmat, double *ea, double *eb, double *ec)
{
	double a, b, c, sinb;

	b = acos(rotmat->zz);
	sinb = sqrt(1.0 - rotmat->zz * rotmat->zz);

	if (fabs(sinb) < 1.0e-7) {
		a = atan2(-rotmat->xy, rotmat->xx);
		c = 0.0;
	}
	else {
		a = atan2(rotmat->xz, -rotmat->yz);
		c = atan2(rotmat->zx, rotmat->zy);
	}

	*ea = a, *eb = b, *ec = c;
}

void points_to_matrix(const double *pts, mat_t *rotmat)
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

	rotmat->xx = r12.x, rotmat->xy = r13.x, rotmat->xz = cross.x;
	rotmat->yx = r12.y, rotmat->yy = r13.y, rotmat->yz = cross.y;
	rotmat->zx = r12.z, rotmat->zy = r13.z, rotmat->zz = cross.z;
}

void torque_to_deriv(int nb, const double *x, double *gx)
{
	for (int i = 0; i < nb; i++, x += 6, gx += 6) {
		double tx = gx[3];
		double ty = gx[4];
		double tz = gx[5];

		double sina = sin(x[3]);
		double cosa = cos(x[3]);
		double sinb = sin(x[4]);
		double cosb = cos(x[4]);

		gx[3] = tz;
		gx[4] = cosa * tx + sina * ty;
		gx[5] = sinb * sina * tx - sinb * cosa * ty + cosb * tz;
	}
}

void move_pt(const vec_t *com, const mat_t *rotmat, const vec_t *pos_int, vec_t *out)
{
	*out = mat_vec(rotmat, pos_int);
	out->x += com->x, out->y += com->y, out->z += com->z;
}

void rotate_t2(const mat_t *rotmat, const double *in, double *out)
{
	for (int i = 0; i < 3 * 3; i++)
		out[i] = 0.0;

	for (int a1 = 0; a1 < 3; a1++)
	for (int b1 = 0; b1 < 3; b1++)
		for (int a2 = 0; a2 < 3; a2++)
		for (int b2 = 0; b2 < 3; b2++)
			out[a2 * 3 + b2] += in[a1 * 3 + b1] *
					mat_get(rotmat, a2, a1) *
					mat_get(rotmat, b2, b1);
}

void rotate_t3(const mat_t *rotmat, const double *in, double *out)
{
	for (int i = 0; i < 3 * 3 * 3; i++)
		out[i] = 0.0;

	for (int a1 = 0; a1 < 3; a1++)
	for (int b1 = 0; b1 < 3; b1++)
	for (int c1 = 0; c1 < 3; c1++)
		for (int a2 = 0; a2 < 3; a2++)
		for (int b2 = 0; b2 < 3; b2++)
		for (int c2 = 0; c2 < 3; c2++)
			out[a2 * 9 + b2 * 3 + c2] += in[a1 * 9 + b1 * 3 + c1] *
					mat_get(rotmat, a2, a1) *
					mat_get(rotmat, b2, b1) *
					mat_get(rotmat, c2, c1);
}
