#ifndef LIBEFP_UTIL_H
#define LIBEFP_UTIL_H

#include "math_util.h"

int check_rotation_matrix(const mat_t *);
int matrix_eigen(const mat_t *, vec_t *, mat_t *);
void matrix_to_euler(const mat_t *, double *, double *, double *);
void euler_to_matrix(double, double, double, mat_t *);
void points_to_matrix(const double *, mat_t *);
void torque_to_deriv(int, const double *, double *);
void move_pt(const vec_t *, const mat_t *, const vec_t *, const vec_t *, vec_t *);
void rotate_t2(const mat_t *, const double *, double *);
void rotate_t3(const mat_t *, const double *, double *);

#endif /* LIBEFP_UTIL_H */
