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

#ifndef LIBEFP_TEST_COMMON_H
#define LIBEFP_TEST_COMMON_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <efp.h>

struct test_data {
	/** Array of paths to EFP data files, NULL-terminated. */
	const char **potential_files;

	/** Array of fragment names, NULL-terminated. */
	const char **fragname;

	/**
	 * If not null geometry will be set up using efp_set_coordinates
	 * with this field as a parameter. */
	const double *geometry_xyzabc;

	/**
	 * If not null geometry will be set up using efp_set_coordinates_2
	 * with this field as a parameter. */
	const double *geometry_points;

	/** Reference energy value. */
	double ref_energy;

	/**
	 * Reference gradient values. Turns on gradient test if not NULL.
	 * The size of this array must be 6N, where N is the number of
	 * fragments. */
	const double *ref_gradient;

	/**
	 * Turns on comparison of analytical and numerical gradients if
	 * nonzero. */
	int test_numerical_gradient;

	/** Simulation settings. */
	struct efp_opts opts;

	/** Callbacks. */
	struct efp_callbacks callbacks;
};

#define ARRAY_SIZE(x) (sizeof(x) / sizeof(x[0]))
#define BOHR(x) ((x) / 0.52917721092)
#define ANGSTROM(x) ((x) * 0.52917721092)

int run_test(const struct test_data *test_data);

enum efp_result st_integrals_from_file(const struct efp_st_block *block,
				       int compute_derivatives,
				       struct efp_st_data *st,
				       void *user_data,
				       int expected_size_i,
				       int expected_size_j,
				       const char *s_path,
				       const char *t_path);

#define DEFINE_TEST(test_data)                                               \
int main(void)                                                               \
{                                                                            \
	return run_test(&(test_data));                                       \
}

#endif /* LIBEFP_TEST_COMMON_H */
