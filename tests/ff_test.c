#include <assert.h>
#include <stdio.h>

#include "ff.h"

#define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))
#define check(x) res = (x); assert(res == FF_OK)

#define GRAD_DELTA 1.0e-3
#define GRAD_ACCURACY 1.0e-3 /* XXX convert to au */

static void test_1(struct ff *ff)
{
	enum ff_res res;

	check(efp_ff_add_atom(ff, "HW"));
	check(efp_ff_add_atom(ff, "OW"));
	check(efp_ff_add_atom(ff, "HW"));

	check(efp_ff_add_bond(ff, 1, 0));
	check(efp_ff_add_bond(ff, 1, 2));

	check(efp_ff_add_angle(ff, 0, 1, 2));

	efp_ff_set_atom_pos(ff, 0, (vec_t) { 0.0,  0.8, -0.5 });
	efp_ff_set_atom_pos(ff, 1, (vec_t) { 0.0,  0.0,  0.1 });
	efp_ff_set_atom_pos(ff, 2, (vec_t) { 0.0, -0.8, -0.5 });
}

static void test_2(struct ff *ff)
{
	enum ff_res res;

	check(efp_ff_add_atom(ff, "CT"));
	check(efp_ff_add_atom(ff, "HC"));
	check(efp_ff_add_atom(ff, "HC"));
	check(efp_ff_add_atom(ff, "HC"));
	check(efp_ff_add_atom(ff, "CT"));
	check(efp_ff_add_atom(ff, "HC"));
	check(efp_ff_add_atom(ff, "HC"));
	check(efp_ff_add_atom(ff, "CT"));
	check(efp_ff_add_atom(ff, "HC"));
	check(efp_ff_add_atom(ff, "HC"));
	check(efp_ff_add_atom(ff, "CT"));
	check(efp_ff_add_atom(ff, "HC"));
	check(efp_ff_add_atom(ff, "HC"));
	check(efp_ff_add_atom(ff, "HC"));

	check(efp_ff_add_bond(ff, 1, 0));
	check(efp_ff_add_bond(ff, 2, 0));
	check(efp_ff_add_bond(ff, 3, 0));
	check(efp_ff_add_bond(ff, 0, 4));
	check(efp_ff_add_bond(ff, 5, 4));
	check(efp_ff_add_bond(ff, 6, 4));
	check(efp_ff_add_bond(ff, 4, 7));
	check(efp_ff_add_bond(ff, 8, 7));
	check(efp_ff_add_bond(ff, 9, 7));
	check(efp_ff_add_bond(ff, 7, 10));
	check(efp_ff_add_bond(ff, 11, 10));
	check(efp_ff_add_bond(ff, 12, 10));
	check(efp_ff_add_bond(ff, 13, 10));

	check(efp_ff_auto_angles(ff));
	check(efp_ff_auto_torsions(ff));

	efp_ff_set_atom_pos(ff,  0, (vec_t) { -6.0492, -0.0134,  0.0000 });
	efp_ff_set_atom_pos(ff,  1, (vec_t) { -6.8305,  0.7465,  0.0000 });
	efp_ff_set_atom_pos(ff,  2, (vec_t) { -6.1478, -0.6352,  0.8899 });
	efp_ff_set_atom_pos(ff,  3, (vec_t) { -6.1478, -0.6352, -0.8899 });
	efp_ff_set_atom_pos(ff,  4, (vec_t) { -4.6686,  0.6688,  0.0000 });
	efp_ff_set_atom_pos(ff,  5, (vec_t) { -4.6279,  1.2966, -0.8901 });
	efp_ff_set_atom_pos(ff,  6, (vec_t) { -4.6279,  1.2966,  0.8901 });
	efp_ff_set_atom_pos(ff,  7, (vec_t) { -3.3875, -0.1858,  0.0000 });
	efp_ff_set_atom_pos(ff,  8, (vec_t) { -3.4281, -0.8135,  0.8901 });
	efp_ff_set_atom_pos(ff,  9, (vec_t) { -3.4281, -0.8135, -0.8901 });
	efp_ff_set_atom_pos(ff, 10, (vec_t) { -2.0069,  0.4965,  0.0000 });
	efp_ff_set_atom_pos(ff, 11, (vec_t) { -1.2255, -0.2635,  0.0000 });
	efp_ff_set_atom_pos(ff, 12, (vec_t) { -1.9083,  1.1182,  0.8899 });
	efp_ff_set_atom_pos(ff, 13, (vec_t) { -1.9083,  1.1182, -0.8899 });
}

static void test_3(struct ff *ff)
{
	enum ff_res res;

	check(efp_ff_add_atom(ff, "H5"));
	check(efp_ff_add_atom(ff, "C"));
	check(efp_ff_add_atom(ff, "O"));
	check(efp_ff_add_atom(ff, "N"));
	check(efp_ff_add_atom(ff, "H"));
	check(efp_ff_add_atom(ff, "H"));

	check(efp_ff_add_bond(ff, 0, 1));
	check(efp_ff_add_bond(ff, 2, 1));
	check(efp_ff_add_bond(ff, 1, 3));
	check(efp_ff_add_bond(ff, 3, 4));
	check(efp_ff_add_bond(ff, 3, 5));

	check(efp_ff_auto_angles(ff));
	check(efp_ff_auto_torsions(ff));

	efp_ff_set_atom_pos(ff, 0, (vec_t) { -7.127, -2.105, -0.006 });
	efp_ff_set_atom_pos(ff, 1, (vec_t) { -6.291, -1.451,  0.000 });
	efp_ff_set_atom_pos(ff, 2, (vec_t) { -6.464, -0.227,  0.008 });
	efp_ff_set_atom_pos(ff, 3, (vec_t) { -4.996, -1.973, -0.000 });
	efp_ff_set_atom_pos(ff, 4, (vec_t) { -4.213, -1.359,  0.005 });
	efp_ff_set_atom_pos(ff, 5, (vec_t) { -4.857, -2.958, -0.007 });
}

static void test_4(struct ff *ff)
{
	enum ff_res res;

	check(efp_ff_add_atom(ff, "CT"));
	check(efp_ff_add_atom(ff, "H1"));
	check(efp_ff_add_atom(ff, "H1"));
	check(efp_ff_add_atom(ff, "H1"));
	check(efp_ff_add_atom(ff, "CT"));
	check(efp_ff_add_atom(ff, "H1"));
	check(efp_ff_add_atom(ff, "H1"));
	check(efp_ff_add_atom(ff, "OH"));
	check(efp_ff_add_atom(ff, "HO"));

	check(efp_ff_add_bond(ff, 0, 1));
	check(efp_ff_add_bond(ff, 0, 2));
	check(efp_ff_add_bond(ff, 0, 3));
	check(efp_ff_add_bond(ff, 0, 4));
	check(efp_ff_add_bond(ff, 4, 5));
	check(efp_ff_add_bond(ff, 4, 6));
	check(efp_ff_add_bond(ff, 4, 7));
	check(efp_ff_add_bond(ff, 7, 8));

	check(efp_ff_auto_angles(ff));
	check(efp_ff_auto_torsions(ff));

	efp_ff_set_atom_pos(ff, 0, (vec_t) { -1.745,  1.106,  1.225 });
	efp_ff_set_atom_pos(ff, 1, (vec_t) { -1.675,  0.626,  2.226 });
	efp_ff_set_atom_pos(ff, 2, (vec_t) { -0.836,  1.733,  1.088 });
	efp_ff_set_atom_pos(ff, 3, (vec_t) { -2.624,  1.789,  1.242 });
	efp_ff_set_atom_pos(ff, 4, (vec_t) { -1.874,  0.067,  0.109 });
	efp_ff_set_atom_pos(ff, 5, (vec_t) { -2.785, -0.555,  0.253 });
	efp_ff_set_atom_pos(ff, 6, (vec_t) { -1.962,  0.576, -0.876 });
	efp_ff_set_atom_pos(ff, 7, (vec_t) { -0.720, -0.754,  0.101 });
	efp_ff_set_atom_pos(ff, 8, (vec_t) { -0.857, -1.462,  0.717 });
}

static void run_test(struct ff *ff)
{
	size_t n_atoms = efp_ff_get_atom_count(ff);
	vec_t grad[n_atoms];

	double energy = efp_ff_compute(ff, grad);
	printf("energy = %.8e au\n", energy);

	for (size_t i = 0; i < n_atoms; i++) {
		vec_t pos = efp_ff_get_atom_pos(ff, i);

		for (size_t j = 0; j < 3; j++) {
			vec_t pos2 = pos;
			double e1, e2, crd;

			crd = vec_get(&pos, j) - GRAD_DELTA;
			vec_set(&pos2, j, crd);
			efp_ff_set_atom_pos(ff, i, pos2);
			e1 = efp_ff_compute(ff, NULL);

			crd = vec_get(&pos, j) + GRAD_DELTA;
			vec_set(&pos2, j, crd);
			efp_ff_set_atom_pos(ff, i, pos2);
			e2 = efp_ff_compute(ff, NULL);

			double numer = (e2 - e1) / (2 * GRAD_DELTA);
			double analyt = vec_get(grad + i, j);

			assert(fabs(numer - analyt) < GRAD_ACCURACY);
		}

		efp_ff_set_atom_pos(ff, i, pos);
	}

	printf("numerical and analytical gradients match\n");
}

int main(void)
{
	enum ff_res res;
	struct ff *ff;

	void (*init[])(struct ff *) = {
		test_1, test_2, test_3, test_4
	};

	for (size_t i = 0; i < ARRAY_SIZE(init); i++) {
		printf("starting test %zu\n", i + 1);
		ff = efp_ff_create();
		assert(ff);
		check(efp_ff_parse(ff, "../fraglib/amber99.ff"));
		init[i](ff);
		run_test(ff);
		efp_ff_free(ff);
	}

	return 0;
}
