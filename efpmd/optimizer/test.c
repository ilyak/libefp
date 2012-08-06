#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "opt.h"

#define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))
#define UNUSED __attribute__((unused))

#define OPT_STEPS 50

static enum opt_result fn_1(size_t n, const double *x, double *fx,
		double *gx, UNUSED void *data)
{
	assert(n == 1);

	*fx = sin(*x);
	*gx = cos(*x);

	return OPT_RESULT_SUCCESS;
}

static enum opt_result fn_2(size_t n, const double *x, double *fx,
		double *gx, UNUSED void *data)
{
	assert(n == 2);

	*fx = x[0] * x[0] + x[1] * x[1];

	gx[0] = 2.0 * x[0];
	gx[1] = 2.0 * x[1];

	return OPT_RESULT_SUCCESS;
}

static enum opt_result fn_3(size_t n, const double *x, double *fx,
		double *gx, UNUSED void *data)
{
	assert(n == 1);

	*fx = (*x) * (*x) - 2.0 * (*x);
	*gx = 2.0 * (*x) - 2.0;

	return OPT_RESULT_SUCCESS;
}

static enum opt_result fn_4(size_t n, const double *x, double *fx,
		double *gx, UNUSED void *data)
{
	assert(n == 2);

	*fx = sin(x[0] * x[1]);

	gx[0] = x[1] * cos(x[0] * x[1]);
	gx[1] = x[0] * cos(x[0] * x[1]);

	return OPT_RESULT_SUCCESS;
}

static const double x_1[] = { 2.0 };
static const double x_2[] = { 4.0, -3.0 };
static const double x_3[] = { 2.3 };
static const double x_4[] = { 0.3, 0.7 };

static const struct test {
	enum opt_result (*fn)(size_t, const double *, double *, double *, void *);
	size_t n_dim;
	const double *x;
} tests[] = {
	{ fn_1, ARRAY_SIZE(x_1), x_1 },
	{ fn_2, ARRAY_SIZE(x_2), x_2 },
	{ fn_3, ARRAY_SIZE(x_3), x_3 },
	{ fn_4, ARRAY_SIZE(x_4), x_4 }
};

static void print_array(size_t n, const double *x)
{
	for (size_t i = 0; i < n; i++, x++)
		printf("%8.3lf", *x);
}

static void do_test(const struct test *test)
{
	struct opt_state *state = opt_create(test->n_dim);
	assert(state);

	enum opt_result res;
	double x[test->n_dim], fx, gx[test->n_dim];

	opt_set_fn(state, test->fn);

	res = opt_init(state, test->x);
	assert(res == OPT_RESULT_SUCCESS);

	for (size_t i = 0; i < OPT_STEPS; i++) {
		res = opt_step(state);
		assert(res == OPT_RESULT_SUCCESS);

		fx = opt_get_fx(state);
		assert(fx != NAN);

		res = opt_get_x(state, test->n_dim, x);
		assert(res == OPT_RESULT_SUCCESS);

		res = opt_get_gx(state, test->n_dim, gx);
		assert(res == OPT_RESULT_SUCCESS);

		printf("step %d:\n", i + 1);
		printf("  fx: %8.3lf\n", fx);
		printf("   x: ");
		print_array(test->n_dim, x);
		printf("\n");
		printf("  gx: ");
		print_array(test->n_dim, gx);
		printf("\n");
	}

	opt_shutdown(state);
}

int main(void)
{
	for (size_t i = 0; i < ARRAY_SIZE(tests); i++) {
		printf("PERFORMING TEST %d\n\n", i + 1);
		do_test(tests + i);
		printf("\n\n");
	}

	return EXIT_SUCCESS;
}
