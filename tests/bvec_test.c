#include <assert.h>

#include "bvec.h"

static void test_1(void)
{
	struct bvec *bvec = efp_bvec_create(64);

	efp_bvec_set(bvec, 7);
	assert(efp_bvec_is_set(bvec, 7));

	efp_bvec_unset(bvec, 7);
	assert(!efp_bvec_is_set(bvec, 7));

	efp_bvec_set(bvec, 63);
	assert(efp_bvec_is_set(bvec, 63));

	efp_bvec_unset(bvec, 63);
	assert(!efp_bvec_is_set(bvec, 63));

	efp_bvec_set_val(bvec, 0, true);
	assert(efp_bvec_is_set(bvec, 0));

	efp_bvec_set_val(bvec, 0, false);
	assert(!efp_bvec_is_set(bvec, 0));

	efp_bvec_free(bvec);
}

static void test_2(void)
{
	struct bvec *bvec = efp_bvec_create(117);

	efp_bvec_set(bvec, 7);
	assert(efp_bvec_is_set(bvec, 7));

	efp_bvec_unset(bvec, 7);
	assert(!efp_bvec_is_set(bvec, 7));

	efp_bvec_set(bvec, 100);
	assert(efp_bvec_is_set(bvec, 100));

	efp_bvec_unset(bvec, 100);
	assert(!efp_bvec_is_set(bvec, 100));

	efp_bvec_set_val(bvec, 64, true);
	assert(efp_bvec_is_set(bvec, 64));

	efp_bvec_set_val(bvec, 64, false);
	assert(!efp_bvec_is_set(bvec, 64));

	efp_bvec_free(bvec);
}

int main(void)
{
	test_1();
	test_2();

	return 0;
}
