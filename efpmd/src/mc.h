#include <stddef.h>

enum mc_result {
	MC_RESULT_SUCCESS = 0,
	MC_RESULT_ERROR
};

struct mc_state;


typedef double (*mc_func_t)(size_t, const double *, void *);

struct mc_state *mc_create(size_t);
enum mc_result mc_init(struct mc_state *, size_t, const double *);
void mc_set_func(struct mc_state *, mc_func_t);
void mc_set_user_data(struct mc_state *, void *);

double mc_get_fx(struct mc_state *);
void mc_get_x(struct mc_state *, size_t, double *);

void mc_shutdown(struct mc_state *);



