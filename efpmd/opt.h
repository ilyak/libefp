typedef double (*opt_func_t)(int, const double *, double *, void *);
double opt_cg_step(opt_func_t, int, double *, double *, void *);
