#ifndef EFPMD_CFG_H
#define EFPMD_CFG_H

#include <efp.h>

enum run_type {
	RUN_TYPE_SP,
	RUN_TYPE_GRAD,
	RUN_TYPE_HESS,
	RUN_TYPE_OPT,
	RUN_TYPE_MD
};

enum ensemble_type {
	ENSEMBLE_TYPE_NVE,
	ENSEMBLE_TYPE_NVT
};

struct frag {
	char *name;
	double coord[9];
	double vel[6];
};

struct config {
	enum run_type run_type;
	enum efp_coord_type coord_type;
	double units_factor;
	unsigned terms;
	enum efp_elec_damp elec_damp;
	enum efp_disp_damp disp_damp;
	enum efp_pol_damp pol_damp;
	double hess_delta;
	int max_steps;
	int print_step;
	double temperature;
	double time_step;
	enum ensemble_type ensemble_type;
	double opt_tol;
	char *fraglib_path;
	char *userlib_path;
	int n_frags;
	struct frag *frags;
};

void parse_config(const char *, struct config *);
void free_config(struct config *);

#endif /* EFPMD_CFG_H */
