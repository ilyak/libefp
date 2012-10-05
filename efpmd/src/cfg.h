#ifndef EFPMD_CFG_H
#define EFPMD_CFG_H

#include <stdbool.h>
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
	double coord[12];
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
	bool enable_pbc;
	double pbc_box[3];
	double swf_cutoff;
	double hess_delta;
	int max_steps;
	int print_step;
	double target_temperature;
	double time_step;
	enum ensemble_type ensemble_type;
	double thermostat_tau;
	double opt_tol;
	char *fraglib_path;
	char *userlib_path;
	int n_frags;
	struct frag *frags;
};

struct config *parse_config(const char *);
void free_config(struct config *);

#endif /* EFPMD_CFG_H */
