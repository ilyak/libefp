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

#include "common.h"

void sim_sp(struct efp *, const struct config *);
void sim_grad(struct efp *, const struct config *);
void sim_hess(struct efp *, const struct config *);
void sim_opt(struct efp *, const struct config *);
void sim_md(struct efp *, const struct config *);

static void run_sim(struct efp *efp, const struct config *config)
{
	static const struct {
		const char *run_type;
		const char *description;
		void (*sim_fn)(struct efp *, const struct config *);
	} sim_list[] = {
		{ "sp",   "SINGLE POINT ENERGY JOB", sim_sp   },
		{ "grad", "ENERGY GRADIENT JOB",     sim_grad },
		{ "hess", "HESSIAN JOB",             sim_hess },
		{ "opt",  "ENERGY MINIMIZATION JOB", sim_opt  },
		{ "md",   "MOLECULAR DYNAMICS JOB",  sim_md   }
	};

	for (size_t i = 0; i < ARRAY_SIZE(sim_list); i++) {
		if (streq(config->run_type, sim_list[i].run_type)) {
			printf("%s\n\n\n", sim_list[i].description);
			sim_list[i].sim_fn(efp, config);
			return;
		}
	}

	error("UNKNOWN RUN TYPE SPECIFIED");
}

static void print_banner(void)
{
	printf("EFPMD ver. " EFPMD_VERSION "\n");
	printf("Copyright (c) 2012 Ilya Kaliman\n\n");
	printf("%s\n", efp_banner());
}

static int string_compare(const void *a, const void *b)
{
	const char *s1 = *(const char *const *)a;
	const char *s2 = *(const char *const *)b;

	return u_strcasecmp(s1, s2);
}

static char *make_name_list(const struct config *config)
{
	size_t len = 0;

	for (int i = 0; i < config->n_frags; i++)
		len += strlen(config->frags[i].name) + 1;

	char *names = xmalloc(len), *ptr = names;

	for (int i = 0; i < config->n_frags; i++) {
		strcpy(ptr, config->frags[i].name);
		ptr += strlen(config->frags[i].name);
		*ptr++ = '\n';
	}

	*(--ptr) = '\0';
	return names;
}

static char *make_potential_file_list(const struct config *config)
{
	/* This function constructs the list of library potential data files.
	 * For each unique fragment if fragment name contains an _l suffix
	 * append fraglib_path prefix and remove _l suffix. Otherwise append
	 * userlib_path prefix. Add .efp extension in both cases. */

	const char *unique[config->n_frags];

	for (int i = 0; i < config->n_frags; i++)
		unique[i] = config->frags[i].name;

	qsort(unique, config->n_frags, sizeof(char *), string_compare);

	int n_unique = 1;

	for (int i = 1; i < config->n_frags; i++) {
		if (!streq(unique[i - 1], unique[i])) {
			unique[n_unique] = unique[i];
			n_unique++;
		}
	}

	size_t len = 0;

	for (int i = 0; i < n_unique; i++) {
		const char *name = unique[i];
		size_t n = strlen(name);

		if (n > 2 && streq(name + n - 2, "_l"))
			len += strlen(config->fraglib_path) + n + 4;
		else
			len += strlen(config->userlib_path) + n + 6;
	}

	char *files = xmalloc(len), *ptr = files;

	for (int i = 0; i < n_unique; i++) {
		const char *name = unique[i];
		size_t n = strlen(name);

		if (n > 2 && streq(name + n - 2, "_l")) {
			strcat(strcpy(ptr, config->fraglib_path), "/");
			strcat(strncat(ptr, name, n - 2), ".efp");
			ptr += strlen(config->fraglib_path) + n + 3;
			*ptr++ = '\n';
		}
		else {
			strcat(strcpy(ptr, config->userlib_path), "/");
			strcat(strcat(ptr, name), ".efp");
			ptr += strlen(config->userlib_path) + n + 5;
			*ptr++ = '\n';
		}
	}

	*(--ptr) = '\0';
	return files;
}

static void init_efp(struct efp **efp, const struct config *config)
{
	enum efp_result res;

	struct efp_opts opts = {
		.terms = config->terms,
		.elec_damp = config->elec_damp,
		.disp_damp = config->disp_damp,
		.pol_damp = config->pol_damp
	};

	char *files = make_potential_file_list(config);
	char *names = make_name_list(config);

	if ((res = efp_init(efp, &opts, NULL, files, names)))
		lib_error(res);

	free(files);
	free(names);
}

static int get_coord_count(enum efp_coord_type coord_type)
{
	switch (coord_type) {
		case EFP_COORD_TYPE_XYZABC: return  6;
		case EFP_COORD_TYPE_POINTS: return  9;
		case EFP_COORD_TYPE_ROTMAT: return 12;
	};
	assert(0);
}

static void set_coord(struct efp *efp, const struct config *config)
{
	int n_coord = get_coord_count(config->coord_type);
	double *coord = xmalloc(n_coord * config->n_frags * sizeof(double));

	for (int i = 0; i < config->n_frags; i++)
		memcpy(coord + n_coord * i, config->frags[i].coord, n_coord * sizeof(double));

	enum efp_result res;

	if ((res = efp_set_coordinates(efp, config->coord_type, coord)))
		lib_error(res);

	free(coord);
}

int main(int argc, char **argv)
{
	if (argc < 2)
		die("usage: efpmd <input>");

	struct efp *efp;
	struct config config;

	print_banner();
	parse_config(argv[1], &config);
	init_efp(&efp, &config);
	set_coord(efp, &config);

	run_sim(efp, &config);
	printf("SIMULATION COMPLETED SUCCESSFULLY\n");

	efp_shutdown(efp);
	free_config(&config);

	return EXIT_SUCCESS;
}
