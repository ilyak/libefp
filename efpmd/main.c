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
#include "parse.h"
#include "sim.h"

typedef void (*sim_fn_t)(struct efp *, const struct config *);

static sim_fn_t get_sim_fn(const char *run_type)
{
	static const struct {
		const char *run_type;
		sim_fn_t sim_fn;
	} sim_list[] = {
		{ "sp",   sim_sp   },
		{ "grad", sim_grad },
		{ "cg",   sim_cg   },
		{ "nve",  sim_nve  },
		{ "nvt",  sim_nvt  }
	};

	for (size_t i = 0; i < ARRAY_SIZE(sim_list); i++)
		if (streq(run_type, sim_list[i].run_type))
			return sim_list[i].sim_fn;

	return NULL;
}

static void print_banner(void)
{
	printf("EFPMD ver. " EFPMDVERSION "\n");
	printf("Copyright (c) 2012 Ilya Kaliman\n\n");
	printf("%s\n\n", efp_banner());
}

static int string_compare(const void *a, const void *b)
{
	const char *s1 = *(const char **)a;
	const char *s2 = *(const char **)b;

	return u_strcasecmp(s1, s2);
}

static char **make_potential_file_list(const struct config *config)
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

	char **list = malloc((n_unique + 1) * sizeof(char *));

	for (int i = 0; i < n_unique; i++) {
		const char *name = unique[i];
		size_t len = strlen(name);
		char *path;

		if (len > 2 && streq(name + len - 2, "_l")) {
			path = malloc(strlen(config->fraglib_path) + len + 4);
			strcat(strcpy(path, config->fraglib_path), "/");
			strcat(strncat(path, name, len - 2), ".efp");
		}
		else {
			path = malloc(strlen(config->userlib_path) + len + 6);
			strcat(strcpy(path, config->userlib_path), "/");
			strcat(strcat(path, name), ".efp");
		}

		list[i] = path;
	}

	list[n_unique] = NULL;
	return list;
}

static void init_efp(struct efp **efp, const struct config *config)
{
	struct efp_opts opts = {
		.terms = config->terms,
		.elec_damp = config->elec_damp,
		.disp_damp = config->disp_damp
	};

	char **files = make_potential_file_list(config);
	char *names[config->n_frags + 1];

	for (int i = 0; i < config->n_frags; i++)
		names[i] = config->frags[i].name;

	names[config->n_frags] = NULL;
	enum efp_result res;

	if ((res = efp_init(efp, &opts, NULL, (const char **)files, (const char **)names)))
		lib_error(res);

	for (int i = 0; files[i]; i++)
		free(files[i]);

	free(files);
}

static int get_coord_count(enum efp_coord_type coord_type)
{
	switch (coord_type) {
		case EFP_COORD_TYPE_XYZABC: return 6;
		case EFP_COORD_TYPE_POINTS: return 9;
	};

	lib_error(EFP_RESULT_INCORRECT_ENUM_VALUE);
}

static void set_coord(struct efp *efp, const struct config *config)
{
	int n_coord = get_coord_count(config->coord_type);
	double *coord = malloc(n_coord * config->n_frags * sizeof(double));

	for (int i = 0; i < config->n_frags; i++)
		memcpy(coord + n_coord * i, config->frags[i].coord,
						n_coord * sizeof(double));

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

	sim_fn_t sim_fn = get_sim_fn(config.run_type);

	if (!sim_fn)
		error("UNKNOWN RUN TYPE SPECIFIED");

	sim_fn(efp, &config);

	efp_shutdown(efp);
	free_config(&config);

	puts("EFP SIMULATION COMPLETED SUCCESSFULLY");
	return EXIT_SUCCESS;
}
