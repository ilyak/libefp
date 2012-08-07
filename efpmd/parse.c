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

#include <assert.h>
#include <ctype.h>
#include <stddef.h>
#include <stdbool.h>

#include "common.h"
#include "parse.h"

struct stream {
	char *buffer;
	char *ptr;
	FILE *in;
};

static char *read_line(FILE *in)
{
	int size = 128;
	int i = 0;
	char *buffer = malloc(size);

	if (!buffer)
		return NULL;

	for(;;) {
		int ch = getc(in);

		switch(ch) {
		case EOF:
			if (i == 0) {
				free(buffer);
				return NULL;
			}
			/* fall through */
		case '\n':
			if (i == size)
				buffer = realloc(buffer, size + 1);
			buffer[i] = '\0';
			return buffer;
		default:
			buffer[i++] = ch;

			if (i == size) {
				size *= 2;
				buffer = realloc(buffer, size);
			}
		}
	}
	return NULL;
}

static void next_line(struct stream *stream)
{
	if (stream->buffer)
		free(stream->buffer);

	stream->buffer = read_line(stream->in);
	stream->ptr = stream->buffer;

	if (stream->buffer) {
		for (char *p = stream->buffer; *p; p++)
			*p = tolower(*p);
	}
}

static void skip_space(struct stream *stream)
{
	while (*stream->ptr && isspace(*stream->ptr))
		stream->ptr++;
}

static bool parse_string(char **str, void *out)
{
	size_t len = 0;
	char *ptr = *str, *start;

	while (*ptr && isspace(*ptr))
		ptr++;

	if (*ptr == '"') {
		start = ++ptr;

		while (*ptr && *ptr != '"')
			ptr++, len++;

		if (!*ptr)
			return false;

		ptr++;
	}
	else {
		start = ptr;

		while (*ptr && !isspace(*ptr))
			ptr++, len++;
	}

	if (*(char **)out)
		free(*(char **)out);

	*(char **)out = u_strndup(start, len);
	*str = ptr;

	return true;
}

static bool parse_int(char **str, void *out)
{
	char *endptr;
	int val = strtol(*str, &endptr, 10);

	if (endptr == *str)
		return false;

	*str = endptr;
	*(int *)out = val;
	return true;
}

static bool parse_double(char **str, void *out)
{
	char *endptr;
	double val = strtod(*str, &endptr);

	if (endptr == *str)
		return false;

	*str = endptr;
	*(double *)out = val;
	return true;
}

static bool parse_coord(char **str, void *out)
{
	static const struct {
		const char *name;
		enum efp_coord_type value;
	} list[] = {
		{ "points", EFP_COORD_TYPE_POINTS },
		{ "xyzabc", EFP_COORD_TYPE_XYZABC }
	};

	for (size_t i = 0; i < ARRAY_SIZE(list); i++) {
		if (strneq(list[i].name, *str, strlen(list[i].name))) {
			*str += strlen(list[i].name);
			*(enum efp_coord_type *)out = list[i].value;
			return true;
		}
	}
	return false;
}

static bool parse_units(char **str, void *out)
{
	static const struct {
		const char *name;
		double value;
	} list[] = {
		{ "bohr", 1.0 },
		{ "angs", 1.0 / BOHR_RADIUS }
	};

	for (size_t i = 0; i < ARRAY_SIZE(list); i++) {
		if (strneq(list[i].name, *str, strlen(list[i].name))) {
			*str += strlen(list[i].name);
			*(double *)out = list[i].value;
			return true;
		}
	}
	return false;
}

static bool parse_terms(char **str, void *out)
{
	static const struct {
		const char *name;
		enum efp_term value;
	} list[] = {
		{ "elec", EFP_TERM_ELEC },
		{ "pol",  EFP_TERM_POL  },
		{ "disp", EFP_TERM_DISP },
		{ "xr",   EFP_TERM_XR   }
	};

	char *ptr = *str;
	unsigned terms = 0;

	while (*ptr) {
		for (size_t i = 0; i < ARRAY_SIZE(list); i++) {
			if (strneq(list[i].name, ptr, strlen(list[i].name))) {
				ptr += strlen(list[i].name);
				terms |= list[i].value;

				while (*ptr && isspace(*ptr))
					ptr++;
			}
			else {
				return false;
			}
		}
	}

	*(unsigned *)out = terms;
	*str = ptr;
	return terms != 0;
}

static bool parse_elec_damp(char **str, void *out)
{
	static const struct {
		const char *name;
		enum efp_elec_damp value;
	} list[] = {
		{ "screen",  EFP_ELEC_DAMP_SCREEN  },
		{ "overlap", EFP_ELEC_DAMP_OVERLAP },
		{ "off",     EFP_ELEC_DAMP_OFF     }
	};

	for (size_t i = 0; i < ARRAY_SIZE(list); i++) {
		if (strneq(list[i].name, *str, strlen(list[i].name))) {
			*str += strlen(list[i].name);
			*(enum efp_elec_damp *)out = list[i].value;
			return true;
		}
	}
	return false;
}

static bool parse_disp_damp(char **str, void *out)
{
	static const struct {
		const char *name;
		enum efp_disp_damp value;
	} list[] = {
		{ "tt",      EFP_DISP_DAMP_TT      },
		{ "overlap", EFP_DISP_DAMP_OVERLAP },
		{ "off",     EFP_DISP_DAMP_OFF     }
	};

	for (size_t i = 0; i < ARRAY_SIZE(list); i++) {
		if (strneq(list[i].name, *str, strlen(list[i].name))) {
			*str += strlen(list[i].name);
			*(enum efp_disp_damp *)out = list[i].value;
			return true;
		}
	}
	return false;
}

static bool int_gt_zero(void *val)
{
	return *(int *)val > 0;
}

static bool double_gt_zero(void *val)
{
	return *(double *)val > 0.0;
}

static const struct {
	const char *name;
	const char *default_value;
	bool (*parse_fn)(char **, void *);
	bool (*check_fn)(void *);
	size_t member_offset;
} config_list[] = {
	{ "run_type",     "sp",               parse_string,    NULL,           offsetof(struct config, run_type)     },
	{ "coord",        "xyzabc",           parse_coord,     NULL,           offsetof(struct config, coord_type)   },
	{ "units",        "angs",             parse_units,     NULL,           offsetof(struct config, units_factor) },
	{ "terms",        "elec pol disp xr", parse_terms,     NULL,           offsetof(struct config, terms)        },
	{ "elec_damp",    "screen",           parse_elec_damp, NULL,           offsetof(struct config, elec_damp)    },
	{ "disp_damp",    "tt",               parse_disp_damp, NULL,           offsetof(struct config, disp_damp)    },
	{ "max_steps",    "100",              parse_int,       int_gt_zero,    offsetof(struct config, max_steps)    },
	{ "print_step",   "1",                parse_int,       int_gt_zero,    offsetof(struct config, print_step)   },
	{ "temperature",  "300.0",            parse_double,    double_gt_zero, offsetof(struct config, temperature)  },
	{ "time_step",    "1.0",              parse_double,    double_gt_zero, offsetof(struct config, time_step)    },
	{ "ls_step_size", "20.0",             parse_double,    double_gt_zero, offsetof(struct config, ls_step_size) },
	{ "opt_tol",      "1.0e-4",           parse_double,    double_gt_zero, offsetof(struct config, opt_tol)      },
	{ "fraglib_path", ".",                parse_string,    NULL,           offsetof(struct config, fraglib_path) },
	{ "userlib_path", ".",                parse_string,    NULL,           offsetof(struct config, userlib_path) }
};

static void parse_field(struct stream *stream, struct config *config)
{
	for (size_t i = 0; i < ARRAY_SIZE(config_list); i++) {
		const char *name = config_list[i].name;
		size_t offset = config_list[i].member_offset;

		if (strneq(name, stream->ptr, strlen(name))) {
			stream->ptr += strlen(name);
			skip_space(stream);

			if (!config_list[i].parse_fn(&stream->ptr, (char *)config + offset))
				error("UNABLE TO READ OPTION %s", name);

			bool (*check_fn)(void *) = config_list[i].check_fn;

			if (check_fn && !check_fn((char *)config + offset))
				error("OPTION %s HAS OUT OF RANGE VALUE", name);

			return;
		}
	}
	error("UNKNOWN OPTION IN INPUT FILE");
}

static void parse_frag(struct stream *stream, enum efp_coord_type coord_type,
		       double units_factor, struct frag *frag)
{
	int n_rows, n_cols;
	frag->name = NULL;

	if (!parse_string(&stream->ptr, &frag->name))
		error("UNABLE TO READ FRAGMENT NAME");

	next_line(stream);

	switch (coord_type) {
		case EFP_COORD_TYPE_XYZABC:
			n_rows = 1;
			n_cols = 6;
			break;
		case EFP_COORD_TYPE_POINTS:
			n_rows = 3;
			n_cols = 3;
			break;
		default:
			lib_error(EFP_RESULT_INCORRECT_ENUM_VALUE);
	}

	for (int i = 0, idx = 0; i < n_rows; i++) {
		for (int j = 0; j < n_cols; j++, idx++) {
			double val;

			if (!parse_double(&stream->ptr, &val))
				error("UNABLE TO PARSE FRAGMENT COORDINATES");

			frag->coord[idx] = val * units_factor;
		}
		next_line(stream);
	}
}

static void set_config_defaults(struct config *config)
{
	memset(config, 0, sizeof(struct config));

	for (size_t i = 0; i < ARRAY_SIZE(config_list); i++) {
		char *str = (char *)config_list[i].default_value;
		size_t offset = config_list[i].member_offset;

		if (!config_list[i].parse_fn(&str, (char *)config + offset))
			error("INCORRECT OPTION DEFAULT VALUE");
	}
}

void parse_config(const char *path, struct config *config)
{
	set_config_defaults(config);

	struct stream stream = {
		.buffer = NULL,
		.ptr = NULL,
		.in = fopen(path, "r")
	};

	if (!stream.in)
		error("UNABLE TO OPEN INPUT FILE");

	next_line(&stream);

	while (stream.buffer) {
		if (*stream.ptr == '#')
			goto next;

		skip_space(&stream);

		if (!*stream.ptr)
			goto next;

		if (strneq(stream.ptr, "fragment", strlen("fragment"))) {
			stream.ptr += strlen("fragment");

			config->n_frags++;
			config->frags = realloc(config->frags,
				config->n_frags * sizeof(struct frag));

			parse_frag(&stream, config->coord_type, config->units_factor,
					config->frags + config->n_frags - 1);
			continue;
		}
		else {
			parse_field(&stream, config);
			skip_space(&stream);

			if (*stream.ptr)
				error("ONLY ONE OPTION PER LINE IS ALLOWED");
		}
next:
		next_line(&stream);
	}

	if (config->n_frags < 1)
		error("AT LEAST ONE FRAGMENT MUST BE SPECIFIED");

	fclose(stream.in);
	free(stream.buffer);
}

void free_config(struct config *config)
{
	free(config->run_type);
	free(config->fraglib_path);
	free(config->userlib_path);

	for (int i = 0; i < config->n_frags; i++)
		free(config->frags[i].name);

	free(config->frags);

	/* don't do free(config) here */
}
