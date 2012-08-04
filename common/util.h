#ifndef LIBEFP_COMMON_UTIL_H
#define LIBEFP_COMMON_UTIL_H

#include <stddef.h>

size_t u_strnlen(const char *, size_t);
char *u_strdup(const char *);
char *u_strndup(const char *, size_t);
int u_strcasecmp(const char *, const char *);
int u_strncasecmp(const char *, const char *, size_t);

#endif /* LIBEFP_COMMON_UTIL_H */
