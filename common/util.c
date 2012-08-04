#include <ctype.h>
#include <stdlib.h>
#include <string.h>

#include "util.h"

size_t u_strnlen(const char *str, size_t n)
{
	size_t len;

	for (len = 0; *str && len < n; len++)
		str++;

	return len;
}

char *u_strdup(const char *str)
{
	size_t len;
	char *copy;

	len = strlen(str) + 1;
	copy = malloc(len);

	if (copy)
		memcpy(copy, str, len);

	return copy;
}

char *u_strndup(const char *str, size_t n)
{
	size_t len;
	char *copy;

	len = u_strnlen(str, n);
	copy = malloc(len + 1);

	if (copy) {
		memcpy(copy, str, len);
		copy[len] = '\0';
	}

	return copy;
}

int u_strcasecmp(const char *s1, const char *s2)
{
	while (tolower(*s1) == tolower(*s2++))
		if (*s1++ == '\0')
			return 0;

	return tolower(*s1) - tolower(*--s2);
}

int u_strncasecmp(const char *s1, const char *s2, size_t n)
{
	if (n != 0) {
		do {
			if (tolower(*s1) != tolower(*s2++))
				return tolower(*s1) - tolower(*--s2);
			if (*s1++ == '\0')
				break;
		} while (--n != 0);
	}
	return 0;
}
