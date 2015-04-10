/*-
 * Copyright (c) 2012-2015 Ilya Kaliman
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

#include "msg.h"

#ifdef WITH_MPI
#include <mpi.h>
#endif

void msg(const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	vfmsg(stdout, fmt, ap);
	va_end(ap);
}

void fmsg(FILE *st, const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);
	vfmsg(st, fmt, ap);
	va_end(ap);
}

void vmsg(const char *fmt, va_list ap)
{
	vfmsg(stdout, fmt, ap);
}

void vfmsg(FILE *st, const char *fmt, va_list ap)
{
#ifdef WITH_MPI
	int rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		vfprintf(st, fmt, ap);
	}
#else
	vfprintf(st, fmt, ap);
#endif
}
