/*-
 * Copyright (c) 2012-2017 Ilya Kaliman
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

#include <stdlib.h>

#include "clapack.h"

void dgemm_(char *,
	    char *,
	    fortranint_t *,
	    fortranint_t *,
	    fortranint_t *,
	    double *,
	    double *,
	    fortranint_t *,
	    double *,
	    fortranint_t *,
	    double *,
	    double *,
	    fortranint_t *);

void dsyev_(char *,
	    char *,
	    fortranint_t *,
	    double *,
	    fortranint_t *,
	    double *,
	    double *,
	    fortranint_t *,
	    fortranint_t *);

void dgesv_(fortranint_t *,
	    fortranint_t *,
	    double *,
	    fortranint_t *,
	    fortranint_t *,
	    double *,
	    fortranint_t *,
	    fortranint_t *);

void
efp_dgemm(char transa, char transb, fortranint_t m, fortranint_t n,
    fortranint_t k, double alpha, double *a, fortranint_t lda, double *b,
    fortranint_t ldb, double beta, double *c, fortranint_t ldc)
{
	dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda,
	    b, &ldb, &beta, c, &ldc);
}

fortranint_t
efp_dsyev(char jobz, char uplo, fortranint_t n, double *a, fortranint_t lda,
    double *w)
{
	fortranint_t info, lwork;
	double *work;

	lwork = n * n;
	work = (double *)malloc(lwork * sizeof *work);

	dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);

	free(work);
	return (info);
}

fortranint_t
efp_dgesv(fortranint_t n, fortranint_t nrhs, double *a, fortranint_t lda,
    fortranint_t *ipiv, double *b, fortranint_t ldb)
{
	fortranint_t info;

	dgesv_(&n, &nrhs, a, &lda, ipiv, b, &ldb, &info);

	return (info);
}
