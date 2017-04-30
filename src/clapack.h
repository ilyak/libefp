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

#ifndef LIBEFP_CLAPACK_H
#define LIBEFP_CLAPACK_H

#ifdef LIBEFP_FORTRAN_INT64
typedef long long int fortranint_t;
#else /* LIBEFP_FORTRAN_INT64 */
typedef int fortranint_t;
#endif /* LIBEFP_FORTRAN_INT64 */

void efp_dgemm(char,
	       char,
	       fortranint_t,
	       fortranint_t,
	       fortranint_t,
	       double,
	       double *,
	       fortranint_t,
	       double *,
	       fortranint_t,
	       double,
	       double *,
	       fortranint_t);

fortranint_t efp_dsyev(char,
		       char,
		       fortranint_t,
		       double *,
		       fortranint_t,
		       double *);

fortranint_t efp_dgesv(fortranint_t,
		       fortranint_t,
		       double *,
		       fortranint_t,
		       fortranint_t *,
		       double *,
		       fortranint_t);

#endif /* LIBEFP_CLAPACK_H */
