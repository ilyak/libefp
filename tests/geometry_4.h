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

static const char *potential_files[] = {
	ABS_TOP_SRCDIR "/fraglib/acetone.efp",
	ABS_TOP_SRCDIR "/fraglib/c2h5oh.efp",
	ABS_TOP_SRCDIR "/fraglib/c6h6.efp",
	ABS_TOP_SRCDIR "/fraglib/ccl4.efp",
	ABS_TOP_SRCDIR "/fraglib/ch3oh.efp",
	ABS_TOP_SRCDIR "/fraglib/ch4.efp",
	ABS_TOP_SRCDIR "/fraglib/cl2.efp",
	ABS_TOP_SRCDIR "/fraglib/dcm.efp",
	ABS_TOP_SRCDIR "/fraglib/dmso.efp",
	ABS_TOP_SRCDIR "/fraglib/h2.efp",
	ABS_TOP_SRCDIR "/fraglib/h2o.efp",
	ABS_TOP_SRCDIR "/fraglib/nh3.efp",
	NULL
};

static const char *fragname[] = {
	"ACETONE_L",
	"C2H5OH_L",
	"C6H6_L",
	"CCL4_L",
	"CH3OH_L",
	"CH4_L",
	"CL2_L",
	"DCM_L",
	"DMSO_L",
	"H2_L",
	"H2O_L",
	"NH3_L",
	 NULL
};

static const double xyzabc[] = { /* some random geometry */
	BOHR( 0.0), BOHR( 0.0), BOHR(0.0), 0.0, 0.2, 0.3,
	BOHR( 7.0), BOHR( 0.0), BOHR(0.0), 0.0, 2.0, 3.7,
	BOHR(14.0), BOHR( 0.0), BOHR(0.0), 3.1, 0.8, 2.0,
	BOHR(21.0), BOHR( 0.0), BOHR(0.0), 0.0, 8.0, 0.0,
	BOHR( 0.0), BOHR( 6.0), BOHR(0.0), 0.7, 2.0, 1.0,
	BOHR( 7.0), BOHR( 6.0), BOHR(0.0), 0.6, 0.0, 4.7,
	BOHR(14.0), BOHR( 6.0), BOHR(0.0), 0.0, 0.0, 0.3,
	BOHR(21.0), BOHR( 6.0), BOHR(0.0), 0.0, 0.4, 0.3,
	BOHR( 0.0), BOHR(12.0), BOHR(0.0), 0.8, 0.0, 0.0,
	BOHR( 7.0), BOHR(12.0), BOHR(0.0), 8.0, 0.7, 0.8,
	BOHR(14.0), BOHR(12.0), BOHR(0.0), 0.0, 0.0, 0.0,
	BOHR(21.0), BOHR(12.0), BOHR(0.0), 0.0, 2.0, 0.0
};
