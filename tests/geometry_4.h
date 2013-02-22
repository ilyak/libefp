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

static const char files[] =
	"./fraglib/acetone.efp\n"
	"./fraglib/c2h5oh.efp\n"
	"./fraglib/c6h6.efp\n"
	"./fraglib/ccl4.efp\n"
	"./fraglib/ch3oh.efp\n"
	"./fraglib/ch4.efp\n"
	"./fraglib/cl2.efp\n"
	"./fraglib/dcm.efp\n"
	"./fraglib/dmso.efp\n"
	"./fraglib/h2.efp\n"
	"./fraglib/h2o.efp\n"
	"./fraglib/nh3.efp";

static const char names[] =
	"ACETONE_L\n"
	"C2H5OH_L\n"
	"C6H6_L\n"
	"CCL4_L\n"
	"CH3OH_L\n"
	"CH4_L\n"
	"CL2_L\n"
	"DCM_L\n"
	"DMSO_L\n"
	"H2_L\n"
	"H2O_L\n"
	"NH3_L";

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
