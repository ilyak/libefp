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

#ifndef EFPMD_PHYS_H
#define EFPMD_PHYS_H

/* Bohr radius in angstroms */
#define BOHR_RADIUS 0.52917721092

/* Boltzmann constant in [Hartree / K] */
#define BOLTZMANN 3.166811429e-6

/* Femtoseconds to atomic units of time conversion */
#define FS_TO_AU (1.0 / 2.41888432650212e-2)

/* AMU to atomic units of mass conversion */
#define AMU_TO_AU (1.0 / 5.485799094622e-4)

/* Hertree energy in Joules */
#define HARTREE 4.35974434e-18

/* Bar to atomic units of pressure */
#define BAR_TO_AU (1.0e-25 * BOHR_RADIUS * BOHR_RADIUS * BOHR_RADIUS / HARTREE)

/* Fine structure constant */
#define FINE_CONST 7.297352569824e-3

#endif /* EFPMD_PHYS_H */
