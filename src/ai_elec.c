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

#include "efp_private.h"

#define DR_A      (dr[a])
#define DR_B      (dr[b])
#define DR_C      (dr[c])

#define D1_A      (pt_i->dipole[a])
#define Q1_AB     (pt_i->quadrupole[t2_idx(a, b)])
#define O1_ABC    (pt_i->octupole[t3_idx(a, b, c)])

#define SUM_A(x)  (a = 0, tmp_a = x, a = 1, tmp_a += x, a = 2, tmp_a += x)
#define SUM_B(x)  (b = 0, tmp_b = x, b = 1, tmp_b += x, b = 2, tmp_b += x)
#define SUM_C(x)  (c = 0, tmp_c = x, c = 1, tmp_c += x, c = 2, tmp_c += x)

static double
compute_ai_pt(struct efp *efp, int frag_idx, int pt_idx, int qm_idx)
{
	struct frag *fr_i = efp->frags + frag_idx;
	struct multipole_pt *pt_i = fr_i->multipole_pts + pt_idx;
	struct efp_qm_atom *qm_atom = efp->qm_data.atoms + qm_idx;

	double dr[3] = {
		qm_atom->x - pt_i->x,
		qm_atom->y - pt_i->y,
		qm_atom->z - pt_i->z
	};

	double r = vec_len(VEC(dr[0]));

	double ri[8];
	powers(1.0 / r, 8, ri);

	double energy = 0.0;

	/* used in SUM_(ABC) */
	int a, b, c;
	double tmp_a, tmp_b, tmp_c;

	energy += ri[1] * qm_atom->charge * pt_i->monopole;
	energy += ri[3] * qm_atom->charge * SUM_A(D1_A * DR_A);
	energy += ri[5] * qm_atom->charge * SUM_A(SUM_B(Q1_AB * DR_A * DR_B));
	energy += ri[7] * qm_atom->charge * SUM_A(SUM_B(SUM_C(
						O1_ABC * DR_A * DR_B * DR_C)));

	if (efp->grad) {
	}
	return energy;
}

static double
compute_ai_frag(struct efp *efp, int frag_idx)
{
	struct frag *fr_i = efp->frags + frag_idx;
	double energy = 0.0;

	for (int i = 0; i < fr_i->n_multipole_pts; i++)
		for (int j = 0; j < efp->qm_data.n_atoms; j++)
			energy += compute_ai_pt(efp, frag_idx, i, j);

	return energy;
}

enum efp_result
efp_compute_ai_elec(struct efp *efp)
{
	if (efp->grad)
		return EFP_RESULT_NOT_IMPLEMENTED;

	double energy = 0.0;

	for (int i = 0; i < efp->n_frag; i++)
		energy += compute_ai_frag(efp, i);

	efp->energy[efp_get_term_index(EFP_TERM_AI_ELEC)] = energy;
	return EFP_RESULT_SUCCESS;
}
