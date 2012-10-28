#include "clapack.h"

void dsyev_(char *,
	    char *,
	    int *,
	    double *,
	    int *,
	    double *,
	    double *,
	    int *,
	    int *);

int u_dsyev(char jobz, char uplo, int n, double *a, int lda,
		double *w, double *work, int lwork)
{
	int info;

	dsyev_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, &info);

	return info;
}
