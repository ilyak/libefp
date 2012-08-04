#include "common.h"
#include "opt.h"

#define MAX_LS_ITER 20

//static double dot(int n, const double *a, const double *b)
//{
//	double sum = 0.0;
//
//	while (n--)
//		sum += (*a++) * (*b++);
//
//	return sum;
//}

double opt_cg_step(opt_func_t f, int n, double *x, double *dfx, void *data)
{
//	int iter;
//
//	double a, fa, xa[n], ga[n];
//	double b, fb, xb[n], gb[n];
//	double c, fc, xc[n], gc[n];
//
//	double init_slope, ftol, fx;
//
//	do {
//		fx = f(n, x, dfx, data);
//
//		if (fb <= fx + init_slope * ftol * b) {
//			memcpy(x, xb, n * sizeof(double));
//			memcpy(dfx, gb, n * sizeof(double));
//			fx = fb;
//		}
//
//		iter++;
//	} while ((fb > fa || fb > fc) && (iter < MAX_LS_ITER));
//
//	if (fb < fx || iter >= MAX_LS_ITER)
//		error("UNABLE TO FIND A POINT WITH LOWER ENERGY");
//
//	if (fa <= fc) {
//		memcpy(x, xa, n * sizeof(double));
//		memcpy(dfx, ga, n * sizeof(double));
//		return fa;
//	}
//	else {
//		memcpy(x, xc, n * sizeof(double));
//		memcpy(dfx, gc, n * sizeof(double));
//		return fc;
//	}

	return 0.0;
}
