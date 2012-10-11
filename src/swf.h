#ifndef LIBEFP_SWF_H
#define LIBEFP_SWF_H

#include "../common/math_util.h"

struct swf {
	double swf;
	vec_t dswf;
	vec_t cell;
};

double get_swf(double, double);
double get_dswf(double, double);

#endif /* LIBEFP_SWF_H */
