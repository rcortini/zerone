#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifndef _PREDICT
#define _PREDICT

#include "zerone.h"

double * extractfeat(ChIP_t *, zerone_t *);
double * zscale(double *);
double   predict(ChIP_t *, zerone_t *);

#endif
