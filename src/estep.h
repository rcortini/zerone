#include <math.h>
#include <float.h>
#include <stdlib.h>

#ifndef _ESTEP_HMMNB
#define _ESTEP_HMMNB
#define EPSILON 2.2e-16
void compute_pratio(int *, int *, int *, double *, double *,
      double *, int *, double *);
void fwdb(int *, int *, double *, double *, double *,
      double *, double *);
#endif
