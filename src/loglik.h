#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include <float.h>
#include "utils.h"

#ifndef _HMMNB_LOGLIK
#define _HMMNB_LOGLIK

typedef struct {
   double alpha;
   double beta;
   double Q[9];
   double *gammas;
} param_set;

void compute_loglik (const int *, const int *, const int *,
      const double *, const double *, const double *, const double *,
      int *, const int *, double *);
void simAnneal (const int *, const int *, const int *, double *);
#endif
