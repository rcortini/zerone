#include <time.h>
#include <stdlib.h>
#include <pthread.h>
#include "utils.h"
#include "mnmultinom.h"
#include "hmm.h"


#ifndef _PSO_HMM
#define _PSO_HMM

void
pso
(
   // input //
   const int *n_states,
   const int *n_obs,
   const int *dim_yz,
   const int *yz,
   // fixed params //
   const double *a,
   const double *t,
   // start conditions //
         double *Q,
         double *p,
         double *q,
   // index //
         int *index,
   // output //
   double *loglikmax
);

#endif
