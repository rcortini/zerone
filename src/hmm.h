#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <err.h>

double
fwd(
   // input //
         int    m,
         int    n,
   const double *Q,
   const double *init,
   // output //
         double *prob
);

void
bwd(
   // input //
         int    m,
         int    n,
   const double *Q,
   // output //
         double *alpha,
         double *phi,
         double *T
);

double fwdb(
   // input //
      int    m,
      int    n,
const double *Q,
const double *init,
   // output //
      double *prob,
      double *phi,
   // not initialized //
      double *trans
);

void block_fwdb(
   // input //
      int *n_states,
      int *nblocks,
      int *size,
   // params //
      double *Q,
      double *init,
   // output //
      double *prob,
      double *phi,
      double *sumtrans,
      double *loglik
);

int *
viterbi(
   // input //
   int m,
   int n,
   const double *log_Q,
   const double *log_init,
   const double *log_prob,
   // output //
   int *path
);

int
block_viterbi(
   // input //
   const int *nstates,
   const int *nblocks,
   const int *size,
   const double *Q,
   const double *init,
   const double *prob,
   // control //
   const int *arglog,
   // output //
   int *path
);
