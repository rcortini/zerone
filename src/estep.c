#include <math.h>
#include <float.h>
#include <pthread.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#define M_LOOK 10000

void spcfwdb(
   // input //
   const int n,
   const double *Q,
   const double init,
   // output //
   double *pratio,
   double *phi,
   // not initialized //
   double *sumtrans
){
// SYNOPSIS:                                                             
//   Special forward-backward algorithm for 2-state chains where only    
//   the ratio of emission probabilities is given.                       
//                                                                       
// ARGUMENTS:                                                            
//   'n': length of the sequence of observations                         
//   'Q': (2) diagonal of the transition matrix                          
//   'init': initial probability of the first state                      
//   'r atio': (n) ratio of emission probabilities                       
//   'phi': (n) probabilities of the first state given observations      
//   'sumtrans': (2,2) sum of conditional transitions probabilties       
//                                                                       
// RETURN:                                                               
//   'void'                                                              
//                                                                       
// SIDE EFFECTS:                                                         
//   Replace 'pratio' by alphas, update 'phi' and 'sumtrans' in place.   

   int k;
   double a;                // Current forward alpha for the first state.
   double tmp[2];           // Computational intermediate.
   double rev[2];           // Diagonal of the reverse kernel.
   double *alpha = pratio;  // For clarity.

   // First iteration of the forward pass.
   a = alpha[0] = (init * pratio[0]) / (init * pratio[0] + 1.0 - init);

   // Next iterations of the forward pass.
   for (k = 1 ; k < n ; k++) {
      tmp[0] = a * Q[0] + (1-a) * (1-Q[1]);
      tmp[1] = a * (1-Q[0]) + (1-a) * Q[0];
      a = alpha[k] = (tmp[0]*pratio[k]) / (tmp[0]*pratio[k]+tmp[1]);
   }

   // NB: 'ratio' now contains the forward alphas for the first state.
   // However, we refer to it as 'alpha' for clarity.
   // First iteration of the backward pass.
   phi[n-1] = alpha[n-1];

   // Next iterations of the backward pass.
   for (k = n-2 ; k >= 0 ; k--) {
      rev[0] = alpha[k] * Q[0];
      rev[0] /= rev[0] + (1 - alpha[k]) * (1 - Q[1]);
      rev[1] = (1 - alpha[k]) * Q[1];
      rev[1] /= rev[1] + alpha[k] * (1 - Q[0]);

      tmp[0] = phi[k+1] * rev[0];
      tmp[1] = (1 - phi[k+1]) * (1 - rev[1]);

      phi[k] = tmp[0] + tmp[1];
      sumtrans[0] += tmp[0];
      sumtrans[1] += tmp[1];
      sumtrans[2] += phi[k+1] - tmp[0];
      sumtrans[3] += (1 - phi[k+1]) - tmp[1];
   }

   return;

}

void compute_pratio(
   // input //
   int *n,
   int *x,
   int *y,
   // params //
   double *alpha,
   double *beta,
   double *gamma,
   // output //
   double *pratio
){

   int i;

   double _r = (1 + 1 / *beta + gamma[1]) / (1 + 1 / *beta + gamma[0]);
   double r = gamma[0] * _r / gamma[1];

   double *lx = (double *) malloc (M_LOOK * sizeof(double));
   double *ly = (double *) malloc (M_LOOK * sizeof(double));
   for (i = 0 ; i < M_LOOK ; i++) lx[i] = ly[i] = -1.0;

   for (i = 0 ; i < *n ; i++) {
      // NAs are passed as negative values to 'int'. Set ratio to
      // 1.0 in case of NA emission (assuming both states have the
      // same probability of producing NAs).
      if (x[i] < 0 || y[i] < 0) {
         pratio[i] = 1.0;
         continue;
      }
      if (lx[x[i]] < 0) lx[x[i]] = pow(r, x[i]);
      if (ly[y[i]] < 0) ly[y[i]] = pow(_r, *alpha + y[i]);
      else pratio[i] = lx[x[i]]*ly[y[i]];
   }

   free(lx);
   free(ly);

   return;
}

void fwdb(
   // input //
   int *nblocks,
   int *lengths,
   // params //
   double *Q,
   double *init,
   // output //
   double *pratio,
   double *phi,
   double *sumtrans
){

   int i;
   int start = 0;

   for (i = 0 ; i < 2*2 ; i++) sumtrans[i] = 0.0;
   for (i = 0 ; i < *nblocks ; i++) {
      spcfwdb(lengths[i], Q, *init, pratio+start, phi+start, sumtrans);
      start += lengths[i];
   }

   return;
}
