#include <math.h>
#include <float.h>
#include <pthread.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

void spcfwdb(
   // input //
   const int n,
   const double *Q,
   const double *init,
   // output //
   double *pratio,
   double *phi,
   // not initialized //
   double *sumtrans
){
// SYNOPSIS:                                                             
//   Special forward-backward algorithm for 3-state chains where only    
//   the ratio of emission probabilities is given.                       
//                                                                       
// ARGUMENTS:                                                            
//   'n': length of the sequence of observations                         
//   'Q': (3,3) diagonal of the transition matrix                        
//   'init': (3) initial probability of the first state                  
//   'pratio': (2,n) ratio of emission probabilities                     
//   'phi': (2,n) smoothed probabilities for the first 2 states          
//   'sumtrans': (3,3) sum of conditional transitions probabilties       
//                                                                       
// ENCODING OF Q.                                                        
//   Q(i,j) is a transition from state i to state j, so Q[0] is the      
//   transition from state 0 to state 0, Q[1] is the transition from     
//   state 1 to state 0, Q[3] is the transition from state 0 to state    
//   1 etc. (the general term Q[i+3*j] is the transition from state i    
//   to state j).                                                        
//                                                                       
// RETURN:                                                               
//   'void'                                                              
//                                                                       
// SIDE EFFECTS:                                                         
//   Replace 'pratio' by alphas, update 'phi' and 'sumtrans' in place.   

   int j;
   int k;
   double sum;              // Computational intermediate.
   double tmp[3];           // Computational intermediate.
   double B[9];             // Diagonal of the reverse kernel.
   double *alpha = pratio;  // For clarity.

   // First iteration of the forward pass.
   sum = init[0]*pratio[0] + init[1]*pratio[1] + init[2];
   alpha[0] = init[0]*pratio[0] / sum;
   alpha[1] = init[1]*pratio[1] / sum;
   // a[2] = alpha[2] = init[2] / tmp;

   // Next iterations of the forward pass.
   for (k = 1 ; k < n ; k++) {
      for (j = 0 ; j < 3 ; j++) {
         tmp[j] = Q[  3*j] * alpha[2*k-2] + 
                  Q[1+3*j] * alpha[2*k-1] +
                  Q[2+3*j] * (1 - alpha[2*k-2] - alpha[2*k-1]);
      }
      sum = tmp[0]*pratio[2*k] + tmp[1]*pratio[2*k+1] + tmp[2];
      if (sum > 0) {
         alpha[2*k  ] = tmp[0]*pratio[2*k  ] / sum;
         alpha[2*k+1] = tmp[1]*pratio[2*k+1] / sum;
      }
      else {
         alpha[2*k+1] = alpha[2*k] = 1.0/3;
      }
   }

   // NB: 'pratio' now contains the forward alphas for the first state.
   // However, we refer to it as 'alpha' for clarity.
   // First iteration of the backward pass.
   phi[2*n-2] = alpha[2*n-2];
   phi[2*n-1] = alpha[2*n-1];

   // Next iterations of the backward pass.
   for (k = n-2 ; k >= 0 ; k--) {
      // The backward recursion kernel B(j,i) gives the probability of
      // a backward transition from state j to state i.
      for (j = 0 ; j < 3 ; j++) {

         sum = Q[  3*j] * alpha[2*k] +
               Q[1+3*j] * alpha[2*k+1] +
               Q[2+3*j] * (1 - alpha[2*k] - alpha[2*k+1]);

         if (sum > 0) {
            B[j  ] = Q[  3*j] * alpha[2*k] / sum;
            B[j+3] = Q[1+3*j] * alpha[2*k+1] / sum;
            B[j+6] = Q[2+3*j] * (1 - alpha[2*k] - alpha[2*k+1]) / sum;
         }
         else {
            B[j] = B[j+3] = B[j+6] = 1.0/3;
         }

         sumtrans[  3*j] = B[j  ] * phi[2*k+2];
         sumtrans[1+3*j] = B[j+3] * phi[2*k+3];
         sumtrans[2+3*j] = B[j+6] * (1 - phi[2*k+2] - phi[2*k+3]);

      }

      phi[2*k]   = B[0] * phi[2*k+2] +
                   B[1] * phi[2*k+3] +
                   B[2] * (1 - phi[2*k+2] - phi[2*k+3]);
      phi[2*k+1] = B[3] * phi[2*k+2] +
                   B[4] * phi[2*k+3] +
                   B[5] * (1 - phi[2*k+2] - phi[2*k+3]);

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
   double *b,
   double *gamma,
   // output //
   double *pratio
){

   int k;
   double beta = *b;

   int max_x = -1;
   int max_y = -1;
   for (k = 0 ; k < *n ; k++) {
      if (x[k] > max_x) max_x = x[k];
      if (y[k] > max_y) max_y = y[k];
   }

   double *l1x = malloc ((max_x+1) * sizeof(double));
   double *l2x = malloc ((max_x+1) * sizeof(double));
   double *l1y = malloc ((max_y+1) * sizeof(double));
   double *l2y = malloc ((max_y+1) * sizeof(double));
   for (k = 0 ; k < max_x+1 ; k++) l1x[k] = l2x[k] = -1.0;
   for (k = 0 ; k < max_y+1 ; k++) l1y[k] = l2y[k] = -1.0;

   // With our parametrization, the negative binomial is written        
   //                                                                   
   //     C * gamma^x * 1/beta^y / (1 + 1/beta + gamma)^(alpha+x+y)     
   //                                                                   
   // where C is a constant that does not depend on (x,y). We directly  
   // take the ratio with the probability of the first state, which     
   // saves computation (the forward probabilities are normalized to    
   // sum out to 1, so we can skip C and we need to compute only m-1    
   // ratios, where m is the number of states).

   // Do the computation in log space (compute log ratios) to avoid
   // numeric overflow.
   double _r0 = log(1 + 1/beta + gamma[2]) - log(1 + 1/beta + gamma[0]);
   double _r1 = log(1 + 1/beta + gamma[2]) - log(1 + 1/beta + gamma[1]);
   double r0 = log(gamma[0]) + _r0 - log(gamma[2]);
   double r1 = log(gamma[1]) + _r1 - log(gamma[2]);
   for (k = 0 ; k < *n ; k++) {
      // NAs are passed as negative values to 'int'. Set ratio to
      // 1.0 in case of NA emission (assuming all states have the
      // same probability of producing NAs).
      if (x[k] < 0 || y[k] < 0) {
         pratio[2*k] = 1.0;
         pratio[2*k+1] = 1.0;
         continue;
      }
      if (l1x[x[k]] < 0) {
         l1x[x[k]] = r0 * x[k];
         l2x[x[k]] = r1 * x[k];
      }
      if (l1y[y[k]] < 0) {
         l1y[y[k]] = _r0 * (*alpha + y[k]);
         l2y[y[k]] = _r1 * (*alpha + y[k]);
      }
      // Take the exponential.
      pratio[2*k] = exp(l1x[x[k]] + l1y[y[k]]);
      pratio[2*k+1] = exp(l2x[x[k]] + l2y[y[k]]);

   }

   free(l1x);
   free(l2x);
   free(l1y);
   free(l2y);

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

   for (i = 0 ; i < 3*3 ; i++) sumtrans[i] = 0.0;
   for (i = 0 ; i < *nblocks ; i++) {
      spcfwdb(lengths[i], Q, init, pratio+start, phi+start, sumtrans);
      start += 2*lengths[i];
   }

   return;
}
