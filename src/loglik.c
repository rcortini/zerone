#include "loglik.h"

int
max3 (
   double array[3]
)
{
   int argmax = 0;
   if (array[1] > array[argmax]) argmax = 1;
   if (array[2] > array[argmax]) argmax = 2;
   return argmax;
}

void
compute_loglik (
   // input //
   int *n_obs,
   int *dim_z,
   int *yz,
   // params //
   double *init,
   double *Q,
   double *a,
   double *b,
   double *gamma,
   // index //
   int *index,
   int *tabulated,
   // output //
   double *loglik
)
// SYNOPSIS:                                                             
//   Compute emission probabilities with a negative multinomial model.   
//   Since those are up to a multiplicative constant in the forward-     
//   backward algorithm, we can compute the ratios of the emission       
//   probabilities relative to a reference state, of which we set the    
//   emission probability to 1. This saves computation for one state.    
//                                                                       
//      Since the negative nultinomial takes discrete values, the time   
//   series takes a limited number of values. This means that we can     
//   cache the results for reuse in order to save computation.           
//                                                                       
//      My parametrization is of the form                                
//                                                                       
//       p_0(i)^alpha * p_1(i)^y * p_2(i)^z_1 * ... * p_r+1(i)^z_r       
//                                                                       
//      By taking the ratio relative to state i=3 the terms to compute   
//   are of the form                                                     
//                                                                       
//       q_0(i)^alpha * q_1(i)^y * q_2(i)^z_1 * ... * q_r+1(i)^z_r       
//                                                                       
//   where q_0(i) = p_0(i)/p_0(3) etc.                                   
//                                                                       
// ARGUMENTS:                                                            
//   'n_obs': (1) length of the sequence of observations                 
//   'dim_z': (1) dimension of z (number of profiles).                   
//   'y': (n_obs) control profile.                                       
//   'z': (n_obs,dim_z) profiles.                                        
//   'a': (1) alias 'alpha', model parameter                             
//   'b': (1) alias 'beta', model parameter                              
//   'gamma': (dim_z,3) model parameter                                  
//                                                                       
// RETURN:                                                               
//   'void'                                                              
{

   int i;
   int j;
   int k;
   int n = *n_obs;
   int r = *dim_z;
   // The following aliases are to avoid the statement '1/*beta'
   // which is interpreted as a comment.
   double alpha = *a;
   double beta = *b;

   double sum;
   double tmp[3];
   double phi[3] = {init[0], init[1], init[2]};

   // Compute p's.
   double logp[(r+2)*3];
   for (i = 0 ; i < 3 ; i++) {
      double denom = 1 + 1/beta;
      for (j = 0 ; j < r ; j++) denom += gamma[j+i*r];
      for (j = 0 ; j < r ; j++) {
         logp[(j+2)+i*(r+2)] = log(gamma[j+i*r] / denom);
      }
      // Fill in p_0 and p_1.
      logp[0+i*(r+2)] = log(1/beta / denom);
      logp[1+i*(r+2)] = log(1 / denom);
   }

   if (*index == -1) {
      // Index the time series, assuming that 'index' has been
      // properly allocated.
      indexts(n, r+1, yz, index);
   }

   double *pem = malloc(3*n * sizeof(double));
   *loglik = 0.0;

   for (k = 0 ; k < n ; k++) {
      tmp[0] = tmp[1] = tmp[2] = 0.0;
      for (i = 0 ; i < 3 ; i++) {
      for (j = 0 ; j < 3 ; j++) {
         tmp[j] += phi[i] * Q[i+3*j];
      }
      }

      // NAs of type 'int' is a large negative value. Set ratio to
      // 1.0 in case of NA emission (assuming all states have the
      // same probability of producing NAs).
      int is_na = 0;
      for (i = 0 ; i < r+1 ; i++) {
         if (yz[i+k*(r+1)] < 0) {
            is_na = 1;
            break;
         }
      }
      if (is_na) {
         for (i = 0 ; i < 3 ; i++) phi[i] = tmp[i];
         continue;
      }

      // Caching by indexing.
      if (index[k] == k) {
         // Compute log-ratio.
         double q0 = alpha * logp[0];
         double q1 = alpha * logp[0+(r+2)*1];
         double q2 = alpha * logp[0+(r+2)*2];
         for (i = 0 ; i < r+1; i++) {
            q0 += yz[i+k*(r+1)] * logp[i+1];
            q1 += yz[i+k*(r+1)] * logp[i+1+(r+2)*1];
            q2 += yz[i+k*(r+1)] * logp[i+1+(r+2)*2];
         }

         // Take exponential.
         pem[3*k  ] = exp(q0);
         pem[3*k+1] = exp(q1);
         pem[3*k+2] = exp(q2);

         // NB: In case of numerical underflow, we store the
         // log probability of emission, which are negative.
         if (!(pem[3*k]+pem[3*k+1]+pem[3*k+2] > 0)) {
            pem[3*k  ] = q0;
            pem[3*k+1] = q1;
            pem[3*k+2] = q2;
         }

      }

      tmp[0] *= pem[3*index[k]  ];
      tmp[1] *= pem[3*index[k]+1];
      tmp[2] *= pem[3*index[k]+2];

      sum = tmp[0] + tmp[1] + tmp[2];
      if (sum > 0) *loglik += log(sum);
      else {
         int argmax = max3(pem+3*index[k]);
         *loglik += pem[3*index[k]+argmax];
      }
      // Break out if computation is unstable.
      if (*loglik != *loglik) break;
      for (i = 0 ; i < 3 ; i++) phi[i] = tmp[i] / sum;

   }

   n = 0;
   for (i = 0 ; tabulated[i] > -1 ; i++) {
      *loglik += tabulated[i] * lgamma(alpha+i);
      n += tabulated[i];
   }
   *loglik -= n * lgamma(alpha);

   free(pem);
   return;
}

