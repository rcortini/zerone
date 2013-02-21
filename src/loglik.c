#include <stdlib.h>
#include <math.h>

void
compute_loglik(
   int *n,
   int *y,
   int *z,
   double *init,
   double *Q,
   double *alpha,
   double *beta,
   double *g,
   double *l
){

// With our parametrization, the negative binomial is written        
//                                                                   
// Gamma(a+y+z)/(Gamma(a)*b^a) * g^z/(1+1/b+g)^{a+y+z} / y!z!
//
// We skip the factorial terms y! and z! because they are constant
// for every model. The term Gamma(a+y+z)/(Gamma(a)*b^a) is computed
// last so only g^z/(1+1/b+g)^{a+y+z} is computed inside the loop.
// The computation is carried out in log space and we cache the terms
// g^z and (1+1/b+g)^{a+y+z}.

   int i;
   int j;
   int k;

   double a = *alpha;
   double b = *beta;

   int max_y = -1;
   int max_z = -1;
   for (k = 0 ; k < *n ; k++) {
      if (y[k] > max_y) max_y = y[k];
      if (z[k] > max_z) max_z = z[k];
   }

   // Caching and counting.
   double *g0 = malloc ((max_z+1) * sizeof(double));
   double *g1 = malloc ((max_z+1) * sizeof(double));
   double *g2 = malloc ((max_z+1) * sizeof(double));
   double *bg0 = malloc ((max_y+max_z+1) * sizeof(double));
   double *bg1 = malloc ((max_y+max_z+1) * sizeof(double));
   double *bg2 = malloc ((max_y+max_z+1) * sizeof(double));
   int *counter = malloc((max_y+max_z+1) * sizeof(int));

   // Forward alphas.
   double phi[3] = {init[0], init[1], init[3]};
   double tmp[3];

   for (k = 0 ; k < max_z+1 ; k++) g0[k] = g1[k] = g2[k] = -1.0;
   for (k = 0 ; k < max_y+max_z+1 ; k++) {
      bg0[k] = bg1[k] = bg2[k] = -1.0;
      counter[k] = 0;
   }

   *l = 0.0;

   for (k = 0 ; k < *n ; k++) {
      // NAs are passed as negative values to 'int'. Set ratio to
      // 1.0 in case of NA emission (assuming all states have the
      // same probability of producing NAs).
      if (y[k] < 0 || z[k] < 0) {
         continue;
      }

      int zk = z[k];
      int yzk = y[k]+z[k];
      counter[yzk]++;

      if (g0[zk] < 0) {
         g0[zk] = zk * log(g[0]);
         g1[zk] = zk * log(g[1]);
         g2[zk] = zk * log(g[2]);
      }
      if (bg0[yzk] < 0) {
         bg0[yzk] = (a+yzk) * log(1+1/b+g[0]);
         bg1[yzk] = (a+yzk) * log(1+1/b+g[1]);
         bg2[yzk] = (a+yzk) * log(1+1/b+g[2]);
      }

      tmp[0] = tmp[1] = tmp[2] = 0.0;
      for (i = 0 ; i < 3 ; i++)
      for (j = 0 ; j < 3 ; j++)
         tmp[j] += phi[i] * Q[i+3*j];

      tmp[0] *= exp(g0[zk] - bg0[zk]);
      tmp[1] *= exp(g1[zk] - bg1[zk]);
      tmp[2] *= exp(g2[zk] - bg2[zk]);

      double sum = tmp[0] + tmp[1] + tmp[2];
      for (i = 0 ; i < 3 ; i++) phi[i] = tmp[i] / sum;

      *l += log(sum);

   }

   // Finish the computation by adding the Gamma terms.
   for (i = 0 ; i < max_y+max_z+1 ; i++) *l += counter[i] * lgamma(a+i);
   *l -= *n * (lgamma(a) + log(b));

   free(g0);
   free(g1);
   free(g2);
   free(bg0);
   free(bg1);
   free(bg2);
   free(counter);

   return;

}
