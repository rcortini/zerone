#include "loglik.h"

// Arguments for `social_search`.
typedef struct {
   int n;
   int r;
   const int *yz;
   int *index;
   const int *tabulated;
   double *globalmax;
   double range;
   double (*cool)(int);
   // control //
   int signal;
} worker_arg;


// Global lock for mutex.
pthread_mutex_t lock;
//FILE *outf;

void
compute_loglik (
   // input //
   const int *n_obs,
   const int *dim_z,
   const int *yz,
   // params //
   const double *Q,
   const double *a,
   const double *b,
   const double *gamma,
   // index //
   int *index,
   const int *tabulated,
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
// ARGUMENTS:                                                            
//   'n_obs': (1) length of the sequence of observations                 
//   'dim_z': (1) dimension of z (number of profiles).                   
//   'yz': (dim_z+1,n_obs) control profile.                              
//   'a': (1) alias 'alpha', model parameter                             
//   'b': (1) alias 'beta', model parameter                              
//   'gamma': (dim_z,2) model parameter                                  
//                                                                       
// RETURN:                                                               
//   'void'                                                              
{

   int i;
   int j;
   int k;
   int m;
   int n = *n_obs;
   int r = *dim_z;
   // The following aliases are to avoid the statement '1/*a'
   // which is interpreted as a comment.
   double alpha = *a;
   double beta = *b;

   double sum;
   double tmp[2];
   double phi[2] = {1.0/2, 1.0/2};

   // Compute p's.
   double logp[(r+2)*2];
   for (i = 0 ; i < 2 ; i++) {
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

   for (m = 0 ; tabulated[m] > -1 ; m++);
   int *sumyz = calloc(n, sizeof(double));
   double *lgamma_alpha_sumyz = malloc(m * sizeof(double));
   double lgamma_alpha = lgamma(alpha);
   for (k = 0 ; k < n ; k++) {
      for (int i = 0 ; i < r+1 ; i++) sumyz[k] += yz[i+k*(r+1)];
   }
   for (i = 0 ; i < m ; i++) {
      if (tabulated[i]) lgamma_alpha_sumyz[i] = lgamma(alpha + i);
   }

   double *pem = malloc(2*n * sizeof(double));
   *loglik = 0.0;

   for (k = 0 ; k < n ; k++) {
      tmp[0] = tmp[1] = 0.0;
      for (i = 0 ; i < 2 ; i++) {
      for (j = 0 ; j < 2 ; j++) {
         tmp[j] += phi[i] * Q[i+2*j];
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
         for (i = 0 ; i < 2 ; i++) phi[i] = tmp[i];
         continue;
      }

      // Caching by indexing.
      if (index[k] == k) {
         // Compute log-ratio.
         double q0 = alpha * logp[0];
         double q1 = alpha * logp[0+(r+2)*1];
         for (i = 0 ; i < r+1; i++) {
            q0 += yz[i+k*(r+1)] * logp[i+1];
            q1 += yz[i+k*(r+1)] * logp[i+1+(r+2)*1];
         }

         // Take exponential.
         pem[2*k  ] = exp(q0);
         pem[2*k+1] = exp(q1);

         // NB: In case of numerical underflow, we store the
         // log probability of emission, which are negative.
         // FIXME make 'pem' a long double??
         if (!(pem[2*k]+pem[2*k+1] > DBL_EPSILON)) {
            pem[2*k  ] = q0;
            pem[2*k+1] = q1;
         }
      }

      tmp[0] *= pem[2*index[k]  ];
      tmp[1] *= pem[2*index[k]+1];

      sum = tmp[0] + tmp[1];
      if (sum > 0) {
         *loglik += log(sum);
         phi[0] = tmp[0] / sum;
         phi[1] = tmp[1] / sum;
      }
      else {
         int argmax = pem[2*index[k]] > pem[2*index[k]+1] ? 0 : 1;
         if (phi[argmax] > 0) {
            *loglik += log(phi[argmax]) + pem[2*index[k]+argmax];
         }
         else {
            *loglik += lgamma_alpha - lgamma_alpha_sumyz[sumyz[k]];
         }
         phi[0] = 0.0;
         phi[1] = 0.0;
         phi[argmax] = 1.0;
      }

      // Break out if computation is unstable.
      if (*loglik != *loglik) break;
      for (i = 0 ; i < 2 ; i++) phi[i] = tmp[i] / sum;

   }

   int n_no_NA = 0;
   for (i = 0 ; tabulated[i] > -1 ; i++) {
      *loglik += tabulated[i] * lgamma_alpha_sumyz[i];
      n_no_NA += tabulated[i];
   }
   *loglik -= n_no_NA * lgamma_alpha;

   free(pem);
   free(sumyz);
   free(lgamma_alpha_sumyz);
   return;
}


double
cool
(int iter)
{
   if (iter > 150) return -1.0;
   return (.15 / (1+iter));
}

param_set *
new_params
(
   param_set *params,
   const param_set *old_params,
   int r,
   double range
)
{
   int i,j;
   params->alpha = old_params->alpha * (1+(range*(drand48()-.5)));
   params->beta = old_params->beta * (1+(range*(drand48()-.5)));
   for (i = 0 ; i < 2*r ; i++) {
      params->gammas[i] = old_params->gammas[i] * 
         (1+(range*(drand48()-.5)));
   }
   for (i = 0 ; i < 2 ; i++) {
      for (j = 0 ; j < 2 ; j++) {
         params->Q[i+2*j] = old_params->Q[i+2*j] *
            (1+(range*(drand48()-.5)));
      }
      double sum = params->Q[i+0] + params->Q[i+2*1];
      for (j = 0 ; j < 2 ; j++) params->Q[i+2*j] /= sum;
   }

   return params;

}

void
destroy_params
(
   param_set *p
)
{
   free(p->gammas);
   free(p);
}

void
update_globalmax
(
   double *globalmax,
   double loglik,
   param_set *params,
   int r
)
{
   globalmax[0] = loglik;
   globalmax[5]= params->alpha;
   globalmax[6] = params->beta;
   memcpy(globalmax + 1, params->Q, 4*sizeof(double));
   memcpy(globalmax + 7, params->gammas, 2*r*sizeof(double));
}


void *
social_search
(void *arg)
{
   // Unpack arguments.
   worker_arg *unpack = (worker_arg *) arg;
   int n = unpack->n;
   int r = unpack->r;
   const int *yz = unpack->yz;
   int *index = unpack->index;
   const int *tabulated = unpack->tabulated;
   double *globalmax = unpack->globalmax;
   double range = unpack->range;
   double (*cool)(int) = unpack->cool;
   int *signal = &unpack->signal;

   // Allocate extra set of parameters.
   param_set *params = malloc(sizeof(param_set));
   params->gammas = malloc(2*r * sizeof(double));

   param_set *old_params = malloc(sizeof(param_set));
   old_params->gammas = malloc(2*r * sizeof(double));
   double loglik;
   double old_loglik = *globalmax;

   // Copy best conditions.
   memcpy(old_params->Q, globalmax+1, 4 * sizeof(double));
   memcpy(old_params->gammas, globalmax+7, 2*r * sizeof(double));
   old_params->alpha = globalmax[5];
   old_params->beta = globalmax[6];

   // Start simulated annealing.
   double T = .15;
   for (int iter = 1 ; T >= 0 ; T = (*cool)(iter++)) {

      if (*signal > 0) {
         // Copy best conditions.
         memcpy(old_params->Q, globalmax+1, 4 * sizeof(double));
         memcpy(old_params->gammas, globalmax+7, 2*r * sizeof(double));
         old_params->alpha = globalmax[5];
         old_params->beta = globalmax[6];
         pthread_mutex_lock(&lock);
         (*signal)--;
         pthread_mutex_unlock(&lock);
      }
      new_params(params, old_params, r, range);
      compute_loglik(&n, &r, yz, params->Q, &params->alpha,
            &params->beta, params->gammas, index, tabulated,
            &loglik);
      loglik /= n;

      // Re-iterate if 'loglik' is NA.
      if (loglik != loglik) {
         iter--;
         continue;
      }

      // Keep the global best hit. 
      if (loglik > *globalmax) {
         fprintf(stderr, "%f (T=%.3f)\n", loglik, T);
         pthread_mutex_lock(&lock);
         *signal = 12;
         update_globalmax(globalmax, loglik, params, r);
         pthread_mutex_unlock(&lock);
      }

      int loglik_increased = (loglik > old_loglik);
      int change_anyway = drand48() < exp((loglik-old_loglik)/T);

      if (loglik_increased || change_anyway) {
         param_set *tmp = old_params;
         old_params = params;
         params = tmp;
         old_loglik = loglik;
      }
      // Save trace in file.
      //pthread_mutex_lock(&lock);
      //fprintf(outf, "%d\t%f\t%f", (int)pthread_self(), old_params->alpha,
      //      old_params->beta);
      //fprintf(outf, "\t%d\n", loglik == *globalmax);
      //for (int j = 0 ; j < 4 ; j++) {
      //   fprintf(outf, "\t%f", old_params->Q[j]);
      //}
      //for (int j = 0 ; j < 2*r ; j++) {
      //   fprintf(outf, "\t%f", old_params->gammas[j]);
      //}
      //fprintf(outf, "\n");
      //pthread_mutex_unlock(&lock);
   }

   destroy_params(params);
   destroy_params(old_params);

   return NULL;

}

void
simAnneal
(
   // input //
   const int *n_obs,
   const int *dim_z,
   const int *yz,
   // output //
   double *globalmax
)
{
   //outf = fopen("tracks.txt", "w");
   int n_threads = 25;

   int i;
   int n = *n_obs;
   int r = *dim_z;

   // Index time series.
   int *tabulated = histsum(yz, n, r+1);
   int *index = malloc(n * sizeof(int));
   indexts(n, r+1, yz, index);

   // Initial parameters.
   param_set *init_params = malloc(sizeof(param_set));
   init_params->gammas = malloc(2*r * sizeof(double));

   for (i = 0 ; i < 4 ; i++) init_params->Q[i] = .01;
   init_params->Q[0] = init_params->Q[2] = .99;
   init_params->alpha = 1.0;
   init_params->beta = mean(yz, n, r+1);
   for (i = 0 ; i < r ; i++) {
      double meanz_by_beta = mean(yz+i+1, n, r+1) / init_params->beta;
      init_params->gammas[i+0*r] = meanz_by_beta;
      init_params->gammas[i+1*r] = meanz_by_beta * 8;
   }

   // Initialize random generator.
   srand48(time(NULL));
   compute_loglik(&n, &r, yz, init_params->Q, &init_params->alpha,
         &init_params->beta, init_params->gammas, index, tabulated,
         globalmax);
   *globalmax /= n;
   update_globalmax(globalmax, *globalmax, init_params, r);
   if (*globalmax != *globalmax) *globalmax = -INFINITY;

   // Initialize mutex.
   int err = pthread_mutex_init(&lock, NULL);
   if (err) {
      fprintf(stderr, "error initializing mutex (%d)\n", err);
      return;
   }

   // Initialize arguments.
   worker_arg arg = {
      .n = n,
      .r = r,
      .yz = yz,
      .index = index,
      .tabulated = (const int *) tabulated,
      .globalmax = globalmax,
      .range = .4,
      .cool = &cool,
   };

   // Allocate 'tid'.
   pthread_t *tid = (pthread_t *) malloc(n_threads * sizeof(pthread_t));

   fprintf(stderr, "starting threads\n");
   arg.signal = 0;
   for (i = 0 ; i < n_threads ; i++) {
      err = pthread_create(&(tid[i]), NULL, &social_search, &arg);
      if (err) {
         fprintf(stderr, "error creating thread (%d)\n", err);
         return;
      }
   }

   // Wait for threads to return.
   for (i = 0 ; i < n_threads ; i++) {
      pthread_join(tid[i], NULL);
   }

   fprintf(stderr, "starting second round of optimization\n");

   // Second round optimization.
   for (int j = 0 ; j < 5 ; j++) {
      memcpy(init_params->Q, globalmax+1, 4 * sizeof(double));
      memcpy(init_params->gammas, globalmax+7, 2*r * sizeof(double));
      init_params->alpha = globalmax[5];
      init_params->beta = globalmax[6];
      arg.range = .05;
      arg.cool = &cool;
      arg.signal = 0;

      // Start search again with new conditions.
      for (i = 0 ; i < n_threads ; i++) {
         err = pthread_create(&(tid[i]), NULL, &social_search, &arg);
         if (err) {
            fprintf(stderr, "error creating thread (%d)\n", err);
            return;
         }
      }

      // Wait for threads to return.
      for (i = 0 ; i < n_threads ; i++) {
         pthread_join(tid[i], NULL);
      }
   }

   fprintf(stderr, "starting final optimization\n");

   // Final optimization.
   memcpy(init_params->Q, globalmax+1, 4 * sizeof(double));
   memcpy(init_params->gammas, globalmax+7, 2*r * sizeof(double));
   init_params->alpha = globalmax[5];
   init_params->beta = globalmax[6];
   arg.range = .01;
   arg.cool = &cool;
   arg.signal = 0;

   // Start search again with fine conditions.
   for (i = 0 ; i < n_threads ; i++) {
      err = pthread_create(&(tid[i]), NULL, &social_search, &arg);
      if (err) {
         fprintf(stderr, "error creating thread (%d)\n", err);
         return;
      }
   }

   // Wait for threads to return.
   for (i = 0 ; i < n_threads ; i++) {
      pthread_join(tid[i], NULL);
   }

   pthread_mutex_destroy(&lock);
   destroy_params(init_params);
   free(index);
   free(tabulated);

   free(tid);
   //fclose(outf);
}
