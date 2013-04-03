#include "loglik.h"

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
   int remaining_restarts;
   int signal;
} worker_arg;


// Global lock for mutex.
pthread_mutex_t lock;

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
   int m;
   int n = *n_obs;
   int r = *dim_z;
   // The following aliases are to avoid the statement '1/*beta'
   // which is interpreted as a comment.
   double alpha = *a;
   double beta = *b;

   double sum;
   double tmp[3];
   double phi[3] = {1.0/3, 1.0/3, 1.0/3};

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

   for (m = 0 ; tabulated[m] > -1 ; m++);
   double *pem = malloc(3*n * sizeof(double));
   int *sumyz = calloc(n, sizeof(double));
   double *lgamma_alpha_sumyz = malloc(m * sizeof(double));
   double lgamma_alpha = lgamma(alpha);
   for (k = 0 ; k < n ; k++) {
      for (int i = 0 ; i < r+1 ; i++) sumyz[k] += yz[i+k*(r+1)];
   }
   for (i = 0 ; i < m ; i++) {
      if (tabulated[i]) lgamma_alpha_sumyz[i] = lgamma(alpha + i);
   }

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
         if (!(pem[3*k]+pem[3*k+1]+pem[3*k+2] > DBL_EPSILON)) {
            pem[3*k  ] = q0;
            pem[3*k+1] = q1;
            pem[3*k+2] = q2;
         }
      }

      tmp[0] *= pem[3*index[k]  ];
      tmp[1] *= pem[3*index[k]+1];
      tmp[2] *= pem[3*index[k]+2];

      sum = tmp[0] + tmp[1] + tmp[2];
      if (sum > 0) {
         *loglik += log(sum);
         phi[0] = tmp[0] / sum;
         phi[1] = tmp[1] / sum;
         phi[2] = tmp[2] / sum;
      }
      else {
         int argmax = max3(pem+3*index[k]);
         if (phi[argmax] > 0) {
            *loglik += log(phi[argmax]) + pem[3*index[k]+argmax];
         }
         else {
            *loglik += lgamma_alpha - lgamma_alpha_sumyz[sumyz[k]];
         }
         phi[0] = 0.0;
         phi[1] = 0.0;
         phi[2] = 0.0;
         phi[argmax] = 1.0;
      }

      // Break out if computation is unstable.
      if (*loglik != *loglik) break;
      for (i = 0 ; i < 3 ; i++) phi[i] = tmp[i] / sum;

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
mean
(
   const int *yz,
   int n,
   int r
)
{
   double value = 0.0;
   int n_obs_no_NA = 0;
   for (int k = 0 ; k < n ; k++) {
      if (yz[(r+1)*k] < 0) continue;
      value += yz[(r+1)*k];
      n_obs_no_NA++;
   }
   return value / n_obs_no_NA;
}

int *
hist
(
   const int *yz,
   int n,
   int r
)
{
   int size = 1024;
   int *counts = calloc(size, sizeof(int));
   if (counts == NULL) {
      fprintf(stderr, "memory error (hist 1)\n");
      return NULL;
   }

   int maxval = 0;
   for (int k = 0 ; k < n ; k++) {
      int sum = 0;
      for (int i = 0 ; i < r+1 ; i++){
         if (yz[i+k*(r+1)] < 0) {
            sum = -1;
            break;
         }
         sum += yz[i+k*(r+1)];
      }
      if (sum < 0) continue;
      if (sum > maxval) maxval = sum;
      if (sum > size-2) {
         int newsize;
         for (newsize = size ; sum > newsize-2 ; newsize *= 2);
         int *newcounts = realloc(counts, newsize * sizeof(int));
         if (newcounts == NULL) {
            fprintf(stderr, "memory error (hist 2)\n");
            return NULL;
         }
         else {
            // Extra memory must be initialized to 0.
            counts = newcounts;
            for (int j = size ; j < newsize ; j++) counts[j] = 0;
            size = newsize;
         }
      }
      counts[sum]++;
   }
   // Add the sentinel.
   counts[maxval+1] = -1;

   return counts;

}



double
cool1
(int iter)
{
   //if (drand48() < T/10) return T - 0.2;
   //if (T < 0.4) return 0.0;
   //return T;
   if (iter > 200) return -1.0;
   if (iter > 150) return 0.0;
   return (5.0 / (1+iter));
}


double
cool2
(int iter)
{
   //if (T > 0.01) T = 0.01;
   //return T - 0.0001;
   //return 0.0;
   if (iter > 200) return -1.0;
   if (iter > 150) return 0.0;
   return (0.5 / (1+iter));
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
   for (i = 0 ; i < 3*r ; i++) {
      params->gammas[i] = old_params->gammas[i] * 
         (1+(range*(drand48()-.5)));
   }
   for (i = 0 ; i < 3 ; i++) {
      for (j = 0 ; j < 3 ; j++) {
         params->Q[i+3*j] = old_params->Q[i+3*j] *
            (1+(range*(drand48()-.5)));
      }
      double sum = params->Q[i+0] + params->Q[i+3*1] + params->Q[i+3*2];
      for (j = 0 ; j < 3 ; j++) {
         params->Q[i+3*j] /= sum;
      }
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
   globalmax[10]= params->alpha;
   globalmax[11] = params->beta;
   memcpy(globalmax + 1, params->Q, 9*sizeof(double));
   memcpy(globalmax + 12, params->gammas, 3*r*sizeof(double));
}


void *
search
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
   int *remaining_restarts = &unpack->remaining_restarts;

   // Allocate extra set of parameters.
   param_set *params = malloc(sizeof(param_set));
   params->gammas = malloc(3*r * sizeof(double));

   param_set *old_params = malloc(sizeof(param_set));
   old_params->gammas = malloc(3*r * sizeof(double));
   double loglik;
   double old_loglik = *globalmax;

   // Restart.
   while (*remaining_restarts > 0)  {
      pthread_mutex_lock(&lock);
      (*remaining_restarts)--;
      pthread_mutex_unlock(&lock);

      // Copy best conditions.
      memcpy(old_params->Q, globalmax+1, 9 * sizeof(double));
      memcpy(old_params->gammas, globalmax+12, 3*r * sizeof(double));
      old_params->alpha = globalmax[10];
      old_params->beta = globalmax[11];

      // Start simulated annealing.
      double T = 1;
      for (int iter = 1 ; T >= 0 ; T = (*cool)(iter++)) {

         if (*signal > 0) {
            // Copy best conditions.
            memcpy(old_params->Q, globalmax+1, 9 * sizeof(double));
            memcpy(old_params->gammas, globalmax+12, 3*r * sizeof(double));
            old_params->alpha = globalmax[10];
            old_params->beta = globalmax[11];
            pthread_mutex_lock(&lock);
            (*signal)--;
            pthread_mutex_unlock(&lock);
         }
         else {
            new_params(params, old_params, r, range);
         }
         compute_loglik(&n, &r, yz, params->Q, &params->alpha,
               &params->beta, params->gammas, index, tabulated,
               &loglik);
         loglik /= n;

         if (loglik != loglik) continue;

         // Keep the global best hit. 
         if (loglik > *globalmax) {
            fprintf(stderr, "%f (T=%.3f)\n", loglik, T);
            pthread_mutex_lock(&lock);
            *signal = 25;
            update_globalmax(globalmax, loglik, params, r);
            pthread_mutex_unlock(&lock);
         }

         int loglik_increased = (loglik > old_loglik);
         int change_anyway = (drand48() < exp(20*(loglik-old_loglik)/T));

         if (loglik_increased || change_anyway) {
            param_set *tmp = old_params;
            old_params = params;
            params = tmp;
            old_loglik = loglik;
         }
      }
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

/*
#ifdef _SC_NPROCESSORS_ONLN
   int n_threads = (int) sysconf(_SC_NPROCESSORS_ONLN) - 1;
   if (n_threads < 1) n_threads = 1;
#else
   int n_threads = 1;
#endif
*/

   int n_threads = 50;

   int i;
   int n = *n_obs;
   int r = *dim_z;

   // Index time series.
   int *tabulated = hist(yz, n, r);
   int *index = malloc(n * sizeof(int));
   indexts(n, r+1, yz, index);

   // Initial parameters.
   param_set *init_params = malloc(sizeof(param_set));
   init_params->gammas = malloc(3*r * sizeof(double));

   for (i = 0 ; i < 9 ; i++) init_params->Q[i] = .025;
   init_params->Q[0] = init_params->Q[4] = init_params->Q[8] = .95;
   init_params->alpha = 1.0;
   init_params->beta = mean(yz, n, r);
   for (i = 0 ; i < r ; i++) {
      double meanz_by_beta = mean(yz+i+1, n, r) / init_params->beta;
      init_params->gammas[i+0*r] = meanz_by_beta / 2;
      init_params->gammas[i+1*r] = meanz_by_beta;
      init_params->gammas[i+2*r] = meanz_by_beta * 8;
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
      .cool = &cool1,
   };

   fprintf(stderr, "starting threads\n");
   // Allocate 'tid'.
   pthread_t *tid = (pthread_t *) malloc(n_threads * sizeof(pthread_t));

   do {
      // Set number of restarts.
      arg.remaining_restarts = 50;
      arg.signal = 0;
      for (i = 0 ; i < n_threads ; i++) {
         err = pthread_create(&(tid[i]), NULL, &search, &arg);
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
   while (*globalmax != *globalmax);

   fprintf(stderr, "starting second round of optimization\n");

   // Second round optimization.
   memcpy(init_params->Q, globalmax+1, 9 * sizeof(double));
   memcpy(init_params->gammas, globalmax+12, 3*r * sizeof(double));
   init_params->alpha = globalmax[10];
   init_params->beta = globalmax[11];
   arg.range = .05;
   arg.cool = &cool1;
   arg.signal = 0;
   arg.remaining_restarts = 50;

   // Start search again with new conditions.
   for (i = 0 ; i < n_threads ; i++) {
      err = pthread_create(&(tid[i]), NULL, &search, &arg);
      if (err) {
         fprintf(stderr, "error creating thread (%d)\n", err);
         return;
      }
   }

   // Wait for threads to return.
   for (i = 0 ; i < n_threads ; i++) {
      pthread_join(tid[i], NULL);
   }

   fprintf(stderr, "starting final optimization\n");

   // Final optimization.
   memcpy(init_params->Q, globalmax+1, 9 * sizeof(double));
   memcpy(init_params->gammas, globalmax+12, 3*r * sizeof(double));
   init_params->alpha = globalmax[10];
   init_params->beta = globalmax[11];
   arg.range = .01;
   arg.cool = &cool2;
   arg.signal = 0;
   arg.remaining_restarts = 50;

   // Start search again with fine conditions.
   for (i = 0 ; i < n_threads ; i++) {
      err = pthread_create(&(tid[i]), NULL, &search, &arg);
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

}
