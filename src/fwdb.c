#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>

// Thread arguments.
typedef struct {
   const int nblocks;
         int *taskQ;
   const int m;
   const double *Q; 
   const double *init;
         double *prob;
         double *phi;
         double *sumtrans;
         double *loglik;
} worker_arg;


// Thread synchronization.
int taskQ_i;
pthread_mutex_t lock;

double fwdb(
   // input //
   const int m,
   const int n,
   const double *Q,
   const double *init,
   // output //
   double *prob,
   double *phi,
   // not initialized //
   double *sumtrans
){
// SYNOPSIS:                                                             
//   Forward-backward algorithm with Markovian backward smooting.        
//                                                                       
// ARGUMENTS:                                                            
//   'm': the number of states                                           
//   'n': the length of the sequence of observations                     
//   'Q': (m,m) transition matrix ('Q[i+j*m]' is a ij transtition).      
//   'init': (m) initial probabilities                                   
//   'prob': (m,n) emission probabilities                                
//   'phi': (m,n) probabilities given observations                       
//   'sumtrans': (m,m) sum of conditional transitions probabilties       
//                                                                       
// RETURN:                                                               
//   The total log-likelihood.                                           
//                                                                       
// SIDE EFFECTS:                                                         
//   Replaces 'prob' by alphas, updates 'phi' in place.                  

   int i;
   int j;
   int k;
   double S;            // Sum used for computation intermediates.
   double tmp[m];
   double alpha[m];
   double rev[m*m];

   // Nomrmalization constants, also used to return log-likelihood.
   double *c = (double *) malloc(n * sizeof(double));

   // 'sumtrans[i+j*m]' is the Sum of conditional transition
   // probabilities from state 'i' to state 'j' (congruent with 'Q').
   // for (j = 0 ; j < m*m ; j++) sumtrans[j] = 0.0;

   // First iteration of the forward pass.
   c[0] = 0.0;
   for (j = 0 ; j < m ; j++) {
      alpha[j] = init[j] * prob[j+0*m];
      c[0] += alpha[j];
   }
   // Sum-normalize 'alpha'.
   for (j = 0 ; j < m ; j++) {
      alpha[j] /= c[0];
      prob[j+0*m] = alpha[j];
   }

   // Next iterations of the forward pass.
   for (k = 1 ; k < n ; k++) {
      for (j = 0 ; j < m ; j++) {
         tmp[j] = 0.0;
         for (i = 0 ; i < m ; i++) tmp[j] += alpha[i] * Q[i+j*m];
      }
      c[k] = 0.0;
      for (j = 0 ; j < m ; j++) {
         alpha[j] = tmp[j] * prob[j+k*m];
         c[k] += alpha[j];
      }
      // Sum-normalize 'alpha'.
      for (j = 0 ; j < m ; j++) {
         alpha[j] /= c[k];
         prob[j+k*m] = alpha[j];
      }
   }

   // NB: 'prob' now contains the forward alphas.
   // First iteration of the backward pass.
   for (j = 0 ; j < m ; j++) phi[j+(n-1)*m] = prob[j+(n-1)*m];

   // Next iterations of the backward pass.
   for (k = n-2 ; k >= 0 ; k--) {
      for (j = 0 ; j < m ; j++) {
         S = 0.0;
         for (i = 0 ; i < m ; i++) S += rev[j+i*m] = prob[i+k*m]*Q[i+j*m];
         for (i = 0 ; i < m ; i++) rev[j+i*m] /= S;
      }
      for (j = 0 ; j < m ; j++) {
         phi[j+k*m] = 0.0;
         for (i = 0 ; i < m ; i++) {
            S = phi[i+(k+1)*m] * rev[i+j*m];
            phi[j+k*m] += S;
            sumtrans[i+j*m] += S;
         }
      }
   }

   S = log(c[0]);
   for (k = 1 ; k < n ; k++) S += log(c[k]);
   free(c);
   return S;

}


void *work(
   void *arg
){
// Just work bitch!!

   worker_arg *myarg = (worker_arg *) arg;
   const int   nblocks = myarg->nblocks;
         int    *taskQ = myarg->taskQ;
   const int         m = myarg->m;
   const double     *Q = myarg->Q;
   const double  *init = myarg->init;
         double  *prob = myarg->prob;
         double   *phi = myarg->phi;

   int i;
   int k;
   int n;
   int job_index;

   double loglik = 0.0;
   double *sumtrans = (double *) malloc(m*m * sizeof(double));
   for (i = 0 ; i < m*m ; i++) sumtrans[i] = 0.0;

   while (1) {
      pthread_mutex_lock(&lock);
         job_index = taskQ_i++;
      pthread_mutex_unlock(&lock);

      // Stop if task queue is empty.
      if (job_index >= nblocks) break;

      k = taskQ[job_index];
      n = taskQ[job_index+1] - k;

      // Update 'prob', 'phi' and local 'sumtrans'.
      loglik += fwdb(m, n, Q, init, prob+k*m, phi+k*m, sumtrans);

   }

   // Request the lock to update global 'sumtrans' and 'loglik'.
   pthread_mutex_lock(&lock);
      myarg->loglik[0] += loglik;
      for (i = 0 ; i < m*m ; i++) myarg->sumtrans[i] += sumtrans[i];
   pthread_mutex_unlock(&lock);

   free(sumtrans);
   return NULL;

}

void R_fwdb(
   // input //
   int *m,
   int *nblocks,
   int *lengths, // The length of the fragments. //
   double *Q,
   double *init,
   // output //
   double *prob,
   double *phi,
   double *sumtrans,
   double *loglik,
   // thread control //
   int *n_threads
){
// SYNOPSIS:
//    Wrapper for the R call using '.C("R_fwdb", ...)'.

   // Get the number of CPUs if 'n_threads' is 0 or lower.
   if (*n_threads < 1) {
   #ifdef _SC_NPROCESSORS_ONLN
         *n_threads = (int) sysconf(_SC_NPROCESSORS_ONLN);
   #else
         *n_threads = 1;
   #endif
   }

   // Do not create more threads than needed.
   *n_threads = *n_threads > *nblocks ? *nblocks : *n_threads;

   int i;
   int err;

   // Initialize output variables updated by increment.
   for (i = 0 ; i < (*m)*(*m) ; i++) sumtrans[i] = 0.0;
   *loglik = 0.0;

   // Create task queue.
   int *taskQ = (int *) malloc((*nblocks+1) * sizeof(int));
   taskQ[0] = 0;
   for (i = 1 ; i < *nblocks+1 ; i++) {
      taskQ[i] = taskQ[i-1] + lengths[i-1];
   }

   err = pthread_mutex_init(&lock, NULL);
   if (err) {
      fprintf(stderr, "error initializing mutex (%d)\n", err);
      free(taskQ);
      return;
   }

   pthread_t *tid = (pthread_t *) malloc(*n_threads * sizeof(pthread_t));

   worker_arg arg = {
      .nblocks = (const int) *nblocks,
      .taskQ = taskQ,
      .m = (const int) *m,
      .Q = (const double *) Q,
      .init = (const double *) init,
      .prob = prob,
      .phi = phi,
      .sumtrans = sumtrans,
      .loglik = loglik,
   };

   taskQ_i = 0;
   // Instantiate threads and start running jobs.
   for (i = 0 ; i < *n_threads ; i++) tid[i] = 0;
   for (i = 0 ; i < *n_threads ; i++) {
      err = pthread_create(&(tid[i]), NULL, &work, &arg);
      if (err) {
         fprintf(stderr, "error creating thread (%d)\n", err);
         free(taskQ);
         free(tid);
         return;
      }
   }

   // Wait for threads to return.
   for (i = 0 ; i < *n_threads ; i++) {
      pthread_join(tid[i], NULL);
   }

   free(taskQ);
   return;

}
