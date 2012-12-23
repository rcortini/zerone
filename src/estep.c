#include <math.h>
#include <float.h>
#include <pthread.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef DBL_EPSILON
#define DBL_EPSILON 2.2e-16
#endif

// There is no need to start more than two threads.
#define MAX_N_THREADS 2
#define M_LOOK 10000

// Thread arguments.
typedef struct {
   const int nthreads;
   const pthread_t *tid;
   const int n;
   const int *x;
   const int *y;
   const double *lx;
   const double *ly;
         double *pratio;
} pratio_worker_arg;

typedef struct {
   const int nthreads;
   const pthread_t *tid;
   const int nblocks;
         int *taskQ;
   const double *Q; 
   const double init;
         double *pratio;
         double *phi;
         double *sumtrans;
} spcfwdb_worker_arg;


// Thread synchronization.
int taskQ_i;
pthread_mutex_t lock;

int nCPU(void) {
#ifdef _SC_NPROCESSORS_ONLN
   return (int) sysconf(_SC_NPROCESSORS_ONLN);
#else
   return = 1;
#endif
}


void spcfwdb(
   // input //
   const int n,
   const double *Q,
   const double init,
   // output //
   double *ratio,
   double *phi,
   // not initialized //
   double *sumtrans
){
// SYNOPSIS:                                                             
//   Forward-backward algorithm with Markovian backward smooting.        
//                                                                       
// ARGUMENTS:                                                            
//   'n': length of the sequence of observations                         
//   'Q': (2) diagonal of the transition matrix                          
//   'init': initial probability of the first state                      
//   'ratio': (n) emission probabilities ratio                           
//   'phi': (n) probabilities given observations                         
//   'sumtrans': (2,2) sum of conditional transitions probabilties       
//                                                                       
// RETURN:                                                               
//   'void'
//                                                                       
// SIDE EFFECTS:                                                         
//   Replace 'ratio' by alphas, update 'phi' in place.                 

   int k;
   double a;
   double tmp[2];
   double rev[2];
   double *alpha = ratio;  // For clarity.

   // First iteration of the forward pass.
   a = alpha[0] = (init * ratio[0]) / (init * ratio[0] + 1.0 - init);

   // Next iterations of the forward pass.
   for (k = 1 ; k < n ; k++) {
      tmp[0] = a * Q[0] + (1-a) * (1-Q[1]);
      tmp[1] = a * (1-Q[0]) + (1-a) * Q[0];
      a = alpha[k] = (tmp[0]*ratio[k]) / (tmp[0]*ratio[k]+tmp[1]);
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

void *pratio_work(
   void *arg
){

   pratio_worker_arg *myarg = (pratio_worker_arg *) arg;
   const int        nthreads = myarg->nthreads;
   const pthread_t      *tid = myarg->tid;
   const int               n = myarg->n;
   const int              *x = myarg->x;
   const int              *y = myarg->y;
   const double          *lx = myarg->lx;
   const double          *ly = myarg->ly;
         double      *pratio = myarg->pratio;

   int j;

   // Identify thread from 'tid'.
   int myrank;
   for (myrank = 0 ; myrank < nthreads ; myrank++) {
      if (pthread_equal(pthread_self(), tid[myrank])) break;
   }

   // Allocate task from thread id.
   int start = myrank * (n / nthreads);
   int end = myrank == nthreads -1 ? n : (myrank + 1) * (n / nthreads);

   // Fill 'pratio' from 'start' to 'end'.
   for (j = start ; j < end ; j++) {
      // NAs are passed as negative values to 'int'. Set ratio to
      // 1.0 in case of NA emission (assuming both states have the
      // same probability of producing NAs).
      if (x[j] < 0 || y[j] < 0) pratio[j] = 1.0;
      else pratio[j] = lx[x[j]]*ly[y[j]];
   }

   return NULL;

}

void *spcfwdb_work(
   void *arg
){

   spcfwdb_worker_arg *myarg = (spcfwdb_worker_arg *) arg;
   const int      nthreads = myarg->nthreads;
   const pthread_t    *tid = myarg->tid;
   const int       nblocks = myarg->nblocks;
         int        *taskQ = myarg->taskQ;
   const double         *Q = myarg->Q;
   const double       init = myarg->init;
         double    *pratio = myarg->pratio;
         double       *phi = myarg->phi;

   int i;
   int j;

   int myrank;
   for (myrank = 0 ; myrank < nthreads ; myrank++) {
      if (pthread_equal(pthread_self(), tid[myrank])) break;
   }

   double *sumtrans = (double *) malloc(2*2 * sizeof(double));
   for (i = 0 ; i < 2*2 ; i++) sumtrans[i] = 0.0;

   for (j = myrank ; j < nblocks ; j += nthreads) {

      int k = taskQ[j];
      int n = taskQ[j+1] - k;

      // Update 'prob', 'phi' and local 'sumtrans'.
      spcfwdb(n, Q, init, pratio+k, phi+k, sumtrans);

   }

   // Request the lock to update global 'sumtrans'.
   pthread_mutex_lock(&lock);
      for (i = 0 ; i < 2*2 ; i++) myarg->sumtrans[i] += sumtrans[i];
   pthread_mutex_unlock(&lock);

   free(sumtrans);
   return NULL;

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
   double *fwd_alpha,
   // thread control //
   int *nthreads
){

   if (*nthreads < 1) *nthreads = nCPU();
   *nthreads = *nthreads > MAX_N_THREADS ? MAX_N_THREADS : *nthreads;

   int i;
   int err;

   double _r = (1 + 1 / *beta + gamma[1]) / (1 + 1 / *beta + gamma[0]);
   double r = gamma[0] * _r / gamma[1];

   double *lx = (double *) malloc (M_LOOK * sizeof(double));
   double *ly = (double *) malloc (M_LOOK * sizeof(double));
   for (i = 0 ; i < M_LOOK ; i++) lx[i] = ly[i] = -1.0;

   for (i = 0 ; i < *n ; i++) {
      if (x[i] >= 0 && lx[x[i]] < 0) lx[x[i]] = pow(r, x[i]);
      if (y[i] >= 0 && ly[y[i]] < 0) ly[y[i]] = pow(_r, *alpha + y[i]);
   }

   pthread_t *tid = (pthread_t *) malloc(*nthreads * sizeof(pthread_t));

   pratio_worker_arg pratio_arg = {
      .nthreads = *nthreads,
      .tid = tid,
      .n = *n,
      .x = x,
      .y = y,
      .lx = lx,
      .ly = ly,
      .pratio = fwd_alpha,     // Store ratios in 'fwd_alpha'.
   };

   // Instantiate threads and start running jobs.
   for (i = 0 ; i < *nthreads ; i++) {
      err = pthread_create(tid+i, NULL, &pratio_work, &pratio_arg);
      if (err) {
         fprintf(stderr, "error creating thread (%d)\n", err);
         free(tid);
         free(lx);
         free(ly);
         return;
      }
   }

   // Wait for threads to return.
   for (i = 0 ; i < *nthreads ; i++) pthread_join(tid[i], NULL);

   free(tid);
   free(lx);
   free(ly);

   return;
}

void fwdb(
   // input //
   int *n,
   int *nblocks,
   int *lengths,
   // params //
   double *Q,
   double *init,
   // output //
   double *fwd_alpha,
   double *phi,
   double *sumtrans,
   // thread control //
   int *nthreads
){

   if (*nthreads < 1) *nthreads = nCPU();
   *nthreads = *nthreads > MAX_N_THREADS ? MAX_N_THREADS : *nthreads;

   int i;
   int err;

   for (i = 0 ; i < 2*2 ; i++) sumtrans[i] = 0.0;

   // Create task queue.
   int *taskQ = (int *) malloc((*nblocks+1) * sizeof(int));
   taskQ[0] = 0;
   for (i = 1 ; i < *nblocks+1 ; i++) taskQ[i] = taskQ[i-1]+lengths[i-1];

   err = pthread_mutex_init(&lock, NULL);
   if (err) {
      fprintf(stderr, "error initializing mutex (%d)\n", err);
      free(taskQ);
      return;
   }

   pthread_t *tid = (pthread_t *) malloc(*nthreads * sizeof(pthread_t));

   spcfwdb_worker_arg spcfwdb_arg = {
      .nthreads = *nthreads,
      .tid = tid,
      .nblocks = (const int) *nblocks,
      .taskQ = taskQ,
      .Q = (const double *) Q,
      .init = *init,
      .pratio = fwd_alpha,    // The ratios were sotred in 'fwd_alpha'.
      .phi = phi,
      .sumtrans = sumtrans,
   };

   taskQ_i = 0;
   // Instantiate threads and start running jobs.
   for (i = 0 ; i < *nthreads ; i++) {
      err = pthread_create(tid+i, NULL, &spcfwdb_work, &spcfwdb_arg);
      if (err) {
         fprintf(stderr, "error creating thread (%d)\n", err);
         free(taskQ);
         free(tid);
         return;
      }
   }

   // Wait for threads to return.
   for (i = 0 ; i < *nthreads ; i++) pthread_join(tid[i], NULL);

   pthread_mutex_destroy(&lock);

   free(taskQ);
   free(tid);

   return;
}

void estep(
   // input //
   int *n,
   int *x,
   int *y,
   int *nblocks,
   int *lengths,
   // params //
   double *alpha,
   double *beta,
   double *gamma,
   double *Q,
   double *init,
   // output //
   double *fwd_alpha,
   double *phi,
   double *sumtrans,
   // thread control //
   int *nthreads
){

   ///  PART I: COMPUTE EMISSION PROBABILITY RATIOS  ///
   compute_pratio(n, x, y, alpha, beta, gamma, fwd_alpha, nthreads);

   ///  PART II: SPECIAL FOWARD-BACKWARD ALGORITHM  ///
   fwdb(n, nblocks, lengths, Q, init, fwd_alpha, phi, sumtrans, nthreads);

   return;


}
