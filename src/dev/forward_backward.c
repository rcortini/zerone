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


