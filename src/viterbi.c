#include <stdlib.h>

void viterbi(
   // input //
   int m,
   int n,
   double *log_init,
   double *log_Q,
   double *log_pem,
   // output //
   int *path
){

   int i;
   int j;
   int k;

   double max;
   double tmp;

   double *maxmat = malloc(m*n * sizeof(double));

   // Forward pass.
   for (j = 0 ; j < m ; j++) maxmat[j] = log_init[j+0*m] + log_pem[j+0*m];
   for (k = 1 ; k < n ; k++) {
      for (j = 0 ; j < m ; j++) {
         max = maxmat[0+(k-1)*m] + log_Q[0+j*m];
         for (i = 1 ; i < m ; i++) {
            tmp = maxmat[i+(k-1)*m] + log_Q[i+j*m];
            if (tmp > max) max = tmp;
         }
         maxmat[j+k*m] = max + log_pem[j+k*m];
      }
   }

   // Backward pass.
   i = 0;
   for (j = 1 ; j < m ; j++)
      if (maxmat[j+(n-1)*m] > maxmat[i+(n-1)*m]) i = j;
   path[n-1] = i;
   for (k = n-2 ; k >= 0 ; k--) {
      max = log_Q[0+path[k+1]*m] + maxmat[0+k*m];
      i = 0;
      for (j = 1 ; j < m ; j++) {
         tmp = log_Q[j+path[k+1]*m] + maxmat[j+k*m];
         if (tmp > max) {
            max = tmp;
            i = j;
         }
      }
      path[k] = i;
   }

   free(maxmat);
   return;

}

void block_viterbi(
   // input //
   int *nblocks,
   int *lengths,
   double *log_init,
   double *log_pem,
   double *log_Q,
   // output //
   int *path
){

   int i;
   int start = 0;

   for (i = 0 ; i < *nblocks ; i++) {
      viterbi(2, lengths[i], log_init, log_Q, log_pem+start, path+start);
      start += 2*lengths[i];
   }

   return;

}
