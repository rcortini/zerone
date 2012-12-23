#include <stdlib.h>

void viterbi(
   // input //
   int *mptr,
   int *nptr,
   double *logini,
   double *logPem,
   double *logQ,
   // output //
   int *path
){

   const int m = *mptr;
   const int n = *nptr;

   int i;
   int j;
   int k;

   double max;
   double tmp;

   double *maxmat = (double *) malloc(m*n * sizeof(double));

   // Forward pass.
   for (j = 0 ; j < m ; j++) maxmat[j] = logini[j+0*m] + logPem[j+0*m];
   for (k = 1 ; k < n ; k++) {
      for (j = 0 ; j < m ; j++) {
         max = maxmat[0+(k-1)*m] + logQ[0+j*m];
         for (i = 1 ; i < m ; i++) {
            tmp = maxmat[i+(k-1)*m] + logQ[i+j*m];
            if (tmp > max) max = tmp;
         }
         maxmat[j+k*m] = max + logPem[j+k*m];
      }
   }

   // Backward pass.
   i = 0;
   for (j = 1 ; j < m ; j++)
      if (maxmat[j+(n-1)*m] > maxmat[i+(n-1)*m]) i = j;
   path[n-1] = i;
   for (k = n-2 ; k >= 0 ; k--) {
      max = logQ[0+path[k+1]*m] + maxmat[0+k*m];
      i = 0;
      for (j = 1 ; j < m ; j++) {
         tmp = logQ[j+path[k+1]*m] + maxmat[j+k*m];
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
