#include <stdlib.h>

void cviterbi(
   // input //
   const int n,
   const int m,
   const double *logini,
   const double *logPem,
   const double *logQ,
   // output //
   int *path
){

   int i;
   int j;
   int k;

   double max;
   double tmp;

   double *maxmat = (double *) malloc(m*n * sizeof(double));

   // 
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
   // 
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

void Rvit(int *n, int *m, double *logini, double *logPem,
      double *logQ, int *path) {
   cviterbi(*n, *m, (const double *) logini, (const double *) logPem,
   (const double *) logQ, path);
   return ;
}
