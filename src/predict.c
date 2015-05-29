#include "predict.h"

double *
readmatrix
(
   char   * fn,
   int      nrow,
   int      ncol
)
{
   FILE   * fp  = fopen(fn, "r");
   if (fp == NULL) return NULL;
   size_t   n   = 512;
   char   * row = malloc(n * sizeof(size_t));
   int      eln = 0; // Element number in matrix.
   double * matrix = malloc(nrow * ncol * sizeof(double));
   if (row == NULL || matrix == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   int slen;
   while ((slen = getline(&row, &n, fp)) > 0) {
      char * el = row;
      for (int i = 0; i < slen; i++) {
         if (row[i] == ',' || row[i] == '\n') {
            row[i] = '\0';
            matrix[eln++] = atof(el);
            el = row + i + 1;
         }
      }
   }

   assert(eln == nrow * ncol);
   free(row);
   fclose(fp);

   return matrix;
}

double *
extractfeats
(
   ChIP_t  * ChIP,
   zerone_t * zerone
)
{
   double * feat = malloc(DIM * sizeof(double));
   double * p = zerone->p;
   unsigned int m = zerone->m;
   unsigned int n = nobs(ChIP);
   unsigned int r = ChIP->r + 1;

   // Add the values of the transition matrix Q to the feature vector...
   for (int i = 0; i < 6; i++) feat[i] = zerone->Q[i];

   // ...the min, max and mean values of the ratios between the values of p...
   double ummin = 1.0;
   double ubmin = 1.0;
   double mbmin = 1.0;
   double ummax = 0.0;
   double ubmax = 0.0;
   double mbmax = 0.0;
   double ummean = 0.0;
   double ubmean = 0.0;
   double mbmean = 0.0;
   for (int i = 2; i < r; i++) {
      double umratio = (p[    i] / (1 - p[0])) / (p[  r + i] / (1 - p[  r]));
      double ubratio = (p[    i] / (1 - p[0])) / (p[2*r + i] / (1 - p[2*r]));
      double mbratio = (p[r + i] / (1 - p[r])) / (p[2*r + i] / (1 - p[2*r]));
      if (umratio < ummin) ummin = umratio;
      if (ubratio < ubmin) ubmin = ubratio;
      if (mbratio < mbmin) mbmin = mbratio;
      if (umratio > ummax) ummax = umratio;
      if (ubratio > ubmax) ubmax = ubratio;
      if (mbratio > mbmax) mbmax = mbratio;
      ummean += umratio;
      ubmean += ubratio;
      mbmean += mbratio;
   }
   feat[ 6] = ummin;
   feat[ 7] = ummean / (r - 2);
   feat[ 8] = ummax;
   feat[ 9] = ubmin;
   feat[10] = ubmean / (r - 2);
   feat[11] = ubmax;
   feat[12] = mbmin;
   feat[13] = mbmean / (r - 2);
   feat[14] = mbmax;

   // ...the mean of the values of phi...
   double meanphi[2] = { 0.0, 0.0 };
   for (int i = 1; i < m; i++) {
      for (int j = 0; j < n; j++) {
         meanphi[i - 1] += zerone->phi[j * m + i];
      }
      feat[14 + i] = meanphi[i - 1] / n;
   }

   // ...and also the mean of the Viterbi path.
   double meanpath = 0.0;
   for (int i = 0; i < n; i++) meanpath += zerone->path[i];
   feat[17] = meanpath / n;

   return feat;
}

double *
zscale
(
   double * feat,
   double * center,
   double * scale
)
{
   double * sfeat = malloc(DIM * sizeof(double));
   // Z-score scaling.
   for (int i = 0; i < DIM; i++) sfeat[i] = (feat[i] - center[i]) / scale[i];
   return sfeat;
}

int
predict
(
   double * feat,
   double * sv,
   double * coefs
)
{
   // Kernel: the exponential of the squared Euclidean distance
   // between the test and support vectors parametrized by gamma.
   double kvals[NSV];
   for (int i = 0; i < NSV; i++) {
      double sum = 0.0;
      for (int j = 0; j < DIM; j++) {
         double dist = feat[j] - sv[i * DIM + j];
         sum += dist * dist;
      }
      kvals[i] = exp(-GAMMA * sum);
   }

   // Inner product between the kernelized test vector
   // and the w hyperplane (here named coefs).
   double label = 0.0;
   for (int i = 0; i < NSV; i++) label += coefs[i] * kvals[i];

   // Add the intercept term of the hyperplane before returning.
  return (label - RHO) + 0.65 > 0 ? 1 : -1;

   // The +0.65 before returning moves the decision boundary
   // towards the negative space by that proportion of the margin width.
   // The classifier is thus more prone to label examples as positive.
}
