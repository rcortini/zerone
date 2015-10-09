#include <stdio.h>
#include <stdlib.h>
#include "predict.h"
#include "svmdata.h"

double *
extractfeat
(
   ChIP_t  * ChIP,
   zerone_t * zerone
)
{
   double * feat = malloc(DIM * sizeof(double));
   if (feat == NULL) {
      fprintf(stderr, "memory error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return NULL;
   }
   double * p = zerone->p;
   unsigned int m = zerone->m;
   unsigned int n = nobs(ChIP);
   unsigned int r = ChIP->r + 1;

   // Add the values of the transition matrix Q to the feature vector...
   int f = m * m - m;
   for (int i = 0; i < f; i++) feat[i] = zerone->Q[i];

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
   feat[f++] = ummin;
   feat[f++] = ummean / (r - 2);
   feat[f++] = ummax;
   feat[f++] = ubmin;
   feat[f++] = ubmean / (r - 2);
   feat[f++] = ubmax;
   feat[f++] = mbmin;
   feat[f++] = mbmean / (r - 2);
   feat[f++] = mbmax;

   // ...the mean of the values of phi...
   if (m <= 1) {
      fprintf(stderr, "you need at least 2 states!");
      return NULL;
   }
   double * meanphi = calloc((size_t) m - 1, sizeof(double));
   if (meanphi == NULL) {
      fprintf(stderr, "memory error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return NULL;
   }
   for (int i = 1; i < m; i++) {
      for (int j = 0; j < n; j++) {
         meanphi[i - 1] += zerone->phi[j * m + i];
      }
      feat[f++] = meanphi[i - 1] / n;
   }
   free(meanphi);

   // ...and also the mean of the Viterbi path.
   double meanpath = 0.0;
   for (int i = 0; i < n; i++) meanpath += zerone->path[i];
   feat[f++] = meanpath / n;

   // TODO: remove this feature!
   feat[f++] = zerone->l;

   return feat;
}

double *
zscale
(
   double * feat
)
{
   double * sfeat = malloc(DIM * sizeof(double));
   if (sfeat == NULL) {
      fprintf(stderr, "memory error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return NULL;
   }
   // Z-score scaling.
   for (int i = 0; i < DIM; i++) sfeat[i] = (feat[i] - center[i]) / scale[i];
   return sfeat;
}

double
predict
(
   ChIP_t  * ChIP,
   zerone_t * zerone
)
{
   double * feat = zscale(extractfeat(ChIP, zerone));

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
   //return (label - RHO) + 0.65 > 0 ? 1 : -1;
   return (label - RHO);// + 0.65;

   // The +0.65 before returning moves the decision boundary
   // towards the negative space by that proportion of the margin width.
   // The classifier is thus more prone to label examples as positive.
}
