#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include "predict.h"
#include "svmdata.h"

// Below are the global constants declared in svmdata.h.
// GAMMA, RHO, CENTER, SCALE, SV, COEFFS, DIM, NSV

double *
extractfeat
(
   zerone_t * zerone,
   double   * feat
)
{
   ChIP_t *ChIP = zerone->ChIP;
   double * p = zerone->p;
   unsigned int m = zerone->m;
   unsigned int n = nobs(ChIP);
   // XXX This bit is confusing. Needs to be discussed. XXX //
   unsigned int r = ChIP->r + 1;

   // 'map' contains the permutation of the states.
   int *map = zerone->map;

   // Add the values of the transition matrix Q to the feature vector...
   feat[0] = zerone->Q[map[0] + map[1]*m];
   feat[1] = zerone->Q[map[1] + map[1]*m];
   feat[2] = zerone->Q[map[2] + map[1]*m];
   feat[3] = zerone->Q[map[0] + map[2]*m];
   feat[4] = zerone->Q[map[1] + map[2]*m];
   feat[5] = zerone->Q[map[2] + map[2]*m];

   int f = m * m - m;

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
//      double umratio = (p[    i] / (1 - p[0])) / (p[  r + i] / (1 - p[  r]));
//      double ubratio = (p[    i] / (1 - p[0])) / (p[2*r + i] / (1 - p[2*r]));
//      double mbratio = (p[r + i] / (1 - p[r])) / (p[2*r + i] / (1 - p[2*r]));
      double umratio = (p[i+map[0]*r] / (1 - p[map[0]*r])) /
         (p[i+map[1]*r] / (1 - p[map[1]*r]));
      double ubratio = (p[i+map[0]*r] / (1 - p[map[0]*r])) /
         (p[i+map[2]*r] / (1 - p[map[2]*r]));
      double mbratio = (p[i+map[1]*r] / (1 - p[map[1]*r])) /
         (p[i+map[2]*r] / (1 - p[map[2]*r]));
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
   // XXX This may be excluded when 'm' is set to 3. XXX //
   if (m <= 1) {
      debug_print("%s", "'m' must be at least 2\n");
      return NULL;
   }
   double * meanphi = calloc((size_t) m - 1, sizeof(double));
   if (meanphi == NULL) {
      debug_print("%s", "memory error\n");
      return NULL;
   }
   for (int i = 1; i < m; i++) {
      for (int j = 0; j < n; j++) {
         meanphi[i - 1] += zerone->phi[j * m + map[i]];
      }
      feat[f++] = meanphi[i - 1] / n;
   }
   free(meanphi);

   // ...and also the mean of the Viterbi path.
   double meanpath = 0.0;
   for (int i = 0; i < n; i++) meanpath += map[zerone->path[i]];
   feat[f++] = meanpath / n;

   for (int i = 0 ; i < DIM ; i++) {
      debug_print("feature %d: %.3f\n", i, feat[i]);
   }

   return feat;

}

double *
zscale
(
   double * feat
)
{
   // Z-score scaling.
   for (int i = 0; i < DIM; i++) {
      feat[i] = (feat[i] - CENTER[i]) / SCALE[i];
      debug_print("scaled feature %d: %.3f\n", i, feat[i]);
   }

   return feat;

}

double
predict
(
   double *feat
)
{

   // Kernel: the exponential of the squared Euclidean distance
   // between the test and support vectors parametrized by gamma.
   double kvals[NSV];
   for (int i = 0; i < NSV; i++) {
      double sum = 0.0;
      for (int j = 0; j < DIM; j++) {
         double dist = feat[j] - SV[i + j * NSV];
         sum += dist * dist;
      }
      kvals[i] = exp(-GAMMA * sum);
   }

   // Inner product between the kernelized test vector
   // and the w hyperplane (here named COEFS).
   double label = 0.0;
   for (int i = 0; i < NSV; i++) label += COEFS[i] * kvals[i];

   // Add the intercept term of the hyperplane before returning.
   //return (label - RHO) + 0.65 > 0 ? 1 : -1;
   return (label - RHO);// + 0.65;

   // The +0.65 before returning moves the decision boundary
   // towards the negative space by that proportion of the margin width.
   // The classifier is thus more prone to label examples as positive.
}

double
zerone_predict
(
   zerone_t * zerone
)
{

   double feat[DIM] = {0};
   zscale(extractfeat(zerone, feat));

   return predict(feat);

}
