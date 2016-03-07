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
   zerone_t * Z,
   double   * feat
)
{
   ChIP_t *ChIP = Z->ChIP;
   double * p = Z->p;
   unsigned int m = Z->m;
   unsigned int n = nobs(ChIP);
   // XXX This bit is confusing. Needs to be discussed. XXX //
   unsigned int r = ChIP->r + 1;

   // 'map' contains the permutation of the states.
   int *map = Z->map;

   // Add the values of the transition matrix Q to the feature vector...
   feat[0] = Z->Q[map[2] + map[0]*m];
   feat[1] = Z->Q[map[2] + map[1]*m];
   feat[2] = Z->Q[map[0] + map[2]*m];

   double ummax = 0.0;
   double ubmax = 0.0;
   double mbmax = 0.0;

   for (int i = 2; i < r; i++) {
      double umratio = (p[i+map[0]*r] / (1 - p[map[0]*r])) /
         (p[i+map[1]*r] / (1 - p[map[1]*r]));
      double ubratio = (p[i+map[0]*r] / (1 - p[map[0]*r])) /
         (p[i+map[2]*r] / (1 - p[map[2]*r]));
      double mbratio = (p[i+map[1]*r] / (1 - p[map[1]*r])) /
         (p[i+map[2]*r] / (1 - p[map[2]*r]));
      if (umratio > ummax) ummax = umratio;
      if (ubratio > ubmax) ubmax = ubratio;
      if (mbratio > mbmax) mbmax = mbratio;
   }

   feat[3] = ummax;
   feat[4] = ubmax;
   feat[5] = mbmax;

   double meanphi = 0;
   double meanpath = 0;
   int ntargets = 0;
   for (int j = 0; j < n; j++) {
      if (Z->path[j] == map[2]) {
         ntargets++;
         meanpath++;
         meanphi += Z->phi[j * m + map[2]];
      }
   }

   feat[6] = (ntargets == 0) ? 0.0 : meanphi / ntargets;
   feat[7] = meanpath / n;

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
   debug_print("SVM score: %.3f\n", label - RHO);
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
