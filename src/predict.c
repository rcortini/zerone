/* Copyright 2015, 2016 Pol Cusco and Guillaume Filion

   This file is part of Zerone.

   Zerone is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Zerone is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Zerone. If not, see <http://www.gnu.org/licenses/>.
*/

#define _GNU_SOURCE
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "debug.h"
#include "predict.h"
#include "svmdata.h"

#define SQ(a) (a)*(a)

// Below are the global constants declared in svmdata.h.
// GAMMA, RHO, CENTER, SCALE, SV, COEFFS, DIM, NSV

double *
extract_features
(
   zerone_t * Z,
   double   * features
)
{
   ChIP_t *ChIP = Z->ChIP;
   double * p = Z->p;
   const unsigned int m = Z->m;
   const unsigned int n = nobs(ChIP);
   const unsigned int r = ChIP->r;

   // Feature 0: tansition from "top" to "mid".
   features[0] = Z->Q[2 + 1*m];

   // Feature 1: smallest "top" to "mid" signal ratio.
   features[1] = 10;
   for (int i = 2; i < r+1; i++) {
      double ratio = (p[i+2*(r+1)] / p[0+2*(r+1)]) /
         (p[i+1*(r+1)] / p[0+1*(r+1)]);
      if (ratio < features[1]) features[1] = ratio;
   }

   double *mean_yes = calloc(r, sizeof(double));
   double *mean_no = calloc(r, sizeof(double));
   double *prod_yes = calloc(r*r, sizeof(double));
   double *var = calloc(r, sizeof(double));

   if (mean_yes == NULL || mean_no == NULL ||
               var == NULL || prod_yes == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      // Fill feature vector with NAs.
      for (int i = 0 ; i < DIM ; i++) features[i] = 0.0/0.0;
      goto clean_and_return;
   }

   double n_yes = 0.0;

   // Scan data and collect intermediate values to
   // compute means, variances and covariances.
   for (int i = 0; i < n; i++) {
      if (Z->path[i] == 2) {
         n_yes++;
         for (int j = 1 ; j < r ; j++) {
            mean_yes[j] += ChIP->y[j+i*r];
            for (int k = j ; k < r ; k++) {
               prod_yes[j+k*r] +=
                  ChIP->y[j+i*r] * ChIP->y[k+i*r];
            }
         }
      }
      else {
         for (int j = 1 ; j < r ; j++) {
            mean_no[j] += ChIP->y[j+i*r];
         }
      }
      for (int j = 1 ; j < r ; j++) {
         var[j] += (ChIP->y[j+i*r]) * (ChIP->y[j+i*r]);
      }
   }

   double n_no = n - n_yes;

   if (n_yes < 1 || n_no < 1) {
      // Limit case: avoid division by 0
      // and set following features to 0.
      features[2] = features[3] = features[4] = 0.0;
      goto clean_and_return;
   }

   for (int j = 1 ; j < r ; j++) {
      var[j] = var[j] / n - SQ((mean_no[j] + mean_yes[j]) / n);
      mean_no[j] /= n_no;
      mean_yes[j] /= n_yes;
      for (int k = j ; k < r ; k++) {
         prod_yes[j+k*r] /= n_yes;
      }
   }

   // Feature 2: average number of targets.
   features[2] = n_yes / n;

   // Feature 3: minimum explained variance.
   features[3] = 1.0;
   
   for (int j = 1 ; j < r ; j++) {
      double v = n_yes*n_no * SQ(mean_yes[j]-mean_no[j]) / (var[j]*SQ(n));
      if (v < features[3]) features[3] = v;
   }

   // Feature 4: minimum correlation on targets.
   // This cannot be computed if only one profile is available.
   if (r < 3) {
      features[4] = 0.0/0.0;
      goto clean_and_return;
   }

   features[4] = 1.0;

   for (int j = 1 ; j < r ; j++) {
   for (int k = j+1 ; k < r ; k++) {
      double c = (prod_yes[j+k*r] - mean_yes[j]*mean_yes[k]) /
         sqrt((prod_yes[j+j*r] - SQ(mean_yes[j])) *
               (prod_yes[k+k*r] - SQ(mean_yes[k])));
      if (c < features[4]) features[4] = c;
   }
   }

clean_and_return:
   free(mean_yes);
   free(mean_no);
   free(prod_yes);
   free(var);

   return features;

}

double *
zscale
(
   double * features
)
{
   // Z-score scaling.
   debug_print("%s", "scaled features:\n");
   for (int i = 0; i < DIM; i++) {
      features[i] = (features[i] - CENTER[i]) / SCALE[i];
      debug_print("| %d: %.3f\n", i, features[i]);
   }

   return features;

}

double
predict
(
   double *scaled_features
)
{

   // Kernel: the exponential of the squared Euclidean distance
   // between the test and support vectors parametrized by gamma.
   double kvals[NSV];
   for (int i = 0; i < NSV; i++) {
      double sum = 0.0;
      for (int j = 0; j < DIM; j++) {
         double dist = scaled_features[j] - SV[i + j * NSV];
         sum += dist * dist;
      }
      kvals[i] = exp(-GAMMA * sum);
   }

   // Inner product between the kernelized test vector
   // and the w hyperplane (here named COEFS).
   double label = 0.0;
   for (int i = 0; i < NSV; i++) label += COEFS[i] * kvals[i];

   // Add the hyperplane intercept.
   return (label - RHO);

}

double
zerone_qc
(
   zerone_t * zerone,
   double   * features
)
{

   // Extract and store unscaled features.
   extract_features(zerone, features);

   // Copy features and scale.
   double scaled_features[DIM] = {0};
   memcpy(scaled_features, features, DIM * sizeof(double));
   zscale(scaled_features);

   return predict(scaled_features);

}
