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

   // Add the values of the transition matrix Q to the feature vector...
   feat[0] = Z->Q[2 + 0*m];
   feat[1] = Z->Q[2 + 1*m];
   feat[2] = Z->Q[0 + 2*m];

   double ummax = 0.0;
   double ubmax = 0.0;
   double mbmax = 0.0;

   for (int i = 2; i < r; i++) {
      double umratio = (p[i+0*r] / (1 - p[0*r])) /
         (p[i+1*r] / (1 - p[1*r]));
      double ubratio = (p[i+0*r] / (1 - p[0*r])) /
         (p[i+2*r] / (1 - p[2*r]));
      double mbratio = (p[i+1*r] / (1 - p[1*r])) /
         (p[i+2*r] / (1 - p[2*r]));
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
      if (Z->path[j] == 2) {
         ntargets++;
         meanpath++;
         meanphi += Z->phi[j * m + 2];
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
   debug_print("SVM score: %.3f\n", label - RHO);
   return (label - RHO);
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
