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

#include <stdio.h>

#include "unittest.h"
#include "predict.c"

void
test_predict
(void)
{

   double case1[5] = {
      0.0, 0.0, 0.0, 0.0, 0.0
   };

   test_assert(fabs(predict(zscale(case1)) + 2.101421) < 1e-4);

   double case2[5] = {
      0.1, 0.2, 0.3, 0.4, 0.5
   };

   test_assert(fabs(predict(zscale(case2)) + 0.164167) < 1e-4);

   double case3[5] = {
      0.04336188,  3.972419, 6.960893e-03, 0.15916909, 0.8062380
   };

   test_assert(fabs(predict(zscale(case3)) - 0.319410) < 1e-4);

}

// Test cases for export.
const test_case_t test_cases_predict[] = {
   {"predict/predict",         test_predict},
   {NULL, NULL},
};

