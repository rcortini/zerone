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

   double case1[8] = {
      0.530507921, 0.043361879, 2.851383e-02, 1.4379088,
      0.4816441, 0.3493823, 0.8808826, 6.960893e-03
   };

   test_assert(fabs(predict(zscale(case1)) - 1.00014284) < 1e-3);

   double case2[8] = {
      0.006321302, 0.006369056, 2.620211e-26, 0.9321181,
      4.2291074, 6.2477326, 0.9995843, 4.922989e-05
   };

   test_assert(fabs(predict(zscale(case2)) + 1.00024599) < 1e-3);

   double case3[8] = {
      0.043084589, 0.478401608, 1.193816e-04, 0.7304060,
      0.5678939, 1.4512983, 0.9649980, 2.195285e-03
   };

   test_assert(fabs(predict(zscale(case3)) - 0.06892158) < 1e-3);

}
