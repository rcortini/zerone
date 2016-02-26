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

#include "unittest.h"
#include "utils.c"

void
test_indexts
(void)
{
   int index[8];

   int ts_1[24] = {
      0, 0, 0,
      0, 0, 1,
      0, 1, 0,
      1, 0, 1,
      0, 0, 1,
      0, 0, 1,
      0, 1, 0,
      0, 0, 0,
   };

   int ts_2[27] = {
      -100, 0, 1,
      -1,   0, 100,
      -1,   0, 1,
      -1,   0, 100,
      -100, 0, 1,
      -1,   0, 1,
      -100, 0, 100,
       0,   0, 0,
      -100, 0, 100,
   };

   int expected_1[8] = {0,1,2,3,1,1,2,0};
   int expected_2[9] = {0,1,2,1,0,2,6,7,6};

   test_assert(indexts(8, 3, ts_1, index) == 0);
   for (int i = 0 ; i < 8 ; i++) {
      test_assert(index[i] == expected_1[i]);
   }

   test_assert(indexts(9, 3, ts_2, index) == 7);
   for (int i = 0 ; i < 8 ; i++) {
      test_assert(index[i] == expected_2[i]);
   }

}
