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
