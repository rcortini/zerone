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

#include "debug.h"
#include "utils.h"

#define U32 uint32_t

// Globals variable.
U32 *xxhash;

int
stblcmp
(
   const void *a,
   const void *b
)
// SYNOPSIS:
//   Comparison function for stable sort on hashes. Used to sort
//   addresses of global pointer 'xxhash'.
{
   U32 A = xxhash[*(int *)a];
   U32 B = xxhash[*(int *)b];
   if (A > B) return 1;
   if (A < B) return -1;
   // Hashes are identical, compare addresses for stable sort.
   if (*(int *)a > *(int *)b)   return  1;
   return -1;
}


int
indexts
(
         int n,
         int r,
   const int *ts,
         int *index
)
// SYNOPSIS:
// Index the time series using xxhash and return the index of the
// first all-0 observation.
{
   int i;
   xxhash = malloc(n*sizeof(U32));
   int *addr = malloc(n * sizeof(int));
   if (xxhash == NULL || addr == NULL) {
      debug_print("%s", "memory error\n");
      return -1;
   }

   for (i = 0 ; i < n ; i++) addr[i] = i;

   // Compute xxhash digests.
   for (i = 0 ; i < n ; i++) {
      xxhash[i] = XXH32(ts+i*r, r*sizeof(int), 0);
   }

   // Compute the xxhash digest of all 0s. We will need it
   // to return the index of the small such observation.
   int *all0 = calloc(r, sizeof(int));
   if (all0 == NULL) {
      debug_print("%s", "memory error\n");
      return -1;
   }
   U32 xxhash0 = XXH32(all0, r*sizeof(int), 0);
   free(all0);

   // Stable sort array indices on digests order. Stability is
   // important because we want the first occurrence of an
   // observation to point to its own index, and all subsequent
   // occurrences to point to it as well. If the reference index
   // was not the first occurence in array order, it would cause
   // difficulties for the purpose of computing emission
   // probabilities.
   qsort(addr, n, sizeof(int), stblcmp);

   int current = 0;
   int index0 = -1;
   index[addr[0]] = addr[0];
   if (xxhash[addr[0]] == xxhash0) index0 = addr[0];
   for (i = 1 ; i < n ; i++) {
      if (xxhash[addr[i]] == xxhash[current]) {
         index[addr[i]] = current;
      }
      else {
         current = index[addr[i]] = addr[i];
         if (xxhash[addr[i]] == xxhash0) index0 = addr[i];
      }
   }

   free(xxhash);
   free(addr);

   return index0;

}
