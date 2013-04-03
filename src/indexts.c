#include "indexts.h"

typedef struct {
   unsigned char sum[20];
} sha1sum;

// Globals.
sha1sum *sha1;

int
sha1cmp
(
   sha1sum A,
   sha1sum B
)
// SYNOPSIS:                                                             
//   Compare two SHA1 digests.                                           
{
   for (int i = 0 ; i < 20 ; i++) {
      if (A.sum[i] == B.sum[i]) continue;
      if (A.sum[i]  > B.sum[i]) return  1;
      else                      return -1;
   }
   return 0;
}

int
stblcmp
(
   const void *a,
   const void *b
)
// SYNOPSIS:                                                             
//   Comparison function for stable sort on SHA1 digests. Used to sort   
//   addresses of global pointer 'sha1'.                                 
{
   sha1sum A = sha1[*(int *)a];
   sha1sum B = sha1[*(int *)b];
   switch(sha1cmp(A,B)) {
      case -1: return -1;
      case  1: return  1;
   }
   // SHA1 sums are identical, compare addresses for stable sort.
   if (*(int *)a > *(int *)b)   return  1;
   else                         return -1;
}


int *
indexts
(
   int n,
   int r,
   const int *ts,
   // output //
   int *index
)
{
   int i;
   sha1 = malloc(n*sizeof(sha1sum));
   int *addr = malloc(n * sizeof(int));
   for (i = 0 ; i < n ; i++) addr[i] = i;

   // Compute SHA1 digests.
   for (i = 0 ; i < n ; i++) {
      compute_sha1((const int *) ts+i*r, r*sizeof(int)/sizeof(char),
            sha1[i].sum);
   }

   // Stable sort on SHA1 digests.
   qsort(addr, n, sizeof(int), stblcmp);

   int current = 0;
   index[addr[0]] = addr[0];
   for (i = 1 ; i < n ; i++) {
      if (sha1cmp(sha1[addr[i]], sha1[current]) == 0) {
         index[addr[i]] = current;
      }
      else {
         current = index[addr[i]] = addr[i];
      }
   }
   free(sha1);
   free(addr);

   return index;

}
