#include "utils.h"
#include <stdint.h>

#define U32 uint32_t


//-----------------------------------------------------------------------
// Standard numeric functions.                                           
//-----------------------------------------------------------------------

double
mean
(
   const int *yz,
   int n,
   int r
)
// SYNOPSIS:                                                             
//   Compute the mean of the first column of an integer array. To com-   
//   pute the mean of another column, call as `mean(yz+1, n, r)`, where  
//   1 is for column 2 etc.                                              
//                                                                       
// PARAMETERS:                                                           
//   'yz': (r,n) the integer array                                       
//   'n': row number of the array                                        
//   'r': column number of the array                                     
//                                                                       
// RETURN:                                                               
//   The mean as a 'double'.                                             
{
   double sum = 0.0;
   int n_obs_no_NA = 0;

   for (int k = 0 ; k < n ; k++) {
      // Casting NA to integer gives -2147483648, which is the 
      // largest negative value stored in 'int'. Here I test for
      // NA by wrapping around.
      if (yz[r*k] == INT_MIN) continue;
      sum += yz[r*k];
      n_obs_no_NA++;
   }
   // The result can be 0/0, which is 'nan'.
   return sum / n_obs_no_NA;
}


int *
histsum
(
   const int *yz,
   int n,
   int r
)
// SYNOPSIS:                                                             
//   Compute the histogram of the row-wise sum of an integer array.      
//                                                                       
// INPUT:                                                                
//   The presence of a negative value makes the whole row ignored.       
//                                                                       
// PARAMETERS:                                                           
//   'yz': (r,n) the integer array                                       
//   'n': row number of the array                                        
//   'r': column number of the array                                     
//                                                                       
// RETURN:                                                               
//   A pointer of 'int' to the histogram.                                
{
   int size = 1024;
   int *counts = calloc(size, sizeof(int));
   if (counts == NULL) {
      fprintf(stderr, "memory error (hist 1)\n");
      return NULL;
   }

   int maxval = 0;
   for (int k = 0 ; k < n ; k++) {
      int sum = 0;
      for (int i = 0 ; i < r ; i++){
         if (yz[i+k*r] < 0) {
            sum = -1;
            break;
         }
         sum += yz[i+k*r];
      }
      if (sum < 0) continue;
      if (sum > maxval) maxval = sum;
      if (sum > size-2) {
         int newsize;
         for (newsize = size ; sum > newsize-2 ; newsize *= 2);
         int *newcounts = realloc(counts, newsize * sizeof(int));
         if (newcounts == NULL) {
            fprintf(stderr, "memory error (hist 2)\n");
            return NULL;
         }
         else {
            // Extra memory must be initialized to 0.
            counts = newcounts;
            for (int j = size ; j < newsize ; j++) counts[j] = 0;
            size = newsize;
         }
      }
      counts[sum]++;
   }
   // Add the sentinel.
   counts[maxval+1] = -1;
   return counts;

}


//-----------------------------------------------------------------------
// Indexing of time series using the SHA1 digest.                        
//-----------------------------------------------------------------------

// Globals.
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
   xxhash = malloc(n*sizeof(U32));
   int *addr = malloc(n * sizeof(int));
   for (i = 0 ; i < n ; i++) addr[i] = i;

   // Compute SHA1 digests.
   for (i = 0 ; i < n ; i++) {
      xxhash[i] = XXH32(ts+i*r, r * sizeof(int), 0);
   }

   // Stable sort on hashes.
   qsort(addr, n, sizeof(int), stblcmp);

   int current = 0;
   index[addr[0]] = addr[0];
   for (i = 1 ; i < n ; i++) {
      if (xxhash[addr[i]] == xxhash[current]) {
         index[addr[i]] = current;
      }
      else {
         current = index[addr[i]] = addr[i];
      }
   }
   free(xxhash);
   xxhash = NULL;
   free(addr);

   return index;

}
