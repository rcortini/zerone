#include <stdlib.h>

int
stablecmp(
   const void *a,
   const void *b
)
{
   int keya = array[*(int *)a];
   int keyb = array[*(int *)b];

   switch((keya > keyb) - (keya < keyb)) {
      case -1: return -1;
      case  1: return  1;
      default: return *(int *)a > *(int *b) - *(int *)a < *(int *b);
   }
}

typedef int (*cmpfn)(const void*, const void*);
cmpfn cmpFactory(int *array) {
   cmpfn fnPtr = &stablecmp;
   return fnPtr;
}

