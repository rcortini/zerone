#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "xxhash.h"

#ifndef _HMMNB_UTILS
#define _HMMNB_UTILS
int *
indexts
(
   int n,
   int r,
   const int *ts,
   // output //
   int *index
);

int *
histsum
(
   const int *yz,
   int n,
   int r
);

double
mean
(
   const int *yz,
   int n,
   int r
);
#endif
