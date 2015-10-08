#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#ifndef _PREDICT
#define _PREDICT

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // Link with -lm
#include <assert.h>
#include "zerone.h"
#include "svmdata.h"

double * extractfeat(ChIP_t *, zerone_t *);
double * zscale(double *);
double   predict(ChIP_t *, zerone_t *);

#endif
