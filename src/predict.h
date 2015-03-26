#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // Link with -lm
#include <assert.h>
#include "jahmm.h"

#define NSV 129
#define DIM 19
#define GAMMA 0.015625
#define RHO 0.276603236917580208587
//#define NTESTS 946

double * readmatrix(char *, int, int);
double * extractfeats(ChIP_t *, jahmm_t *);
double * zscale(double *, double *, double *);
int      predict(double *, double *, double *);
