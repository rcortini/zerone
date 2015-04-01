#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> // Link with -lm
#include <assert.h>
#include "jahmm.h"

#define NSV 200
#define DIM 18
#define GAMMA 0.0625
#define RHO 0.7399845397264748214639
//#define NTESTS 946

double * readmatrix(char *, int, int);
double * extractfeats(ChIP_t *, jahmm_t *);
double * zscale(double *, double *, double *);
int      predict(double *, double *, double *);
