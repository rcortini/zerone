#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

#ifndef _PUBLIC_ZINM_HEADER
#define _PUBLIC_ZINM_HEADER

#define HISTO_INIT_SIZE 128
#define ZINM_MAXITER 32
#define ZINM_TOL 1e-6

#define sq(x) ((x)*(x))

struct zinm_part_t;

typedef struct zinm_par_t zinm_par_t;

struct zinm_par_t {
   size_t r;
   double pi;
   double alpha;
   double p[];
};

double       eval_nb_dfda(double, const tab_t *);
double       eval_nb_f(double, const tab_t *);
double 	     eval_zinm_dfda(double, double, unsigned int);
double 	     eval_zinm_dfdp(double, double, unsigned int, double);
double 	     eval_zinm_dgda(double, double, const tab_t *);
double 	     eval_zinm_f(double, double, unsigned int, double);
double 	     eval_zinm_g(double, double, const tab_t *);
double	     ll_zinm(double, double, double, const tab_t *);
zinm_par_t * mle_nm(size_t *, size_t, size_t);
zinm_par_t * mle_zinm(size_t *, size_t, size_t);
double       nb_est_alpha(tab_t *);
zinm_par_t * new_zinm_par(size_t);
void         zinm_prob (const int *, const int *, const int *, const int *,
                const double *, const double *, const double *,
                int *, const int *, double *);
#endif
