#include <emmintrin.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef _ZINB_HEADER
#define _ZINB_HEADER

// Control parameters.
#define HISTO_INIT_SIZE 128
#define ZINB_MAXITER 32
#define ZINB_TOL 1e-6


typedef struct zinb_par_t zinb_par_t;
typedef struct histo_t histo_t;
typedef struct tab_t tab_t;
typedef void (*zinb_err_handler_t) (const char *, const char *, int);

struct zinb_par_t {
   double   a;
   double   p;
   double   pi;
};

struct histo_t {
   size_t size;
   unsigned int num[];
};

struct tab_t {
   size_t size;
   unsigned int *num;
   unsigned int val[];
};

histo_t * new_histo (void);
int       histo_push (histo_t **, size_t);
tab_t   * compress_histo (histo_t *);
// Compute the maximum likelihood estimates for NB and
// ZINB distributions, and return the paramters.
zinb_par_t * mle_nb   (int *, size_t);
zinb_par_t * mle_zinb (int *, size_t);
// Change the default error handler. Pass a NULL argument
// to rest to default behavior.
void         set_zinb_err_handler(zinb_err_handler_t);
#endif
