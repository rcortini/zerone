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
