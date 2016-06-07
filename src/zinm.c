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

#include "zinm.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define sq(x) ((x)*(x))

#define handle_error() (*zinb_err_handler)(__FILE__, __func__, __LINE__)


// Declaration of private types. //
//struct histo_t;
//struct tab_t;
//
//typedef struct histo_t histo_t;
//typedef struct tab_t tab_t;
//
//struct histo_t {
//   size_t size;
//   unsigned int num[];
//};
//
//struct tab_t {
//   size_t size;
//   unsigned int *num;
//   unsigned int val[];
//};


// Declaration of private functions. //
//tab_t   * compress_histo (histo_t *);
void      default_warning (const char *, const char *, int);
double    digamma (double);
double    eval_nb_dfda (double, const tab_t *);
double    eval_nb_f (double, const tab_t *);
double    eval_zinb_dfda (double, double, unsigned int);
double    eval_zinb_dfdp (double, double, unsigned int, double);
double    eval_zinb_dgda (double, double, const tab_t *);
double    eval_zinb_f (double, double, unsigned int, double);
double    eval_zinb_g (double, double, const tab_t *);
//int       histo_push (histo_t **, size_t);
double    ll_zinb (double, double, double, const tab_t *);
double    nb_est_alpha (tab_t *);
//histo_t * new_histo (void);
tab_t   * tabulate (int *, unsigned int);
double    trigamma (double);

// Global variables. //
zinb_err_handler_t zinb_err_handler = default_warning;


// Error handling functions. //
void
default_warning
(
   const char * file,
   const char * func,
         int    line
)
{
   fprintf(stderr, "memory error %s:%d (%s)\n", file, line, func);
}


void
set_zinb_err_handler
(
   zinb_err_handler_t h
)
{
   zinb_err_handler = (h == NULL) ? default_warning : h;
}


// High level exported functions 'mle_nb()' and 'mle_zinm()'.

zinb_par_t *
mle_nb
(
   int *x,
   size_t nobs
)
// SYNOPSIS:
//   Return the maximum likelihood estimate of the Negative Binomial
//   distribution on the sample. Only positive integer values are
//   counted, negative values (used for NA) are ignored.
//
// PARAMETERS:
//   x: observed sample
//   nobs: sample size
//
// RETURN:
//   A 'zinb_par_t' pointer with estimated paraters. Note that the
//   mixture parameter is always set to 1.0.
//
// SIDE EFFECTS:
//   None.
{

   // NB: Negative values are excluded from tabulation.
   tab_t *tab = tabulate(x, nobs);
   if (tab == NULL) return NULL;

   double a = nb_est_alpha(tab);

   // Return NULL if failed.
   if (a < 0 || a != a) {
      free(tab);
      return NULL;
   }

   // Allocate return struct.
   zinb_par_t *par = calloc(1, sizeof(zinb_par_t));
   if (par == NULL) {
      free(tab);
      handle_error();
      return NULL;
   }

   double sum = 0.0;
   unsigned int nona = 0;
   for (size_t i = 0 ; i < tab->size ; i++) {
      nona += tab->num[i];
      sum += tab->val[i]*tab->num[i];
   }

   free(tab);

   par->a = a;
   par->p = a / (a + (sum/nona));
   par->pi = 1.0;

   return par;

}

zinb_par_t *
mle_zinb
(
   int *x,
   size_t nobs
)
// SYNOPSIS:
//   Return the maximum likelihood estimate of the Zero-Inflated
//   Negative Binomial distribution on the sample. Only positive
//   integer values are counted, negative values (used for NA)
//   are ignored.
//
// PARAMETERS:
//   x: observed sample
//   nobs: sample size
//
// RETURN:
//   A 'zinb_par_t' pointer with estimated paraters.
//
// SIDE EFFECTS:
//   None.
{

   tab_t *tab = tabulate(x, nobs);
   if (tab == NULL) return NULL;

   double sum = 0.0;
   unsigned int nona = 0;
   for (size_t i = 0 ; i < tab->size ; i++) {
      sum += tab->val[i]*tab->num[i];
      nona += tab->num[i];
   }

   // Extract the number of zero observaions.
   const unsigned int z0 = tab->val[0] == 0 ? tab->num[0] : 0;

   double deficit[11] = {0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1};
   double init_a[12] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1};
   double init_p[12] = {0,0,0,0,0,0,0,0,0,0,0,.5};

   // Deplete some 0s from the observations, compute alpha
   // and p0 with standard negative binomial estimation
   // and keep the values to be used as initial conditions.
   for (size_t i = 0 ; i < 11 ; i++) {
      if (tab->val[0] == 0) tab->num[0] = z0 * (1-deficit[i]);
      double newmean = sum / nona / (1.0 - z0*deficit[i]/nobs);
      double alpha = nb_est_alpha(tab);
      init_a[i] = alpha;
      init_p[i] = alpha / (alpha + newmean);
   }

   // Reset 'tab'.
   if (tab->val[0] == 0) tab->num[0] = z0;

   zinb_par_t *par = calloc(1, sizeof(zinb_par_t));
   if (par == NULL) {
      free(tab);
      handle_error();
      return NULL;
   }

   // Try initial conditions. Number 12 is a safety in case
   // all the rest failed during the first phase.
   double max_loglik = -1.0/0.0;
   for (size_t i = 0 ; i < 12 ; i++) {

      // Skip failures.
      if (init_a[i] < 0 || init_a[i] != init_a[i]) continue;

      double a = init_a[i];
      double p = init_p[i];

      double grad;
      unsigned int iter = 0;

      double f = eval_zinb_f(a, p, nobs-z0, sum);
      double g = eval_zinb_g(a, p, tab);

      // Newton-Raphson iterations.
      while ((grad = f*f+g*g) > sq(ZINB_TOL) && iter++ < ZINB_MAXITER) {

         double dfda, dfdp, dgda, dgdp;
         dfda = dgdp = eval_zinb_dfda(a, p, nobs-z0);
         dfdp = eval_zinb_dfdp(a, p, nobs-z0, sum);
         dgda = eval_zinb_dgda(a, p, tab);

         double denom = dfdp*dgda - dfda*dgdp;
         double da = (f*dgdp - g*dfdp) / denom;
         double dp = (g*dfda - f*dgda) / denom;
         // Maintain 'a' and 'p' in their domain of definition.
         while (a+da < 0 || p+dp < 0 || p+dp > 1) {
            da /= 2;
            dp /= 2;
         }
         f = eval_zinb_f(a+da, p+dp, nobs-z0, sum);
         g = eval_zinb_g(a+da, p+dp, tab);
         // Backtrack if necessary.
         for (int j = 0 ; j < ZINB_MAXITER && f*f+g*g > grad ; j++) {
            da /= 2;
            dp /= 2;
            f = eval_zinb_f(a+da, p+dp, nobs-z0, sum);
            g = eval_zinb_g(a+da, p+dp, tab);
         }

         a = a+da;
         p = p+dp;

      }

      double pi = (nobs-z0) / (1-pow(p,a)) / nobs;
      if (pi > 1) pi = 1.0;
      if (pi < 0) pi = 0.0;
      // Compute the log-likelihood. The update of 'pi' above
      // may yield a non optimal solution.
      double loglik = ll_zinb(a, p, pi, tab);
      if (loglik > max_loglik) {
         max_loglik = loglik;
         par->a = a;
         par->p = p;
         par->pi = pi;
      }

   }

   free(tab);

   return par;

}



// ---- Private functions (used only internally) ---- //


// The functions below are computational intermidates required
// for the Newton-Raphson iterations performed in 'mle_nb()' and
// 'mle_zinb()'.

double
eval_nb_f
(
         double a,
   const tab_t *tab
)
{

   double retval;
   double prev;
   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;

   size_t nobs = num[0];
   double mean = num[0] * val[0];

   prev = digamma(a + val[0]);
   retval = num[0] * prev;
   // Iterate over the occurrences and compute the new value
   // of digamma either by the recurrence relation, or by
   // a new call to 'digamma()', whichever is faster.
   for (size_t i = 1 ; i < tab->size ; i++) {
      nobs += num[i];
      mean += num[i] * val[i];
      prev = (val[i] - val[i-1] == 1) ?
         prev + 1.0 / (a-1 + val[i]) :
         digamma(a + val[i]);
      retval += num[i] * prev;
   }

   mean /= nobs;
   retval += nobs*(log(a) - digamma(a) - log(a + mean));

   return retval;

}

double
eval_nb_dfda
(
   double a,
   const tab_t *tab
)
{

   double retval;
   double prev;
   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;

   size_t nobs = num[0];
   double mean = num[0] * val[0];

   prev = trigamma(a + val[0]);
   retval = num[0] * prev;
   // Iterate over the occurrences and compute the new value
   // of trigamma either by the recurrence relation, or by
   // a new call to 'trigamma()', whichever is faster.
   for (size_t i = 1 ; i < tab->size ; i++) {
      nobs += num[i];
      mean += num[i] * val[i];
      prev = (val[i] - val[i-1] == 1) ?
         prev - 1.0 / sq(a-1 +val[i]) :
         trigamma(a + val[i]);
      retval += num[i] * prev;
   }

   mean /= nobs;
   retval += nobs*(mean/(a*(a+mean)) - trigamma(a));

   return retval;

}

double
eval_zinb_f
(
   double a,
   double p,
   unsigned int   nz,
   double sum
)
{
   return a*nz / (p*(1-pow(p,a))) - sum/(1-p);
}

double
eval_zinb_g
(
         double   a,
         double   p,
   const tab_t  * tab
)
{

   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;

   unsigned int nz = 0;
   double retval = 0.0;
   double prev = digamma(a + val[0]);

   if (val[0] > 0) {
      retval += num[0] * prev;
      nz += num[0];
   }

   // Iterate over the occurrences and compute the new value
   // of digamma either by the recurrence relation, or by
   // a new call to 'digamma()', whichever is faster.
   const size_t imin = val[0] == 0 ? 1 : 0;
   for (size_t i = imin ; i < tab->size ; i++) {
      nz += num[i];
      prev = (val[i] - val[i-1] == 1) ?
         prev + 1.0 / (a-1 + val[i]) :
         digamma(a + val[i]);
      retval += num[i] * prev;
   }

   retval += nz*(log(p) / (1-pow(p,a)) - digamma(a));
   return retval;

}

double
eval_zinb_dfda
(
   double a,
   double p,
   unsigned int   nz
)
{
   const double ppa = pow(p,a);
   return nz*(1-ppa+a*ppa*log(p)) / (p*sq(1-ppa));
}

double
eval_zinb_dfdp
(
   double a,
   double p,
   unsigned int   nz,
   double m
)
{
   const double ppa = pow(p,a);
   return -(nz*a*(1-(a+1)*ppa) / sq(p*(1-ppa)) + m/sq(1-p));
}

double
eval_zinb_dgda
(
         double   a,
         double   p,
   const tab_t  * tab
)
{

   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;
   const double ppa = pow(p,a);

   unsigned int nz = 0;
   double retval = 0.0;
   double prev = trigamma(a + val[0]);

   if (val[0] > 0) {
      retval += num[0] * prev;
      nz += num[0];
   }

   // Iterate over the occurrences and compute the new value
   // of digamma either by the recurrence relation, or by
   // a new call to 'trigamma()', whichever is faster.
   const size_t imin = val[0] == 0 ? 1 : 0;
   for (size_t i = imin ; i < tab->size ; i++) {
      nz += num[i];
      prev = (val[i] - val[i-1] == 1) ?
         prev - 1.0 / sq(a-1 + val[i]) :
         trigamma(a + val[i]);
      retval += num[i] * prev;
   }

   retval += nz*(sq(log(p))*ppa / sq(1-ppa) - trigamma(a));
   return retval;

}

// This function compute the log-likelihood of a ZINB distribution
// on a tabulated sample. It is used in 'mle_zinb()' to make sure
// that the solution is indeed the maximum likelihood estimate.

double
ll_zinb
(
         double   a,
         double   p,
         double   pi,
   const tab_t  * tab
)
{

   // Convenience variables.
   const unsigned int *val = tab->val;
   const unsigned int *num = tab->num;
   const unsigned int z0 = val[0] == 0 ? num[0] : 0;
   const double logp_ = log(1-p);

   unsigned int nobs = z0;
   double retval = z0*log(pi*pow(p,a) + 1-pi);

   const size_t imin = val[0] == 0 ? 1 : 0;
   for (size_t i = imin ; i < tab->size ; i++) {
      retval += num[i] * (lgamma(a+val[i]) + val[i]*logp_);
      nobs += num[i];
   }
   // Ignore constant factorial terms.
   retval += (nobs-z0) * (a*log(p) - lgamma(a) + log(pi));

   return retval;

}


// This function finds the optimum value of alpha (the shape parameter
// of the Negative Binomial distribution) by the Newton-Raphson method.

double
nb_est_alpha
(
   tab_t *tab
)
{

   // Find upper and lower bouds for a(lpha).
   double a = 1.0;
   double a_lo;
   double a_hi;
   if (eval_nb_f(a, tab) < 0) {
      a /= 2;
      for (int i = 0 ; i < ZINB_MAXITER ; i++) {
         if (eval_nb_f(a, tab) > 0) break;
         a /= 2;
      }
      a_lo = a;
      a_hi = a*2;
   }
   else {
      a *= 2;
      for (int i = 0 ; i < ZINB_MAXITER ; i++) {
         if (eval_nb_f(a, tab) < 0) break;
         a *= 2;
      }
      a_lo = a/2;
      a_hi = a;
   }

   // Input is pathological.
   if (a_lo > 128) return -1.0;

   // Newton-Raphson iterations.
   double new_a = (a_lo + a_hi) / 2;
   for (int i = 0 ; i < ZINB_MAXITER ; i++) {
      a = (new_a < a_lo || new_a > a_hi) ?
         (a_lo + a_hi) / 2 :
         new_a;
      double f = eval_nb_f(a, tab);
      if   (f > 0) a_lo = a;
      else         a_hi = a;
      if ((a_hi - a_lo) < ZINB_TOL)  return a;
      double dfda = eval_nb_dfda(a, tab);
      new_a = a - f / dfda;
   }

   return -1.0;

}

// The functions below are used to manipulate histograms, which are
// used to speed up the computations.

tab_t *
tabulate
(
   int *x,
   unsigned int nobs
)
{

   histo_t *histo = new_histo();
   if (histo == NULL) {
      return NULL;
   }
   for (size_t i = 0 ; i < nobs ; i++) {
      // Skip negative values (used for NA).
      if (x[i] < 0) continue;
      if (histo_push(&histo, x[i])) {
         free(histo);
         return NULL;
      }
   }

   tab_t *tab = compress_histo(histo);
   free(histo);
   return tab;

}

histo_t *
new_histo
(void)
{

   size_t initsize = sizeof(histo_t) + HISTO_INIT_SIZE * sizeof(int);
   histo_t *new = calloc(1, initsize);
   if (new == NULL) {
      handle_error();
      return NULL;
   }
   new->size = HISTO_INIT_SIZE;

   return new;

}

int
histo_push
(
   histo_t **histo_addr,
   size_t val
)
{

   // Convenience variable.
   histo_t *histo = *histo_addr;
   if (val >= histo->size) {
      size_t newsize = 2*val * (sizeof(int));
      histo_t *new = realloc(histo, sizeof(histo_t) + newsize);
      if (new == NULL) {
         handle_error();
         return 1;
      }
      *histo_addr = histo = new;
      size_t added_size = (2*val - histo->size) * sizeof(int);
      memset(histo->num + histo->size, 0, added_size);
      histo->size = 2*val;
   }

   histo->num[val]++;
   return 0;

}

tab_t *
compress_histo
(
   histo_t *histo
)
{

   size_t size = 0;
   for (size_t i = 0 ; i < histo->size ; i++) {
      size += (histo->num[i] > 0);
   }
   tab_t *new = malloc(sizeof(tab_t) + 2*size * sizeof(int));
   if (new == NULL) {
      handle_error();
      return NULL;
   }
   new->size = size;
   new->num = new->val + size;

   size_t j = 0;
   for (size_t i = 0 ; i < histo->size ; i++) {
      if (histo->num[i] > 0) {
         new->val[j] = i;
         new->num[j] = histo->num[i];
         j++;
      }
   }

   return new;

}

// The code below was copied from the following link
// http://pmtksupport.googlecode.com/svn/trunk/lightspeed2.3/util.c
// Written by Tom Minka (unless otherwise noted).

/* The digamma function is the derivative of gammaln.

   Reference:
    J Bernardo,
    Psi ( Digamma ) Function,
    Algorithm AS 103,
    Applied Statistics,
    Volume 25, Number 3, pages 315-317, 1976.

    From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
    (with modifications for negative numbers and extra precision)
*/


double
digamma
(
   double x
)
{
  double neginf = -INFINITY;
  static const double c = 12,
    digamma1 = -0.57721566490153286,
    trigamma1 = 1.6449340668482264365, /* pi^2/6 */
    s = 1e-6,
    s3 = 1./12,
    s4 = 1./120,
    s5 = 1./252,
    s6 = 1./240,
    s7 = 1./132;
    //s8 = 691./32760,
    //s9 = 1./12,
    //s10 = 3617./8160;
  double result;
  /* Illegal arguments */
  if((x == neginf) || isnan(x)) {
    return NAN;
  }
  /* Singularities */
  if((x <= 0) && (floor(x) == x)) {
    return neginf;
  }
  /* Negative values */
  /* Use the reflection formula (Jeffrey 11.1.6):
   * digamma(-x) = digamma(x+1) + pi*cot(pi*x)
   *
   * This is related to the identity
   * digamma(-x) = digamma(x+1) - digamma(z) + digamma(1-z)
   * where z is the fractional part of x
   * For example:
   * digamma(-3.1) = 1/3.1 + 1/2.1 + 1/1.1 + 1/0.1 + digamma(1-0.1)
   *               = digamma(4.1) - digamma(0.1) + digamma(1-0.1)
   * Then we use
   * digamma(1-z) - digamma(z) = pi*cot(pi*z)
   */
  if(x < 0) {
    return digamma(1-x) + M_PI / tan(-M_PI*x);
  }
  /* Use Taylor series if argument <= S */
  if(x <= s) return digamma1 - 1/x + trigamma1*x;
  /* Reduce to digamma(X + N) where (X + N) >= C */
  result = 0;
  while(x < c) {
    result -= 1/x;
    x++;
  }
  /* Use de Moivre's expansion if argument >= C */
  /* This expansion can be computed in Maple via asympt(Psi(x),x) */
  if(x >= c) {
    double r = 1/x;
    result += log(x) - 0.5*r;
    r *= r;
    result -= r * (s3 - r * (s4 - r * (s5 - r * (s6 - r * s7))));
  }
  return result;
}

/* The trigamma function is the derivative of the digamma function.

   Reference:

    B Schneider,
    Trigamma Function,
    Algorithm AS 121,
    Applied Statistics,
    Volume 27, Number 1, page 97-99, 1978.

    From http://www.psc.edu/~burkardt/src/dirichlet/dirichlet.f
    (with modification for negative arguments and extra precision)
*/


double trigamma(double x)
{
  double neginf = -INFINITY,
    small = 1e-4,
    large = 8,
    trigamma1 = 1.6449340668482264365, /* pi^2/6 = Zeta(2) */
    tetragamma1 = -2.404113806319188570799476,  /* -2 Zeta(3) */
    b2 =  1./6,  /* B_2 */
    b4 = -1./30, /* B_4 */
    b6 =  1./42, /* B_6 */
    b8 = -1./30, /* B_8 */
    b10 = 5./66; /* B_10 */
  double result;
  /* Illegal arguments */
  if((x == neginf) || isnan(x)) {
    return NAN;
  }
  /* Singularities */
  if((x <= 0) && (floor(x) == x)) {
    return neginf;
  }
  /* Negative values */
  /* Use the derivative of the digamma reflection formula:
   * -trigamma(-x) = trigamma(x+1) - (pi*csc(pi*x))^2
   */
  if(x < 0) {
    result = M_PI/sin(-M_PI*x);
    return -trigamma(1-x) + result*result;
  }
  /* Use Taylor series if argument <= small */
  if(x <= small) {
    return 1/(x*x) + trigamma1 + tetragamma1*x;
  }
  result = 0;
  /* Reduce to trigamma(x+n) where ( X + N ) >= B */
  while(x < large) {
    result += 1/(x*x);
    x++;
  }
  /* Apply asymptotic formula when X >= B */
  /* This expansion can be computed in Maple via asympt(Psi(1,x),x) */
  if(x >= large) {
    double r = 1/(x*x);
    result += 0.5*r + (1 + r*(b2 + r*(b4 + r*(b6 + r*(b8 + r*b10)))))/x;
  }
  return result;
}
