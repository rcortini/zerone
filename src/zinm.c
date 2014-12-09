#include "zinm.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


int
is_invalid
(
    const int *y,
    int k,
    int r
)
// SYNOPSIS:                                                              
//   Helper function for `zinm_prob`. NAs of type 'int' is the largest
//   negative value. More generally, any negative value in 'y' is
//   invalid.                                     
{
   for (int i = 0 ; i < r ; i++) if (y[i + k*r] < 0) return 1;
   return 0;
}

int
is_all_zero
(
   const int *y,
   int k,
   int r
)
// SYNOPSIS:                                                              
//   Helper function for `zinm_prob`. Returns 1 if and only if all
//   the observations are 0.
{
   for (int i = 0 ; i < r ; i++) if (y[i + k*r] != 0) return 0;
   return 1;
}

void
zinm_prob
(
   // input //
   const int *n_states,
   const int *n_obs,
   const int *dim_y,
   const int * restrict y,
   // params //
   const double * restrict a,
   const double * restrict p,
   const double * restrict pi,
   // index //
         int * restrict index,
   // control //
   const int * restrict output,
   // output //
   double * restrict pem
)
// SYNOPSIS:                                                             
//   Compute emission probabilities with a mixture negative multinomial  
//   model. Since those are up to a multiplicative constant in the       
//   forward-backward algorithm, we can drop the multiplicative terms    
//   that do not depend on the state of the HMM.                         
//   Since the negative nultinomial takes discrete values, we can cache  
//   the results for reuse in order to save computation. This is done    
//   by indexing the series.                                             
//                                                                       
//   My parametrization is of the form:                                  
//                                                                       
//        p_0(i)^a * p_1(i)^y_1 * p_2(i)^y_2 * ... * p_r+1(i)^y_r        
//                                                                       
//   And in the case that all emissions are 0                            
//                                                                       
//                      pi * p_0(i)^a + (1-pi)                           
//                                                                       
// NUMERICAL STABILITY:                                                  
//   Each term of the sum above is computed in log space, the result is  
//   the computed as the sum of two exponentials. NA emissions are       
//   allowed and yield NA for the whole line of emissions.               
//                                                                       
// ARGUMENTS:                                                            
//   'n_states': (1) number of states in the HMM (alias 'm')             
//   'n_obs': (1) length of the sequence of observations (alias 'n')     
//   'dim_y': (1) number of columns of 'y' (alias 'r')                 
//   'y': (n_obs,dim_y) profiles                                       
//   'pi': (1) model parameter                                           
//   'a': (1) model parameter                                            
//   'p': (dim_y,m) model parameters                                    
//   'output': the type of output to produce (see below)                 
//   'pem': (n_obs,n_states) emission probability                        
//                                                                       
// RETURN:                                                               
//   'void'                                                              
//                                                                       
// SIDE EFFECTS:                                                         
//   Update 'pem' in place.                                              
//                                                                       
// OUTPUT:                                                               
//   The output type for 'pem' can be the emission probability in        
//   log space (1), the same emission probability in linear space (2),   
//   or in linear by default and in log space in case of underflow (0).  
//   'output' also controls the verbosity. If the third bit is set,      
//   i.e. the value is set to 4, 5 or 6, the function will suppress      
//   warnings. Setting the fourth bit of 'output' forces to compute      
//   the constant terms in emission probabilities.                       
{

   char *depends   = "compute in lin space, log space if underflow";
   char *log_space = "always compute in log space";
   char *lin_space = "always compute in linear space";
   char *cases[3] = {depends, log_space, lin_space};
   char *output_type = cases[*output & 3];

   int compute_constant_terms = (*output >> 3) & 1;

   int n = *n_obs;
   int m = *n_states;
   int r = *dim_y;

   if (*index < 0) indexts(n, r, y, index);

   double *logp = malloc((r+1)*m * sizeof(double));
   if (logp == NULL) {
      fprintf(stderr, "memory error (%s:%d)\n", __FILE__, __LINE__);
      return;
   }

   // If the third bit of 'output' is set, suppress warnings
   // by setting 'warned' to 1.
   int warned = (*output >> 2) & 1;

   // Make sure that 'p' defines a probability.
   for (int i = 0 ; i < m ; i++) {
      double sump = 0.0;
      for (int j = 0 ; j < r+1 ; j++) {
         // Cannot normalize negative values. Sorry folks.
         if (p[j+i*(r+1)] < 0) {
            fprintf(stderr, "error: 'p' negative\n");
            return;
         }
         sump += p[j+i*(r+1)];
      }
      int p_normalized_no = fabs(sump - 1.0) > DBL_EPSILON;
      if (!warned && p_normalized_no) {
         fprintf(stderr, "warning: renormalizing 'p'\n");
         warned = 1;
      }
      for (int j = 0 ; j < r+1 ; j++) {
         logp[j+i*(r+1)] = log(p[j+i*(r+1)] / sump);
      }
   }

   // The following variable 'row_of_na' comes in handy to write
   // full lines of NAs in the emissions.
   double *row_of_na = malloc(m * sizeof(double));
   if (row_of_na == NULL) {
      fprintf(stderr, "memory error (%s:%d)\n", __FILE__, __LINE__);
      return;
   }
   for (int i = 0 ; i < m ; i++) row_of_na[i] = NAN;

   for (int k = 0 ; k < n ; k++) {
      // Indexing allows to compute the terms only once. If the term
      // has been computed before, copy the value and move on.
      if (index[k] < k) {
         memcpy(pem + k*m, pem + index[k]*m, m * sizeof(double));
         continue;
      }

      // This is the firt occurrence of the emission in the times
      // series. We need to compute the emission probability.
      // Test the presence of invalid/NA emissions in the row.
      // If so, fill the row with NAs and move on.
      if (is_invalid(y, k, r)) {
         memcpy(pem + k*m, row_of_na, m * sizeof(double));
         continue;
      }

      if (is_all_zero(y, k, r)) {
         // Emissions are all zeros, use the zero-inflated
         // term from the zinm model.
         for (int i = 0 ; i < m ; i++) {
            pem[i+k*m] = log(*pi*exp(*a*logp[0+i*(r+1)]) + (1.0-*pi));
         }
      }
      else {
         // Otherwise use the standard probability.
         for (int i = 0 ; i < m ; i++) {
            pem[i+k*m] = (*a) * logp[0+i*(r+1)];
            for (int j = 0 ; j < r ; j++) {
               pem[i+k*m] += y[j+k*r] * logp[(j+1)+i*(r+1)];
            }
         }
      }

      if (compute_constant_terms) {
         double c_term = -lgamma(*a);
         double sum = *a;
         for (int j = 0 ; j < r ; j++) {
            int term = y[j+k*r];
            sum += term;
            c_term -= lgamma(term+1);
         }
         c_term += lgamma(sum);
         for (int i = 0 ; i < m ; i++) {
            pem[i+k*m] += c_term;
         }
      }

      if (output_type == log_space) continue;

      double sum = 0.0;
      double lin[m];
      for (int i = 0 ; i < m ; i++) sum += lin[i] = exp(pem[i+k*m]);
      if (sum > 0 || output_type == lin_space) {
         memcpy(pem+k*m, lin, m * sizeof(double));
      }

   }

   free(logp);
   free(row_of_na);
   return;

}


zinm_par_t *
new_zinm_par
(
   size_t r
)
{

   zinm_par_t *new = calloc(1, sizeof(zinm_par_t) + (r+1)*sizeof(double));
   if (new == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }
   new->r = r;
   new->pi = 1.0;

   return new;

}


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
eval_zinm_f
(
   double a,
   double p,
   unsigned int nz,
   double sum
)
{
   return a*nz / (p*(1-pow(p,a))) - sum/(1-p);
}


double
eval_zinm_g
(
         double a,
         double p,
   const tab_t *tab
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
eval_zinm_dfda
(
   double a,
   double p,
   unsigned int nz
)
{
   const double ppa = pow(p,a);
   return nz*(1-ppa+a*ppa*log(p)) / (p*sq(1-ppa));
}

double
eval_zinm_dfdp
(
   double a,
   double p,
   unsigned int nz,
   double m
)
{
   const double ppa = pow(p,a);
   return -(nz*a*(1-(a+1)*ppa) / sq(p*(1-ppa)) + m/sq(1-p));
}


double
eval_zinm_dgda
(
         double a,
         double p,
   const tab_t *tab
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


double
ll_zinm
(
         double a,
         double p,
         double pi,
   const tab_t *tab
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


zinm_par_t *
mle_zinm
(
   size_t *x,
   size_t dim,
   size_t nobs
)
{

   tab_t *tab = tabulate(x, dim, nobs);
   double *means = malloc(dim * sizeof(double));
   if (means == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   // Compute the means in all dimensions and the grand mean.
   compute_means(x, dim, nobs, means);
   double mean = 0.0;
   for (size_t i = 0 ; i < dim ; i++) mean += means[i];
   double sum = mean * nobs;

   // Extract the number of all-zero observaions.
   const unsigned int z0 = tab->val[0] == 0 ? tab->num[0] : 0;

   double deficit[11] = {0,.1,.2,.3,.4,.5,.6,.7,.8,.9,1};
   double init_a[12] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,1};
   double init_p[12] = {0,0,0,0,0,0,0,0,0,0,0,.5};

   // Deplete some 0s from the observations, compute alpha
   // and p0 with standard negative binomial estimation
   // and keep the values to be used as initial conditions.
   for (size_t i = 0 ; i < 11 ; i++) {
      double newmean = mean;
      if (tab->val[0] == 0) {
         tab->num[0] = z0 * (1-deficit[i]);
         newmean /= (1.0 - z0*deficit[i]/nobs);
      }
      double alpha = nb_est_alpha(tab);
      init_a[i] = alpha;
      init_p[i] = alpha / (alpha + newmean);
   }

   // Reset 'tab'.
   if (tab->val[0] == 0) tab->num[0] = z0;

   // Try initial conditions. Number 12 is a safety in case
   // all the rest failed during the first phase.
   double max_loglik = -1.0/0.0; // -inf
   zinm_par_t *par = new_zinm_par(dim);
   if (par == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   for (size_t i = 0 ; i < 12 ; i++) {

      if (init_a[i] < 0) continue;  // Skip failures.

      double a = init_a[i];
      double p = init_p[i];

      double grad;
      unsigned int iter = 0;

      double f = eval_zinm_f(a, p, nobs-z0, sum);
      double g = eval_zinm_g(a, p, tab);

      // Newton-Raphson iterations.
      while ((grad = f*f+g*g) > sq(ZINM_TOL) && iter++ < ZINM_MAXITER) {

         double dfda, dfdp, dgda, dgdp;
         dfda = dgdp = eval_zinm_dfda(a, p, nobs-z0);
         dfdp = eval_zinm_dfdp(a, p, nobs-z0, sum);
         dgda = eval_zinm_dgda(a, p, tab);

         double denom = dfdp*dgda - dfda*dgdp;
         double da = (f*dgdp - g*dfdp) / denom;
         double dp = (g*dfda - f*dgda) / denom;
         f = eval_zinm_f(a+da, p+dp, nobs-z0, sum);
         g = eval_zinm_g(a+da, p+dp, tab);

         for (int j = 0 ; j < ZINM_MAXITER && f*f+g*g > grad ; j++) {
            da /= 2;
            dp /= 2;
            f = eval_zinm_f(a+da, p+dp, nobs-z0, sum);
            g = eval_zinm_g(a+da, p+dp, tab);
         }

         a = a+da;
         p = p+dp;

      }

      double pi = (nobs-z0) / (1-pow(p,a)) / nobs;
      if (pi > 1) pi = 1.0;
      if (pi < 0) pi = 0.0;
      double loglik = ll_zinm(a, p, pi, tab);
      if (loglik > max_loglik) {
         max_loglik = loglik;
         par->alpha = a;
         par->pi = pi;
         par->p[0] = p;
      }
            
   }

   free(means);
   return par;

}


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
      while (eval_nb_f(a, tab) < 0) a /= 2;
      a_lo = a;
      a_hi = a*2;
   }
   else {
      a *= 2;
      while (eval_nb_f(a, tab) > 0) a *= 2;
      a_lo = a/2;
      a_hi = a;
   }

   // Input is pathological.
   if (a_lo > 128) return -1.0;

   double new_a = (a_lo + a_hi) / 2;
   for (int i = 0 ; i < ZINM_MAXITER ; i++) {
      a = (new_a < a_lo || new_a > a_hi) ?
         (a_lo + a_hi) / 2 :
         new_a;
      double f = eval_nb_f(a, tab);
      if (f < 0) a_hi = a; else a_lo = a;
      if ((a_hi - a_lo) < ZINM_TOL) break;
      double dfda = eval_nb_dfda(a, tab);
      new_a = a - f / dfda;
   }

   return a;

}

zinm_par_t *
mle_nm
(
   size_t *x,
   size_t dim,
   size_t nobs
)
{

   // Estimate alpha.
   tab_t *tab = tabulate(x, dim, nobs);
   double alpha = nb_est_alpha(tab);
   free(tab);

   // Return NULL if failed.
   if (alpha < 0) return NULL;

   // Compute the means in all dimensions.
   double *means = malloc(dim * sizeof(double));
   if (means == NULL) {
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   compute_means(x, dim, nobs, means);

   zinm_par_t *par = new_zinm_par(dim);
   if (par == NULL) {
      free(means);
      fprintf(stderr, "memory error: %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   double mean = 0;
   for (size_t i = 0 ; i < dim ; i++) mean += means[i];
   par->alpha = alpha;
   par->p[0] = alpha / (alpha + mean);
   for (size_t i = 1 ; i < dim+1 ; i++) {
      par->p[i] = par->p[0] / alpha * means[i-1];
   }

   free(means);

   return par;

}
