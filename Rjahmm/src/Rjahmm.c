#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "jahmm.h"


// Register the method 'jahmm_R_call()'.
SEXP jahmm_R_call (SEXP, SEXP);

R_CallMethodDef callMethods[] = {
   {"jahmm_R_call", (DL_FUNC) &jahmm_R_call, 2},
   {NULL, NULL, 0},
};

void R_init_jahmm(DllInfo *info) {
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

SEXP
jahmm_R_call
(
   SEXP nstates,
   SEXP Y
)
{

   const unsigned int r = length(Y) - 1;
   const unsigned int n = length(VECTOR_ELT(Y, 0));

   // Get the values.
   int *y = malloc(n*r * sizeof(int));
   if (y == NULL) {
      Rprintf("Rjahmm memory error %s:%d\n", __FILE__, __LINE__);
      return R_NilValue;
   }

   for (size_t i = 1 ; i < r+1 ; i++) {
      int *v = INTEGER(coerceVector(VECTOR_ELT(Y, i), INTSXP));
      for (size_t j = 0 ; j < n ; j++) y[i-1+j*r] = v[j];
   }

   // Get the blocks.
   histo_t *histo = new_histo();
   if (histo == NULL) {
      Rprintf("Rjahmm memory error %s:%d\n", __FILE__, __LINE__);
      return R_NilValue;
   }

   int *id;
   char prev[256] = {0};

   SEXP S = VECTOR_ELT(Y, 0);
   switch (TYPEOF(S)) {
      case STRSXP:
         *id = -1;
         for (size_t i = 0 ; i < n ; i++) {
            const char *block = CHAR(STRING_ELT(S, i));
            if (strcmp(prev, block) != 0) {
               strncpy(prev, block, 256);
               (*id)++;
            }
            histo_push(&histo, *id);
         }
         break;
      case INTSXP:
         id = INTEGER(S);
         for (size_t i = 0 ; i < n ; i++) {
            histo_push(&histo, id[i]);
         }
   }

   tab_t *tab = compress_histo(histo);
   if (tab == NULL) {
      Rprintf("Rjahmm memory error %s:%d\n", __FILE__, __LINE__);
      return R_NilValue;
   }
   ChIP_t *ChIP = new_ChIP(r, tab->size, y, tab->num);
   free(histo);
   free(tab);

   unsigned int m = INTEGER(coerceVector(nstates, INTSXP))[0];
   jahmm_t * jahmm = do_jahmm(m, ChIP);

   if (jahmm == NULL) {
      Rprintf("Rjahmm error\n");
      return R_NilValue;
   }

   SEXP PHI;
   PROTECT(PHI = allocVector(REALSXP, m*n));
   double *phi = REAL(PHI);
   memcpy(phi, jahmm->phi, m*n * sizeof(double));

   SEXP PATH;
   PROTECT(PATH = allocVector(INTSXP, n));
   int *path = INTEGER(PATH);
   memcpy(path, jahmm->path, n * sizeof(int));

   SEXP dimPHI;
   PROTECT(dimPHI = allocVector(INTSXP, 2));
   INTEGER(dimPHI)[0] = n;
   INTEGER(dimPHI)[1] = m;
   setAttrib(PHI, R_DimSymbol, dimPHI);

   SEXP RETLIST;
   PROTECT(RETLIST = allocVector(VECSXP, 2));
   SET_VECTOR_ELT(RETLIST, 0, PHI);
   SET_VECTOR_ELT(RETLIST, 1, PATH);
   
   UNPROTECT(4);

   return RETLIST;


}
