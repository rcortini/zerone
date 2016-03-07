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

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "../../src/zerone.h"


// Register the method 'zerone_R_call()'.
SEXP zerone_R_call (SEXP, SEXP, SEXP);

R_CallMethodDef callMethods[] = {
   {"zerone_R_call", (DL_FUNC) &zerone_R_call, 3},
   {NULL, NULL, 0},
};

void R_init_zerone(DllInfo *info) {
   R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

SEXP
zerone_R_call
(
   SEXP RNAMES,
   SEXP RSIZE,
   SEXP RY
)
{

   const unsigned int nb = length(RNAMES);
   const unsigned int r  = length(RY) - 1;
   const unsigned int n  = length(VECTOR_ELT(RY, 0));

   // Get/copy the values from R objects.
   const char **name = malloc(nb * sizeof(char *));
   if (name == NULL) {
      Rprintf("Rzerone memory error %s:%d\n", __FILE__, __LINE__);
      return R_NilValue;
   }

   for (int i = 0 ; i < nb ; i++) {
      name[i] = CHAR(STRING_ELT(RNAMES, i));
   }

   int *size = INTEGER(RSIZE);

   int *y = malloc(n*r * sizeof(int));
   if (y == NULL) {
      Rprintf("Rzerone memory error %s:%d\n", __FILE__, __LINE__);
      return R_NilValue;
   }

   for (size_t i = 1 ; i < r+1 ; i++) {
      int *v = INTEGER(coerceVector(VECTOR_ELT(RY, i), INTSXP));
      for (size_t j = 0 ; j < n ; j++) y[i-1+j*r] = v[j];
   }

   ChIP_t *ChIP = new_ChIP(r, nb, y, name, size);
   zerone_t * zerone = do_zerone(ChIP);

   if (zerone == NULL) {
      Rprintf("Rzerone error\n");
      return R_NilValue;
   }

   // Zerone uses 3 states.
   const unsigned int m = 3;

   SEXP Q;
   PROTECT(Q = allocVector(REALSXP, m*m));
   memcpy(REAL(Q), zerone->Q, m*m * sizeof(double));

   SEXP A;
   PROTECT(A = allocVector(REALSXP, 1));
   *REAL(A) =  zerone->a;

   SEXP PI_;
   PROTECT(PI_ = allocVector(REALSXP, 1));
   *REAL(PI_) =  zerone->pi;

   // The parameters 'p' are coded "row-wise". Since R
   // objects are coded "column-wise" we neeed to
   // disentangle them.
   SEXP P;
   PROTECT(P = allocVector(REALSXP, m*(r+1)));
   for (size_t i = 0 ; i < r+1 ; i++) {
   for (size_t j = 0 ; j < m ; j++) {
      REAL(P)[j+i*m] = zerone->p[i+j*(r+1)];
   }
   }

   // For reason of cache-friendlyness, 'phi' and 'pem'
   // are also coded "row-wise".
   SEXP PHI;
   SEXP PEM;
   PROTECT(PHI = allocVector(REALSXP, m*n));
   PROTECT(PEM = allocVector(REALSXP, m*n));
   for (size_t i = 0 ; i < n ; i++) {
   for (size_t j = 0 ; j < m ; j++) {
      REAL(PEM)[i+j*n] = zerone->pem[j+i*m];
      REAL(PHI)[i+j*n] = zerone->phi[j+i*m];
   }
   }

   SEXP PATH;
   PROTECT(PATH = allocVector(INTSXP, n));
   memcpy(INTEGER(PATH), zerone->path, n * sizeof(int));

   SEXP L;
   PROTECT(L = allocVector(REALSXP, 1));
   *REAL(L) = zerone->l;

   SEXP dimQ;
   PROTECT(dimQ = allocVector(INTSXP, 2));
   INTEGER(dimQ)[0] = m;
   INTEGER(dimQ)[1] = m;
   setAttrib(Q, R_DimSymbol, dimQ);

   SEXP dimP;
   PROTECT(dimP = allocVector(INTSXP, 2));
   INTEGER(dimP)[0] = m;
   INTEGER(dimP)[1] = r+1;
   setAttrib(P, R_DimSymbol, dimP);

   SEXP dimPHI;
   PROTECT(dimPHI = allocVector(INTSXP, 2));
   INTEGER(dimPHI)[0] = n;
   INTEGER(dimPHI)[1] = m;
   setAttrib(PHI, R_DimSymbol, dimPHI);

   SEXP dimPEM;
   PROTECT(dimPEM = allocVector(INTSXP, 2));
   INTEGER(dimPEM)[0] = n;
   INTEGER(dimPEM)[1] = m;
   setAttrib(PEM, R_DimSymbol, dimPEM);

   SEXP RETLIST;
   PROTECT(RETLIST = allocVector(VECSXP, 8));
   SET_VECTOR_ELT(RETLIST, 0, Q);
   SET_VECTOR_ELT(RETLIST, 1, A);
   SET_VECTOR_ELT(RETLIST, 2, PI_);
   SET_VECTOR_ELT(RETLIST, 3, P);
   SET_VECTOR_ELT(RETLIST, 4, PHI);
   SET_VECTOR_ELT(RETLIST, 5, PEM);
   SET_VECTOR_ELT(RETLIST, 6, PATH);
   SET_VECTOR_ELT(RETLIST, 7, L);

   UNPROTECT(13);

   destroy_zerone_all(zerone);

   return RETLIST;

}
