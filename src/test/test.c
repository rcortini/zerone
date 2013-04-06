#define _GNU_SOURCE
#include <time.h>
#include <glib.h>
#include <stdio.h>
#include <string.h>
#include "estep.h"
#include "loglik.h"

typedef struct {
   double *Q;
   // double *init;
   double *alpha;
   double *beta;
   double *gamma;
   int *n;
   int *r;
   int *yz;
   int *index;
   int *tabulated;
} loglik_fixture;


// Misc functions.

int *
readinput
(void)
{
   //FILE *f = fopen("H3K9ac.txt", "r");
   FILE *f = fopen("H3K9ac.txt", "r");
   char *line = NULL;
   size_t len = 0;
   ssize_t read;

   int n = 1031884;
   int i = 0;

   char seqname[20];
   int *yz = malloc (3*n * sizeof(int));
   memset(yz, (int) -1, 3*n * sizeof(int));

   // Discard header.
   getline(&line, &len, f);
   while ((read = getline(&line, &len, f)) != -1) {
      sscanf(line, "%s\t%d\t%d\t%d\n", seqname, 
         yz+i, yz+i+1, yz+i+2);
      i += 3;
   }
   free(line);
   fclose(f);
   return yz;
}

// Setup functions.

void
setup_1
(
   loglik_fixture *f,
   gconstpointer test_data
)
{
   static double Q[9] = {
      1.0/3, 1.0/3, 1.0/3,
      1.0/3, 1.0/3, 1.0/3,
      1.0/3, 1.0/3, 1.0/3
   };

   // static double init[3] = {1.0/3, 1.0/3, 1.0/3};
   static double alpha = 1.0;
   static double beta = 1.0;
   static double gamma[3] = {1.0, 1.0, 1.0};
   static int n = 3;
   static int r = 1;
   static int yz[6] = {1,1,1,1,1,1};
   static int index[3] = {-1,-1,-1};
   static int tabulated[4] = {0,0,3,-1};

   // f->init = init;
   f->Q = Q;
   f->alpha = &alpha;
   f->beta = &beta;
   f->gamma = gamma;
   f->n = &n;
   f->r = &r;
   f->yz = yz;
   f->index = index;
   f->tabulated = tabulated;
}

void
setup_underflow
(
   loglik_fixture *f,
   gconstpointer test_data
)
{
   static double Q[9] = {
      1.0, 1.0, 1.0,
      0, 0, 0,
      0, 0, 0
   };

   //static double init[3] = {1.0, 0, 0};
   static double alpha = 1.0;
   static double beta = 1.0;
   static double gamma[3] = {1.0, 2.0, 1.0};
   static int n = 2;
   static int r = 1;
   static int yz[4] = {1,999,1,999};
   static int index[2] = {-1,-1};
   int *tabulated = calloc(1002, sizeof(int));
   tabulated[1000] = 2;
   tabulated[1001] = -1;

   //f->init = init;
   f->Q = Q;
   f->alpha = &alpha;
   f->beta = &beta;
   f->gamma = gamma;
   f->n = &n;
   f->r = &r;
   f->yz = yz;
   f->index = index;
   f->tabulated = tabulated;
}

void
setup_NA
(
   loglik_fixture *f,
   gconstpointer test_data
)
{
   static double Q[9] = {
      1.0, 1.0, 1.0,
      0, 0, 0,
      0, 0, 0
   };

   //static double init[3] = {1.0, 0, 0};
   static double alpha = 1.0;
   static double beta = 1.0;
   static double gamma[3] = {1.0, 1.0, 1.0};
   static int n = 2;
   static int r = 1;
   static int yz[4] = {-1,1,1,1};
   static int index[2] = {-1,-1};
   static int tabulated[4] = {0,0,1,-1};

   //f->init = init;
   f->Q = Q;
   f->alpha = &alpha;
   f->beta = &beta;
   f->gamma = gamma;
   f->n = &n;
   f->r = &r;
   f->yz = yz;
   f->index = index;
   f->tabulated = tabulated;
}

// Test functions.

void
test_indexts
(void)
{
   int index[8];

   int ts1[24] = {
      0, 0, 0,
      0, 0, 1,
      0, 1, 0,
      1, 0, 1,
      0, 0, 1,
      0, 0, 1,
      0, 1, 0,
      0, 0, 0,
   };

   int ts2[24] = {
      -100, 0, 1,
      -1,   0, 100,
      -1,   0, 1,
      -1,   0, 100,
      -100, 0, 1,
      -1,   0, 1,
      -100, 0, 100,
      -100, 0, 100,
   };

   int expected1[8] = {0,1,2,3,1,1,2,0};
   int expected2[8] = {0,1,2,3,0,2,6,6};
   
   indexts(8, 3, ts1, index);
   for (int i = 0 ; i < 8 ; i++) g_assert(index[i] == expected1[i]);

   indexts(8, 3, ts2, index);
   for (int i = 0 ; i < 8 ; i++) g_assert(index[i] == expected2[i]);

}

// Teardown functions.

void
teardown_underflow (
   loglik_fixture *f,
   gconstpointer test_data
)
{
   free(f->tabulated);
}


void
test_1(
   loglik_fixture *f,
   gconstpointer ignore
)
{
   double loglik;
   compute_loglik(
      (const int *) f->n,
      (const int *) f->r,
      (const int *) f->yz,
      (const double *) f->Q,
      (const double *) f->alpha,
      (const double *) f->beta,
      (const double *) f->gamma,
      f->index,
      (const int *) f->tabulated,
      &loglik
   );
   g_assert(loglik == 3*(lgamma(3)+3*log(1.0/3)));
}

void
test_underflow(
   loglik_fixture *f,
   gconstpointer ignore
)
{
   double loglik;
   compute_loglik(
      (const int *) f->n,
      (const int *) f->r,
      (const int *) f->yz,
      (const double *) f->Q,
      (const double *) f->alpha,
      (const double *) f->beta,
      (const double *) f->gamma,
      f->index,
      (const int *) f->tabulated,
      &loglik
   );
   g_assert(loglik == 
         log(1.0/3)+999*log(2.0/4)+2*log(1.0/4)+lgamma(1001));
}

void
test_NA(
   loglik_fixture *f,
   gconstpointer ignore
)
{
   double loglik;
   compute_loglik(
      (const int *) f->n,
      (const int *) f->r,
      (const int *) f->yz,
      (const double *) f->Q,
      (const double *) f->alpha,
      (const double *) f->beta,
      (const double *) f->gamma,
      f->index,
      (const int *) f->tabulated,
      &loglik
   );
   g_assert(loglik == lgamma(3)+3*log(1.0/3));
}

void
test_simAnneal
(void)
{
   int *yz = readinput();
   int n = 1031884;
   int r = 2;
   double *globalmax = malloc(18 * sizeof(double));
   double time_taken[10];
   double perf[10];
   for (int i = 0 ; i < 10 ; i++) {
      memset(globalmax, 0.0, 18 * sizeof(double));
      time_t start = time(NULL);
      simAnneal(&n, &r, yz, globalmax);
      time_t end = time(NULL);
      time_taken[i] = difftime(end, start);
      perf[i] = *globalmax;
   }
   free(yz);
   free(globalmax);
   for (int i = 0 ; i < 10 ; i++) {
      fprintf(stderr, "[%.2f] %.3f\n", time_taken[i], perf[i]);
   }
}

int
main(
   int argc,
   char **argv
)
{
   g_test_init(&argc, &argv, NULL);
   g_test_add_func("/indexts/test", test_indexts);
   //g_test_add("/loglik/test_1", loglik_fixture, NULL, setup_1,
   //      test_1, NULL);
   //g_test_add("/loglik/test_underflow", loglik_fixture, NULL,
   //      setup_underflow, test_underflow, teardown_underflow);
   //g_test_add("/loglik/test_NA", loglik_fixture, NULL, setup_NA,
   //      test_NA, NULL);
   g_test_add_func("/simAnneal/run", test_simAnneal);
   return g_test_run();
}
