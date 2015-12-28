#include "unittest.h"
#include "utils.h"
#include "zerone.c"

void
test_zinm_prob
(void)
{

   const uint m = 3;
   const uint n = 7;
   const uint r = 3;

   const double pi = 0.8;
   const double a = 1.2;
   // The constant C1 facilitates the definition of 'p'.
   const double C1 = 2.5;
   const double p[12] = {
      C1, 1.0, 1.0, 1.0, 
      C1, 1.0, 2.0, 1.0,
      C1, 1.0, 0.2, 0.1,
   };

   int y[21] = {
         1, 2, 2,
         0, 4, 2,
         1, 2, 2,
         1, 2, 0,
         0, 0, 0,
      // Underflow
      1500, 1, 2,
      // Negative values
        -1, 4, 2,
   };

   ChIP_t *ChIP = new_ChIP(m, 1, y, NULL, &n);

   double dummy[9] = {0};
   zerone_t *zerone = new_zerone(r, ChIP);
   set_zerone_par(zerone, dummy, a, pi, p);

   int index[7] = {-1,-1,-1,-1,-1,-1,-1};
   double pem[21];
   indexts(n, r, y, index);

   //--               Test output type 0               --//
   int out = 0;
   redirect_stderr();
   zinm_prob(zerone, index, out, pem);
   unredirect_stderr();

   double expected_pem_1[21] = {
      // Checked manually.
          7.713996e-05,      0.0001095280,      3.054426e-07,
          1.402544e-05,      6.740185e-05,      3.215185e-09,
          7.713996e-05,      0.0001095280,      3.054426e-07,
          0.0023334833,      0.0046275583,      0.0004410592,
          0.5105866396,      0.4541686425,      0.6840360201,
      -2563.1825314667,  -2813.7721384365,  -2013.2236637989,
          NAN,               NAN,               NAN,
   };

   for (int i = 0 ; i < 18 ; i++) {
      test_assert(fabs(expected_pem_1[i] - pem[i]) < 1e-9);
   }
   for (int i = 18 ; i < 21 ; i++) {
      test_assert(pem[i] != pem[i]);
   }

   // Test the warning message.
   test_assert_stderr("warning: renormalizing 'p'\n");

   //--               Test output type 1               --//
   out = 1;
   redirect_stderr();
   zinm_prob(zerone, index, out, pem);
   unredirect_stderr();

   double expected_pem_2[21] = {
      // Checked manually.
      log(7.713996e-05), log(0.0001095280), log(3.054426e-07),
      log(1.402544e-05), log(6.740185e-05), log(3.215185e-09),
      log(7.713996e-05), log(0.0001095280), log(3.054426e-07),
      log(0.0023334833), log(0.0046275583), log(0.0004410592),
      log(0.5105866396), log(0.4541686425), log(0.6840360201),
      -2563.1825314667,  -2813.7721384365,  -2013.2236637989,
      NAN,          NAN,          NAN,
   };

   for (int i = 0 ; i < 18 ; i++) {
      test_assert(fabs(expected_pem_2[i] - pem[i]) < 1e-6);
   }
   for (int i = 18 ; i < 21 ; i++) {
      test_assert(pem[i] != pem[i]);
   }
   // Test the warning message.
   test_assert_stderr("warning: renormalizing 'p'\n");

   //--               Test output type 2               --//
   out = 2;
   redirect_stderr();
   zinm_prob(zerone, index, out, pem);
   unredirect_stderr();

   double expected_pem_3[21] = {
      // Checked manually.
      7.713996e-05,      0.0001095280,      3.054426e-07,
      1.402544e-05,      6.740185e-05,      3.215185e-09,
      7.713996e-05,      0.0001095280,      3.054426e-07,
      0.0023334833,      0.0046275583,      0.0004410592,
      0.5105866396,      0.4541686425,      0.6840360201,
      0.0000000000,      0.0000000000,      0.0000000000,
      NAN,               NAN,               NAN,
   };

   for (int i = 0 ; i < 15 ; i++) {
      test_assert(fabs(expected_pem_3[i] - pem[i]) < 1e-8);
   }
   for (int i = 15 ; i < 18 ; i++) {
      test_assert(pem[i] == 0.0);
   }
   for (int i = 18 ; i < 21 ; i++) {
      test_assert(pem[i] != pem[i]);
   }
   // Test the warning message.
   test_assert_stderr("warning: renormalizing 'p'\n");

   free(ChIP);
   zerone->ChIP = NULL;
   destroy_zerone_all(zerone);

   return;

}


void
test_bw_zinm
(void)
{
   int y[90] = {
      2,2,2,
      5,0,2,
      3,3,2,
      2,1,1,
      8,0,2,
      2,2,0,
      0,1,0,
      1,2,0,
      5,2,1,
      3,0,0,
      2,1,0,
      9,1,1,
      2,2,0,
      2,0,2,
      2,0,1,
      2,2,3,
      4,2,3,
      6,12,9,
      2,1,2,
      4,3,10,
      6,8,8,
      1,7,5,
      5,8,3,
      5,5,9,
      4,5,3,
      4,2,0,
      6,12,7,
      2,11,6,
      4,5,6,
      4,3,5,
   };

   unsigned size[2] = {15,15};
   ChIP_t *ChIP = new_ChIP(3, 2, y, NULL, size);

   double p[8] = {
      .3448276, .4137931, .1379310, .1034483,
      .1612903, .1935484, .3225806, .3225806,
   };
   double Q[4] = { .8, .2, .2, .8 };
   zerone_t *zerone = new_zerone(2, ChIP);
   set_zerone_par(zerone, Q, 3.4, 1.0, p);
   bw_zinm(zerone);

   double expected_Q[4] = {
      // transpose //
      1.000000, 0.000001,
      0.000000, 0.999999, 
   };
   for (size_t i = 0 ; i < 4 ; i++) {
      test_assert(fabs(zerone->Q[i] - expected_Q[i]) < 1e-6);
   }

   free(ChIP);
   zerone->ChIP = NULL;
   destroy_zerone_all(zerone);

   return;

}



void
test_update_trans
(void)
{
   double Q[9] = {0};
   double trans[9] = {
      // transpose //
      14, 56, 44,
      14, 57, 43,
      14, 58, 42,
   };

   update_trans(3, Q, trans);

   double expected_Q[9] = {
      // transpose //
      .3333, 0.3274, 0.3410,
      .3333, 0.3333, 0.3333,
      .3333, 0.3391, 0.3255,
   };

   for (size_t i = 0 ; i < 9 ; i++) {
      test_assert(fabs(Q[i] - expected_Q[i]) < 1e-3);
   }

   return;

}

#if 0
void
test_read_file
(void)
{

   FILE *inputf = fopen("sample_file.txt", "r");
   test_assert_critical(inputf != NULL);
   ChIP_t *ChIP = read_file(inputf);
   fclose(inputf);

   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 3);

   test_assert(ChIP->nb == 1);
   test_assert(ChIP->sz[0] == 5);

   int *y = ChIP->y;
   int expected_y[15] = {
      -1,-1,-1,
      -1,-1,-1,
      -1,-1,-1,
       0, 0, 0,
       2, 1, 0,
   };
   for (size_t i = 0 ; i < 15 ; i++) {
      test_assert(y[i] == expected_y[i]);
   }

   free(ChIP->y);
   free(ChIP);

   return;

}
#endif


void
test_eval_bw_f
(void)
//   double a,
//   double pi,
//   double p0,
//   double A,
//   double B,
//   double C,
//   double D,
//   double E
{
   return;
//   double term1 = (D + a*A) / p0;
//   double term2 = B * pi*a*pow(p0,a-1) / (pi*pow(p0,a)+1-pi);
//   return p0 + E/(term1 + term2) - 1.0 / C;
}
