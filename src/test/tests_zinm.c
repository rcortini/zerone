#include "unittest.h"
#include "zinm.c"

void
test_new_histo
(void)
{

   histo_t *histo;
   histo = new_histo();
   test_assert_critical(histo != NULL);
   for (size_t i = 0 ; i < HISTO_INIT_SIZE ; i++) {
      test_assert(histo->num[i] == 0);
   }
   free(histo);

   redirect_stderr();
   set_alloc_failure_rate_to(1.0);
   histo = new_histo();
   reset_alloc();
   unredirect_stderr();
   test_assert(histo == NULL);
   test_assert(strcmp(caught_in_stderr(), "") != 0);

   return;

}

void
test_histo_push
(void)
{

   histo_t *histo;
   histo = new_histo();
   test_assert_critical(histo != NULL);
   for (size_t i = 0 ; i < (1 << 16) ; i++) {
      test_assert(histo_push(&histo, i) == 0);
   }
   test_assert(histo->size == (1 << 16));
   for (size_t i = 0 ; i < (1 << 16) ; i++) {
      test_assert(histo->num[i] == 1);
   }
   free(histo);

   histo = new_histo();
   test_assert_critical(histo != NULL);
   for (size_t i = (1 << 16) ; i > 0 ; i--) {
      test_assert(histo_push(&histo, i-1) == 0);
   }
   test_assert(histo->size == 2*((1 << 16)-1));
   for (size_t i = 0 ; i < (1 << 16) ; i++) {
      test_assert(histo->num[i] == 1);
   }
   free(histo);

   histo = new_histo();
   test_assert_critical(histo != NULL);
   redirect_stderr();
   set_alloc_failure_rate_to(1.0);
   test_assert(histo_push(&histo, HISTO_INIT_SIZE) == 1);
   reset_alloc();
   unredirect_stderr();
   test_assert(strcmp(caught_in_stderr(), "") != 0);
   free(histo);

   return;

}

void
test_compress_histo
(void)
{

   histo_t *histo = new_histo();
   test_assert_critical(histo != NULL);
   for (size_t i = 0 ; i < 4096 ; i += 2) {
      test_assert(histo_push(&histo, i) == 0);
   }

   tab_t *tab;
   tab = compress_histo(histo);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 2048);
   for (size_t i = 0 ; i < tab->size ; i++) {
      test_assert(tab->val[i] == 2*i);
      test_assert(tab->num[i] == 1);
   }

   free(tab);

   redirect_stderr();
   set_alloc_failure_rate_to(1.0);
   tab = compress_histo(histo);
   reset_alloc();
   unredirect_stderr();
   test_assert(tab == NULL);
   test_assert(strcmp(caught_in_stderr(), "") != 0);

   free(histo);
   return;

}

void
test_tabulate
(void)
{

   int x1[12] = {1,2,3,4,1,2,3,4,1,2,3,4};
   tab_t *tab;
   tab = tabulate(x1, 12);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 4);
   test_assert(tab->val[0] == 1);
   test_assert(tab->val[1] == 2);
   test_assert(tab->val[2] == 3);
   test_assert(tab->val[3] == 4);
   for (int i = 0 ; i < 4 ; i++) {
      test_assert(tab->num[i] == 3);
   }
   free(tab);

   int x2[4] = {4096,2048,1024,512};
   tab = tabulate(x2, 4);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 4);
   test_assert(tab->val[0] == 512);
   test_assert(tab->val[1] == 1024);
   test_assert(tab->val[2] == 2048);
   test_assert(tab->val[3] == 4096);
   for (int i = 0 ; i < 4 ; i++) {
      test_assert(tab->num[i] == 1);
   }
   free(tab);

   int x3[5] = {-1,2,3,-1,2};
   tab = tabulate(x3, 5);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 2);
   test_assert(tab->val[0] == 2);
   test_assert(tab->val[1] == 3);
   test_assert(tab->num[0] == 2);
   test_assert(tab->num[1] == 1);
   free(tab);

   return;

}


void
test_eval_nb_f
(void)
{

   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     1,1,1,1,1,2,2,2,2,3,5 };
   tab_t *tab = tabulate(x, 25);
   test_assert(fabs(eval_nb_f(1.0, tab)+0.12747262) < 1e-6);
   test_assert(fabs(eval_nb_f(1.1, tab)+0.24215981) < 1e-6);
   test_assert(fabs(eval_nb_f(1.2, tab)+0.31636395) < 1e-6);
   test_assert(fabs(eval_nb_f(1.3, tab)+0.36350700) < 1e-6);
   test_assert(fabs(eval_nb_f(1.4, tab)+0.39225466) < 1e-6);
   test_assert(fabs(eval_nb_f(1.5, tab)+0.40834322) < 1e-6);
   test_assert(fabs(eval_nb_f(2.0, tab)+0.39975512) < 1e-6);

   free(tab);

   return;

}


void
test_eval_nb_dfda
(void)
{

   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    1,1,1,1,1,2,2,2,2,3,5 };
   tab_t *tab = tabulate(x, 25);
   test_assert(fabs(eval_nb_dfda(1.0, tab)+1.41167874) < 1e-6);
   test_assert(fabs(eval_nb_dfda(1.2, tab)+0.58911102) < 1e-6);
   test_assert(fabs(eval_nb_dfda(1.3, tab)+0.36790287) < 1e-6);
   test_assert(fabs(eval_nb_dfda(1.4, tab)+0.21643981) < 1e-6);
   test_assert(fabs(eval_nb_dfda(1.5, tab)+0.11168877) < 1e-6);
   test_assert(fabs(eval_nb_dfda(2.0, tab)-0.08773865) < 1e-6);

   free(tab);

   return;

}


void
test_eval_zinb_f
(void)
{

   test_assert(fabs(eval_zinb_f(1.0, .5, 1, 2.0)) < 1e-6);
   test_assert(fabs(eval_zinb_f(1.0, .5, 2, 4.0)) < 1e-6);
   test_assert(fabs(eval_zinb_f(1.3, .7, 141, 8.3)-678.0838351) < 1e-6);

   return;

}


void
test_eval_zinb_g
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(eval_zinb_g(1, .5, tab)+0.1325713) < 1e-6);
   free(tab);

   return;

}


void
test_eval_zinb_dfda
(void)
{

   test_assert(fabs(eval_zinb_dfda(1, .5, 1)-1.2274112) < 1e-6);
   test_assert(fabs(eval_zinb_dfda(1, .5, 9)-11.0467015) < 1e-6);
   test_assert(fabs(eval_zinb_dfda(2, .5, 9)-12.9096451) < 1e-6);
   test_assert(fabs(eval_zinb_dfda(2, .3, 9)-25.1159846) < 1e-6);
   test_assert(fabs(eval_zinb_dfda(2.4, .3, 9)-26.3620864) < 1e-6);

   return;

}


void
test_eval_zinb_dfdp
(void)
{

   test_assert(eval_zinb_dfdp(1, .5, 0, 0) == 0.0);
   test_assert(eval_zinb_dfdp(1, .5, 1, 0) == 0.0);
   test_assert(fabs(eval_zinb_dfdp(2, .5, 1, 0)+3.5555555) < 1e-6);
   test_assert(fabs(eval_zinb_dfdp(2, .3, 1, 0)+19.5896899) < 1e-6);
   test_assert(fabs(eval_zinb_dfdp(2, .3, 9, 0)+176.3072092) < 1e-6);
   test_assert(fabs(eval_zinb_dfdp(2, .3, 9, 1.7)+179.7765970) < 1e-6);

   return;

}


void
test_eval_zinb_dgda
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(eval_zinb_dgda(1, .5, tab)+2.2547559) < 1e-6);
   test_assert(fabs(eval_zinb_dgda(2, .5, tab)+1.2605630) < 1e-6);
   test_assert(fabs(eval_zinb_dgda(2, .3, tab)+1.8764955) < 1e-6);
   free(tab);

   return;

}


void
test_ll_zinb
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(ll_zinb(1, .5, 1, tab)+22.5329303) < 1e-6);
   test_assert(fabs(ll_zinb(1, .5, .7, tab)+22.7832550) < 1e-6);
   test_assert(fabs(ll_zinb(2, .5, .7, tab)+23.7608409) < 1e-6);
   test_assert(fabs(ll_zinb(2, .3, .7, tab)+31.6978553) < 1e-6);
   free(tab);

   return;

}


void
test_nb_est_alpha
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-0.9237) < 1e-3);
   free(tab);

   // 0:27, 1:12, 2:8, 3:1, 4:1, 5:1
   int x2[50] = {3,0,1,2,0,0,1,0,0,0,0,1,1,0,0,1,2,2,0,0,0,1,2,
      0, 0,0,0,0,4,0,0,0,1,5,1,0,1,2,1,2,2,2,0,0,0,1,0,1,0,0};
   tab = tabulate(x2, 50);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-1.3436) < 1e-3);
   free(tab);

   // 0:12, 1:7, 2:13, 3:4, 4:6, 5:2, 6:1, 7:3, 8:1, 9:1
   int x3[50] = {4,5,2,1,2,4,2,2,0,4,2,1,3,6,0,0,7,3,0,8,4,2,0,
      0,2,3,2,3,7,9,2,4,0,4,2,0,0,2,5,1,1,2,1,0,0,0,1,2,1,7};
   tab = tabulate(x3, 50);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-1.7969) < 1e-3);
   free(tab);

   // 0:39, 1:8, 2:2, 3:1
   int x4[50] = {1,0,0,0,0,0,3,1,0,1,0,0,0,0,0,0,0,1,0,0,0,2,0,
      2,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0};
   tab = tabulate(x4, 50);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-0.7073) < 1e-3);
   free(tab);

   // 0:59, 1:83, 2:99, 3:67, 4:67, 5:49, 6:27, 7:22, 8:11, 9:6
   // 10:6, 11:3, 12:2, 13:3
   int x5[500] = {1,0,0,1,1,2,1,0,5,7,1,3,3,1,6,0,2,5,7,0,5,2,1,
       10,5,3,4,5,7,0,8,6,3,0,2,1,1,0,2,3,7,2,3,2,2,1,0,4,4,2,4,2,
       0,6,3,2,5,2,1,4,3,4,2,2,5,3,2,0,2,8,1,3,1,7,5,1,4,1,1,0,2,
       2,4,1,1,1,4,1,3,4,4,10,5,2,0,7,1,6,1,3,6,4,0,2,4,1,12,2,5,
       6,5,4,1,11,0,1,3,2,4,2,0,2,3,4,0,2,9,9,7,4,2,1,3,3,3,4,2,9,
       2,4,3,2,2,4,2,5,3,0,1,3,2,0,3,3,4,1,3,3,5,7,3,3,2,1,5,5,4,
       6,1,1,1,2,9,5,1,2,4,0,2,1,0,3,2,4,3,1,4,2,1,4,1,6,0,6,5,3,
       5,2,0,1,2,1,0,5,3,2,7,6,4,3,2,5,7,5,5,1,1,3,10,2,0,5,0,1,2,
       0,5,1,2,3,6,4,0,3,1,2,2,4,3,0,3,2,5,4,10,1,2,4,4,2,13,4,3,
       1,5,4,8,5,6,2,3,4,3,1,5,5,1,8,2,0,5,7,3,2,2,4,2,3,1,5,3,7,
       13,1,4,7,5,5,0,3,0,4,2,3,1,2,4,2,8,1,2,5,6,1,1,0,7,2,2,3,5,
       12,2,2,2,0,3,3,4,0,2,5,1,10,0,7,6,5,0,11,2,3,7,3,5,4,2,1,2,
       4,0,2,2,2,0,6,2,3,4,2,3,7,3,5,2,5,0,4,4,6,3,1,2,7,3,0,2,5,
       7,2,2,0,0,0,6,3,0,1,1,5,5,2,6,2,4,6,0,1,2,3,2,2,2,3,4,1,1,
       4,0,2,0,1,3,4,1,2,2,3,1,4,4,3,4,4,1,5,2,13,4,10,5,6,1,0,5,
       0,0,5,6,0,1,8,5,1,3,1,8,1,8,1,6,7,2,8,2,2,3,3,0,4,2,1,9,6,
       0,6,7,1,8,2,2,1,11,3,0,4,2,5,1,6,8,3,4,7,0,4,2,4,1,1,1,6,0,
       4,4,6,2,1,3,1,0,4,9,3,1,4,2,2,0,1};
   tab = tabulate(x5, 500);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-3.0057) < 1e-3);
   free(tab);
}


void
test_mle_nb
(void)
{

   // These test cases have been verified with R.
   zinb_par_t *par;
   
   // 0:14, 1:5, 2:4, 3:1, 5:1
   int x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   par = mle_nb(x1, 25);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a - 0.9237) < 1e-3);
   test_assert(par->pi == 1.0);
   test_assert(fabs(par->p - 0.5237) < 1e-3);
   free(par);

   // 0:27, 1:12, 2:8, 3:1, 4:1, 5:1
   int x2[50] = {3,0,1,2,0,0,1,0,0,0,0,1,1,0,0,1,2,2,0,0,0,1,2,
      0, 0,0,0,0,4,0,0,0,1,5,1,0,1,2,1,2,2,2,0,0,0,1,0,1,0,0};
   par = mle_nb(x2, 50);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a - 1.3436) < 1e-3);
   test_assert(par->pi == 1.0);
   test_assert(fabs(par->p - 0.6267) < 1e-3);
   free(par);

   // 0:12, 1:7, 2:13, 3:4, 4:6, 5:2, 6:1, 7:3, 8:1, 9:1
   int x3[50] = {4,5,2,1,2,4,2,2,0,4,2,1,3,6,0,0,7,3,0,8,4,2,0,
      0,2,3,2,3,7,9,2,4,0,4,2,0,0,2,5,1,1,2,1,0,0,0,1,2,1,7};
   par = mle_nb(x3, 50);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a - 1.7969) < 1e-3);
   test_assert(par->pi == 1.0);
   test_assert(fabs(par->p - 0.4221) < 1e-3);
   free(par);

   // 0:39, 1:8, 2:2, 3:1
   int x4[50] = {1,0,0,0,0,0,3,1,0,1,0,0,0,0,0,0,0,1,0,0,0,2,0,
      2,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0};
   par = mle_nb(x4, 50);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a - 0.7073) < 1e-3);
   test_assert(par->pi == 1.0);
   test_assert(fabs(par->p - 0.7021) < 1e-3);
   free(par);

   // 0:59, 1:83, 2:99, 3:67, 4:67, 5:49, 6:27, 7:22, 8:11, 9:6
   // 10:6, 11:3, 12:2, 13:3
   int x5[500] = {1,0,0,1,1,2,1,0,5,7,1,3,3,1,6,0,2,5,7,0,5,2,1,
       10,5,3,4,5,7,0,8,6,3,0,2,1,1,0,2,3,7,2,3,2,2,1,0,4,4,2,4,2,
       0,6,3,2,5,2,1,4,3,4,2,2,5,3,2,0,2,8,1,3,1,7,5,1,4,1,1,0,2,
       2,4,1,1,1,4,1,3,4,4,10,5,2,0,7,1,6,1,3,6,4,0,2,4,1,12,2,5,
       6,5,4,1,11,0,1,3,2,4,2,0,2,3,4,0,2,9,9,7,4,2,1,3,3,3,4,2,9,
       2,4,3,2,2,4,2,5,3,0,1,3,2,0,3,3,4,1,3,3,5,7,3,3,2,1,5,5,4,
       6,1,1,1,2,9,5,1,2,4,0,2,1,0,3,2,4,3,1,4,2,1,4,1,6,0,6,5,3,
       5,2,0,1,2,1,0,5,3,2,7,6,4,3,2,5,7,5,5,1,1,3,10,2,0,5,0,1,2,
       0,5,1,2,3,6,4,0,3,1,2,2,4,3,0,3,2,5,4,10,1,2,4,4,2,13,4,3,
       1,5,4,8,5,6,2,3,4,3,1,5,5,1,8,2,0,5,7,3,2,2,4,2,3,1,5,3,7,
       13,1,4,7,5,5,0,3,0,4,2,3,1,2,4,2,8,1,2,5,6,1,1,0,7,2,2,3,5,
       12,2,2,2,0,3,3,4,0,2,5,1,10,0,7,6,5,0,11,2,3,7,3,5,4,2,1,2,
       4,0,2,2,2,0,6,2,3,4,2,3,7,3,5,2,5,0,4,4,6,3,1,2,7,3,0,2,5,
       7,2,2,0,0,0,6,3,0,1,1,5,5,2,6,2,4,6,0,1,2,3,2,2,2,3,4,1,1,
       4,0,2,0,1,3,4,1,2,2,3,1,4,4,3,4,4,1,5,2,13,4,10,5,6,1,0,5,
       0,0,5,6,0,1,8,5,1,3,1,8,1,8,1,6,7,2,8,2,2,3,3,0,4,2,1,9,6,
       0,6,7,1,8,2,2,1,11,3,0,4,2,5,1,6,8,3,4,7,0,4,2,4,1,1,1,6,0,
       4,4,6,2,1,3,1,0,4,9,3,1,4,2,2,0,1};
   par = mle_nb(x5, 500);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a-3.0057) < 1e-3);
   test_assert(par->pi == 1.0);
   test_assert(fabs(par->p-0.4854) < 1e-3);
   free(par);

   return;

}


void
test_mle_zinb
(void)
{

   // Cases checked by simulated annealing.
   
   zinb_par_t *par;
   
   // 0:53, 1:8, 2:9, 3:8, 4:4, 5:7, 6:3, 7:3, 8:1, 9:3, 10:1
   int x1[100] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,3,7,4,5,3,5,1,5,7,3,6,1,2,9,1,10,6,2,2,2,3,2,1,1,5,0,2,3,
      9,4,2,9,3,5,7,3,5,2,1,0,6,1,4,2,3,4,5,8,1 };

   par = mle_zinb(x1, 100);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a - 3.6855) < 1e-3);
   test_assert(fabs(par->pi - 0.5110) < 1e-3);
   test_assert(fabs(par->p - 0.5044) < 1e-3);
   free(par);

   // 0:73, 1:7, 2:7, 3:3, 4:4, 5:2, 7:1, 8:2, 9:1
   int x2[100] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
      1,0,0,8,1,1,3,0,0,8,2,4,1,2,0,3,2,0,0,4,0,0,3,1,0,2,0,0,5,
      7,0,0,2,4,0,2,1,0,0,0,0,0,0,0,1,2,9,0,4 };

   par = mle_zinb(x2, 100);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->a - 1.8251) < 1e-3);
   test_assert(fabs(par->pi - 0.3363) < 1e-3);
   test_assert(fabs(par->p - 0.4109) < 1e-3);
   free(par);

   return;

}

void
test_fail_mle_nb
(void)
{

   int obs[] = {6,2,2,3,2,0,1,2,5,2};
   zinb_par_t *par = NULL;

   redirect_stderr();
   set_alloc_failure_countdown_to(0);
   par = mle_nb(obs, 10);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strncmp(caught_in_stderr(), "memory error", 12) == 0);

   redirect_stderr();
   set_alloc_failure_countdown_to(1);
   par = mle_nb(obs, 10);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strncmp(caught_in_stderr(), "memory error", 12) == 0);

   redirect_stderr();
   set_alloc_failure_countdown_to(2);
   par = mle_nb(obs, 10);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strncmp(caught_in_stderr(), "memory error", 12) == 0);
   
   redirect_stderr();
   set_alloc_failure_countdown_to(3);
   par = mle_nb(obs, 10);
   reset_alloc();
   unredirect_stderr();
   test_assert_critical(par != NULL);
   test_assert(strncmp(caught_in_stderr(), "$", 1) == 0);

   free(par);

   return;

}

void
test_fail_mle_zinb
(void)
{

   int obs[] = {6,2,2,3,2,0,1,2,5,2,0,0,0};
   zinb_par_t *par = NULL;

   // Test memory allocation errors.
   redirect_stderr();
   set_alloc_failure_countdown_to(0);
   par = mle_zinb(obs, 13);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strncmp(caught_in_stderr(), "memory error", 12) == 0);

   redirect_stderr();
   set_alloc_failure_countdown_to(0);
   par = mle_nb(obs, 13);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strncmp(caught_in_stderr(), "memory error", 12) == 0);

   redirect_stderr();
   set_alloc_failure_countdown_to(1);
   par = mle_zinb(obs, 13);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strncmp(caught_in_stderr(), "memory error", 12) == 0);

   redirect_stderr();
   set_alloc_failure_countdown_to(1);
   par = mle_nb(obs, 13);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strncmp(caught_in_stderr(), "memory error", 12) == 0);

   redirect_stderr();
   set_alloc_failure_countdown_to(2);
   par = mle_zinb(obs, 13);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strncmp(caught_in_stderr(), "memory error", 12) == 0);
   
   redirect_stderr();
   set_alloc_failure_countdown_to(2);
   par = mle_nb(obs, 13);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strncmp(caught_in_stderr(), "memory error", 12) == 0);
   
   redirect_stderr();
   set_alloc_failure_countdown_to(3);
   par = mle_zinb(obs, 13);
   reset_alloc();
   unredirect_stderr();
   test_assert_critical(par != NULL);
   test_assert(strncmp(caught_in_stderr(), "$", 1) == 0);
   free(par);

   redirect_stderr();
   set_alloc_failure_countdown_to(3);
   par = mle_nb(obs, 13);
   reset_alloc();
   unredirect_stderr();
   test_assert_critical(par != NULL);
   test_assert(strncmp(caught_in_stderr(), "$", 1) == 0);
   free(par);

   // Test fit failure.
   int fail[] = {0,0,0,0,0,0};
   redirect_stderr();
   par = mle_nb(fail, 6);
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strcmp(caught_in_stderr(), "$") == 0);

   free(par);
   return;

}


void
dummy_handler
(
   const char * a,
   const char * b,
          int   c
)
{
   fprintf(stderr, "dummy");
}


void
test_err_handler
(void)
{

   int obs[] = {6,2,2,3,2,0,1,2,5,2,0,0,0};
   zinb_par_t *par = NULL;

   set_zinb_err_handler(dummy_handler);

   redirect_stderr();
   set_alloc_failure_countdown_to(0);
   par = mle_zinb(obs, 13);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strcmp(caught_in_stderr(), "dummy") == 0);

   redirect_stderr();
   set_alloc_failure_countdown_to(0);
   par = mle_zinb(obs, 13);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strcmp(caught_in_stderr(), "dummy") == 0);

   // Reset error handler.
   set_zinb_err_handler(NULL);

   redirect_stderr();
   set_alloc_failure_countdown_to(0);
   par = mle_zinb(obs, 13);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strncmp(caught_in_stderr(), "memory error", 12) == 0);

   redirect_stderr();
   set_alloc_failure_countdown_to(0);
   par = mle_nb(obs, 13);
   reset_alloc();
   unredirect_stderr();
   test_assert(par == NULL);
   test_assert(strncmp(caught_in_stderr(), "memory error", 12) == 0);

   free(par);

}
