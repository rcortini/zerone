#include "hmm.h"
#include "unittest.h"
#include "utils.h"
#include "zinm.h"

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

   size_t x[] = {1,2,3,4,5,6,7,8,9,10,11,12};
   tab_t *tab;
   tab = tabulate(x, 3, 4);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 4);
   test_assert(tab->val[0] == 15);
   test_assert(tab->val[1] == 18);
   test_assert(tab->val[2] == 21);
   test_assert(tab->val[3] == 24);
   for (int i = 0 ; i < 4 ; i++) {
      test_assert(tab->num[i] == 1);
   }
   free(tab);

   tab = tabulate(x, 4, 3);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 3);
   test_assert(tab->val[0] == 22);
   test_assert(tab->val[1] == 26);
   test_assert(tab->val[2] == 30);
   for (int i = 0 ; i < 3 ; i++) {
      test_assert(tab->num[i] == 1);
   }
   free(tab);

   size_t y[] = {4096,2048,1024,512};
   tab = tabulate(y, 1, 4);
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

   tab = tabulate(y, 2, 2);
   test_assert_critical(tab != NULL);
   test_assert(tab->size == 2);
   test_assert(tab->val[0] == 2560);
   test_assert(tab->val[1] == 5120);
   for (int i = 0 ; i < 2 ; i++) {
      test_assert(tab->num[i] == 1);
   }
   free(tab);

   return;

}


void
test_compute_means
(void)
{

   double mean;

   size_t x1[1] = {2595};
   compute_means(x1, 1, 1, &mean);
   test_assert(mean == 2595.0);

   size_t x2[7] = {0,0,0,1,1,2,5};
   compute_means(x2, 1, 7, &mean);
   test_assert(fabs(mean-1.28571428) < 1e-6);

   size_t x3[5] = {0,89,231,55,309};
   compute_means(x3, 1, 5, &mean);
   test_assert(fabs(mean-136.8) < 1e-6);

   // 0:112, 1:94, 2:28, 3:12, 4:3, 7:1
   size_t x4[250] = {0,0,0,3,0,0,1,1,1,1,1,2,0,2,0,0,1,0,0,0,1,1,0,1,
      1,0,1,1,0,2,1,0,2,1,1,0,2,1,1,1,1,1,0,0,2,0,2,1,1,1,2,1,0,0,
      1,0,1,0,0,1,0,0,3,2,0,0,0,0,0,2,1,1,1,0,0,1,0,0,1,0,0,1,0,1,
      0,1,2,1,2,1,0,0,0,2,0,0,0,1,2,1,0,1,1,1,2,0,0,0,0,0,2,1,3,0,
      2,3,0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,1,0,3,1,0,0,0,1,1,0,0,0,0,
      0,1,0,0,0,0,1,2,1,0,2,4,0,1,0,1,0,1,0,1,1,1,0,1,1,1,1,0,0,1,
      0,1,1,3,1,1,1,1,0,0,0,0,3,1,3,0,1,1,0,0,0,1,1,0,1,2,4,2,0,0,
      4,0,2,1,0,0,2,1,2,1,7,1,2,3,0,0,1,1,0,3,1,1,1,3,1,1,0,0,0,0,
      1,2,2,0,1,1,0,1,1,1,0,2,3,0,0,0};
   compute_means(x4, 1, 250, &mean);
   test_assert(fabs(mean-0.82) < 1e-6);

   return;

}


void
test_eval_nb_f
(void)
{

   // 0:14, 1:5, 2:4, 3:1, 5:1
   size_t x[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                     1,1,1,1,1,2,2,2,2,3,5 };
   tab_t *tab = tabulate(x, 1, 25);
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
   size_t x[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                    1,1,1,1,1,2,2,2,2,3,5 };
   tab_t *tab = tabulate(x, 1, 25);
   double mean;
   compute_means(x, 1, 25, &mean);
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
test_eval_zinm_f
(void)
{

   test_assert(fabs(eval_zinm_f(1.0, .5, 1, 2.0)) < 1e-6);
   test_assert(fabs(eval_zinm_f(1.0, .5, 2, 4.0)) < 1e-6);
   test_assert(fabs(eval_zinm_f(1.3, .7, 141, 8.3)-678.0838351) < 1e-6);

   return;

}


void
test_eval_zinm_g
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   size_t x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(eval_zinm_g(1, .5, tab)+0.1325713) < 1e-6);
   free(tab);

   return;

}


void
test_eval_zinm_dfda
(void)
{

   test_assert(fabs(eval_zinm_dfda(1, .5, 1)-1.2274112) < 1e-6);
   test_assert(fabs(eval_zinm_dfda(1, .5, 9)-11.0467015) < 1e-6);
   test_assert(fabs(eval_zinm_dfda(2, .5, 9)-12.9096451) < 1e-6);
   test_assert(fabs(eval_zinm_dfda(2, .3, 9)-25.1159846) < 1e-6);
   test_assert(fabs(eval_zinm_dfda(2.4, .3, 9)-26.3620864) < 1e-6);

   return;

}


void
test_eval_zinm_dfdp
(void)
{

   test_assert(eval_zinm_dfdp(1, .5, 0, 0) == 0.0);
   test_assert(eval_zinm_dfdp(1, .5, 1, 0) == 0.0);
   test_assert(fabs(eval_zinm_dfdp(2, .5, 1, 0)+3.5555555) < 1e-6);
   test_assert(fabs(eval_zinm_dfdp(2, .3, 1, 0)+19.5896899) < 1e-6);
   test_assert(fabs(eval_zinm_dfdp(2, .3, 9, 0)+176.3072092) < 1e-6);
   test_assert(fabs(eval_zinm_dfdp(2, .3, 9, 1.7)+179.7765970) < 1e-6);

   return;

}


void
test_eval_zinm_dgda
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   size_t x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(eval_zinm_dgda(1, .5, tab)+2.2547559) < 1e-6);
   test_assert(fabs(eval_zinm_dgda(2, .5, tab)+1.2605630) < 1e-6);
   test_assert(fabs(eval_zinm_dgda(2, .3, tab)+1.8764955) < 1e-6);
   free(tab);

   return;

}


void
test_ll_zinm
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   size_t x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(ll_zinm(1, .5, 1, tab)+22.5329303) < 1e-6);
   test_assert(fabs(ll_zinm(1, .5, .7, tab)+22.7832550) < 1e-6);
   test_assert(fabs(ll_zinm(2, .5, .7, tab)+23.7608409) < 1e-6);
   test_assert(fabs(ll_zinm(2, .3, .7, tab)+31.6978553) < 1e-6);
   free(tab);

   return;

}


void
test_nb_est_alpha
(void)
{

   tab_t *tab;
   // 0:14, 1:5, 2:4, 3:1, 5:1
   size_t x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   tab = tabulate(x1, 1, 25);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-0.9237) < 1e-3);
   free(tab);

   // 0:27, 1:12, 2:8, 3:1, 4:1, 5:1
   size_t x2[50] = {3,0,1,2,0,0,1,0,0,0,0,1,1,0,0,1,2,2,0,0,0,1,2,
      0, 0,0,0,0,4,0,0,0,1,5,1,0,1,2,1,2,2,2,0,0,0,1,0,1,0,0};
   tab = tabulate(x2, 1, 50);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-1.3436) < 1e-3);
   free(tab);

   // 0:12, 1:7, 2:13, 3:4, 4:6, 5:2, 6:1, 7:3, 8:1, 9:1
   size_t x3[50] = {4,5,2,1,2,4,2,2,0,4,2,1,3,6,0,0,7,3,0,8,4,2,0,
      0,2,3,2,3,7,9,2,4,0,4,2,0,0,2,5,1,1,2,1,0,0,0,1,2,1,7};
   tab = tabulate(x3, 1, 50);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-1.7969) < 1e-3);
   free(tab);

   // 0:39, 1:8, 2:2, 3:1
   size_t x4[50] = {1,0,0,0,0,0,3,1,0,1,0,0,0,0,0,0,0,1,0,0,0,2,0,
      2,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0};
   tab = tabulate(x4, 1, 50);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-0.7073) < 1e-3);
   free(tab);

   // 0:59, 1:83, 2:99, 3:67, 4:67, 5:49, 6:27, 7:22, 8:11, 9:6
   // 10:6, 11:3, 12:2, 13:3
   size_t x5[500] = {1,0,0,1,1,2,1,0,5,7,1,3,3,1,6,0,2,5,7,0,5,2,1,
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
   tab = tabulate(x5, 1, 500);
   test_assert_critical(tab != NULL);
   test_assert(fabs(nb_est_alpha(tab)-3.0057) < 1e-3);
   free(tab);
}


void
test_mle_nm
(void)
{

   // These test cases have been verified with R.
   zinm_par_t *par;
   
   // 0:14, 1:5, 2:4, 3:1, 5:1
   size_t x1[25] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                  1,1,1,1,1,2,2,2,2,3,5 };
   par = mle_nm(x1, 1, 25);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->alpha-0.9237) < 1e-3);
   test_assert(fabs(par->p[0]-0.5237) < 1e-3);
   test_assert(fabs(par->p[1]-0.4763) < 1e-3);
   free(par);

   // 0:27, 1:12, 2:8, 3:1, 4:1, 5:1
   size_t x2[50] = {3,0,1,2,0,0,1,0,0,0,0,1,1,0,0,1,2,2,0,0,0,1,2,
      0, 0,0,0,0,4,0,0,0,1,5,1,0,1,2,1,2,2,2,0,0,0,1,0,1,0,0};
   par = mle_nm(x2, 1, 50);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->alpha-1.3436) < 1e-3);
   test_assert(fabs(par->p[0]-0.6267) < 1e-3);
   test_assert(fabs(par->p[1]-0.3732) < 1e-3);
   free(par);

   // 0:12, 1:7, 2:13, 3:4, 4:6, 5:2, 6:1, 7:3, 8:1, 9:1
   size_t x3[50] = {4,5,2,1,2,4,2,2,0,4,2,1,3,6,0,0,7,3,0,8,4,2,0,
      0,2,3,2,3,7,9,2,4,0,4,2,0,0,2,5,1,1,2,1,0,0,0,1,2,1,7};
   par = mle_nm(x3, 1, 50);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->alpha-1.7969) < 1e-3);
   test_assert(fabs(par->p[0]-0.4221) < 1e-3);
   test_assert(fabs(par->p[1]-0.5779) < 1e-3);
   free(par);

   // 0:39, 1:8, 2:2, 3:1
   size_t x4[50] = {1,0,0,0,0,0,3,1,0,1,0,0,0,0,0,0,0,1,0,0,0,2,0,
      2,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0};
   par = mle_nm(x4, 1, 50);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->alpha-0.7073) < 1e-3);
   test_assert(fabs(par->p[0]-0.7021) < 1e-3);
   test_assert(fabs(par->p[1]-0.2978) < 1e-3);
   free(par);

   // 0:59, 1:83, 2:99, 3:67, 4:67, 5:49, 6:27, 7:22, 8:11, 9:6
   // 10:6, 11:3, 12:2, 13:3
   size_t x5[500] = {1,0,0,1,1,2,1,0,5,7,1,3,3,1,6,0,2,5,7,0,5,2,1,
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
   par = mle_nm(x5, 1, 500);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->alpha-3.0057) < 1e-3);
   test_assert(fabs(par->p[0]-0.4854) < 1e-3);
   test_assert(fabs(par->p[1]-0.5145) < 1e-3);
   free(par);

   return;

}


void
test_mle_zinm
(void)
{

   // Cases checked by simulated annealing.
   
   zinm_par_t *par;
   
   // 0:53, 1:8, 2:9, 3:8, 4:4, 5:7, 6:3, 7:3, 8:1, 9:3, 10:1
   size_t x1[100] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,3,7,4,5,3,5,1,5,7,3,6,1,2,9,1,10,6,2,2,2,3,2,1,1,5,0,2,3,
      9,4,2,9,3,5,7,3,5,2,1,0,6,1,4,2,3,4,5,8,1 };

   par = mle_zinm(x1, 1, 100);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->alpha-3.6855) < 1e-3);
   test_assert(fabs(par->p[0]-0.5044) < 1e-3);
   test_assert(fabs(par->pi-0.5110) < 1e-3);
   free(par);

   // 0:73, 1:7, 2:7, 3:3, 4:4, 5:2, 7:1, 8:2, 9:1
   size_t x2[100] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,
      1,0,0,8,1,1,3,0,0,8,2,4,1,2,0,3,2,0,0,4,0,0,3,1,0,2,0,0,5,
      7,0,0,2,4,0,2,1,0,0,0,0,0,0,0,1,2,9,0,4 };

   par = mle_zinm(x2, 1, 100);
   test_assert_critical(par != NULL);
   test_assert(fabs(par->alpha-1.8251) < 1e-3);
   test_assert(fabs(par->p[0]-0.4109) < 1e-3);
   test_assert(fabs(par->pi-0.3363) < 1e-3);
   free(par);

   return;

}

void
test_indexts
(void)
{
   int index[8];

   int ts_1[24] = {
      0, 0, 0,
      0, 0, 1,
      0, 1, 0,
      1, 0, 1,
      0, 0, 1,
      0, 0, 1,
      0, 1, 0,
      0, 0, 0,
   };

   int ts_2[24] = {
      -100, 0, 1,
      -1,   0, 100,
      -1,   0, 1,
      -1,   0, 100,
      -100, 0, 1,
      -1,   0, 1,
      -100, 0, 100,
      -100, 0, 100,
   };

   int expected_1[8] = {0,1,2,3,1,1,2,0};
   int expected_2[8] = {0,1,2,1,0,2,6,6};

   indexts(8, 3, ts_1, index);
   for (int i = 0 ; i < 8 ; i++) {
      test_assert(index[i] == expected_1[i]);
   }

   indexts(8, 3, ts_2, index);
   for (int i = 0 ; i < 8 ; i++) {
      test_assert(index[i] == expected_2[i]);
   }

}


void
test_fwdb
(void)
{
   int n = 3;
   int m = 2;
   const double Q[4] = {
      // transpose //
      0.8, 0.1,
      0.2, 0.9,
   };
   double init[2] = {.8, .2};
   double prob[6] = {
      // First series.
      0.1,      0.3,
      0.4,      0.8,
      // `fwd` can deal with emission probabilities in log space.
      log(0.2), log(0.7),
   };

   double *phi = malloc(m*n * sizeof(double));
   double *trans = malloc(m*m * sizeof(double));
   test_assert_critical(phi != NULL && trans != NULL);

   double ll = fwdb(m, n, Q, init, prob, phi, trans);

   double expected_loglik = -3.105547;
   double expected_alpha[6] = {
      0.5714286, 0.4285714,
      0.3333333, 0.6666667,
      0.1250000, 0.8750000,
   };
   double expected_phi[6] = {
      0.3571429, 0.6428571,
      0.1875000, 0.8125000,
      0.1250000, 0.8750000,
   };
   double expected_trans[4] = {
      // transpose //
      0.2714285, 0.0410714,
      0.2732143, 1.4142857,
   };

   test_assert(fabs(expected_loglik - ll) < 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      test_assert(fabs(expected_alpha[i] - prob[i]) < 1e-6);
      test_assert(fabs(expected_phi[i] - phi[i]) < 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      test_assert(fabs(expected_trans[i] - trans[i]) < 1e-6);
   }

   free(phi);
   free(trans);

   return;

}


void
test_fwdb_NA
(void)
{

   int n = 3;
   int m = 2;
   const double Q[4] = {
      // transpose //
      0.8, 0.1,
      0.2, 0.9,
   };
   double init[2] = {.8, .2};
   double prob[6] = {
      // First series (NA in first position).
      log(-1.0),     0.3,
            0.4,     0.8,
      // `fwd` can deal log emission probabilities.
      log(0.2), log(0.7),
   };

   double *phi = malloc(m*n * sizeof(double));
   double *trans = malloc(m*m * sizeof(double));
   test_assert_critical(phi != NULL && trans != NULL);

   double ll = fwdb(m, n, Q, init, prob, phi, trans);

   double expected_loglik = -1.362578;
   double expected_alpha[12] = {
      0.8000000, 0.2000000,
      0.4925373, 0.5074627,
      0.1862500, 0.8137500,
   };
   double expected_phi[12] = {
      0.6250000, 0.3750000,
      0.3093750, 0.6906250,
      0.1862500, 0.8137500,
   };
   double expected_trans[4] = {
      // transpose //
      0.4650000, 0.0306250,
      0.4693750, 1.0350000,
   };

   test_assert(fabs(expected_loglik - ll) < 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      test_assert(fabs(expected_alpha[i] - prob[i]) < 1e-6);
      test_assert(fabs(expected_phi[i] - phi[i]) < 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      test_assert(fabs(expected_trans[i] - trans[i]) < 1e-6);
   }

   free(phi);
   free(trans);

   return;

}


void
test_underflow
(void)
{
   int n = 3;
   int m = 2;
   const double Q[4] = {
      // transpose //
      0.8, 0.1,
      0.2, 0.9,
   };
   double init[2] = {.8, .2};
   double prob_1[6] = {
      // First series (underflow at first step).
      0.0,      0.0,
      0.4,      0.8,
      // `fwd` can deal log emission probabilities.
      log(0.2), log(0.7),
   };

   double *phi = malloc(m*n * sizeof(double));
   double *trans = malloc(m*m * sizeof(double));
   test_assert_critical(phi != NULL && trans != NULL);

   double ll_1 = fwdb(m, n, Q, init, prob_1, phi, trans);

   double expected_loglik_1 = -1.362578;
   double expected_alpha_1[6] = {
      0.8000000, 0.2000000,
      0.4925373, 0.5074627,
      0.1862500, 0.8137500,
   };
   double expected_phi_1[6] = {
      0.6250000, 0.3750000,
      0.3093750, 0.6906250,
      0.1862500, 0.8137500,
   };
   double expected_trans_1[4] = {
      // transpose //
      0.4650000, 0.0306250,
      0.4693750, 1.0350000,
   };

   test_assert(fabs(expected_loglik_1 - ll_1) < 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      test_assert(fabs(expected_alpha_1[i] - prob_1[i]) < 1e-6);
      test_assert(fabs(expected_phi_1[i] - phi[i]) < 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      test_assert(fabs(expected_trans_1[i] - trans[i]) < 1e-6);
   }

   // Create the underflow at step 2.
   double prob_2[6] = {
      0.1,      0.3,
      0.0,      0.0,
      // `fwd` can deal log emission probabilities.
      log(0.2), log(0.7),
   };

   double ll_2 = fwdb(m, n, Q, init, prob_2, phi, trans);

   double expected_loglik_2 = -2.710553;
   double expected_alpha_2[6] = {
      0.5714286, 0.4285714,
      0.5000000, 0.5000000,
      0.1894737, 0.8105263,
   };
   double expected_phi_2[6] = {
      0.4451128, 0.5548872,
      0.3157895, 0.6842105,
      0.1894737, 0.8105263,
   };
   double expected_trans_2[4] = {
      // transpose //
      0.4571429, 0.0481203,
      0.3037594, 1.1909774,
   };

   test_assert(fabs(expected_loglik_2 - ll_2) < 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      test_assert(fabs(expected_alpha_2[i] - prob_2[i]) < 1e-6);
      test_assert(fabs(expected_phi_2[i] - phi[i]) < 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      test_assert(fabs(expected_trans_2[i] - trans[i]) < 1e-6);
   }

   free(phi);
   free(trans);

   return;

}

void
test_block_fwdb
(void)
{

   int m = 2;
   int n = 6;
   double Q[4] = {
      // transpose //
      0.8, 0.1,
      0.2, 0.9,
   };
   double init[2] = {.8, .2};
   double prob[12] = {
      // First series.
      0.1,      0.3,
      0.4,      0.8,
      // `fwd` can deal log emission probabilities.
      log(0.2), log(0.7),
      // Second series (identical).
      0.1,      0.3,
      0.4,      0.8,
      0.2,      0.7,
   };

   double *phi = malloc(m*n * sizeof(double));
   double *trans = malloc(m*m * sizeof(double));
   test_assert_critical(phi != NULL && trans != NULL);

   int nblocks = 2;
   int size[2] = {3,3};

   double ll;
   block_fwdb(&m, &nblocks, size, Q, init, prob, phi, trans, &ll);

   double expected_loglik = 2 * -3.105547;
   double expected_phi[12] = {
      0.3571429, 0.6428571,
      0.1875000, 0.8125000,
      0.1250000, 0.8750000,
      0.3571429, 0.6428571,
      0.1875000, 0.8125000,
      0.1250000, 0.8750000,
   };
   double expected_trans[4] = {
      // transpose //
      2 * 0.2714285, 2 * 0.0410714,
      2 * 0.2732143, 2 * 1.4142857,
   };

   test_assert(fabs(expected_loglik - ll) < 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      test_assert(fabs(expected_phi[i] - phi[i]) < 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      test_assert(fabs(expected_trans[i] - trans[i]) < 1e-6);
   }

   free(phi);
   free(trans);

   return;

}


void
test_block_fwdb_NA
(void)
{
   int m = 2;
   int n = 6;
   double Q[4] = {
      // transpose //
      0.8, 0.1,
      0.2, 0.9,
   };
   double init[2] = {.8, .2};
   double prob[12] = {
      // First series.
      log(-1.0), 0.3,
      0.4,       0.8,
      // `fwd` can deal log emission probabilities.
      log(0.2),  log(0.7),
      // Second series (identical).
      0.1,       0.3,
      log(-1.0), 0.8,
      0.2,       0.7,
   };

   double *phi = malloc(m*n * sizeof(double));
   double *trans = malloc(m*m * sizeof(double));
   test_assert_critical(phi != NULL && trans != NULL);

   int nblocks = 2;
   int size[2] = {3,3};

   double ll;
   block_fwdb(&m, &nblocks, size, Q, init, prob, phi, trans, &ll);

   double expected_loglik = -1.362578 -2.710553;
   double expected_alpha[12] = {
      0.8000000, 0.2000000,
      0.4925373, 0.5074627,
      0.1862500, 0.8137500,
      0.5714286, 0.4285714,
      0.5000000, 0.5000000,
      0.1894737, 0.8105263,
   };
   double expected_phi[12] = {
      0.6250000, 0.3750000,
      0.3093750, 0.6906250,
      0.1862500, 0.8137500,
      0.4451128, 0.5548872,
      0.3157895, 0.6842105,
      0.1894737, 0.8105263,
   };
   double expected_trans[4] = {
      // transpose //
      0.4650000+0.4571429, 0.0306250+0.0481203,
      0.4693750+0.3037594, 1.0350000+1.1909774,
   };

   test_assert(fabs(expected_loglik - ll) < 1e-6);
   for (int i = 0 ; i < m*n ; i++) {
      test_assert(fabs(expected_alpha[i] - prob[i]) < 1e-6);
      test_assert(fabs(expected_phi[i] - phi[i]) < 1e-6);
   }
   for (int i = 0 ; i < m*m ; i++) {
      test_assert(fabs(expected_trans[i] - trans[i]) < 1e-6);
   }

   free(phi);
   free(trans);

   return;

}


void
test_viterbi
(void)
{
   // -- Test 1, this is a test case I worked out manually with R -- //
   int m = 4;
   int n = 5;
   double Q_1[16] = {
      // transpose //
      log(0.80), log(0.05), log(0.05), log(0.05),
      log(0.10), log(0.90), log(0.30), log(0.35),
      log(0.05), log(0.05), log(0.50), log(0.45),
      log(0.05), log(0.00), log(0.25), log(0.15),
   };

   double init_1[4] = {log(0.05), log(0.10), log(0.05), log(0.80)};
   double prob_1[20] = {
      log(0.1), log(0.3), log(0.2), log(0.3),
      log(0.2), log(0.4), log(0.8), log(0.1),
      log(0.2), log(0.7), log(0.1), log(0.2),
      log(0.8), log(0.1), log(0.3), log(0.4),
      log(0.9), log(0.1), log(0.2), log(0.1),
   };

   int *path = malloc(n * sizeof(int));
   test_assert_critical(path != NULL);

   viterbi(
                    m,
                    n,
   (const double *) Q_1,
   (const double *) init_1,
   (const double *) prob_1,
                    path
   );

   //   state 0    state 1    state 2    state 3
   // -5.298317  -3.506558  -4.605170  -1.427116
   // -6.032287  -3.393229  -2.448768  -5.626821
   // -7.053938  -3.855265  -5.444500  -5.444500
   // -7.074140  -6.263210  -7.341620  -7.747085
   // -7.402645  -8.671156  -9.644205 -11.030499

   double expected_path_1[5] = {
      3,1,1,0,0,
   };

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path_1[i] == path[i]);
   }

   free(path);
   path = NULL;

   //-- Test 2, this is a test using Wikipedia's Python example --//
   m = 2;
   n = 23;
   double Q_2[4] = {
      // transpose //
      log(0.7), log(0.3),
      log(0.4), log(0.6),
   };

   double init_2[2] = {log(0.6), log(0.4)};
   double prob_2[46] = {
      log(0.5), log(0.1), // normal
      log(0.4), log(0.3), // cold  
      log(0.1), log(0.6), // dizzy 
      log(0.1), log(0.6), // dizzy 
      log(0.5), log(0.1), // normal
      log(0.5), log(0.1), // normal
      log(0.5), log(0.1), // normal
      log(0.1), log(0.6), // dizzy 
      log(0.1), log(0.6), // dizzy 
      log(0.1), log(0.6), // dizzy 
      log(0.4), log(0.3), // cold  
      log(0.4), log(0.3), // cold  
      log(0.4), log(0.3), // cold  
      log(0.1), log(0.6), // dizzy 
      log(0.5), log(0.1), // normal
      log(0.5), log(0.1), // normal
      log(0.4), log(0.3), // cold  
      log(0.5), log(0.1), // normal
      log(0.4), log(0.3), // cold  
      log(0.1), log(0.6), // dizzy 
      log(0.1), log(0.6), // dizzy 
      log(0.5), log(0.1), // normal
      log(0.5), log(0.1), // normal
   };

   path = malloc(n * sizeof(int));
   test_assert_critical(path != NULL);

   viterbi(
                    m,
                    n,
   (const double *) Q_2,
   (const double *) init_2,
   (const double *) prob_2,
                    path
   );

   double expected_path_2[23] = {
      0,0,1,1,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,1,1,0,0,
   };

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path_2[i] == path[i]);
   }

   free(path);

   return;

}


void
test_block_viterbi
(void)
{

   int m = 4;
   int n = 10;
   const double Q[16] = {
      // transpose //
      0.80, 0.05, 0.05, 0.05,
      0.10, 0.90, 0.30, 0.35,
      0.05, 0.05, 0.50, 0.45,
      0.05, 0.00, 0.25, 0.15,
   };
   double init[4] = {0.05, 0.10, 0.05, 0.80};
   double prob[40] = {
      // First series.
      0.1, 0.3, 0.2, 0.3,
      0.2, 0.4, 0.8, 0.1,
      0.2, 0.7, 0.1, 0.2,
      0.8, 0.1, 0.3, 0.4,
      0.9, 0.1, 0.2, 0.1,
      // Second series (identical).
      0.1, 0.3, 0.2, 0.3,
      0.2, 0.4, 0.8, 0.1,
      0.2, 0.7, 0.1, 0.2,
      0.8, 0.1, 0.3, 0.4,
      0.9, 0.1, 0.2, 0.1,
   };

   int *path = malloc(n * sizeof(int));
   test_assert_critical(path != NULL);

   int nblocks = 2;
   static int size[2] = {5,5};
   int nolog = 0;

   block_viterbi(
         &m,
         &nblocks,
         size,
         Q,
         init,
         prob,
         &nolog,
         path
   );

   double expected_path[10] = {
      3,1,1,0,0,3,1,1,0,0,
   };

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path[i] == path[i]);
   }


   int yeslog = 1;
   const double log_Q[16] = {
      // transpose //
      log(0.80), log(0.05), log(0.05), log(0.05),
      log(0.10), log(0.90), log(0.30), log(0.35),
      log(0.05), log(0.05), log(0.50), log(0.45),
      log(0.05), log(0.00), log(0.25), log(0.15),
   };
   double log_init[4] = {log(0.05), log(0.10), log(0.05), log(0.80)};
   double log_prob[40] = {
      // First series.
      log(0.1), log(0.3), log(0.2), log(0.3),
      log(0.2), log(0.4), log(0.8), log(0.1),
      log(0.2), log(0.7), log(0.1), log(0.2),
      log(0.8), log(0.1), log(0.3), log(0.4),
      log(0.9), log(0.1), log(0.2), log(0.1),
      // Second series (identical).
      log(0.1), log(0.3), log(0.2), log(0.3),
      log(0.2), log(0.4), log(0.8), log(0.1),
      log(0.2), log(0.7), log(0.1), log(0.2),
      log(0.8), log(0.1), log(0.3), log(0.4),
      log(0.9), log(0.1), log(0.2), log(0.1),
   };

   block_viterbi(
         &m,
         &nblocks,
         size,
         log_Q,
         log_init,
         log_prob,
         &yeslog,
         path
   );

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path[i] == path[i]);
   }

   return;

}


void
test_block_viterbi_NA
(void)
{

   int m = 4;
   int n = 10;
   const double Q[16] = {
      // transpose //
      0.80, 0.05, 0.05, 0.05,
      0.10, 0.90, 0.30, 0.35,
      0.05, 0.05, 0.50, 0.45,
      0.05, 0.00, 0.25, 0.15,
   };
   double init[4] = {0.05, 0.10, 0.05, 0.80};
   double prob[40] = {
      // First series.
     -1.0, 0.3, 0.2, 0.3,
      0.2, 0.4, 0.8, 0.1,
      0.2, 0.7, 0.1, 0.2,
      0.8, 0.1, 0.3, 0.4,
      0.9, 0.1, 0.2, 0.1,
      // Second series.
      0.1, 0.3, 0.2, 0.3,
     -1.0, 0.4, 0.8, 0.1,
      0.0, 0.0, 0.0, 0.0,
      0.8, 0.1, 0.3, 0.4,
      0.9, 0.1, 0.2, 0.1,
   };

   int *path = malloc(n * sizeof(int));
   test_assert_critical(path != NULL);

   int nblocks = 2;
   static int size[2] = {5,5};
   int nolog = 0;

   block_viterbi(
         &m,
         &nblocks,
         size,
         Q,
         init,
         prob,
         &nolog,
         path
   );

   double expected_path[10] = {
      3,1,1,0,0,3,0,0,0,0,
   };

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path[i] == path[i]);
   }

   // Set invalid initial probabilities.
   static double invalid_init[4] = {-1.0, 0.10, 0.05, 0.80};

   // Catch stderr.
   redirect_stderr();
   int exit_status = block_viterbi(
         &m,
         &nblocks,
         size,
         Q,
         invalid_init,
         prob,
         &nolog,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("invalid 'init' parameter in 'block_viterbi'\n");

   // Set invalid transition probabilities.
   static double invalid_Q[16] = {
      // transpose //
      -1.0, 0.05, 0.05, 0.05,
      0.10, 0.90, 0.30, 0.35,
      0.05, 0.05, 0.50, 0.45,
      0.05, 0.00, 0.25, 0.15,
   };

   // Catch stderr.
   redirect_stderr();
   exit_status = block_viterbi(
         &m,
         &nblocks,
         size,
         invalid_Q,
         init,
         prob,
         &nolog,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("invalid 'Q' parameter in 'block_viterbi'\n");

   free(path);

   return;

}


void
test_zinm_prob
(void)
{

   int m = 3;
   int n = 7;
   int r = 3;

   const double pi = 0.8;
   const double a = 1.2;
   // The constant C1 facilitates the definition of 'p'.
   const double C1 = 2.5;
   const double p[12] = {
      C1, 1.0, 1.0, 1.0, 
      C1, 1.0, 2.0, 1.0,
      C1, 1.0, 0.2, 0.1,
   };

   const int yz[21] = {
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

   int index[7] = {-1,-1,-1,-1,-1,-1,-1};
   double pem[21];
   indexts(n, r, (const int *) yz, index);

   //--               Test output type 0               --//
   int out = 0;
   redirect_stderr();
   zinm_prob(&m, &n, &r, yz, &a, p, &pi, index, &out, pem);
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
   zinm_prob(&m, &n, &r, yz, &a, p, &pi, index, &out, pem);
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
   zinm_prob(&m, &n, &r, yz, &a, p, &pi, index, &out, pem);
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

   return;

}


int
main(
   int argc,
   char **argv
)
{

   // Register test cases //
   const static test_case_t test_cases[] = {
      {"utils/new_histo", test_new_histo},
      {"utils/histo_push", test_histo_push},
      {"utils/compress_histo", test_compress_histo},
      {"utils/tabulate", test_tabulate},
      {"utils/compute_means", test_compute_means},
      {"utils/indexts", test_indexts},
      {"eval_nb_f", test_eval_nb_f},
      {"eval_nb_dfda", test_eval_nb_dfda},
      {"eval_zinm_f", test_eval_zinm_f},
      {"eval_zinm_g", test_eval_zinm_g},
      {"eval_zinm_dfda", test_eval_zinm_dfda},
      {"eval_zinm_dfdp", test_eval_zinm_dfdp},
      {"eval_zinm_dgda", test_eval_zinm_dgda},
      {"ll_zinm", test_ll_zinm},
      {"nb_est_alpha", test_nb_est_alpha},
      {"mle_nm", test_mle_nm},
      {"mle_zinm", test_mle_zinm},
      {"HMM/fwdb", test_fwdb},
      {"HMM/fwdb (NAs)", test_fwdb_NA},
      {"HMM/fwdb (underflow)", test_underflow},
      {"HMM/block_fwdb", test_block_fwdb},
      {"HMM/block_fwdb (NAs)", test_block_fwdb_NA},
      {"HMM/viterbi", test_viterbi},
      {"HMM/block_viterbi", test_block_viterbi},
      {"HMM/block_viterbi (NAs)", test_block_viterbi_NA},
      {"zinm_prob", test_zinm_prob},
      {NULL, NULL}
   };

   return run_unittest(argc, argv, test_cases);

}

#if 0










void
test_mean
(void)
{

   // -- First test case -- //
   int yz_1[8] = {
      1, 2,
      3, 4,
      5, 6,
      7, 8,
   };

   g_assert_cmpfloat(mean(yz_1, 4, 2), ==, 4.0);
   g_assert_cmpfloat(mean(yz_1+1, 4, 2), ==, 5.0);

   // -- Second test case -- //
   int NA = (int) log(-1);
   int yz_2[8] = {
       1, NA,
       3, NA,
      NA, NA,
       8, NA,
   };

   g_assert_cmpfloat(mean(yz_2, 4, 2), ==, 4.0);
   g_assert_cmpfloat(mean(yz_2+1, 4, 2), !=, mean(yz_2+1, 4, 2));

}

void
test_histsum
(void)
{

   // -- First test case -- //
   int yz_1[8] = {
      1, 3,
      2, 2,
      3, 1,
      4, 0,
   };

   int expected_hist_1[5] = {0,0,0,0,4};
   int *hist_1 = histsum(yz_1, 4, 2);

   for (int i = 0 ; i < 5 ; i++) {
      g_assert_cmpint(hist_1[i], ==, expected_hist_1[i]);
   }

   free(hist_1);

   // -- Second test case -- //
   int yz_2[8] = {
       1, -3,
       2,  2,
      -3, -1,
       4,  0,
   };

   int expected_hist_2[5] = {0,0,0,0,2};
   int *hist_2 = histsum(yz_2, 4, 2);

   for (int i = 0 ; i < 5 ; i++) {
      g_assert_cmpint(hist_2[i], ==, expected_hist_2[i]);
   }

   free(hist_2);

   // -- Third test case -- //
   int yz_3[8] = {
       1, -3,
       2, -2,
      -3, -1,
       0,  1024,
   };
   int *hist_3 = histsum(yz_3, 4, 2);

   g_assert_cmpint(hist_3[1024], ==, 1);
   for (int i = 0 ; i < 1024 ; i++) {
      g_assert_cmpint(hist_3[i], ==, 0);
   }

   free(hist_3);

}



// -------------------------  mnmultinom.c ------------------------- //



void
test_perf_mnmultinom_prob
(void)
{
   // -- Open input file (the code is hopelessly not portable) -- //
   FILE *f = fopen("RXRA.txt", "r");

   int n = 1031884;
   int i = 0;

   char seqname[32];
   int *yz = malloc (3*n * sizeof(int));
   memset(yz, (int) -1, 3*n * sizeof(int));

   // `getline` is available because we define _GNU_SOURCE.
   char *line = NULL;
   size_t len = 0;
   ssize_t read;
   // Discard header (the 'if' turns off unused variable warnings).
   if (getline(&line, &len, f));
   while ((read = getline(&line, &len, f)) != -1) {
      sscanf(line, "%s\t%d\t%d\t%d\n", seqname, 
         yz+i, yz+i+1, yz+i+2);
      i += 3;
   }
   free(line);
   fclose(f);

   // --               Run the performance test               -- //
   int m = 3;
   int r = 3;

   double t = 0.91;
   double a = 2.844083;

   double p[12] = {
      // transpose //
      0.05155538, 0.75760266, 0.11489932, 0.07594263,
      0.05066482, 0.74451597, 0.15030383, 0.05451538,
      0.04539102, 0.66701779, 0.15180765, 0.13578353,
   };
   double q[12] = {
      // transpose //
      0.54836008, 0.38167944, 0.04426865, 0.02569183,
      0.52542413, 0.36571515, 0.07734838, 0.03151234,
      0.47888700, 0.33332350, 0.16109890, 0.02669060,
   };
   double Q[9] = {
      // transpose //
      0.8712810672, 0.0004510562, 0.4013737724,
      0.0016887750, 0.9669194640, 0.3386870900,
      0.1270301600, 0.0326294800, 0.2599391400,
   };
   double init[3] = {.33333, .33333, .33333};
   double trans[9];

   int *index = malloc(n * sizeof(int));
   double *pem = malloc(n*m * sizeof(double));
   double *phi = malloc(n*m * sizeof(double));
   indexts(n, r, (const int *) yz, index);

   int out = 0;
   double ll;
   int nblocks = 24;
   int size[24] = {
      83083, 45178, 45002, 44617, 38389, 35783, 34177,
      30118, 27065, 26025, 19709, 81066, 21008, 16043,
      17101, 66007, 63718, 60305, 57038, 53046, 48788,
      47071, 51756, 19791,
   };
   redirect_stderr_to(error_buffer);
   mnmultinom_prob(&m, &n, &r, yz, &t, &a, p, q, index, &out, pem);
   block_fwdb(&m, &nblocks, size, Q, init, pem, phi, trans, &ll);
   unredirect_sderr();

   free(yz);
   free(index);
   free(pem);
   free(phi);

}

void
test_params
(void)
{

   const int m = 3;
   const int r = 3;

   double Q[9] = {
      0.970078, 0.007204, 0.025482,
      0.010285, 0.961420, 0.017606,
      0.019637, 0.031376, 0.956912,
   };
   double p[12] = {
      .0548, .8062, .1304, .0086,
      .0459, .6758, .1569, .1214,
      .0399, .5865, .2199, .1537,
   };
   double q[12] = {
      .5331, .3711, .0611, .0347,
      .5387, .3750, .0245, .0618,
      .4807, .3345, .0933, .0915,
   };

   // Test `params_new`.
   params *par = params_new(m, r);
   params *try = params_new(m, r);
   g_assert(par != NULL);
   g_assert(try != NULL);
   g_assert_cmpint(par->m, ==, m);
   g_assert_cmpint(par->r, ==, r);
   g_assert_cmpint(try->m, ==, m);
   g_assert_cmpint(try->r, ==, r);

   // Test `params_set`.
   params_set(par, 0.9172, 2.8440, Q, p, q);
   g_assert_cmpfloat(par->t, ==, 0.9172);
   g_assert_cmpfloat(par->a, ==, 2.8440);
   for (int i = 0 ; i < m*m ; i++) g_assert(par->Q[i] == Q[i]);
   for (int i = 0 ; i < m*(r+1) ; i++) {
      g_assert_cmpfloat(par->p[i], ==, p[i]);
      g_assert_cmpfloat(par->q[i], ==, q[i]);
   }

   // Test `params_cpy`.
   params_cpy(try, par);
   g_assert_cmpfloat(try->t, ==, 0.9172);
   g_assert_cmpfloat(try->a, ==, 2.8440);
   for (int i = 0 ; i < m*m ; i++) {
      g_assert_cmpfloat(try->Q[i], ==, Q[i]);
   }
   for (int i = 0 ; i < m*(r+1) ; i++) {
      g_assert_cmpfloat(try->p[i], ==, p[i]);
      g_assert_cmpfloat(try->q[i], ==, q[i]);
   }
   
   // Test `params_change`.
   for (int i = 0 ; i < 262144 ; i++) {
      params_change(try, par, i % 3);
      for (int j = 0 ; j < m ; j++) {
         double sumQ = 0.0;
         for (int k = 0 ; k < m ; k++) {
            sumQ += try->Q[j+k*m];
            g_assert_cmpfloat(try->Q[j+k*m], <, 1);
            g_assert_cmpfloat(try->Q[j+k*m], >, 0);
         }
         g_assert_cmpfloat(fabs(sumQ - 1.0), <, 1e-6);
         double sump = 0.0;
         double sumq = 0.0;
         for (int k = 0 ; k < r+1 ; k++) {
            sump += try->p[k+j*(r+1)];
            sumq += try->q[k+j*(r+1)];
            g_assert_cmpfloat(try->p[k+j*(r+1)], <, 1);
            g_assert_cmpfloat(try->p[k+j*(r+1)], >, 0);
            g_assert_cmpfloat(try->q[k+j*(r+1)], <, 1);
            g_assert_cmpfloat(try->q[k+j*(r+1)], >, 0);
         }
         g_assert_cmpfloat(fabs(sump - 1.0), <, 1e-6);
         g_assert_cmpfloat(fabs(sumq - 1.0), <, 1e-6);
      }
   }

   params_destroy(try);
   params_destroy(par);

   return;

}

void
test_loglik
(void)
{

   int m = 3;
   int n = 6;
   int r = 3;

   double t = 0.8;
   double a = 1.2;

   double Q[9] = {
      // transpose //
      0.7, 0.1, 0.9,
      0.2, 0.5, 0.1,
      0.1, 0.4, 0.0,
   };

   // 'C1' and 'C2' simplify the definition of 'p' and 'q'.
   double C1 = 0.7692308;
   double C2 = 2.5000000;
   double p[12] = {
      C1/(C1+3.0), 1.0/(C1+3.0), 1.0/(C1+3.0), 1.0/(C1+3.0), 
      C1/(C1+4.0), 1.0/(C1+4.0), 2.0/(C1+4.0), 1.0/(C1+4.0),
      C1/(C1+1.3), 1.0/(C1+1.3), 0.2/(C1+1.3), 0.1/(C1+1.3),
   };
   double q[12] = {
      C2/(C2+3.0),   1.0/(C2+3.0),   1.0/(C2+3.0),   1.0/(C2+3.0), 
      C2/(C2+4.0),   1.0/(C2+4.0),   2.0/(C2+4.0),   1.0/(C2+4.0),
      // This line is not properly normalized. It should
      // trigger a warning but not cause failure.
      C2/(C2+1.3)*2, 1.0/(C2+1.3)*2, 0.2/(C2+1.3)*2, 0.1/(C2+1.3)*2,
   };

   params *par = params_new(m, r);
   params_set(par, t, a, Q, p, q);

   int yz[18] = {
         1, 2, 2,
         0, 4, 2,
         1, 2, 2,
         1, 2, 0,
         // Underflow (should give 0.0).
      1500, 1, 2,
         // Negative values (should give NA).
        -1, 4, 2,
   };

   int index[6] = {-1,-1,-1,-1,-1,-1};
   double pem[18];
   indexts(n, r, (const int *) yz, index);

   redirect_stderr_to(error_buffer);
   double ll = loglik(n, par, yz, index, pem);
   unredirect_sderr();

   double expected_ll = -1135.57472829;
   double expected_alpha[18] = {
      // Checked manually.
          0.502683584781,  0.489600586793,  0.007715828424,
          0.278067313112,  0.721741573153,  0.000191113733,
          0.394077474900,  0.598752057680,  0.007170467419,
          0.322075687622,  0.561610874783,  0.116313437595,
          0.000000000000,  0.000000000000,  1.000000000000,
          0.900000000000,  0.100000000000,  0.000000000000,
   };   

   char expected_warning[] = "fill the buffer";

   g_assert_cmpfloat(abs(ll - expected_ll), <, 1e-6);

   for (int i = 0 ; i < 18 ; i++) {
      g_assert_cmpfloat(fabs(expected_alpha[i] - pem[i]), <, 1e-6);
   }

   // Test the warning message.
   g_assert_cmpstr(error_buffer, ==, expected_warning);

   params_destroy(par);

   return;
}
#endif
