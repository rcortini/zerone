#include "hmm.c"

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

   size_t m = 2;
   size_t n = 6;
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

   unsigned int nblocks = 2;
   unsigned int size[2] = {3,3};

   double l = block_fwdb(m, nblocks, size, Q, init, prob, phi, trans);

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

   test_assert(fabs(expected_loglik - l) < 1e-6);
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
   size_t m = 2;
   size_t n = 6;
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

   unsigned int nblocks = 2;
   unsigned int size[2] = {3,3};

   double l = block_fwdb(m, nblocks, size, Q, init, prob, phi, trans);

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

   test_assert(fabs(expected_loglik - l) < 1e-6);
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
      0.05, 0.05, 0.40, 0.45,
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
   if (path == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   int nblocks = 2;
   uint size[2] = {5,5};

   block_viterbi(
         m,
         nblocks,
         size,
         Q,
         init,
         prob,
         path
   );

   double expected_path[10] = {
      3,1,1,0,0,3,1,1,0,0,
   };

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path[i] == path[i]);
   }


   const double log_Q[16] = {
      // transpose //
      log(0.80), log(0.05), log(0.05), log(0.05),
      log(0.10), log(0.90), log(0.30), log(0.35),
      log(0.05), log(0.05), log(0.40), log(0.45),
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
         m,
         nblocks,
         size,
         log_Q,
         log_init,
         log_prob,
         path
   );

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path[i] == path[i]);
   }

   free(path);

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
      0.05, 0.05, 0.40, 0.45,
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
   uint size[2] = {5,5};

   block_viterbi(
         m,
         nblocks,
         size,
         Q,
         init,
         prob,
         path
   );

   double expected_path[10] = {
      3,1,1,0,0,3,0,0,0,0,
   };

   for (int i = 0 ; i < n ; i++) {
      test_assert(expected_path[i] == path[i]);
   }

   // Set invalid initial probabilities.
   static double invalid_init_1[4] = {-1.0, 0.10, 0.05, 0.80};

   // Catch stderr.
   redirect_stderr();
   int exit_status = block_viterbi(
         m,
         nblocks,
         size,
         Q,
         invalid_init_1,
         prob,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("mixed log/lin arguments in 'block_viterbi()'\n");

   static double invalid_init_2[4] = {.3, .3, .6, 0.0/0.0};

   // Catch stderr.
   redirect_stderr();
   exit_status = block_viterbi(
         m,
         nblocks,
         size,
         Q,
         invalid_init_2,
         prob,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("invalid 'init' argument in 'block_viterbi()'\n");

   static double invalid_init_3[4] = {.3, .3, .6, .9};

   // Catch stderr.
   redirect_stderr();
   exit_status = block_viterbi(
         m,
         nblocks,
         size,
         Q,
         invalid_init_3,
         prob,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr(
         "'init' is not a probability in 'block_viterbi()'\n");

   // Set invalid transition probabilities.
   static double invalid_Q_1[16] = {
      // transpose //
      -1.0, 0.05, 0.05, 0.05,
      0.10, 0.90, 0.30, 0.35,
      0.05, 0.05, 0.40, 0.45,
      0.05, 0.00, 0.25, 0.15,
   };

   // Catch stderr.
   redirect_stderr();
   exit_status = block_viterbi(
         m,
         nblocks,
         size,
         invalid_Q_1,
         init,
         prob,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("mixed log/lin arguments in 'block_viterbi()'\n");

   static double invalid_Q_2[16] = {
      // transpose //
      0.80, 0.05, 0.05, 0.05,
      0.10, 0.90, 0.30, 0.35,
      0.05, 0.05, 0.40, 0.45,
      0.05, 0.00, 0.25, 0.0/0.0,
   };

   // Catch stderr.
   redirect_stderr();
   exit_status = block_viterbi(
         m,
         nblocks,
         size,
         invalid_Q_2,
         init,
         prob,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("invalid 'Q' argument in 'block_viterbi()'\n");

   static double invalid_Q_3[16] = {
      // transpose //
      0.80, 0.05, 0.05, 0.05,
      0.10, 0.90, 0.30, 0.35,
      0.05, 0.05, 0.40, 0.45,
      0.05, 0.00, 0.25, 0.10,
   };

   // Catch stderr.
   redirect_stderr();
   exit_status = block_viterbi(
         m,
         nblocks,
         size,
         invalid_Q_3,
         init,
         prob,
         path
   );
   unredirect_stderr();

   test_assert(exit_status == -1);
   test_assert_stderr("'Q' is not stochastic in 'block_viterbi()'\n");

   free(path);

   return;

}
