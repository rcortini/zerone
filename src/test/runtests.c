#include "unittest.h"

// Test cases form tests_hmm.c
void test_fwdb (void);
void test_fwdb_NA (void);
void test_underflow (void);
void test_block_fwdb (void);
void test_block_fwdb_NA (void);
void test_viterbi (void);
void test_block_viterbi (void);
void test_block_viterbi_NA (void);
void test_update_trans (void);
void test_bw_zinm (void);


// Test cases from tests_utils.c
void test_new_histo (void);
void test_histo_push (void);
void test_compress_histo (void);
void test_tabulate (void);
void test_indexts (void);


// Test cases from tests_zinb.c
void test_zinm_prob (void);
void test_eval_nb_f (void);
void test_eval_nb_dfda (void);
void test_eval_zinb_f (void);
void test_eval_zinb_g (void);
void test_eval_zinb_dfda (void);
void test_eval_zinb_dfdp (void);
void test_eval_zinb_dgda (void);
void test_ll_zinb (void);
void test_nb_est_alpha (void);
void test_mle_nb (void);
void test_mle_zinb (void);
void test_fail_mle_nb (void);
void test_fail_mle_zinb (void);
void test_err_handler (void);


// Test cases from tests_parse.c
void test_bloom (void);
void test_add_to_rod (void);
void test_djb2 (void);
void test_lookup_or_insert (void);
void test_stress_hash (void);
void test_merge_hashes (void);
void test_choose_iterator (void);
void test_parse_gem (void);
void test_parse_sam (void);
void test_parse_bed (void);
void test_parse_wig (void);
void test_autoparse (void);
void test_getgzline (void);
void test_getgzline_err (void);
void test_parse_input_files (void);


// Test cases from test_predict.c
void test_predict (void);


// test cases from test_zerone.c
void test_eval_bw_f (void);


int
main(
   int argc,
   char **argv
)
{

   // Register test cases //
   const static test_case_t test_cases[] = {
      {"utils/new_histo",         test_new_histo},
      {"utils/histo_push",        test_histo_push},
      {"utils/compress_histo",    test_compress_histo},
      {"utils/tabulate",          test_tabulate},
      {"utils/indexts",           test_indexts},
      {"zinm/zinm_prob",          test_zinm_prob},
      {"zinm/eval_nb_f",          test_eval_nb_f},
      {"zinm/eval_nb_dfda",       test_eval_nb_dfda},
      {"zinm/eval_zinb_f",        test_eval_zinb_f},
      {"zinm/eval_zinb_g",        test_eval_zinb_g},
      {"zinm/eval_zinb_dfda",     test_eval_zinb_dfda},
      {"zinm/eval_zinb_dfdp",     test_eval_zinb_dfdp},
      {"zinm/eval_zinb_dgda",     test_eval_zinb_dgda},
      {"zinm/ll_zinb",            test_ll_zinb},
      {"zinm/nb_est_alpha",       test_nb_est_alpha},
      {"zinm/mle_nb",             test_mle_nb},
      {"zinm/mle_zinb",           test_mle_zinb},
      {"zinm/fail_mle_nb",        test_fail_mle_nb},
      {"zinm/fail_mle_zinb",      test_fail_mle_zinb},
      {"zinm/err_handler",        test_err_handler},
      {"hmm/fwdb",                test_fwdb},
      {"hmm/fwdb (NAs)",          test_fwdb_NA},
      {"hmm/fwdb (underflow)",    test_underflow},
      {"hmm/block_fwdb",          test_block_fwdb},
      {"hmm/block_fwdb (NAs)",    test_block_fwdb_NA},
      {"hmm/viterbi",             test_viterbi},
      {"hmm/block_viterbi",       test_block_viterbi},
      {"hmm/block_viterbi (NAs)", test_block_viterbi_NA},
      {"hmm/update_trans",        test_update_trans},
      {"hmm/bw_zinm",             test_bw_zinm},
      {"parse/bloom",             test_bloom},
      {"parse/add_to_rod",        test_add_to_rod},
      {"parse/djb2",              test_djb2},
      {"parse/lookup_or_insert",  test_lookup_or_insert},
      {"parse/stress_hash",       test_stress_hash},
      {"parse/merge_hashes",      test_merge_hashes},
      {"parse/choose_iterator",   test_choose_iterator},
      {"parse/parse_gem",         test_parse_gem},
      {"parse/parse_sam",         test_parse_sam},
      {"parse/parse_bed",         test_parse_bed},
      {"parse/parse_wig",         test_parse_wig},
      {"parse/autoparse",         test_autoparse},
      {"parse/getgzline",         test_getgzline},
      {"parse/getgzline_err",     test_getgzline_err},
      {"parse/parse_input_files", test_parse_input_files},
      {"predict/predict",         test_predict},
      {"zerone/eval_bw_f",        test_eval_bw_f},
      {NULL, NULL}
   };

   // Run the tests.
   return run_unittest(argc, argv, test_cases);

}
