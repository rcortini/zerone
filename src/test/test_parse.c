#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "parse.c"


void
test_bloom
(void)
{

   bloom_t bloom = calloc(1, BSIZE);
   if (bloom == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   test_assert(bloom_query_and_set("text", 33, bloom) == 0);
   test_assert(bloom_query_and_set("text", 33, bloom) == 1);

   test_assert(bloom_query_and_set("aaa", 159753, bloom) == 0);
   test_assert(bloom_query_and_set("aaa", 159753, bloom) == 1);

   free(bloom);

}


void
test_add_to_rod
(void)
{
   rod_t *rod = malloc(sizeof(rod_t) + 32*sizeof(uint32_t));
   if (rod == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   rod->sz= 32;
   rod->mx = 0;
   memset(rod->array, 0, 32*sizeof(uint32_t));

   test_assert(add_to_rod(&rod, 0));
   test_assert(rod->mx == 0);
   test_assert(rod->array[0] == 1);

   test_assert(add_to_rod(&rod, 1));
   test_assert(rod->mx == 1);
   test_assert(rod->array[1] == 1);

   test_assert(add_to_rod(&rod, 32));
   test_assert_critical(rod->sz == 64);
   test_assert(rod->mx == 32);
   test_assert(rod->array[32] == 1);

   test_assert(add_to_rod(&rod, 1));
   test_assert(rod->mx == 32);
   test_assert(rod->array[1] == 2);

   test_assert(add_to_rod(&rod, 63));
   test_assert(rod->mx == 63);
   test_assert(rod->array[63] == 1);

   test_assert(add_to_rod(&rod, 65));
   test_assert_critical(rod->sz == 130);
   test_assert(rod->mx == 65);
   test_assert(rod->array[65] == 1);

   free(rod);

}


void
test_djb2
(void)
{
   test_assert(djb2("") == 5381);
   test_assert(djb2("\1") == 177572);
}


void
test_lookup_or_insert
(void)
{

   hash_t * hashtab = calloc(HSIZE, sizeof(link_t *));
   if (hashtab == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;   
   }

   link_t *lnk = NULL;
   lnk = lookup_or_insert("chr1", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(strcmp(lnk->seqname, "chr1") == 0);
   test_assert(lnk->next == NULL);
   // Test that counts have been allocated and initialized.
   test_assert_critical(lnk->counts != NULL);
   test_assert(lnk->counts->sz == 32);
   test_assert(lnk->counts->mx == 0);

   // Test that the pointer is the same.
   test_assert_critical(lnk == lookup_or_insert("chr1", hashtab));

   test_assert(add_to_rod(&lnk->counts, 157));
   test_assert(lnk->counts->sz == 314);
   test_assert(lnk->counts->mx == 157);

   lnk = lookup_or_insert("chr2", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(strcmp(lnk->seqname, "chr2") == 0);
   test_assert(lnk->next == NULL);
   test_assert_critical(lnk->counts != NULL);
   test_assert(lnk->counts->sz == 32);
   test_assert(lnk->counts->mx == 0);

   test_assert(add_to_rod(&lnk->counts, 39));
   test_assert(lnk->counts->sz == 78);
   test_assert(lnk->counts->mx == 39);

   // Test systematic update.
   for (int i = 0 ; i < 10000 ; i++) {
      lnk = lookup_or_insert("chr2", hashtab);;
      test_assert_critical(lnk != NULL);
      test_assert(add_to_rod(&lnk->counts, i));
   }

   test_assert(lnk->counts->sz == 19968);
   test_assert(lnk->counts->mx == 9999);

   destroy_hash(hashtab);

}

void
test_stress_hash
(void)
{

   hash_t * hashtab = calloc(HSIZE, sizeof(link_t *));
   if (hashtab == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;   
   }

   srand(123);

   // Insert way too many items in the hash table (there are
   // HSIZE = 997 buckets). Make sure insertion works until
   // the end.
   for (int i = 0 ; i < 100000 ; i++) {
      char s[32] = {0};
      for (int j = 0 ; j < 32 ; j++) s[j] = rand();
      test_assert_critical(lookup_or_insert(s, hashtab));
   }

   destroy_hash(hashtab);

}

void
test_merge_hashes
(void)
{

   hash_t * hashes[2] = {0};
   hashes[0] = calloc(HSIZE, sizeof(link_t *));
   hashes[1] = calloc(HSIZE, sizeof(link_t *));
   if (hashes[0] == NULL || hashes[1] == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;   
   }

   link_t *lnk = NULL;

   lnk = lookup_or_insert("chr1", hashes[0]);
   test_assert_critical(lnk != NULL);
   test_assert(add_to_rod(&lnk->counts, 5));
   lnk = lookup_or_insert("chr2", hashes[0]);
   test_assert_critical(lnk != NULL);
   test_assert(add_to_rod(&lnk->counts, 76));

   lnk = lookup_or_insert("chr1", hashes[1]);
   test_assert_critical(lnk != NULL);
   test_assert(add_to_rod(&lnk->counts, 81));
   lnk = lookup_or_insert("chr3", hashes[1]);
   test_assert_critical(lnk != NULL);
   test_assert(add_to_rod(&lnk->counts, 8));

   ChIP_t *ChIP = merge_hashes(hashes, 2);
   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 2);
   test_assert(ChIP->nb == 3);
   test_assert(ChIP->sz[0] == 81);
   test_assert(ChIP->sz[1] == 76);
   test_assert(ChIP->sz[2] == 8);

   destroy_hash(hashes[0]);
   destroy_hash(hashes[1]);

   // Manually destroy 'ChIP'.
   free(ChIP->y);
   free(ChIP);

}


void
test_choose_iterator
(void)
{

   generic_state_t * g_state = NULL;
   bgzf_state_t    * b_state = NULL;
   iter_t            iterate = NULL;
   
   iterate = choose_iterator("test_file_good.map");
   test_assert(iterate == generic_iterator);

   g_state = (generic_state_t *) STATE;
   test_assert_critical(g_state != NULL);
   test_assert(g_state->parser == parse_gem);
   test_assert(g_state->reader == getline);
   test_assert(g_state->file != NULL);
   test_assert(g_state->buff != NULL);
   test_assert(g_state->bsz == 32);

   // Clean.
   iterate(NULL);
   test_assert(STATE == NULL);

   iterate = choose_iterator("test_file_good.map.gz");
   test_assert(iterate == generic_iterator);

   g_state = (generic_state_t *) STATE;
   test_assert_critical(g_state != NULL);
   test_assert(g_state->parser == parse_gem);
   test_assert(g_state->reader == getgzline);
   test_assert(g_state->file != NULL);
   test_assert(g_state->buff != NULL);
   test_assert(g_state->bsz == 32);

   // Clean.
   iterate(NULL);
   test_assert(STATE == NULL);

   iterate = choose_iterator("test_file_good.bam");
   test_assert(iterate == bgzf_iterator);

   b_state = (bgzf_state_t *) STATE;
   test_assert_critical(b_state != NULL);
   test_assert(b_state->file != NULL);
   test_assert(b_state->hdr != NULL);
   test_assert(b_state->bam != NULL);

   // Clean.
   iterate(NULL);
   test_assert(STATE == NULL);

   redirect_stderr();
   iterate = choose_iterator("no_such_file.map");
   unredirect_stderr();
   test_assert(iterate == NULL);
   test_assert_stderr("cannot open file no_such_file.map\n");
   test_assert(STATE == NULL);

   return;

}


void
test_parse_gem
(void)
{

      loc_t loc;

      // NOTE that we cannot pass constant strings to
      // 'parse_gem()' because it modifies them in general.
      char line1[] = "a\tb\tc\t0:0:0\t-";
      test_assert(parse_gem(&loc, line1));
      test_assert(loc.name == NULL);
      test_assert(loc.pos == 0);

      char line2[] = "a\tb\tc\t0\tchr18:-:16507402:A35";
      test_assert(parse_gem(&loc, line2));
      test_assert(strcmp(loc.name, "chr18") == 0);
      test_assert(loc.pos == 16507402);

      char line3[] = "abc0chr18:-:16507402:A35";
      test_assert(!parse_gem(&loc, line3));

      char line4[] = "a\tb\tc\t0\tchr1816507402:A35";
      test_assert(!parse_gem(&loc, line4));

      char line5[] = "a\tb\tc\t0\tchr18:-:wrong:A35";
      test_assert(!parse_gem(&loc, line5));

      return;

}


void
test_parse_sam
(void)
{

      loc_t loc;

      // NOTE that we cannot pass constant strings to
      // 'parse_sam()' because it modifies them in general.
      char line1[] = "a\tb\tchr1\t456";
      test_assert(parse_sam(&loc, line1));
      test_assert(strcmp(loc.name, "chr1") == 0);
      test_assert(loc.pos == 456);

      char line2[] = "a\tb\tchr2\t345";
      test_assert(parse_sam(&loc, line2));
      test_assert(strcmp(loc.name, "chr2") == 0);
      test_assert(loc.pos == 345);

      char line3[] = "abc0chr18:-:16507402:A35";
      test_assert(!parse_sam(&loc, line3));

      char line4[] = "a\tb\tc\twrong";
      test_assert(!parse_sam(&loc, line4));

      char line5[] = "a\tb\tc\t0\t...";
      test_assert(!parse_sam(&loc, line5));

      return;

}


void
test_parse_bed
(void)
{

      loc_t loc;

      // NOTE that we cannot pass constant strings to
      // 'parse_sam()' because it modifies them in general.
      char line1[] = "chr1\t1\t3";
      test_assert(parse_bed(&loc, line1));
      test_assert(strcmp(loc.name, "chr1") == 0);
      test_assert(loc.pos == 2);

      char line2[] = "chr2\t345\t345";
      test_assert(parse_bed(&loc, line2));
      test_assert(strcmp(loc.name, "chr2") == 0);
      test_assert(loc.pos == 345);

      char line3[] = "abc0chr18:-:16507402:A35";
      test_assert(!parse_bed(&loc, line3));

      char line4[] = "a\twrong\t45";
      test_assert(!parse_bed(&loc, line4));

      char line5[] = "a\t45\twrong";
      test_assert(!parse_bed(&loc, line5));

      return;

}

void
test_autoparse
(void)
{

   hash_t *hashtab = NULL;
   link_t *lnk = NULL;

   hashtab = calloc(HSIZE, sizeof(link_t *));

   if (hashtab == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }
      
   test_assert(autoparse("test_file_good.map", hashtab));
   test_assert_critical(hashtab != NULL);

   lnk = lookup_or_insert("chr6", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[317504] == 1);

   lnk = lookup_or_insert("chr18", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[55024] == 1);

   lnk = lookup_or_insert("chrX", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[445385] == 1);

   // Try parsing a non gem file.
   redirect_stderr();
   test_assert(!autoparse("test_parse.c", hashtab));
   unredirect_stderr();
   test_assert(hashtab != NULL);
   test_assert(strncmp("unknown format for file test_parse.c",
            caught_in_stderr(), 25) == 0);

   // Try parsing bad gem files.
   redirect_stderr();
   test_assert(!autoparse("test_file_bad1.map", hashtab));
   unredirect_stderr();
   test_assert(hashtab != NULL);
   test_assert(strncmp("format conflict in line:\n",
            caught_in_stderr(), 25) == 0);

   redirect_stderr();
   test_assert(!autoparse("test_file_bad2.map", hashtab));
   unredirect_stderr();
   test_assert(hashtab != NULL);
   test_assert(strncmp("format conflict in line:\n",
            caught_in_stderr(), 25) == 0);

   // Now test .bam format. Note that there is a sequence
   // at position 0 in this file, so the parser will choke
   // on it. But there are 8 sequences before, that should
   // fill the hash as tested below.
   test_assert(!autoparse("test_file_good.bam", hashtab));
   test_assert_critical(hashtab != NULL);

   lnk = lookup_or_insert("insert", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[0] == 2);

   lnk = lookup_or_insert("ref1", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[0] == 5);

   destroy_hash(hashtab);

}


void
test_getgzline
(void)
{

   FILE *fp = fopen("test_file_bad1.map.gz", "r");

   if (fp == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   size_t bsz = 32;
   char *buff = malloc(bsz * sizeof(char));

   if (buff == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   test_assert_critical(getgzline(&buff, &bsz, fp) > 0);
   test_assert(strncmp((char *) buff,
            "SRR574805.4 ROCKFORD:3:1:1728:956/1", 35) == 0);

   test_assert_critical(getgzline(&buff, &bsz, fp) > 0);
   test_assert(strncmp((char *) buff,
            "SRR574805.6 ROCKFORD:3:1:2065:964/1", 35) == 0);

   test_assert(getgzline(&buff, &bsz, fp) == -1);

   fclose(fp);
   free(buff);

}

void
test_getgzline_err
(void)
{

   FILE *fp = NULL;
      
   fp = fopen("test_file_bad1.map", "r");

   if (fp == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   size_t bsz = 32;
   char *buff = malloc(bsz * sizeof(char));

   if (buff == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   test_assert(getgzline(&buff, &bsz, fp) == -2);

   fclose(fp);

   // Check that the error did not mess up
   // the internal state of 'getgzline'.
   fp = fopen("test_file_bad1.map.gz", "r");

   if (fp == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   test_assert_critical(getgzline(&buff, &bsz, fp) > 0);
   test_assert(strncmp((char *) buff,
            "SRR574805.4 ROCKFORD:3:1:1728:956/1", 35) == 0);

   // Interrupt inflation (to prevent memory leak).
   test_assert_critical(getgzline(NULL, &bsz, fp) == -1);

   fclose(fp);

   fp = fopen("test_file_corrupt.map.gz", "r");

   if (fp == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   test_assert(getgzline(&buff, &bsz, fp) == -2);

   fclose(fp);
   free(buff);

}

void
test_parse_input_files
(void)
{

   ChIP_t *ChIP = NULL;

   char *mock_fnames_1[] = { "test_file_good.map", NULL };
   char *ChIP_fnames_1[] = { "test_file_good.map", NULL };

   ChIP = parse_input_files(mock_fnames_1, ChIP_fnames_1);
   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 2);
   test_assert(ChIP->nb == 3);
   test_assert(ChIP->sz[0] == 445385);
   test_assert(ChIP->sz[1] == 317504);
   test_assert(ChIP->sz[2] == 55024);

   // Manually destroy 'ChIP'.
   free(ChIP->y);
   free(ChIP);
   ChIP = NULL;

   // Same test with gzipped files.
   char *mock_fnames_2[] = { "test_file_good.map.gz", NULL };
   char *ChIP_fnames_2[] = { "test_file_good.map.gz", NULL };

   ChIP = parse_input_files(mock_fnames_2, ChIP_fnames_2);
   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 2);
   test_assert(ChIP->nb == 3);
   test_assert(ChIP->sz[0] == 445385);
   test_assert(ChIP->sz[1] == 317504);
   test_assert(ChIP->sz[2] == 55024);

   // Manually destroy 'ChIP'.
   free(ChIP->y);
   free(ChIP);

}
