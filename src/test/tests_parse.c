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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "unittest.h"
#include "parse.c"


void
test_bitf
(void)
{

   hash_t *hashtab = calloc(HSIZE, sizeof(link_t *));
   link_t *lnk = NULL;

   // Insert an entry for chromosome 1.
   lnk = lookup_or_insert("chr1", hashtab);
   test_assert_critical(lnk != NULL);

   // Add a read at position 123 and check it is added.
   test_assert(bitf_query_and_set(123, lnk) == 0);
   test_assert(bitf_query_and_set(123, lnk) == 1);

   // Add a read at position 456 and check it is added.
   test_assert(bitf_query_and_set(456, lnk) == 0);
   test_assert(bitf_query_and_set(456, lnk) == 1);

   // Insert an entry for chromosome 2.
   lnk = lookup_or_insert("chr2", hashtab);
   test_assert_critical(lnk != NULL);

   // Add a read at position 789 and check it is added.
   test_assert(bitf_query_and_set(789, lnk) == 0);
   test_assert(bitf_query_and_set(789, lnk) == 1);

   // Add a read at position 159 and check it is added.
   test_assert(bitf_query_and_set(159, lnk) == 0);
   test_assert(bitf_query_and_set(159, lnk) == 1);

   // Query the entry for chromosome 1.
   lnk = lookup_or_insert("chr1", hashtab);
   test_assert_critical(lnk != NULL);

   // Check that the reads are still there.
   test_assert(bitf_query_and_set(123, lnk) == 1);
   test_assert(bitf_query_and_set(456, lnk) == 1);

   // Clean.
   destroy_bitfields(hashtab);
   destroy_hash(hashtab);

}


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
   bzero(rod->array, 32*sizeof(uint32_t));

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

   // Merge hashes assuming using the first as mock.
   // Note that this modifies the first hash.
   ChIP_t *ChIP = merge_hashes(hashes, 2, 0);
   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 2);
   test_assert(ChIP->nb == 3);
   test_assert(ChIP->sz[0] == 82);
   test_assert(ChIP->sz[1] == 77);
   test_assert(ChIP->sz[2] == 9);

   // Manually destroy 'ChIP'.
   free(ChIP->y);
   free(ChIP);

   // Further update hash.
   lnk = lookup_or_insert("chr17", hashes[1]);
   test_assert_critical(lnk != NULL);
   test_assert(add_to_rod(&lnk->counts, 5));

   // Merge hash without using the first as mock.
   ChIP = merge_hashes(hashes, 2, 1);
   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 2);
   test_assert(ChIP->nb == 4);
   test_assert(ChIP->sz[0] == 82);
   test_assert(ChIP->sz[1] == 77);
   test_assert(ChIP->sz[2] == 9);
   test_assert(ChIP->sz[3] == 6);
   // This should fill the arrays with 1s.
   for (int i = 0 ; i < nobs(ChIP) ; i++) {
      test_assert(ChIP->y[0+i*ChIP->r] == 1);
   }

   // Manually destroy 'ChIP'.
   free(ChIP->y);
   free(ChIP);

   destroy_hash(hashes[0]);
   destroy_hash(hashes[1]);

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
      char line1[] = "a\tb\tchr1\t456\t60";
      test_assert(parse_sam(&loc, line1));
      test_assert(strcmp(loc.name, "chr1") == 0);
      test_assert(loc.pos == 456);

      char line2[] = "a\tb\tchr2\t345\t60";
      test_assert(parse_sam(&loc, line2));
      test_assert(strcmp(loc.name, "chr2") == 0);
      test_assert(loc.pos == 345);

      char line3[] = "abc0chr18:-:16507402:A35";
      test_assert(!parse_sam(&loc, line3));

      char line4[] = "a\tb\tc\twrong\twrong";
      test_assert(!parse_sam(&loc, line4));

      char line5[] = "a\tb\tc\t0\t...";
      test_assert(!parse_sam(&loc, line5));

      char line6[] = "a\tb\t*\t0\t60";
      test_assert(parse_sam(&loc, line6));
      test_assert(loc.name == NULL);

      char line7[] = "#\tb\tc\t1\t60";
      test_assert(parse_sam(&loc, line7));
      test_assert(strcmp(loc.name, "c") == 0);

      return;

}


void
test_parse_bed
(void)
{

      loc_t loc;

      // NOTE that we cannot pass constant strings to
      // 'parse_bed()' because it modifies them in general.
      char line1[] = "chr1\t1\t3";
      test_assert(parse_bed(&loc, line1));
      test_assert(strcmp(loc.name, "chr1") == 0);
      test_assert(loc.pos == 3);

      char line2[] = "chr2\t345\t345";
      test_assert(parse_bed(&loc, line2));
      test_assert(strcmp(loc.name, "chr2") == 0);
      test_assert(loc.pos == 346);

      char line3[] = "abc0chr18:-:16507402:A35";
      test_assert(!parse_bed(&loc, line3));

      char line4[] = "a\twrong\t45";
      test_assert(!parse_bed(&loc, line4));

      char line5[] = "a\t45\twrong";
      test_assert(!parse_bed(&loc, line5));

      char line6[] = "a\t45\t";
      test_assert(!parse_bed(&loc, line6));

      char line7[] = "a\t45";
      test_assert(!parse_bed(&loc, line7));

      char line8[] = "a\t45\123y";
      test_assert(!parse_bed(&loc, line8));

      char line9[] = "a\t\t45";
      test_assert(!parse_bed(&loc, line9));

      char line10[] = "a\t123y\t45";
      test_assert(!parse_bed(&loc, line10));

      char line11[] = "a\t123\t45\n";
      test_assert(parse_bed(&loc, line11));

      char line12[] = "a\t123 \t45";
      test_assert(parse_bed(&loc, line12));

      return;

}


void
test_parse_wig
(void)
{

      loc_t loc = {0};

      char line1[] = "1";
      test_assert(!parse_wig(&loc, line1));

      char line2[] = "track type=wiggle_0";
      test_assert(parse_wig(&loc, line2));

      char line3[] = "variableStep\tchrom=chr1";
      test_assert(parse_wig(&loc, line3));
      test_assert(strcmp(loc.name, "chr1") == 0);

      char line4[] = "1";
      test_assert(!parse_wig(&loc, line4));

      char line5[] = "1\twrong";
      test_assert(!parse_wig(&loc, line5));

      char line6[] = "wrong\t0";
      test_assert(!parse_wig(&loc, line6));

      char line7[] = "1\t0";
      test_assert(parse_wig(&loc, line7));
      test_assert(strcmp(loc.name, "chr1") == 0);
      test_assert(loc.pos == 1);

      char line8[] = "variableStep\tchrom=chr2\tspan=wrong";
      test_assert(!parse_wig(&loc, line8));
      test_assert(strcmp(loc.name, "chr2") == 0);

      char line9[] = "variableStep\tchrom=chr2\tspan=3";
      test_assert(parse_wig(&loc, line9));
      test_assert(strcmp(loc.name, "chr2") == 0);

      char line10[] = "1\t10";
      test_assert(parse_wig(&loc, line10));
      test_assert(strcmp(loc.name, "chr2") == 0);
      test_assert(loc.pos == 2);

      char line11[] = "fixedStep\tchrom=chr1\tstart=0\tstep=300";
      test_assert(!parse_wig(&loc, line11));
      test_assert(strcmp(loc.name, "chr1") == 0);

      char line12[] = "fixedStep\tchrom=chr1\tstart=1\tstep=wrong";
      test_assert(!parse_wig(&loc, line12));
      test_assert(strcmp(loc.name, "chr1") == 0);

      char line13[] = "fixedStep\tchrom=chr1\tstart=1\tstep=300";
      test_assert(parse_wig(&loc, line13));
      test_assert(strcmp(loc.name, "chr1") == 0);

      char line14[] = "10";
      test_assert(parse_wig(&loc, line14));
      test_assert(strcmp(loc.name, "chr1") == 0);
      test_assert(loc.pos == 1);

      char line15[] = "10";
      test_assert(parse_wig(&loc, line15));
      test_assert(strcmp(loc.name, "chr1") == 0);
      test_assert(loc.pos == 301);

      char line16[] = "fixedStep\tchrom=chr2\tstart=1\tstep=300\tspan=-1";
      test_assert(!parse_wig(&loc, line16));
      test_assert(strcmp(loc.name, "chr2") == 0);

      char line17[] = "fixedStep\tchrom=chr2\tstart=1\tstep=300\tspan=300";
      test_assert(parse_wig(&loc, line17));
      test_assert(strcmp(loc.name, "chr2") == 0);

      char line18[] = "10";
      test_assert(parse_wig(&loc, line18));
      test_assert(strcmp(loc.name, "chr2") == 0);
      test_assert(loc.pos == 150);

      char line19[] = "10";
      test_assert(parse_wig(&loc, line19));
      test_assert(strcmp(loc.name, "chr2") == 0);
      test_assert(loc.pos == 450);

      char line20[] = "wrong";
      test_assert(!parse_wig(&loc, line20));

      char line21[] = "1\t0";
      test_assert(!parse_wig(&loc, line21));

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

   test_assert(autoparse("test_file_good.map", hashtab, 300));
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

   // Try a different window size.
   destroy_hash(hashtab);
   hashtab = calloc(HSIZE, sizeof(link_t *));

   if (hashtab == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   test_assert(autoparse("test_file_good.map", hashtab, 200));
   test_assert_critical(hashtab != NULL);

   lnk = lookup_or_insert("chr6", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[476256] == 1);

   lnk = lookup_or_insert("chr18", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[82537] == 1);

   lnk = lookup_or_insert("chrX", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[668077] == 1);

   // Try parsing a non gem file.
   redirect_stderr();
   test_assert(!autoparse("tests_parse.c", hashtab, 300));
   unredirect_stderr();
   test_assert(hashtab != NULL);
   test_assert(strncmp("unknown format for file tests_parse.c",
            caught_in_stderr(), 25) == 0);

   // Try parsing bad gem files.
   redirect_stderr();
   test_assert(!autoparse("test_file_bad1.map", hashtab, 300));
   unredirect_stderr();
   test_assert(hashtab != NULL);
   test_assert(strncmp("format conflict in line 2:\n",
            caught_in_stderr(), 27) == 0);

   redirect_stderr();
   test_assert(!autoparse("test_file_bad2.map", hashtab, 300));
   unredirect_stderr();
   test_assert(hashtab != NULL);
   test_assert(strncmp("format conflict in line 2:\n",
            caught_in_stderr(), 27) == 0);

   // Now test .bam format.
   test_assert(autoparse("test_file_good.bam", hashtab, 300));
   test_assert_critical(hashtab != NULL);

   lnk = lookup_or_insert("insert", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[0] == 1);
   test_assert(lnk->counts->array[1] == 1);

   lnk = lookup_or_insert("ref1", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[0] == 5);

   lnk = lookup_or_insert("ref2", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[0] == 7);

   lnk = lookup_or_insert("ref3", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[0] == 1);

   // Now test .bed format. Also check the error message.
   redirect_stderr();
   test_assert(!autoparse("test_file_bad1.bed", hashtab, 300));
   unredirect_stderr();
   test_assert(hashtab != NULL);
   test_assert(strcmp("format conflict in line 1:\nchr1\t34\twrong\n",
            caught_in_stderr()) == 0);

   redirect_stderr();
   test_assert(!autoparse("test_file_bad2.bed", hashtab, 300));
   unredirect_stderr();
   test_assert(hashtab != NULL);
   test_assert(strcmp("format conflict in line 1:\n"
               "a_very_long_chromosome_name\t1\t"
               "some_shit_that_will_cause_a_failu...\n",
            caught_in_stderr()) == 0);


   // Now test .sam format.
   destroy_hash(hashtab);
   hashtab = calloc(HSIZE, sizeof(link_t *));

   if (hashtab == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   // This is a paired-end sam file.
   // Half of the lines are ignored.
   test_assert(autoparse("test_file_good.sam", hashtab, 300));
   test_assert_critical(hashtab != NULL);

   lnk = lookup_or_insert("chrM", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[19] == 0); // Low quality.
   test_assert(lnk->counts->array[32] == 0); // Low quality.
   test_assert(lnk->counts->array[55] == 1);

   lnk = lookup_or_insert("chr1", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[830835] == 1);

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

   zerone_parser_args_t args;
   args.window = 300;
   args.minmapq = 20;

   ChIP = parse_input_files(mock_fnames_1, ChIP_fnames_1, args);
   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 2);
   test_assert(ChIP->nb == 3);
   test_assert(ChIP->sz[0] == 445386);
   test_assert(ChIP->sz[1] == 317505);
   test_assert(ChIP->sz[2] == 55025);

   // Manually destroy 'ChIP'.
   free(ChIP->y);
   free(ChIP);
   ChIP = NULL;

   // Try again with a different window size.
   args.window = 200;
   ChIP = parse_input_files(mock_fnames_1, ChIP_fnames_1, args);
   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 2);
   test_assert(ChIP->nb == 3);
   test_assert(ChIP->sz[0] == 668078);
   test_assert(ChIP->sz[1] == 476257);
   test_assert(ChIP->sz[2] == 82538);

   // Manually destroy 'ChIP'.
   free(ChIP->y);
   free(ChIP);
   ChIP = NULL;

   // Same test with gzipped files.
   char *mock_fnames_2[] = { "test_file_good.map.gz", NULL };
   char *ChIP_fnames_2[] = { "test_file_good.map.gz", NULL };

   args.window = 300;
   ChIP = parse_input_files(mock_fnames_2, ChIP_fnames_2, args);
   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 2);
   test_assert(ChIP->nb == 3);
   test_assert(ChIP->sz[0] == 445386);
   test_assert(ChIP->sz[1] == 317505);
   test_assert(ChIP->sz[2] == 55025);

   // Manually destroy 'ChIP'.
   free(ChIP->y);
   free(ChIP);

}

// Test cases for export.
const test_case_t test_cases_parse[] = {
   {"parse/bitf",              test_bitf},
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
   {NULL, NULL},
};

