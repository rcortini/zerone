#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "parse.c"

void
test_add
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

   test_assert(add(&rod, 0));
   test_assert(rod->mx == 0);
   test_assert(rod->array[0] == 1);

   test_assert(add(&rod, 1));
   test_assert(rod->mx == 1);
   test_assert(rod->array[1] == 1);

   test_assert(add(&rod, 32));
   test_assert_critical(rod->sz == 64);
   test_assert(rod->mx == 32);
   test_assert(rod->array[32] == 1);

   test_assert(add(&rod, 1));
   test_assert(rod->mx == 32);
   test_assert(rod->array[1] == 2);

   test_assert(add(&rod, 63));
   test_assert(rod->mx == 63);
   test_assert(rod->array[63] == 1);

   test_assert(add(&rod, 65));
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

   test_assert(add(&lnk->counts, 157));
   test_assert(lnk->counts->sz == 314);
   test_assert(lnk->counts->mx == 157);

   lnk = lookup_or_insert("chr2", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(strcmp(lnk->seqname, "chr2") == 0);
   test_assert(lnk->next == NULL);
   test_assert_critical(lnk->counts != NULL);
   test_assert(lnk->counts->sz == 32);
   test_assert(lnk->counts->mx == 0);

   test_assert(add(&lnk->counts, 39));
   test_assert(lnk->counts->sz == 78);
   test_assert(lnk->counts->mx == 39);

   // Test systematic update.
   for (int i = 0 ; i < 10000 ; i++) {
      lnk = lookup_or_insert("chr2", hashtab);;
      test_assert_critical(lnk != NULL);
      test_assert(add(&lnk->counts, i));
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
      char s[16] = {0};
      for (int j = 0 ; j < 16 ; j++) s[j] = rand();
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
   test_assert(add(&lnk->counts, 5));
   lnk = lookup_or_insert("chr2", hashes[0]);
   test_assert_critical(lnk != NULL);
   test_assert(add(&lnk->counts, 76));

   lnk = lookup_or_insert("chr1", hashes[1]);
   test_assert_critical(lnk != NULL);
   test_assert(add(&lnk->counts, 81));
   lnk = lookup_or_insert("chr3", hashes[1]);
   test_assert_critical(lnk != NULL);
   test_assert(add(&lnk->counts, 8));

   ChIP_t *ChIP = merge_hashes(hashes, 2);
   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 2);
   test_assert(ChIP->nb == 3);
   test_assert(ChIP->size[0] == 81);
   test_assert(ChIP->size[1] == 76);
   test_assert(ChIP->size[2] == 8);

   destroy_hash(hashes[0]);
   destroy_hash(hashes[1]);

   // Manually destroy 'ChIP'.
   free(ChIP->y);
   free(ChIP);

}

void
test_read_and_count
(void)
{

   hash_t *hashtab = NULL;
   link_t *lnk = NULL;
      
   hashtab = read_and_count("test_file_good.map");
   test_assert_critical(hashtab != NULL);

   lnk = lookup_or_insert("chr6", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[317504] == 2);

   lnk = lookup_or_insert("chr18", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[55024] == 1);

   lnk = lookup_or_insert("chrX", hashtab);
   test_assert_critical(lnk != NULL);
   test_assert(lnk->counts->array[445385] == 1);

   destroy_hash(hashtab);
   hashtab = NULL;

   // Try parsing a non gem file.
   redirect_stderr();
   hashtab = read_and_count("test_parse.c");
   unredirect_stderr();
   test_assert(hashtab == NULL);
   test_assert(strncmp("format conflict in line:\n",
            caught_in_stderr(), 25) == 0);

   // Try parsing bad gem files.
   redirect_stderr();
   hashtab = read_and_count("test_file_bad1.map");
   unredirect_stderr();
   test_assert(hashtab == NULL);
   test_assert(strncmp("format conflict in line:\n",
            caught_in_stderr(), 25) == 0);

   redirect_stderr();
   hashtab = read_and_count("test_file_bad2.map");
   unredirect_stderr();
   test_assert(hashtab == NULL);
   test_assert(strncmp("format conflict in line:\n",
            caught_in_stderr(), 25) == 0);

}

void
test_readgzline
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

   test_assert_critical(readgzline(&buff, &bsz, fp) > 0);
   test_assert(strncmp((char *) buff,
            "SRR574805.4 ROCKFORD:3:1:1728:956/1", 35) == 0);

   test_assert_critical(readgzline(&buff, &bsz, fp) > 0);
   test_assert(strncmp((char *) buff,
            "SRR574805.6 ROCKFORD:3:1:2065:964/1", 35) == 0);

   test_assert(readgzline(&buff, &bsz, fp) == -1);

   fclose(fp);
   free(buff);

}

void
test_readgzline_err
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

   test_assert(readgzline(&buff, &bsz, fp) == -2);

   fclose(fp);

   // Check that the error did not mess up
   // the internal state of 'readgzline'.
   fp = fopen("test_file_bad1.map.gz", "r");

   if (fp == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   test_assert_critical(readgzline(&buff, &bsz, fp) > 0);
   test_assert(strncmp((char *) buff,
            "SRR574805.4 ROCKFORD:3:1:1728:956/1", 35) == 0);

   // Interrupt inflation (to prevent memory leak).
   test_assert_critical(readgzline(NULL, &bsz, fp) == -1);

   fclose(fp);

   fp = fopen("test_file_corrupt.map.gz", "r");

   if (fp == NULL) {
      fprintf(stderr, "error in test function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return;
   }

   test_assert(readgzline(&buff, &bsz, fp) == -2);

   fclose(fp);
   free(buff);

}


void
test_read_gem
(void)
{

   ChIP_t *ChIP = NULL;

   const char *fnames[] = {"test_file_good.map", "test_file_good.map"};

   ChIP = read_gem(fnames, 2);
   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 2);
   test_assert(ChIP->nb == 3);
   test_assert(ChIP->size[0] == 445385);
   test_assert(ChIP->size[1] == 317504);
   test_assert(ChIP->size[2] == 55024);

   // Manually destroy 'ChIP'.
   free(ChIP->y);
   free(ChIP);
   ChIP = NULL;

   // Same test with gzipped files.
   const char *fnames_gz[] = {"test_file_good.map.gz",
      "test_file_good.map.gz"};

   ChIP = read_gem(fnames_gz, 2);
   test_assert_critical(ChIP != NULL);
   test_assert(ChIP->r == 2);
   test_assert(ChIP->nb == 3);
   test_assert(ChIP->size[0] == 445385);
   test_assert(ChIP->size[1] == 317504);
   test_assert(ChIP->size[2] == 55024);

   // Manually destroy 'ChIP'.
   free(ChIP->y);
   free(ChIP);

}
