#include <string.h>
#include "bgzf.h"
#include "sam.h"
#include "parsam.h"

#define SUCCESS 1
#define FAILURE 0

#define BIN_SIZE 300

struct counter_t {
   int32_t    n_ref;    // # of chromosomes
   int32_t  * n_bins;   // # of bins
   int32_t ** bins;     // read count per bin
};

typedef struct counter_t counter_t;

// Declaration of local functions.
counter_t * read_count (const char *);
void destroy_counter (counter_t *);


ChIP_t *
read_sam
(
            char * fn[],
   unsigned int    nfiles
)
{
   ChIP_t * ChIP = NULL;       // Return pointer.
   unsigned int * size = NULL; // Size of the blocks.
   int * y = NULL;             // Data counts.
   
   // Build the counter_t's from which to build the final ChIP_t.
   counter_t * counters[nfiles];
   for (int i = 0; i < nfiles; i++) {
      counters[i] = read_count(fn[i]);
      if (counters[i] == NULL) {
         fprintf(stderr, "error in function '%s()' %s:%d\n",
               __func__, __FILE__, __LINE__);
         goto clean_and_return;
      }
   }

   // Define the number of blocks.
   unsigned int nb = counters[0]->n_ref;
   for (int i = 1; i < nfiles; i++) {
      if (counters[i]->n_ref != nb) {
         fprintf(stderr, "%s\n", "files have discordant dimensions");
         goto clean_and_return;
      }
   }

   // Define the size of each block.
   size = malloc(nb * sizeof(unsigned int));
   if (size == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      goto clean_and_return;
   }
   for (int i = 0; i < nb; i++) {
      size[i] = counters[0]->n_bins[i];
      for (int j = 1; j < nfiles; j++) {
         if (counters[j]->n_bins[i] != size[i]) {
            fprintf(stderr, "%s\n", "files have discordant dimensions");
            goto clean_and_return;
         }
      }
   }

   // Build the observations vector.
   int nobs = 0;
   for (int i = 0; i < nb; i++) nobs += size[i];
   y = malloc(nfiles * nobs * sizeof(int));
   if (y == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      goto clean_and_return;
   }
   unsigned int offset = 0;
   for (int i = 0; i < nb; i++) {
      for (int j = 0; j < size[i]; j++) {
         for (int k = 0; k < nfiles; k++) {
            y[offset++] = counters[k]->bins[i][j];
         }
      }
   }

   ChIP = new_ChIP((unsigned int)nfiles, nb, y, size);
   if (ChIP == NULL) {
      fprintf(stderr, "error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
   }

clean_and_return:
   if (size != NULL) free(size);
   for (int i = 0; i < nfiles; i++) {
      if (counters[i] != NULL) destroy_counter(counters[i]);
   }
   // Clean 'y' only if something went wrong.
   if (ChIP == NULL) free(y);

   return ChIP;

}

counter_t *
read_count
(
   const char * fn
)
{

   bam1_t * b = NULL;
   bam_hdr_t * h = NULL;
   BGZF * fp = NULL;
   counter_t * counter = NULL;

   // Open .bam file and get header.
   fp = bgzf_open(fn, "r");
   if (fp == NULL) {
      fprintf(stderr, "cannot open file %s\n", fn);
      goto clean_and_return;
   }

   h = bam_hdr_read(fp);
   if (h == NULL) {
      fprintf(stderr, "cannot read header of file %s\n", fn);
      goto clean_and_return;
   }

   b = bam_init1();
   if (b == NULL) {
      fprintf(stderr, "bam! error %s:%d\n", __FILE__, __LINE__);
      goto clean_and_return;
   }

   // Initialize counter_t.
   counter = malloc(sizeof(counter_t));
   if (counter == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      goto clean_and_return;
   }

   int32_t n_ref = h->n_targets;
   counter->n_ref  = n_ref;
   counter->n_bins = calloc(n_ref, sizeof(int32_t));
   counter->bins   = calloc(n_ref, sizeof(int32_t *));

   if (counter->n_bins == NULL || counter->bins == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      destroy_counter(counter);
      counter = NULL;
      goto clean_and_return;
   }

   // Initialize bins of counter_t.
   for (int i = 0; i < n_ref; i++) {
      uint32_t tlen = h->target_len[i];
      tlen = tlen / BIN_SIZE + (tlen % BIN_SIZE > 0);
      counter->n_bins[i] = tlen;
      counter->bins[i] = calloc(tlen, sizeof(int32_t));
      if (counter->bins[i] == NULL) {
         fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
         destroy_counter(counter);
         counter = NULL;
         goto clean_and_return;
      }
   }

   // Count the number of alignments per bin.
   int bytesread;
   while((bytesread = bam_read1(fp, b)) >= 0) {
      if (bytesread > 0 && b->core.tid != -1) {
         int32_t pos = b->core.pos / BIN_SIZE;
         counter->bins[b->core.tid][pos]++;
      }
   }

clean_and_return:
   if (b != NULL) bam_destroy1(b);
   if (h != NULL) bam_hdr_destroy(h);
   if (fp != NULL) bgzf_close(fp);

   return counter;

}

void
destroy_counter
(
   counter_t * counter
)
{
   if (counter->bins != NULL) {
      for (int i = 0; i < counter->n_ref; i++) {
         if (counter->bins[i] != NULL) free(counter->bins[i]);
      }
      free(counter->bins);
   }
   if (counter->n_bins != NULL) free(counter->n_bins);
   free(counter);
}
