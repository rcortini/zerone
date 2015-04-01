#include "samread.h"

int
is_sam
(
   const char * fn
)
{
   return strcmp(".sam", fn + strlen(fn) - 4) == 0 ||
          strcmp(".bam", fn + strlen(fn) - 4) == 0;
}

ChIP_t *
read_sam
(
   char * fn[],
   int    nfiles
)
{
   // Build the counter_t's from which to build the final ChIP_t.
   counter_t * counters[nfiles];
   for (int i = 0; i < nfiles; i++) counters[i] = read_count(fn[i]);

   // Define the number of blocks.
   unsigned int nb = counters[0]->n_ref;
   for (int i = 1; i < nfiles; i++) {
      if (counters[i]->n_ref != nb) {
         fprintf(stderr, "%s\n", "files have discordant dimensions");
         return NULL;
      }
   }

   // Define the size of each block.
   unsigned int * size = malloc(nb * sizeof(unsigned int));
   if (size == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }
   for (int i = 0; i < nb; i++) {
      size[i] = counters[0]->n_bins[i];
      for (int j = 1; j < nfiles; j++) {
         if (counters[j]->n_bins[i] != size[i]) {
            fprintf(stderr, "%s\n", "files have discordant dimensions");
            return NULL;
         }
      }
   }

   // Build the observations vector.
   int nobs;
   for (int i = 0; i < nb; i++) nobs += size[i];
   int * y = malloc(nfiles * nobs * sizeof(int));
   if (y == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }
   unsigned int offset = 0;
   for (int i = 0; i < nb; i++) {
      for (int j = 0; j < size[i]; j++) {
         for (int k = 0; k < nfiles; k++) {
            y[offset++] = counters[k]->bins[i][j];
         }
      }
   }

fprintf(stderr, "%s\n", "ok!");
   ChIP_t * ChIP = new_ChIP(nfiles, nb, y, size);
   free(size);

   // Create ChIP_t and define the number of blocks.
   //unsigned int nb = counters[0]->n_ref;
   //ChIP_t * ChIP = malloc(sizeof(ChIP_t) + nb * sizeof(unsigned int));
   //if (ChIP == NULL) {
   //   fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
   //   return NULL;
   //}
   //ChIP->r = nfiles;
   //ChIP->nb = nb;
   //for (int i = 1; i < nfiles; i++) {
   //   if (counters[i]->n_ref != ChIP->nb) {
   //      fprintf(stderr, "%s\n", "files have discordant dimensions");
   //      return NULL;
   //   }
   //}

   // Define the size of each block.
   //for (int i = 0; i < ChIP->nb; i++) {
   //   ChIP->size[i] = counters[0]->n_bins[i];
   //   for (int j = 1; j < nfiles; j++) {
   //      if (counters[j]->n_bins[i] != ChIP->size[i]) {
   //         fprintf(stderr, "%s\n", "files have discordant dimensions");
   //         return NULL;
   //      }
   //   }
   //}

   // Build the observations vector.
   //ChIP->y = malloc(ChIP->r * nobs(ChIP) * sizeof(int));
   //if (ChIP->y == NULL) {
   //   fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
   //   return NULL;
   //}
   //unsigned int offset = 0;
   //for (int i = 0; i < ChIP->nb; i++) {
   //   for (int j = 0; j < ChIP->size[i]; j++) {
   //      for (int k = 0; k < nfiles; k++) {
   //         ChIP->y[offset++] = counters[k]->bins[i][j];
   //      }
   //   }
   //}

fprintf(stderr, "%s\n", "5/5 finish read_sam");
   for (int i = 0; i < nfiles; i++) destroy_counter(counters[i]);

   return ChIP;
}

counter_t *
read_count
(
   const char * fn
)
{
   samfile_t * fp = samopen(fn, "r", NULL);

   counter_t * counter = malloc(sizeof(counter_t));
   int32_t n_ref = fp->header->n_targets;
   counter->n_ref  = n_ref;
   counter->n_bins = malloc(n_ref * sizeof(int32_t));
   counter->bins   = malloc(n_ref * sizeof(int32_t *));
   if (counter == NULL || counter->n_bins == NULL || counter->bins == NULL) {
      fprintf(stderr, "memory error %s:%d\n", __FILE__, __LINE__);
      return NULL;
   }

   for (int i = 0; i < n_ref; i++) {
      uint32_t tlen = fp->header->target_len[i];
      tlen = tlen / BIN_SIZE + (tlen % BIN_SIZE > 0);
      counter->n_bins[i] = tlen;
      counter->bins[i] = calloc(tlen, sizeof(int32_t));
   }

fprintf(stderr, "%s\n", "1/2 read_count");
   bam1_t * b;
   int bytesread;
   while((bytesread = samread(fp, (b = bam_init1()))) >= 0) {
//fprintf(stderr, "%s:%d\n", "br", bytesread);
      if (bytesread > 0) {
         int32_t pos = b->core.pos / BIN_SIZE + (b->core.pos % BIN_SIZE > 0);
//fprintf(stderr, "%s:%d\n", "pos", pos);
         counter->bins[b->core.tid][pos]++;
      }
      bam_destroy1(b);
   }

//int bytesread;
//do {
//   bam1_t * b = bam_init1();
//   bytesread = samread(fp, b);
//   int32_t pos = b->core.pos / BIN_SIZE + (b->core.pos % BIN_SIZE > 0);
//   counter->bins[b->core.tid][pos]++;
//   bam_destroy1(b);
//} while (bytesread >= 0);

fprintf(stderr, "%s\n", "2/2 printing counter values");
for (int i = 0; i < counter->n_ref; i++) {
   for (int j = 0; j < counter->n_bins[i]; j++) {
      fprintf(stderr, "%s:%d\n", "count", counter->bins[i][j]);
   }
}

   samclose(fp);
   return counter;
}

void
destroy_counter
(
   counter_t * counter
)
{
   for (int i = 0; i < counter->n_ref; i++) free(counter->bins[i]);
   free(counter->n_bins);
   free(counter->bins);
   free(counter);
}
