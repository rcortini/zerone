#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include "bgzf.h"
#include "debug.h"
#include "parse.h"
#include "sam.h"
#include "zerone.h"


//  ----- Globals ----- //
void * STATE;
int    ERR;


#define SUCCESS 1
#define FAILURE 0

// XXX This has to become a function parameter XXX //
#define BIN_SIZE 300

// Prime number for the size of hash table.
#define HSIZE 997
#define BSIZE 100000000

// Shortcut to acces position of hash table.
#define hv(a) djb2(a) % HSIZE

// Size of the buffer for gzip decompression.
#define CHUNK 16384

// Macro to check all applicable zlib errors.
#define is_zerr(a) ((a) == Z_NEED_DICT || (a) == Z_DATA_ERROR \
      || (a) == Z_MEM_ERROR)

// Type declarations.
struct hash_t;
struct link_t;
struct rod_t;
struct loc_t;

struct bgzf_state_t;
struct generic_state_t;

// Type definitions.
typedef struct link_t link_t;
typedef struct loc_t loc_t;
typedef struct rod_t rod_t;

typedef struct bgzf_state_t bgzf_state_t;
typedef struct generic_state_t generic_state_t;

// Shortcuts.
typedef link_t * hash_t;
typedef char * bloom_t;

// Special functions.
typedef int     (*iter_t)   (loc_t *);
typedef int     (*parser_t) (loc_t *, char *);
typedef ssize_t (*reader_t) (char **, size_t *, FILE *);


// Type definitions.
struct link_t {
   char     seqname[32]; // Chromosome or sequence name.
   rod_t  * counts;      // Counts in chromosome bins.
   link_t * next;        // Next node on the link list.
};

struct rod_t {
   size_t sz;
   size_t mx;
   uint32_t array[];
};

struct loc_t {
   char * name;
   int    pos;
};

// Iterator states.
struct bgzf_state_t {
   BGZF      * file;
   bam_hdr_t * hdr;
   bam1_t    * bam;
};

struct generic_state_t {
   FILE      * file;
   reader_t    reader;
   parser_t    parser;
   char      * buff;
   size_t      bsz;
};


//  ---- Declaration of local functions  ---- //
int      autoparse (const char *, hash_t *);
iter_t   choose_iterator (const char *);

// Iterators.
int      bgzf_iterator (loc_t *);
int      generic_iterator (loc_t *);

// Raders.
ssize_t  getgzline (char **, size_t *, FILE *);

// Parsers.
int      parse_gem (loc_t *, char *);
int      parse_sam (loc_t *, char *);
int      parse_bed (loc_t *, char *);
int      parse_wig (loc_t *, char *);

// Hash handling functions.
int      add_to_rod (rod_t **, uint32_t);
int      bloom_query_and_set(const char *, int, bloom_t);
void     destroy_hash(hash_t *);
link_t * lookup_or_insert (const char *, hash_t *);
ChIP_t * merge_hashes (hash_t **, int);

// Convenience functions.
uint32_t djb2 (const char *);
int      is_gzipped (FILE *);




//  -- Definitions of exported functions  --- //

ChIP_t *
parse_input_files
(
   char * mock_fnames[],
   char * ChIP_fnames[]
)
{

   ChIP_t * ChIP = NULL;         // Return value.
   hash_t * hashtab[512] = {0};  // Array of hash tables.
   int      nhashes = 0;         // Number of hashes.

   // Create a unique hash for all mock files.
   hashtab[0] = calloc(HSIZE, sizeof(link_t *));

   if (hashtab[0] == NULL) {
      debug_print("%s", "memory error\n");
      goto clean_and_return;
   }

   nhashes = 1;

   // Fill in the hash with mock files.
   for (int i = 0; mock_fnames[i] != NULL; i++) {

      debug_print("%s %s\n", "autoparsing mock file", mock_fnames[i]);

      if(!autoparse(mock_fnames[i], hashtab[0])) {
         debug_print("%s", "autoparse failed\n");
         goto clean_and_return;
      }

   }


   // Create distinct hash for each ChIP file.
   for (int i = 0; ChIP_fnames[i] != NULL; i++) {

      debug_print("%s %s\n", "autoparsing ChIP file", ChIP_fnames[i]);

      if (nhashes >= 512) {
         // XXX non debug error XXX //
         fprintf(stderr, "too many files\n");
         goto clean_and_return;
      }

      hashtab[nhashes] = calloc(HSIZE, sizeof(link_t *));
      if (hashtab[nhashes] == NULL) {
         debug_print("%s", "memory error\n");
         goto clean_and_return;
      }

      nhashes++;

      if (!autoparse(ChIP_fnames[i], hashtab[nhashes-1])) {
         debug_print("%s", "autoparse failed\n");
         goto clean_and_return;
      }

   }

   // Merge hash tables in to a 'ChIP_t'.
   ChIP = merge_hashes(hashtab, nhashes);
   if (ChIP == NULL) {
      debug_print("%s", "memory error\n");
      goto clean_and_return;
   }

clean_and_return:
   for (int i = 0 ; i < nhashes ; i++) destroy_hash(hashtab[i]);
   return ChIP;

}



//  ---- Definitions of local functions  ---- //


iter_t
choose_iterator
(
   const char * fname
)
{

   // First check if format is .bam.
   if (strcmp(".bam", fname + strlen(fname) - 4) == 0) {

      bgzf_state_t * state = calloc(1, sizeof(bgzf_state_t));
      if (state == NULL) {
         debug_print("%s", "memory error\n");
         goto exit_bgzf_error;
      }

      state->file = bgzf_open(fname, "r");
      if (state->file == NULL) {
         // XXX non debug error XXX //
         fprintf(stderr, "cannot open file %s\n", fname);
         goto exit_bgzf_error;
      }

      state->hdr = bam_hdr_read(state->file);
      if (state->hdr == NULL) {
         // XXX non debug error XXX //
         fprintf(stderr, "cannot read header from file %s\n", fname);
         goto exit_bgzf_error;
      }

      state->bam = bam_init1();
      if (state->bam == NULL) {
         // XXX non debug error XXX //
         fprintf(stderr, "bam! error (sorry)\n");
         goto exit_bgzf_error;
      }

      STATE = state;
      return bgzf_iterator;

exit_bgzf_error:
      if (state != NULL) {
         if (state->file != NULL) bgzf_close(state->file);
         if (state->hdr != NULL) bam_hdr_destroy(state->hdr);
         if (state->bam != NULL) bam_destroy1(state->bam);
      }
      free(state);
      return NULL;

   }

   // Format is not .bam, so it is generic.
   generic_state_t * state = calloc(1, sizeof(generic_state_t));

   if (state == NULL) {
      debug_print("%s", "memory error\n");
      goto exit_generic_error;
   }

   state->file = fopen(fname, "r");

   if (state->file == NULL) {
      // XXX non debug error XXX //
      fprintf(stderr, "cannot open file %s\n", fname);
      goto exit_generic_error;
   }

   state->bsz = 32;
   state->buff = malloc(state->bsz * sizeof(char));

   if (state->buff == NULL) {
      debug_print("%s", "memory error\n");
      goto exit_generic_error;
   }

   if (strcmp(".map", fname + strlen(fname) - 4) == 0) {
      state->reader = getline;
      state->parser = parse_gem;
   }
   else if (strcmp(".sam", fname + strlen(fname) - 4) == 0) {
      state->reader = getline;
      state->parser = parse_sam;
   }
   else if (strcmp(".bed", fname + strlen(fname) - 4) == 0) {
      state->reader = getline;
      state->parser = parse_bed;
   }
   else if (strcmp(".wig", fname + strlen(fname) - 4) == 0) {
      state->reader = getline;
      state->parser = parse_wig;
   }
   else if (strcmp(".map.gz", fname + strlen(fname) - 7) == 0) {
      state->reader = getgzline;
      state->parser = parse_gem;
   }
   else if (strcmp(".sam.gz", fname + strlen(fname) - 7) == 0) {
      state->reader = getgzline;
      state->parser = parse_sam;
   }
   else if (strcmp(".bed.gz", fname + strlen(fname) - 7) == 0) {
      state->reader = getgzline;
      state->parser = parse_bed;
   }
   else if (strcmp(".wig.gz", fname + strlen(fname) - 7) == 0) {
      state->reader = getgzline;
      state->parser = parse_wig;
   }
   else {
      // XXX non debug error XXX //
      fprintf(stderr, "unknown format for file %s\n", fname);
      goto exit_generic_error;
   }

   STATE = state;
   return generic_iterator;

exit_generic_error:
   if (state != NULL) {
      if (state->file != NULL) fclose(state->file);
      free(state->buff);
   }
   free(state);
   return NULL;

}


int
autoparse
(
   const char   * fname,
         hash_t * hashtab
)
{

   int status = SUCCESS;

   bloom_t bloom = NULL;
   loc_t loc = {0};

   // Find an iterator for the given file type.
   iter_t iterate = choose_iterator(fname);

   if (iterate == NULL) {
      debug_print("%s", "choosing iterator failed\n");
      status = FAILURE;
      goto clean_and_return;
   }

   bloom = calloc(1, BSIZE);
   if (bloom == NULL) {
      debug_print("%s", "memory error\n");
      status = FAILURE;
      goto clean_and_return;
   }

   ERR = 0;
   while (iterate(&loc) > 0) {

      // 'loc.name' is set to NULL for unmapped reads.
      if (loc.name == NULL) continue;

      // Check in Bloom filter whether the read was seen before.
      if (bloom_query_and_set(loc.name, loc.pos, bloom)) continue;

      link_t *lnk = lookup_or_insert(loc.name, hashtab);

      if (lnk == NULL) {
         debug_print("%s", "hash query failed\n");
         status = FAILURE;
         goto clean_and_return;
      }

      // Add read to counts.
      if (!add_to_rod(&lnk->counts, loc.pos / BIN_SIZE)) {
         debug_print("%s", "adding read failed\n");
         status = FAILURE;
         goto clean_and_return;
      }

   }

   if (ERR) {
      // Unexpected exit from iteration.
      status = FAILURE;
      goto clean_and_return;
   }

clean_and_return:
   free(bloom);
   return status;

}


int
generic_iterator
(
   loc_t *loc
)
{
   // Cast as state for generic iterator.
   generic_state_t *state = (generic_state_t *) STATE;
   parser_t parse = state->parser;
   reader_t read = state->reader;

   if (loc == NULL) {
      // Interruption.
      goto clean_and_return;
   }

   // Read one line.
   int nbytes = read(&state->buff, &state->bsz, state->file);

   if (nbytes < -1) {
      ERR = __LINE__;
      goto clean_and_return;
   }

   if (nbytes == -1) {
      // End of file.
      goto clean_and_return;
   }

   if (!parse(loc, state->buff)) {
      // XXX non debug error XXX //
      fprintf(stderr, "format conflict in line:\n%s", state->buff);
      ERR = __LINE__;
      goto clean_and_return;
   }

   return nbytes;

clean_and_return:
   free(state->buff);
   fclose(state->file);
   free(state);

   STATE = NULL;
   return -1;

}

int
bgzf_iterator
(
   loc_t *loc
)
{

   // Cast as state for bgzf iterator.
   bgzf_state_t *state = (bgzf_state_t *) STATE;

   if (loc == NULL) {
      // Interruption.
      goto clean_and_return;
   }

   int bytesread = bam_read1(state->file, state->bam);

   if (bytesread < -1 || state->bam->core.tid < 0) {
      ERR = __LINE__;
      goto clean_and_return;
   }

   if (bytesread == -1) {
      // End of file.
      goto clean_and_return;
   }


   loc->name = state->hdr->target_name[state->bam->core.tid];
   loc->pos = state->bam->core.pos;

   return bytesread;

clean_and_return:
   // This bit is executed when the iterator needs to be
   // cleaned (end of file, interruption or error).
   bam_hdr_destroy(state->hdr);
   bgzf_close(state->file);
   bam_destroy1(state->bam);
   free(state);

   STATE = NULL;
   return -1;

}


int
parse_gem
(
   loc_t *loc,
   char  *line
)
{

   // Get last field of the line.
   char *map = strrchr(line, '\t');

   // No tab character found in line.
   if (map++ == NULL) return FAILURE;

   if (map[0] == '-') {
      loc->name = NULL;
      loc->pos = 0;
      return SUCCESS;
   }

   // Parse mapping data.
   char *chrom = strsep(&map, ":");
                 strsep(&map, ":"); // Discard strand.
   char *tmp =   strsep(&map, ":");

   // Cannot find chromosome or position.
   if (chrom == NULL || tmp == NULL) return FAILURE;

   // Positions in the genome cannot be 0, so we can identify
   // failures of 'atoi' to convert numbers.
   int pos = atoi(tmp);

   // Field position is not a number.
   if (pos == 0) return FAILURE;

   loc->name = chrom;
   loc->pos = pos;

   return SUCCESS;

}

int
parse_sam
(
   loc_t *loc,
   char  *line
)
{
   // Ignore header.
   if (line[0] == '@') {
      loc->name = NULL;
      return SUCCESS;
   }

                 strsep(&line, "\t"); // Discard QNAME.
                 strsep(&line, "\t"); // Discard FLAG.
   char *chrom = strsep(&line, "\t");
   char *tmp   = strsep(&line, "\t");

   // Cannot find chromosome or position.
   if (chrom == NULL || tmp == NULL) return FAILURE;

   // Unmapped read.
   if (strcmp(chrom, "*") == 0) {
      loc->name = NULL;
      return SUCCESS;
   }

   // Positions in the genome cannot be 0, so we can identify
   // failures of 'atoi' to convert numbers.
   int pos = atoi(tmp);

   // Field position is not a number.
   if (pos == 0) return FAILURE;

   loc->name = chrom;
   loc->pos = pos;

   return SUCCESS;
}


int
parse_bed
(
   loc_t *loc,
   char  *line
)
{

   char *chrom = strsep(&line, "\t");
   char *tmp1  = strsep(&line, "\t");
   char *tmp2  = strsep(&line, "\t");

   // Cannot find chromosome or position.
   if (chrom == NULL || tmp1 == NULL || tmp2 == NULL) return FAILURE;

   // Positions in the genome cannot be 0, so we can identify
   // failures of 'atoi' to convert numbers.
   int start = atoi(tmp1);
   int end = atoi(tmp2);

   // Field position is not a number.
   if (start == 0 || end == 0) return FAILURE;

   loc->name = chrom;
   loc->pos = (start + end) / 2;

   return SUCCESS;
}


int
parse_wig
(
   loc_t *loc,
   char  *line
)
{

   // XXX CAUTION: this function is non-reentrant. XXX //
   // XXX Use multithreading at your own risk.     XXX //

   static int   fixedstep = 0;
   static char *chrom     = NULL;
   static int   fstart    = 1;
   static int   step      = 1;
   static int   span      = 1;
   static int   iter      = 0;

   // Ignore track definition lines.
   if (strncmp(line, "track", 5) == 0) return SUCCESS;

   // Detect type of format.
   if  (strncmp(line, "variableStep", 12) == 0) fixedstep = 0;
   else if (strncmp(line, "fixedStep", 9) == 0) fixedstep = 1;

   // Parse data line.
   else {
      int start, end, reads;

      if (fixedstep) start = fstart + iter++ * step;
      else {
         start = atoi(strsep(&line, "\t"));
         if (line == NULL) return FAILURE;
      }

      char *endptr;
      // Not using read counts for now.
      reads = strtol(line, &endptr, 10);
      if (line == endptr) return FAILURE;

      strsep(&line, "\t");
      if (line != NULL) return FAILURE;

      if (chrom == NULL || start <= 0 || reads < 0) return FAILURE;

      end = start + span - 1;
      loc->pos = (start + end) / 2;

      return SUCCESS;
   }

   // Parse data definition lines.
   free(chrom);
                  strsep(&line, "\t");
   chrom = strdup(strsep(&line, "\t") + 6);
   loc->name = chrom;

   if (fixedstep) {
      fstart = atoi(strsep(&line, "\t") + 6);
      step   = atoi(strsep(&line, "\t") + 5);
      iter   = 0;
   }

   if (line == NULL) span = 1;
   else              span = atoi(line + 5);

   if (fstart <= 0 || step <= 0 || span <= 0) return FAILURE;

   return SUCCESS;
}


//  ---------  Utility functions  ----------  //

uint32_t
djb2
(
   const char * s
)
// The magic djb2 (http://www.cse.yorku.ca/~oz/hash.html).
{
   uint32_t val = 5381;
   for (int i = 0 ; i < 31 && s[i] != '\0' ; i++)
      val = val * 33 ^ s[i];
   return val;
}


void
destroy_hash
(
   hash_t * hashtab
)
{
   for (int i = 0 ; i < HSIZE ; i++) {
      for (link_t *lnk = hashtab[i] ; lnk != NULL ; ) {
         link_t *tmp = lnk->next;
         free(lnk->counts);
         free(lnk);
         lnk = tmp;
      }
   }
   free(hashtab);
}


int
is_gzipped
(
   FILE *fp
)
{
   // Get first two bytes.
   char c1 = getc(fp);
   char c2 = getc(fp);
   // Put them back.
   ungetc(c2, fp);
   ungetc(c1, fp);

   return c1 == '\x1f' && c2 == '\x8b';

}


ssize_t
getgzline
(
   char  ** buff,
   size_t * bsz,
   FILE   * gzfile
)
// The function is absolutely not re-entrant.
{

   static z_stream strm;
   static unsigned char *start;

   // 'sentinel' is always 'out + CHUNK - strm.avail_out'.
   static unsigned char *sentinel;
   static unsigned char out[CHUNK];
   static unsigned char in[CHUNK];

   int zstat;

   // Set to 0 upon first call.
   static int is_initialized;

   // Set 'buff' to NULL to interrupt inflation.
   if (is_initialized && buff == NULL) {
      (void) inflateEnd(&strm);
      is_initialized = 0;
      return -1;
   }

   if (!is_initialized) {

      // Allocate inflate state.
      strm.zalloc = Z_NULL;
      strm.zfree = Z_NULL;
      strm.opaque = Z_NULL;
      strm.next_in = Z_NULL;
      strm.avail_in = 0;
      strm.avail_out = CHUNK;

      start = out;
      sentinel = out; // 'out + CHUNK - strm.avail_out'.

      zstat = inflateInit2(&strm, 16+MAX_WBITS);
      if (zstat != Z_OK) goto exit_io_error;

      is_initialized = 1;

   }

   // Try getting newline character from buffer.
   // If not found, 'end' will point to the sentinel.
   unsigned char *end = memchr(start, '\n', sentinel - start);

   if (end == NULL) {

      // No hit: we need to feed the buffer.
      // First we shift 'out' as much as possible.
      memmove(out, start, sentinel - start);

      strm.avail_out += (start - out);
      strm.next_out = out + (sentinel - start);

      // Then refill 'out' with data from 'in'.
      zstat = inflate(&strm, Z_SYNC_FLUSH);
      if (is_zerr(zstat)) goto exit_io_error;

      if (strm.avail_out > 0) {

         // If 'out' is not full, read in more data.
         strm.avail_in = fread(in, 1, CHUNK, gzfile);
         if (ferror(gzfile)) goto exit_io_error;
         strm.next_in = in;

         // Inflate again to fill 'out'.
         zstat = inflate(&strm, Z_SYNC_FLUSH);
         if (is_zerr(zstat)) goto exit_io_error;

      }

      // Unless errors occurred, 'start' now points
      // to the beginning of 'out' and 'sentinel'
      // points to the end of inflated data.

      start = out;
      sentinel = out + CHUNK - strm.avail_out;

      // We can run the search again.
      end = memchr(start, '\n', sentinel - start);

   }

   if (end == NULL) {
      // No newline character. Check if input file is
      // read in full. If not, something went wrong.
      if (zstat == Z_STREAM_END) {
         (void) inflateEnd(&strm);
         is_initialized = 0;
         return -1;
      }
      else goto exit_io_error;
   }
   else {
      // The newline character was found. Copy line
      // to write buffer and return read bytes.
      int nbytes = end - start + 1;

      if (*bsz < nbytes + 1) {
         size_t newbsz = 2*nbytes;
         char *tmp = realloc(*buff, newbsz);
         if (tmp == NULL) abort();
         *buff = tmp;
         *bsz = newbsz;
      }

      memcpy(*buff, start, nbytes);
      (*buff)[nbytes] = '\0';

      start = end+1;
      return nbytes;
   }

exit_io_error:
   (void) inflateEnd(&strm);
   is_initialized = 0;
   ERR = __LINE__;
   return -2;

}



int
add_to_rod
(
   rod_t   ** addr,
   uint32_t   pos
)
{

   // Convenience variable.
   rod_t *rod = *addr;

   if (pos >= rod->sz) {
      // Double size of the rod.
      size_t oldsz = rod->sz;
      size_t newsz = 2*pos;
      rod_t *tmp = realloc(rod, sizeof(rod_t) + newsz*sizeof(uint32_t));
      if (tmp == NULL) {
         debug_print("%s", "memory error\n");
         return FAILURE;
      }
      // Update all.
      *addr = rod = tmp;
      memset(rod->array + oldsz, 0, (newsz-oldsz)*sizeof(uint32_t));
      rod->sz = newsz;
   }

   rod->array[pos]++;
   if (pos > rod->mx) rod->mx = pos;

   return SUCCESS;

}


ChIP_t *
merge_hashes
(
   hash_t ** hashes,
   int       nhashes
)
{
   // Add keys and update a reference hash table.
   hash_t *refhash = hashes[0];

   for (int i = 1 ; i < nhashes ; i++) {
      hash_t *hashtab = hashes[i];
      for (int j = 0 ; j < HSIZE ; j++) {
         for(link_t *lnk = hashtab[j] ; lnk != NULL ; lnk = lnk->next) {
            // Get chromosome from reference hash (or create it).
            link_t *reflnk = lookup_or_insert(lnk->seqname, refhash);
            if (reflnk == NULL) {
               debug_print("%s", "hash query failed\n");
               return NULL;
            }
            // Update the max in reference hash. Note that
            // the size of this 'link_t' may be smaller.
            if (lnk->counts->mx > reflnk->counts->mx)
               reflnk->counts->mx = lnk->counts->mx;
         }
      }
   }

   // Now use the data in reference hash to creat 'ChIP_t'.
   int nkeys = 0;
   int nbins = 0;
   for (int j = 0 ; j < HSIZE ; j++) {
      for(link_t *lnk = refhash[j] ; lnk != NULL ; lnk = lnk->next) {
         nkeys++;
         nbins += lnk->counts->mx;
      }
   }

   unsigned int *size = malloc(nkeys * sizeof(unsigned int));
   char *name = malloc(nkeys * 32);
   char **nptr = malloc(nkeys * sizeof(char *));
   int *y = calloc(nhashes * nbins, sizeof(int));

   if (size == NULL || name == NULL || y == NULL || nptr == NULL) {
      debug_print("%s", "memory error\n");
      return NULL;
   }

   // Fill in the observations.
   int m = 0;
   size_t offset = 0;
   for (int j = 0 ; j < HSIZE ; j++) {
   for (link_t *rlnk = refhash[j] ; rlnk != NULL ; rlnk = rlnk->next) {

      // Get seqname and sequence length.
      char *key = rlnk->seqname;
      nptr[m] = name + 32*m;
      strncpy(nptr[m], key, 32);
      size_t blksz = size[m++] = rlnk->counts->mx;

      // Go through all the hashes to get the data.
      for (int i = 0 ; i < nhashes ; i++) {
         hash_t *hashtab = hashes[i];
         link_t *lnk = lookup_or_insert(key, hashtab);
         if (lnk == NULL) {
            debug_print("%s", "hash query failed\n");
            return NULL;
         }
         rod_t *counts = lnk->counts;
         size_t kmax = counts->sz < blksz ? counts->sz : blksz;
         for (int k = 0 ; k < kmax ; k++) {
            y[offset + (nhashes * k) + i] = counts->array[k];
         }
      }

      // Update offset.
      offset += nhashes * blksz;

   }
   }

   ChIP_t *ChIP = new_ChIP(nhashes, nkeys, y, (const char **) nptr, size);

   free(size);
   free(name);
   free(nptr);

   return ChIP;

}


link_t *
lookup_or_insert
(
   const char   * s,
         hash_t * htab
)
// Famous K&R simple hash lookup.
{

   for (link_t *lnk = htab[hv(s)] ; lnk != NULL ; lnk = lnk->next) {
      if (strcmp(s, lnk->seqname) == 0) return lnk;
   }

   // Entry was not found.
   link_t *new = malloc(sizeof(link_t));

   if (new == NULL) {
      debug_print("%s", "memory error\n");
      return NULL;
   }

   rod_t *rod = malloc(sizeof(rod_t) + 32*sizeof(uint32_t));
   if (rod == NULL) {
      debug_print("%s", "memory error\n");
      free(new);
      return NULL;
   }

   // Initiatile 'rod_t'.
   memset(rod->array, 0, 32*sizeof(uint32_t));
   rod->sz = 32;
   rod->mx = 0;

   // Update 'link_t'.
   strncpy(new->seqname, s, 31);
   new->counts = rod;

   // Add 'link_t' to table hash table.
   new->next = htab[hv(s)];
   htab[hv(s)] = new;

   return new;

}

int
bloom_query_and_set
(
   const char  * name,
         int     pos,
         bloom_t bloom
)
{

   uint32_t a = (djb2(name) + pos) % (8*BSIZE);
   uint32_t b = (a + pos) % (8*BSIZE);

   int abyte = a/8, abit = a%8;
   int bbyte = a/8, bbit = b%8;

   // Get bits.
   int a_yes = (bloom[abyte] >> (abit)) & 1;
   int b_yes = (bloom[bbyte] >> (bbit)) & 1;

   // Set bits.
   bloom[abyte] |= (1 << abit);
   bloom[bbyte] |= (1 << bbit);

   return a_yes && b_yes;

}
