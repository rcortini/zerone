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

#include <errno.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <zlib.h>
#include "bgzf.h"
#include "ctype.h"
#include "debug.h"
#include "parse.h"
#include "sam.h"
#include "zerone.h"


//  ----- Globals ----- //
void * STATE;
int    ERR;


#define SUCCESS 1
#define FAILURE 0

// Macro function.
#define min(a,b) ((a) < (b) ? (a) : (b))
#define max(a,b) ((a) > (b) ? (a) : (b))
#define abs(a)   ((a) < 0 ? -(a) : (a))

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
struct bitf_t;
struct loc_t;

struct bgzf_state_t;
struct generic_state_t;

// Type definitions.
typedef struct link_t link_t;
typedef struct loc_t loc_t;
typedef struct rod_t rod_t;
typedef struct bitf_t bitf_t;
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
   bitf_t * repeats;     // Read bitfield.
   link_t * next;        // Next node on the link list.
};

struct rod_t {
   size_t sz;            // Size (in bytes) of the array.
   size_t mx;            // Highest index with nonzero count.
   uint32_t array[];
};

struct bitf_t {
   size_t sz;            // Size (in bytes) of the array.
   uint8_t array[];
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
int      autoparse (const char *, hash_t *, const int);
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
int      bitf_query_and_set (int, link_t *);
void     destroy_hash(hash_t *);
void     destroy_bitfields(hash_t *);
link_t * lookup_or_insert (const char *, hash_t *);
ChIP_t * merge_hashes (hash_t **, int, int);

// Convenience functions.
uint32_t djb2 (const char *);
int      is_gzipped (FILE *);




//  -- Definitions of exported functions  --- //

ChIP_t *
parse_input_files
(
         char * mock_fnames[],
         char * ChIP_fnames[],
   const int    window
)
{

   ChIP_t * ChIP = NULL;         // Return value.
   hash_t * hashtab[512] = {0};  // Array of hash tables.
   int      nhashes = 0;         // Number of hashes.

   // Create a unique hash table for all mock files.
   hashtab[0] = calloc(HSIZE, sizeof(link_t *));

   if (hashtab[0] == NULL) {
      debug_print("%s", "memory error\n");
      goto clean_and_return;
   }

   nhashes = 1;

   // Fill in the hash with mock files.
   // Because the same hash table is used every time,
   // the reads in the same window are summed even if
   // they are from different files.
   for (int i = 0; mock_fnames[i] != NULL; i++) {

      debug_print("%s %s\n", "autoparsing mock file", mock_fnames[i]);

      if(!autoparse(mock_fnames[i], hashtab[0], window)) {
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

      if (!autoparse(ChIP_fnames[i], hashtab[nhashes-1], window)) {
         debug_print("%s", "autoparse failed\n");
         goto clean_and_return;
      }

   }

   // Merge hash tables in to a 'ChIP_t'. The last argument
   // says whether any mock file was provided.
   ChIP = merge_hashes(hashtab, nhashes, mock_fnames[0] == NULL);
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
         hash_t * hashtab,
   const int      window
)
{

   int status = SUCCESS;

   loc_t loc = {0};

   // Find an iterator for the given file type.
   iter_t iterate = choose_iterator(fname);

   if (iterate == NULL) {
      debug_print("%s", "choosing iterator failed\n");
      status = FAILURE;
      goto clean_and_return;
   }

   ERR = 0;
   while (iterate(&loc) > 0) {

      // 'loc.name' is set to NULL for unmapped reads.
      if (loc.name == NULL) continue;

      // Check in Bloom filter whether the read was seen before.
      // if (bloom_query_and_set(loc.name, loc.pos, bloom)) continue;

      link_t *lnk = lookup_or_insert(loc.name, hashtab);

      if (lnk == NULL) {
         debug_print("%s", "hash query failed\n");
         status = FAILURE;
         goto clean_and_return;
      }

      // Check in bit field whether the read was seen before.
      if (bitf_query_and_set(loc.pos, lnk)) continue;

      // Add read to counts.
      if (!add_to_rod(&lnk->counts, loc.pos / window)) {
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
   destroy_bitfields(hashtab);
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
   reader_t doread = state->reader;

   if (loc == NULL) {
      // Interruption.
      goto clean_and_return;
   }

   // Read one line.
   int nbytes = doread(&state->buff, &state->bsz, state->file);

   if (nbytes < -1) {
      ERR = __LINE__;
      goto clean_and_return;
   }

   if (nbytes == -1) {
      // End of file.
      goto clean_and_return;
   }

   // Many parsers modify the line in place. For purposes
   // of reporting errors we make a copy of it beforehand.
   size_t sz = strlen(state->buff);
   char buffer[64] = {0};
   strncpy(buffer, state->buff, 63);

   if (!parse(loc, state->buff)) {
      // XXX non debug error XXX //
      fprintf(stderr, "format conflict in line:\n%s%s", buffer,
            sz > 63 ? "...\n" : "");
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
   bam1_core_t core = state->bam->core;

   if (bytesread < -1) {
      ERR = __LINE__;
      goto clean_and_return;
   }

   if (bytesread == -1) {
      // End of file.
      goto clean_and_return;
   }

   // Discard unmapped and 2nd align of PE files.
   if (core.tid < 0 || core.flag & BAM_FREAD2)
      loc->name = NULL;
   else 
      loc->name = state->hdr->target_name[core.tid];

   // Note that the bam format is 0-based, so we add 1 to the
   // position because genomic positions are 1-based.
   
   // Compute middle point for PE intervals.
   if (core.flag & BAM_FPAIRED) {
      if (core.mtid != core.tid) loc->name = NULL;
      loc->pos = 1 + min(core.pos, core.mpos) + abs(core.isize)/2;
   }
   else
      loc->pos = 1 + core.pos;

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
   char *flag  = strsep(&line, "\t");
   char *chrom = strsep(&line, "\t");
   char *tmp   = strsep(&line, "\t");

   // Cannot find chromosome or position.
   if (chrom == NULL || tmp == NULL) return FAILURE;

   // Unmapped read.
   if (strcmp(chrom, "*") == 0) {
      loc->name = NULL;
      return SUCCESS;
   }

   // Binary sam flags.
   int blag = atoi(flag);

   // Always skip read 2.
   if (blag & BAM_FREAD2) {
      loc->name = NULL;
      return SUCCESS;
   }
   
   // Note that the sam format is 1-based, so if 'atoi()'
   // has returned 0, something is wrong with the format.
   int pos = atoi(tmp);
   if (pos == 0) return FAILURE;

   // If read is paired use mid-point of the mapping.
   if (blag & BAM_FPAIRED) {
                     strsep(&line, "\t"); // Discard MAPQ.
                     strsep(&line, "\t"); // Discard CIGAR.
      char *pchrom = strsep(&line, "\t");
      char *mpos   = strsep(&line, "\t");
      char *isize  = strsep(&line, "\t");

      // If pairing failed, skip read.
      if (strcmp(pchrom, "=") != 0) {
         loc->name = NULL;
         return SUCCESS;
      }

      pos = min(pos, atoi(mpos)) + abs(atoi(isize))/2;

   }

   loc->name = chrom;
   loc->pos = pos;

   return SUCCESS;

}


int
strtoul_check
(
  char *endptr,
  char *nptr
)
{
   if (endptr == NULL)    return 0;
   if (endptr == nptr)    return 0;
   if (*endptr == '\0')   return 1;
   if (isspace(*endptr)) return 1;
   return 0;
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
   if (chrom == NULL || tmp1 == NULL || tmp2 == NULL)
      return FAILURE;

   // The bed format is 0-based. We cannot use 'atoi()' to
   // parse the coodrinates.
   char *endptr = NULL;

   errno = 0;
   int start = 1 + strtoul(tmp1, &endptr, 10);
   if (!strtoul_check(endptr, tmp1))
      return FAILURE;

   int end = 1 + strtoul(tmp2, &endptr, 10);
   if (!strtoul_check(endptr, tmp2))
      return FAILURE;

   // 'strtoul' may set 'errno' in case of overflow.
   if (errno)
      return FAILURE;

   loc->name = chrom;
   loc->pos = (start + end) / 2;

   return SUCCESS;

}


int
is_wig_defline
(
   char *line,
   int  *fixedstep
)
// Check whether line is a .wig definition line. If so,
// updates 'fixedstep' accordingly.
{

   if (strncmp(line, "variableStep", 12) == 0) {
      *fixedstep = 0;
      return 1;
   }

   else if (strncmp(line, "fixedStep", 9) == 0) {
      *fixedstep = 1;
      return 1;
   }

   // Not a definition line.
   else return 0;

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

   static char chrom[32] = {0};
   static int  fixedstep =  0;
   static int  fstart    =  1;
   static int  step      =  1;
   static int  span      =  1;
   static int  iter      =  0;

   // Ignore track definition lines.
   if (strncmp(line, "track", 5) == 0) return SUCCESS;

   // Check if line is definition.
   if (is_wig_defline(line, &fixedstep)) {

      // Definition lines follow either of the following formats.
      // [A] variableStep\tchrom=chr\t[span=x]
      // [B] fixedStep\tchrom=chr\tstart=x\tstep=y\t[span=z]

      // Skip first token (either "fixedStep" or "variableStep").
      strsep(&line, "\t");

      // Copy chromosome name (and remove the 6 characters of "chrom=").
      strncpy(chrom, strsep(&line, "\t") + 6, 31);
      loc->name = chrom;

      if (fixedstep) {
         // Format [B] has more fields.
         // Remove the characters of "start=" and "step=".
         fstart = atoi(strsep(&line, "\t") + 6);
         step   = atoi(strsep(&line, "\t") + 5);
         iter   = 0;
      }

      // Check if there is optional "span" field at end of line.
      if (line == NULL) span = 1;
      // Remove the 5 characters of "span=".
      else              span = atoi(line + 5);

      // Final sanity check (return SUCCESS if all pass).
      return (fstart > 0 && step > 0 && span > 0);

   }

   // Line is data.
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

      // Final sanity check.
      if (line != NULL || chrom[0] == '\0' || start <= 0 || reads < 0)
         return FAILURE;

      end = start + span - 1;
      loc->pos = (start + end) / 2;

      return SUCCESS;

   }

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
         if (lnk->repeats != NULL) free(lnk->repeats);
         free(lnk);
         lnk = tmp;
      }
   }
   free(hashtab);
}

void
destroy_bitfields
(
 hash_t * hashtab
)
{
   for (int i = 0 ; i < HSIZE ; i++) {
      for (link_t *lnk = hashtab[i] ; lnk != NULL ; lnk = lnk->next) {
         if (lnk->repeats != NULL) free(lnk->repeats);
         lnk->repeats = NULL;
      }
   }
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
         int       nhashes,
   const int       no_mock
)
{
   // Add keys and update a reference hash table (the
   // first, i.e. the hash table of mock controls).
   hash_t *refhash = hashes[0];

   for (int i = 1 ; i < nhashes ; i++) {
      // Run over all other hash tables.
      hash_t *hashtab = hashes[i];
      for (int j = 0 ; j < HSIZE ; j++) {
         // Visit every cell of the hash table and
         // through every link of the list.
         for (link_t *lnk = hashtab[j] ; lnk != NULL ; lnk = lnk->next) {
            // Get chromosome from reference hash (or create it).
            link_t *reflnk = lookup_or_insert(lnk->seqname, refhash);
            if (reflnk == NULL) {
               debug_print("%s", "hash query failed\n");
               return NULL;
            }
            // Update the max in reference hash. Note that
            // the size of this 'link_t' may be smaller.
            if (lnk->counts->mx > reflnk->counts->mx) {
               reflnk->counts->mx = lnk->counts->mx;
            }
         }
      }
   }

   // The reference table now contains all the chromosomes and the
   // max value was updated as the upper bound of all the tables.
   // Now use the data in reference hash to create 'ChIP_t'.
   int nkeys = 0;
   int nbins = 0;
   for (int j = 0 ; j < HSIZE ; j++) {
      for (link_t *lnk = refhash[j] ; lnk != NULL ; lnk = lnk->next) {
         nkeys++;
         nbins += lnk->counts->mx + 1;
      }
   }

   // Allocate all.
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

   // Iterate over the blocks (chromosomes).
   for (int j = 0 ; j < HSIZE ; j++) {
   for (link_t *rlnk = refhash[j] ; rlnk != NULL ; rlnk = rlnk->next) {

      // Get seqname and sequence length.
      char *key = rlnk->seqname;
      nptr[m] = name + 32*m;
      strncpy(nptr[m], key, 32);
      size_t blksz = size[m++] = rlnk->counts->mx + 1;

      // Go through all the hashes to get the data.
      for (int i = 0 ; i < nhashes ; i++) {
         // In case no mock file was provided, the first
         // column of the data is replaced by 1s all along.
         if (no_mock && i == 0) {
            for (int k = 0 ; k < blksz ; k++) {
               y[offset + (nhashes * k)] = 1;
            }
            continue;
         }
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

   bitf_t *rep = malloc(sizeof(rod_t) + 32*sizeof(uint8_t));
   if (rep == NULL) {
      debug_print("%s", "memory error\n");
      free(new);
      return NULL;
   }

   // Initiatile count 'rod_t'.
   memset(rod->array, 0, 32*sizeof(uint32_t));
   rod->sz = 32;
   rod->mx = 0;

   // Initiatile repeat 'rod_t'.
   memset(rep->array, 0, 32*sizeof(uint8_t));
   rep->sz = 32*8;

   // Update 'link_t'.
   strncpy(new->seqname, s, 31);
   new->counts  = rod;
   new->repeats = rep;

   // Add 'link_t' to table hash table.
   new->next = htab[hv(s)];
   htab[hv(s)] = new;

   return new;
}

int
bitf_query_and_set
(
 int      pos,
 link_t * lnk
)
{
   // Extend bf if necessary.
   if (pos >= lnk->repeats->sz) {
      uint64_t newbeg = lnk->repeats->sz/8+1;
      uint64_t newend = pos/8+1;
      lnk->repeats = realloc(lnk->repeats, sizeof(rod_t)+newend*sizeof(uint8_t));
      if (lnk->repeats == NULL) {
         fprintf(stderr, "memory error\n");
         return -1;
      }
      memset(lnk->repeats->array + newbeg, 0, newend-newbeg+1);
      lnk->repeats->sz = newend*8;
   }
   // Check bit.
   int byte = pos/8;
   int bit  = pos%8;
   int set = (lnk->repeats->array[byte] >> bit) & 1;
   // Set bit.
   lnk->repeats->array[byte] |= (1 << bit);

   return set;
}

int
bloom_query_and_set
(
   const char  * name,
         int     pos,
         bloom_t bloom
)
// SYNOPSIS:                                                             
//   The Bloom filter stores the position of mapped reads in order to    
//   remove duplicates. When parsing the mapping data, a read is first   
//   looked up in the Bloom filter before beind added. If a match is     
//   found, the read is skipped. The Bloom filter has 800,000,000 bits   
//   and uses two hash functions. These days, a ChIP-seq lane is         
//   typically 200 million reads, filling a maximum of 400,000,000 or    
//   half of the Bloom filter, but at that stage the false positive      
//   rate is already close to 1/4 (meaning that only 3/4 of the reads    
//   will be included by the end of the parsing). Until approximately    
//   41,000,000 unique reads haveebeen queried, the false positive rate  
//   is under 1%.                                                        
//                                                                       
//   For a Bloom filter of size N, after inserting n distinct objects,   
//   the number of set bits is approximately equal to                    
//                                                                       
//                           N(1-exp(-2n/N)),                            
//                                                                       
//   which is accurate within 1%. The total number of false positives    
//   (here the number of lost reads) is approximately equal to           
//                                                                       
//              n - N(3/4 - exp(-2n/N)(1-exp(-2n/N)/4))                  
//                                                                       
//   which understimates the true value by about 4%. For a lane of 200   
//   million distinct reads, the amount of loss is approximately 6%.     
//                                                                       
//   The first hash function uses the efficient 'djb2()' to produce a    
//   32 bit integer from the chomosome name and adds the position of     
//   the read. The second hash function is computed from the first by    
//   adding again the position of the read. The scheme can be schema-    
//   tized as shown below for a read mapping on chr1 at position 8.      
//                                                                       
//       .|.|X|.|.|.|.|.|.|.|1|.|.|.|.|.|.|.|1|.|.|.|.|.|.|.|.|.|.|.     
//           |               |               |                           
//       djb2(chr1)         +8              +16                          
//                                                                       
{

   char text_a[48];
   char text_b[48];
   int len_a = sprintf(text_a, "%s%d", name, pos);
   int len_b = sprintf(text_b, "%d%s", pos, name);
   uint32_t a = XXH32(text_a, len_a, 0) % (8*BSIZE);
   uint32_t b = XXH32(text_b, len_b, 1) % (8*BSIZE);

   int abyte = a/8, abit = a%8;
   int bbyte = b/8, bbit = b%8;

   // Get bits.
   int a_yes = (bloom[abyte] >> (abit)) & 1;
   int b_yes = (bloom[bbyte] >> (bbit)) & 1;

   // Set bits.
   bloom[abyte] |= (1 << abit);
   bloom[bbyte] |= (1 << bbit);

   return a_yes && b_yes;

}
