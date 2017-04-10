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
#include <stdio.h>
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

// Snippets.
int check_strtoX (char *, char *);

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
   int    count;
   int    mapq;
};

// Iterator states.
struct bgzf_state_t {
   BGZF      * file;
   bam_hdr_t * hdr;
   bam1_t    * bam;
   int         minmapq;
};

struct generic_state_t {
   FILE      * file;
   reader_t    reader;
   parser_t    parser;
   char      * buff;
   size_t      bsz;
   int         minmapq;
};


//  ---- Declaration of local functions  ---- //
int      autoparse (const char *, hash_t *, zerone_parser_args_t);
iter_t   choose_iterator (const char *);

// Iterators.
int      bgzf_iterator (loc_t *);
int      generic_iterator (loc_t *);

// Readers.
ssize_t  getgzline (char **, size_t *, FILE *);
ssize_t  getline   (char **, size_t *, FILE *);
ssize_t  getdelim  (char **, size_t *, int, FILE *);

// Parsers.
int      parse_gem (loc_t *, char *);
int      parse_sam (loc_t *, char *);
int      parse_bed (loc_t *, char *);
int      parse_wig (loc_t *, char *);

// Hash handling functions.
int      add_to_rod (rod_t **, uint32_t, int);
int      bloom_query_and_set(const char *, int, bloom_t);
int      bitf_query_and_set (int, link_t *);
void     destroy_hash(hash_t *);
void     destroy_bitfields(hash_t *);
void     reset_bitfields(hash_t *);
link_t * lookup_or_insert (const char *, hash_t *);
ChIP_t * merge_hashes (hash_t **, int, int);

// Convenience functions.
uint32_t djb2 (const char *);
int      is_gzipped (FILE *);




//  -- Definitions of exported functions  --- //

ChIP_t *
parse_input_files
(
   char                 * mock_fnames[],
   char                 * ChIP_fnames[],
   zerone_parser_args_t   args
)
{

   // debug info //
   {
      for (int i = 0; mock_fnames[i] != NULL; i++) {
         debug_print("| mock file: %s\n", mock_fnames[i]);
      }
      for (int i = 0; ChIP_fnames[i] != NULL; i++) {
         debug_print("| ChIP file: %s\n", ChIP_fnames[i]);
      }
   }

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

      if(!autoparse(mock_fnames[i], hashtab[0], args)) {
         debug_print("%s", "autoparse failed\n");
         goto clean_and_return;
      }

      reset_bitfields(hashtab[0]);
   }

   destroy_bitfields(hashtab[0]);


   // Create distinct hash for each ChIP file.
   for (int i = 0; ChIP_fnames[i] != NULL; i++) {

      debug_print("%s %s\n", "autoparsing ChIP file", ChIP_fnames[i]);

      if (nhashes >= 512) {
         fprintf(stderr, "too many files\n");
         goto clean_and_return;
      }

      hashtab[nhashes] = calloc(HSIZE, sizeof(link_t *));
      if (hashtab[nhashes] == NULL) {
         debug_print("%s", "memory error\n");
         goto clean_and_return;
      }

      nhashes++;
      if (!autoparse(ChIP_fnames[i], hashtab[nhashes-1], args)) {
         debug_print("%s", "autoparse failed\n");
         goto clean_and_return;
      }

      destroy_bitfields(hashtab[nhashes-1]);

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
         fprintf(stderr, "cannot open file %s\n", fname);
         goto exit_bgzf_error;
      }

      state->hdr = bam_hdr_read(state->file);
      if (state->hdr == NULL) {
         fprintf(stderr, "cannot read header from file %s\n", fname);
         goto exit_bgzf_error;
      }

      state->bam = bam_init1();
      if (state->bam == NULL) {
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
   const char                 * fname,
         hash_t               * hashtab,
         zerone_parser_args_t   args
)
{

   // Unpack arguments
   int window = args.window;
   int minmapq = args.minmapq;

   // debug info //
   {
      debug_print("%s", "arguments:'\n");
      debug_print("| window: %d\n", window);
      debug_print("| minmapq: %d\n", minmapq);
   }

   int status = SUCCESS;

   loc_t loc = {0};

   // Find an iterator for the given file type.
   iter_t iterate = choose_iterator(fname);

   if (iterate == NULL) {
      debug_print("%s", "unrecognized iterator\n");
      status = FAILURE;
      goto clean_and_return;
   }

   ERR = 0;

   while (iterate(&loc) > 0) {
      // 'loc.name' is set to NULL for unmapped reads or
      // any read that needs to be ignored. Also ignore
      // reads with low quality.
      if (loc.name == NULL || loc.mapq < minmapq) continue;

      link_t *lnk = lookup_or_insert(loc.name, hashtab);

      if (lnk == NULL) {
         debug_print("%s", "hash query failed\n");
         status = FAILURE;
         goto clean_and_return;
      }

      // Check in bit field whether the read was seen before.
      if (bitf_query_and_set(loc.pos, lnk)) continue;

      // Add read to counts.
      if (!add_to_rod(&lnk->counts, loc.pos / window, loc.count)) {
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
   return status;

}


int
generic_iterator
(
   loc_t *loc
)
{

   static unsigned int lineno;

   // Cast as state for generic iterator.
   generic_state_t *state = (generic_state_t *) STATE;
   parser_t parse = state->parser;
   reader_t doread = state->reader;

   if (loc == NULL) {
      // Interruption by the caller.
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

   lineno++;
   if (!parse(loc, state->buff)) {
      // XXX non debug error XXX //
      fprintf(stderr, "format conflict in line %d:\n%s%s",
            lineno, buffer, sz > 63 ? "...\n" : "");
      ERR = __LINE__;
      goto clean_and_return;
   }

   return nbytes;

clean_and_return:
   free(state->buff);
   fclose(state->file);
   free(state);

   lineno = 0;

   STATE = NULL;
   return -1;

}

int
bgzf_iterator
(
   loc_t *loc
)
{

   static int n_parsed_header_targets;

   // Cast as state for bgzf iterator.
   bgzf_state_t *state = (bgzf_state_t *) STATE;
   bam_hdr_t *hdr = state->hdr;

   if (loc == NULL) {
      // Interruption by the caller.
      goto clean_and_return;
   }


   // First parse the header one item at a time.
   if (n_parsed_header_targets < hdr->n_targets) {
      loc->name = hdr->target_name[n_parsed_header_targets];
      loc->pos = hdr->target_len[n_parsed_header_targets];
      // Make sure this is not ignored when checking 'mapq'.
      loc->mapq = 2147483647;
      // Add 0 read. This will still set assocaited
      // 'rod_t' to the right length.
      loc->count = 0;
      n_parsed_header_targets++;
      return SUCCESS;
   }

   // Parse alignments.

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

#ifdef ATAC_SEQ
      // Discard unmapped and 2nd align of PE files.
      if (core.tid < 0)
         loc->name = NULL;
      else
         loc->name = hdr->target_name[core.tid];
      // Compute middle point for PE intervals.
      loc->pos = 1 + core.pos;
      loc->count = 1;
#else
      // Discard unmapped reads and 2nd align of PE files.
      // The read will be ignored by setting 'name' to NULL.
      if (core.tid < 0 || core.flag & BAM_FREAD2)
         loc->name = NULL;
      else
         loc->name = hdr->target_name[core.tid];

      // Note that the bam format is 0-based, so we add 1 to the
      // position because genomic positions are 1-based.
      // Compute middle point for PE intervals.
      if (core.flag & BAM_FPAIRED) {
         if (core.mtid != core.tid) loc->name = NULL;
         loc->pos = 1 + min(core.pos, core.mpos) + abs(core.isize)/2;
      }
      else {
         loc->pos = 1 + core.pos;
      }
      loc->mapq = core.qual;
      loc->count = 1;
#endif

   return bytesread;

clean_and_return:
   // This bit is executed when the iterator needs to be
   // cleaned (end of file, interruption or error).
   bam_hdr_destroy(hdr);
   bgzf_close(state->file);
   bam_destroy1(state->bam);
   free(state);

   n_parsed_header_targets = 0;

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
   // GEM has no mapping quality.
   loc->mapq = 2147483647;
   loc->count = 1;

   return SUCCESS;

}


int
parse_sam
(
   loc_t *loc,
   char  *line
)
{

   // Parse header.
   if (line[0] == '@') {
      if (strncmp(line, "@SQ", 3) != 0) {
         // Ignore non "@SQ" lines of the header.
         loc->name = NULL;
         return SUCCESS;
      }

      // Get chromosome lengths from "@SQ" lines.
                    strsep(&line, "\t"); // Discard "@SQ".
      char *chrom = strsep(&line, "\t");
      char *len   = strsep(&line, "\t");

      if (strncmp(chrom, "SN:", 3) != 0 ||
                  strncmp(len, "LN:", 3) != 0) {
         // Cannot parse. Skip and ignore.
         loc->name = NULL;
         return SUCCESS;
      }

      // Remove "SN:" and "LN:" field identifiers and
      // create a pseudo count at the last position of
      // the chromosome.
      loc->name = chrom + 3;
      loc->pos = atoi(len+3);
      // Make sure this is not ignored when checking 'mapq'.
      loc->mapq = 2147483647;

      return SUCCESS;

   }

   // Parse alignment.

                 strsep(&line, "\t"); // Discard QNAME.
   char *flag  = strsep(&line, "\t");
   char *chrom = strsep(&line, "\t");
   char *Xpos  = strsep(&line, "\t");
   char *Xmapq = strsep(&line, "\t");

   // Cannot find chromosome or position.
   if (chrom == NULL || Xpos == NULL || Xmapq == NULL)
      return FAILURE;

   // Unmapped read.
   if (strcmp(chrom, "*") == 0) {
      loc->name = NULL;
      return SUCCESS;
   }

   // The mapping quality can be 0. Use 'strtoul()' to
   // convert the input to an integer because it allows
   // error-checking.
   char *endptr = NULL;
   errno = 0;
   int mapq = strtoul(Xmapq, &endptr, 10);
   if (!check_strtoX(Xmapq, endptr))
      return FAILURE;

   // Note that the sam format is 1-based, so if 'atoi()'
   // has returned 0, something is wrong with the format.
   int pos = atoi(Xpos);
   if (pos == 0) return FAILURE;

#ifndef ATAC_SEQ
   // Binary sam flags.
   int blag = atoi(flag);

   // Always skip read 2.
   if (blag & BAM_FREAD2) {
      loc->name = NULL;
      return SUCCESS;
   }

   // If read is paired use mid-point of the mapping.
   if (blag & BAM_FPAIRED) {
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
#endif

   loc->name = chrom;
   loc->pos = pos;
   loc->mapq = mapq;
   loc->count = 1;

   return SUCCESS;

}


int
parse_bed
(
   loc_t *loc,
   char  *line
)
{

   char *chrom  = strsep(&line, "\t"); // chrom
   char *tmp1   = strsep(&line, "\t"); // chromStart
   char *tmp2   = strsep(&line, "\t"); // chromEnd
                  strsep(&line, "\t"); // name
   char *Xcount = strsep(&line, "\t"); // score

   // Cannot find chromosome, position or count.
   if (chrom == NULL || tmp1 == NULL || tmp2 == NULL || Xcount == NULL)
      return FAILURE;

   // The bed format is 0-based. We cannot use 'atoi()' to
   // parse the coodrinates.
   char *endptr = NULL;

   // 'strtoul' may set 'errno' in case of overflow.
   errno = 0;

   int start = 1 + strtoul(tmp1, &endptr, 10);
   if (!check_strtoX(tmp1, endptr) || errno)
      return FAILURE;

   int end = 1 + strtoul(tmp2, &endptr, 10);
   if (!check_strtoX(tmp2, endptr) || errno)
      return FAILURE;

   int count = strtoul(Xcount, &endptr, 10);
   if (!check_strtoX(Xcount, endptr) || errno)
      return FAILURE;

   loc->name = chrom;
   loc->pos = (start + end) / 2;
   // BED has no mapping quality.
   loc->mapq = 2147483647;
   loc->count = count;

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

      int start;
      char *endptr = NULL;

      if (fixedstep) {
         // Line is just a count.
         start = fstart + iter++ * step;
      }
      else {
         // Line is a position and a count.
         char *Xstart = strsep(&line, "\t");

         // Cannot parse if position cannot be found
         // or if it is the last token of the line.
         if (Xstart == NULL || line == NULL)
            return FAILURE;

         // 'strtoul' may set 'errno' in case of overflow.
         errno = 0;

         start = strtoul(Xstart, &endptr, 10);
         if (!check_strtoX(Xstart, endptr) || errno)
            return FAILURE;
      }

      // Now only count is left on 'line'.
      // 'strtoul' may set 'errno' in case of overflow.
      errno = 0;

      int count = strtoul(line, &endptr, 10);
      if (!check_strtoX(line, endptr) || errno)
         return FAILURE;

      loc->pos = (2*start + span - 1) / 2;
      // WIG has no mapping quality.
      loc->mapq = 2147483647;
      loc->count = count;

      return SUCCESS;

   }

}


//  ---------  Utility functions  ----------  //

uint32_t
djb2
(
   const char * s
)
// The reference for djb2 (http://www.cse.yorku.ca/~oz/hash.html).
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

void
reset_bitfields
(
 hash_t * hashtab
)
{
   for (int i = 0 ; i < HSIZE ; i++) {
      for (link_t *lnk = hashtab[i] ; lnk != NULL ; lnk = lnk->next) {
         if (lnk->repeats == NULL) {
            lnk->repeats = malloc(sizeof(bitf_t)+32*sizeof(uint8_t));
            lnk->repeats->sz = 32*8;
         }
         memset(lnk->repeats->array, 0, lnk->repeats->sz/8);
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


ssize_t
getline
(
   char   ** lineptr,
   size_t *  n,
   FILE   *  stream
)
{
   return getdelim(lineptr, n, '\n', stream);
}

/* Default value for line length.  */
static const int line_size = 128;

ssize_t
getdelim
(
   char   ** lineptr,
   size_t *  n,
   int       delim,
   FILE   *  stream
)
{
  int indx = 0;
  int c;

  /* Sanity checks.  */
  if (lineptr == NULL || n == NULL || stream == NULL)
    return -1;

  /* Allocate the line the first time.  */
  if (*lineptr == NULL)
    {
      *lineptr = malloc (line_size);
      if (*lineptr == NULL)
       return -1;
      *n = line_size;
    }

  while ((c = getc (stream)) != EOF)
    {
      /* Check if more memory is needed.  */
      if (indx >= *n)
       {
         *lineptr = realloc (*lineptr, *n + line_size);
         if (*lineptr == NULL)
           return -1;
         *n += line_size;
      }

      /* Push the result in the line.  */
      (*lineptr)[indx++] = c;

      /* Bail out.  */
      if (c == delim)
       break;
   }

  /* Make room for the null character.  */
  if (indx >= *n)
    {
      *lineptr = realloc (*lineptr, *n + line_size);
      if (*lineptr == NULL)
       return -1;
      *n += line_size;
      }

  /* Null terminate the buffer.  */
  (*lineptr)[indx++] = 0;

  /* The last line may not have the delimiter, we have to
   * return what we got and the error will be seen on the
   * next iteration.  */
  return (c == EOF && (indx - 1) == 0) ? -1 : indx - 1;
}

int
add_to_rod
(
   rod_t   ** addr,
   uint32_t   pos,
   int        count
)
{

   // Convenience variable.
   rod_t *rod = *addr;

   if (pos >= rod->sz) {
      // At least double the size of the rod.
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

   rod->array[pos] += count;
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

   bitf_t *rep = malloc(sizeof(bitf_t) + 32*sizeof(uint8_t));
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
      uint64_t newbeg = lnk->repeats->sz/8;
      uint64_t newend = pos/8+1;
      lnk->repeats = realloc(lnk->repeats, sizeof(bitf_t)+newend*sizeof(uint8_t));
      if (lnk->repeats == NULL) {
         fprintf(stderr, "memory error\n");
         return -1;
      }
      memset(lnk->repeats->array + newbeg, 0, newend-newbeg);
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
