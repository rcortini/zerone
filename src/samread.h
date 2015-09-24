#include "zerone.h"

#ifndef _PARSE_H
#define _PARSE_H
#define _GNU_SOURCE

#define BIN_SIZE 300

struct counter_t {
   int32_t    n_ref;    // # of chromosomes
   int32_t  * n_bins;   // # of bins
   int32_t ** bins;     // read count per bin
};

typedef struct counter_t counter_t;

int         is_sam(const char * fn);
ChIP_t    * read_sam(char *fn[], unsigned int nfiles);
counter_t * read_count(const char * fn);
void        destroy_counter(counter_t * counter);

#endif
