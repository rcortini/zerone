#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "zerone.h"

#define SUCCESS 1
#define FAILURE 0

// XXX This has to become a function parameter XXX //
#define BIN_SIZE 300

// Prime number for the size of hash table.
#define HSIZE 997

// Shortcut to acces position of hash table.
#define h(a) djb2(a) % HSIZE

// Type declarations.
struct hash_t;
struct link_t;
struct rod_t;

typedef struct link_t link_t;
typedef struct rod_t rod_t;
typedef link_t * hash_t;


// Type definitions.
struct link_t {
   char     seqname[16]; // Chromosome or sequence name.
   rod_t  * counts;      // Counts in chromosome bins.
   link_t * next;        // Next node on the link list.
};

struct rod_t {
   size_t sz; 
   size_t mx;
   uint32_t array[];
};


//  ---- Declaration of local functions  ---- //

int add (rod_t **, uint32_t);
void destroy_hash(hash_t *);
uint32_t djb2 (const char *);
link_t * lookup_or_insert (const char *, hash_t *);
ChIP_t * merge_hashes(hash_t **, int);
hash_t * read_and_count (const char *);


//  -- Definitions of exported functions  --- //

ChIP_t *
read_gem
(
   const    char * fn[],
   unsigned int    nfiles
)
{

   ChIP_t *ChIP = NULL; // Return value.
   // Array of hash tables (one per file).
   hash_t *hashtab[nfiles];
   memset(hashtab, 0, nfiles * sizeof(hash_t *));

   for (int i = 0; i < nfiles; i++) {
      hashtab[i] = read_and_count(fn[i]);
      if (hashtab[i] == NULL) {
         fprintf(stderr, "error in function '%s()' %s:%d\n",
               __func__, __FILE__, __LINE__);
         goto clean_and_return;
      }
   }

   // Merge hash tables. 
   ChIP = merge_hashes(hashtab, nfiles);
   if (ChIP == NULL) {
      fprintf(stderr, "error in function '%s()' %s:%d\n",
               __func__, __FILE__, __LINE__);
      goto clean_and_return;
   }

clean_and_return:
   for (int i = 0 ; i < nfiles ; i++) {
      if (hashtab[i] != NULL) destroy_hash(hashtab[i]);
   }

   return ChIP;

}



//  ---- Definitions of local functions  ---- //

int
add
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
         fprintf(stderr, "memory error in function '%s()' %s:%d\n",
               __func__, __FILE__, __LINE__);
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
   int       n
)
{
   // Add keys and update a reference hash table.
   hash_t *refhash = hashes[0];

   for (int i = 1 ; i < n ; i++) {
      hash_t *hashtab = hashes[i];
      for (int j = 0 ; j < HSIZE ; j++) {
         for(link_t *lnk = hashtab[j] ; lnk != NULL ; lnk = lnk->next) {
            // Get chromosome from reference hash (or create it).
            link_t *reflnk = lookup_or_insert(lnk->seqname, refhash);
            if (reflnk == NULL) {
               fprintf(stderr, "error in function '%s()' %s:%d\n",
                     __func__, __FILE__, __LINE__);
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
   int *y = calloc(n*nbins, sizeof(int));
   if (size == NULL || y == NULL) {
      fprintf(stderr, "memory error in funtion '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return NULL;
   }

   // Fill in the observations.
   int m = 0;
   size_t offset = 0;
   for (int j = 0 ; j < HSIZE ; j++) {
   for (link_t *rlnk = refhash[j] ; rlnk != NULL ; rlnk = rlnk->next) {

      // Get seqname and sequence length.
      char *key = rlnk->seqname;
      size_t blksz = size[m++] = rlnk->counts->mx;

      // Go through all the hashes to get the data.
      for (int i = 0 ; i < n ; i++) {
         hash_t *hashtab = hashes[i];
         link_t *lnk = lookup_or_insert(key, hashtab);
         if (lnk == NULL) {
            fprintf(stderr, "error in function '%s()' %s:%d\n",
                  __func__, __FILE__, __LINE__);
            return NULL;
         }
         rod_t *counts = lnk->counts;
         size_t kmax = counts->sz < blksz ? counts->sz : blksz;
         for (int k = 0 ; k < kmax ; k++) {
            y[offset + (n*k) + i] = counts->array[k];
         }
      }

      // Update offset.
      offset += n*blksz;

   }
   }

   ChIP_t *ChIP = new_ChIP(n, nkeys, y, size);
   free(size);

   return ChIP;

}

hash_t *
read_and_count
(
   const char * fn
)
{

   int status = SUCCESS;
   FILE *fp = NULL;

   // Create new hash for each file (it is the return value).
   hash_t * hashtab = calloc(HSIZE, sizeof(link_t *));
   if (hashtab == NULL) {
      fprintf(stderr, "memory error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      status = FAILURE;
      goto clean_and_return;
   }

   fp = fopen(fn, "r");
   if (fp == NULL) {
      fprintf(stderr, "cannot open file %s\n", fn);
      status = FAILURE;
      goto clean_and_return;
   }

   size_t bsz = 64;
   char *buff = malloc(bsz * sizeof(char));
   if (buff == NULL) {
      fprintf(stderr, "memory error\n");
      status = FAILURE;
      goto clean_and_return;
   }

   ssize_t bytesread;
   while ((bytesread = getline(&buff, &bsz, fp)) > -1) {
      // Get last field of the line.
      char *map = strrchr(buff, '\t');
      if (map++ == NULL) {
         fprintf(stderr, "format conflict in line:\n%s", buff);
         status = FAILURE;
         goto clean_and_return;
      }

      if (map[0] == '-') continue;

      // Parse mapping data.
      char *chrom = strsep(&map, ":");
      strsep(&map, ":"); // Discard strand.
      char *tmp = strsep(&map, ":");

      if (chrom == NULL || tmp == NULL) {
         fprintf(stderr, "format conflict in line:\n%s", buff);
         status = FAILURE;
         goto clean_and_return;
      }

      // Positions in the genome cannot be 0, so we can identify
      // failures of 'atoi' to convert numbers.
      int pos = atoi(tmp);
      if (pos == 0) {
         fprintf(stderr, "format conflict in line:\n%s", buff);
         status = FAILURE;
         goto clean_and_return;
      }

      link_t *lnk = lookup_or_insert(chrom, hashtab);
      if (lnk == NULL) {
         fprintf(stderr, "error in function '%s()' %s:%d\n",
               __func__, __FILE__, __LINE__);
         status = FAILURE;
         goto clean_and_return;
      }

      // Add read to counts.
      if (!add(&lnk->counts, pos / BIN_SIZE)) {
         fprintf(stderr, "error in function '%s()' %s:%d\n",
               __func__, __FILE__, __LINE__);
         status = FAILURE;
         goto clean_and_return;
      }

   }

clean_and_return:
   if (buff != NULL) free(buff);
   if (fp != NULL) fclose(fp);

   if (status == FAILURE) {
      destroy_hash(hashtab);
      hashtab = NULL;
   }

   return hashtab;

}


link_t *
lookup_or_insert
(
   const char   * s,
         hash_t * hashtab
)
// Famous K&R simple hash lookup.
{

   for (link_t *lnk = hashtab[h(s)] ; lnk != NULL ; lnk = lnk->next) {
      if (strcmp(s, lnk->seqname) == 0)  return lnk;
   }

   // Entry was not found.
   link_t *new = malloc(sizeof(link_t));

   if (new == NULL) {
      fprintf(stderr, "memory error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      return NULL;
   }

   rod_t  *rod = malloc(sizeof(rod_t) + 32*sizeof(uint32_t));
   if (rod == NULL) {
      fprintf(stderr, "memory error in function '%s()' %s:%d\n",
            __func__, __FILE__, __LINE__);
      free(new);
      return NULL;
   }

   // Initiatile 'rod_t'.
   memset(rod->array, 0, 32*sizeof(uint32_t));
   rod->sz = 32;
   rod->mx = 0;

   // Update 'link_t'.
   strncpy(new->seqname, s, 15);
   new->counts = rod;

   // Add 'link_t' to table hash table.
   new->next = hashtab[h(s)];
   hashtab[h(s)] = new;

   return new;

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
   for (int i = 0 ; i < 15 && s[i] != '\0' ; i++) 
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
