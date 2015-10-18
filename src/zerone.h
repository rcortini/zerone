#include <stdint.h>
#include "hmm.h"
#include "utils.h"
#include "zinb.h"

#ifndef _JAHMM_HEADER
#define _JAHMM_HEADER

#define BW_MAXITER 100     // BW iterations //
#define BT_MAXITER 20      // Backtrack iterations //
#define TOLERANCE 1e-6

struct ChIP_t;
struct zerone_t;
struct block_t;

typedef unsigned int uint;
typedef struct ChIP_t ChIP_t;
typedef struct zerone_t zerone_t;
typedef struct block_t block_t;


struct ChIP_t {
   size_t  r;       // number of dimensions of 'y' //
   int     nb;      // number of blocks //
   int   * y;       // observations //
   char  * nm;      // block names //
   uint    sz[];    // block sizes //
};

struct zerone_t {
   int      m;      // number of states (3) //
   int      r;      // number of dimensions of 'y' //
   int      map[3]; // state map (map[2] is target) //
   ChIP_t * ChIP;   // observations //
   double * Q;      // transitions //
   double   a;      // emission par //
   double   pi;     // emission par //
   double * p;      // emission par //
   double * phi;    // posterior probs //
   double * pem;    // emission probs //
   double   l;      // log-likelihood //
   int    * path;   // Viterbi path //
   int      iter;   // number of BW iterations //
};

void       bw_zinm(zerone_t *);
void       destroy_zerone_all(zerone_t *);
zerone_t * do_zerone(ChIP_t *);
ChIP_t   * new_ChIP(uint, uint, int *, char **, const uint *);
zerone_t * new_zerone(uint, ChIP_t *);
uint       nobs(const ChIP_t *);
ChIP_t   * read_file(FILE *);
void       set_zerone_par(zerone_t *, const double *, double, double, const double *);
void       update_trans(size_t, double *, const double *);
void       zinm_prob(zerone_t *, const int *, int, double *);

#endif
