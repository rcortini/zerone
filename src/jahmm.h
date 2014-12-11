#include "hmm.h"
#include "utils.h"
#include "zinm.h"

#define MAXITER 500
#define JAHMM_MAXITER 25
#define TOLERANCE 1e-6

struct ChIP_t;
struct jahmm_t;

typedef struct ChIP_t ChIP_t;
typedef struct jahmm_t jahmm_t;
typedef unsigned int uint;


struct ChIP_t {
   size_t   r;      // number of dimensions of 'y' //
   uint     nb;     // number of blocks //
   int    * y;      // observations //
   uint     size[]; // size of the blocks //
};

struct jahmm_t {
   uint     m;    // number of states //
   ChIP_t * ChIP; // observations //
   double * Q;    // transitions //
   double   a;
   double   pi;   // emissions //
   double * p;
};

void      bw_zinm(jahmm_t *);
void      destroy_ChIP(ChIP_t *);
ChIP_t  * new_ChIP(uint, uint, int *, const uint *);
jahmm_t * new_jahmm(uint, ChIP_t *);
ChIP_t  * read_file(FILE *);
void      set_jahmm_par(jahmm_t *, const double *, double,
            double, const double *);
void      update_trans(size_t, double *, const double *);
void      zinm_prob(jahmm_t *, const int *, int, double *);
