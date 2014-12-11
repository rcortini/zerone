#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <err.h>

#ifndef _HMM_HEADER_
#define _HMM_HEADER_

struct hmm_p_t;

typedef unsigned int uint;
typedef struct hmm_p_t hmm_p_t;

struct hmm_p_t {
   uint     m;     // number of states //
   void   * par;   // emission parameters //
   double   Q[];   // transition parameters //
};

#define D double
#define H hmm_p_t
#define I int
#define U uint
#define V void
#define cD const double
#define cU const uint

//  function      ( 1    2   3    4    5    6   7   8   9  )
V   block_fwdb    (  U,  U, cU*,  D*,  D*,  D*, D*, D*, D* );
I   block_viterbi ( cU, cU, cU*, cD*, cD*, cD*, U,  U*     );
V   bwd           (  U,  U, cD*,  D*,  D*,  D*             );
V   destroy_hmm_p (  H*                                    );
D   fwd           (  U,  U, cD*, cD*,  D*                  );
D   fwdb          (  U,  U, cD*, cD*,  D*,  D*, D*         );
H * new_hmm_p     (  U                                     );
U * viterbi       (  U,  U, cD*, cD*,  cD*, U*             );

#undef D
#undef H
#undef U
#undef V
#undef cD
#undef cU

#endif
