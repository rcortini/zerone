#include <stdio.h>
#include "jahmm.h"

int main(int argc, char **argv) {

   FILE *inputf = fopen(argv[1], "r");
   if (inputf == NULL) {
      fprintf(stderr, "file not found: %s\n", argv[1]);
      return 1;
   }

   ChIP_t *ChIP = read_file(inputf);
   fclose(inputf);

   const unsigned int m = 2; // number of states.
   jahmm_t *jahmm = do_jahmm(m, ChIP);
   if (jahmm == NULL) {
      fprintf(stderr, "now work...\n");
      return 1;
   }

   for (size_t i = 0 ; i < nobs(ChIP) ; i++) {
      fprintf(stdout, "%d\t%f\t%f\n", jahmm->path[i],
            jahmm->phi[0+i*m], jahmm->phi[1+i*m]);
   }

   return 0;

}
