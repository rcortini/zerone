#include "predict.h"
#include <stdio.h>
#include "samread.h"
#include "jahmm.h"

int main(int argc, char **argv) {

   int any_sam = 0, all_sam = 0;
   for (int i = 1; i < argc; i++) {
      if (is_sam(argv[i])) {
         any_sam = 1;
         all_sam += 1;
      }
   }

   ChIP_t *ChIP = NULL;
   if (any_sam && !(all_sam == argc-1)) {
      fprintf(stderr, "%s\n", "different file formats.");
      return 1;

   } else if (any_sam && (all_sam == argc-1)) {
      int nfiles = 0;
      char * samfiles[argc-1];
      for (int i = 1; i < argc; i++) samfiles[nfiles++] = argv[i];
      ChIP = read_sam(nfiles, samfiles);

   } else if (argc == 2 && !any_sam) {
      FILE *inputf = fopen(argv[1], "r");
      if (inputf == NULL) {
         fprintf(stderr, "file not found: %s\n", argv[1]);
         return 1;
      }
      ChIP = read_file(inputf);
      fclose(inputf);
   }

   const unsigned int m = 3; // number of states.
   jahmm_t *jahmm = do_jahmm(m, ChIP);
   if (jahmm == NULL) {
      fprintf(stderr, "now work...\n"); // why now work?
      return 1;
   }

   for (size_t i = 0 ; i < nobs(ChIP) ; i++) {
      fprintf(stdout, "%d\t%f\t%f\n", jahmm->path[i],
            jahmm->phi[0+i*m], jahmm->phi[1+i*m]);
   }

   char * centerfn = "/home/pcusco/jahmm/classifier/SVM_19x1_center.csv";
   char * scalefn  = "/home/pcusco/jahmm/classifier/SVM_19x1_scale.csv";
   char * svfn     = "/home/pcusco/jahmm/classifier/SVM_129x19_sv.csv";
   char * coefsfn  = "/home/pcusco/jahmm/classifier/SVM_129x1_coefs.csv";

   double * coefs  = readmatrix(coefsfn, NSV, 1);
   double * center = readmatrix(centerfn, DIM, 1);
   double * scale  = readmatrix(scalefn, DIM, 1);
   double * sv     = readmatrix(svfn, NSV, DIM);

   double * feat = extractfeats(ChIP, jahmm);
   double * sfeat = zscale(feat, center, scale);

   fprintf(stdout, "%d\n", predict(sfeat, sv, coefs));

   free(center);
   free(scale);
   free(sv);
   free(coefs);
   free(feat);
   free(sfeat);

   return 0;

}
