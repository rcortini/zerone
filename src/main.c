#include "predict.h"
#include <stdio.h>
#include "samread.h"
#include "jahmm.h"

int main(int argc, char **argv) {

   // Read files and build ChIP structure.
   int any_sam = 0;
   int all_sam = 1;
   for (int i = 1; i < argc; i++) {
      if (is_sam(argv[i])) any_sam = 1;
      else all_sam = 0;
   }

   ChIP_t *ChIP = NULL;
   if (any_sam && !all_sam) {
      fprintf(stderr, "%s\n", "different file formats.");
      return 1;

   } else if (any_sam && all_sam) {
      char * samfiles[argc-1];
      unsigned int nfiles = 0;
      for (int i = 1; i < argc; i++) samfiles[nfiles++] = argv[i];
      ChIP = read_sam(samfiles, nfiles);

   } else if (argc == 2 && !any_sam) {
      FILE *inputf = fopen(argv[1], "r");
      if (inputf == NULL) {
         fprintf(stderr, "file not found: %s\n", argv[1]);
         return 1;
      }
      ChIP = read_file(inputf);
      fclose(inputf);
   }

   // Do jahmm.
   const unsigned int m = 3; // number of states.
   jahmm_t *jahmm = do_jahmm(m, ChIP);
   if (jahmm == NULL) return 1;

    Print results (Viretbi path and phi matrix).
   for (size_t i = 0 ; i < nobs(ChIP) ; i++) {
      fprintf(stdout, "%d\t%f\t%f\t%f\n", jahmm->path[i],
            jahmm->phi[0+i*m], jahmm->phi[1+i*m], jahmm->phi[2+i*m]);
   }

   // Quality control.
   char * centerfn = "/home/pcusco/jahmm/classifier/SVM_18x1_center.csv";
   char * scalefn  = "/home/pcusco/jahmm/classifier/SVM_18x1_scale.csv";
   char * svfn     = "/home/pcusco/jahmm/classifier/SVM_200x18_sv.csv";
   char * coefsfn  = "/home/pcusco/jahmm/classifier/SVM_200x1_coefs.csv";

   double * coefs  = readmatrix(coefsfn, NSV, 1);
   double * center = readmatrix(centerfn, DIM, 1);
   double * scale  = readmatrix(scalefn, DIM, 1);
   double * sv     = readmatrix(svfn, NSV, DIM);

   double * feat = extractfeats(ChIP, jahmm);
   double * sfeat = zscale(feat, center, scale);

   char * advice = predict(sfeat, sv, coefs) >= 0 ? "accept" : "reject";
   fprintf(stderr, "advice: %s discretization.\n", advice);

//char * featsfn = "/home/pcusco/jahmm/classifier/SVM_946x18_features.csv";
//char * labelsfn = "/home/pcusco/jahmm/classifier/SVM_946x1_labels.csv";
//double * feats = readmatrix(featsfn, 946, DIM);
//double * labels = readmatrix(labelsfn, 946, 1);
//
//int sum = 0;
//for (int i = 0; i < 946; i++) {
//   double * sfeat = zscale(&feats[i*DIM], center, scale);
//   int p = predict(sfeat, sv, coefs);
//   free(sfeat);
//   if (p == -1) p = 0;
//   if (p != (int)labels[i]) {
//      sum++;
//      fprintf(stderr, "%d:%d\n", (int)labels[i], p);
//   }
//}
//fprintf(stderr, "sum: %d\n", sum);

   free(center);
   free(scale);
   free(sv);
   free(coefs);
   free(feat);
   free(sfeat);
   
   destroy_jahmm_all(jahmm); // Also frees ChIP.

   return 0;
}
