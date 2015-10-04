#include "predict.h"
#include <stdio.h>
#include "parsam.h"
#include "parse.h"
#include "zerone.h"

#define has_map(a) strcmp(".map", (a) + strlen(a) - 4) == 0 || \
   strcmp(".map.gz", (a) + strlen(a) - 7) == 0 ||
#define has_sam(a) strcmp(".sam", (a) + strlen(a) - 4) == 0
#define has_bam(a) strcmp(".bam", (a) + strlen(a) - 4) == 0

int main(int argc, char **argv) {

   // Read files and build ChIP structure.
   int any_sam = 0;
   int all_sam = 1;
   // Are input files sam/bam?
   for (int i = 1; i < argc; i++) {
      if (has_sam(argv[i]) || has_bam(argv[i])) any_sam = 1;
      else all_sam = 0;
   }

   int any_gem = 0;
   int all_gem = 1;
   // Are input files gem?
   for (int i = 1; i < argc; i++) {
      if (has_map(argv[i])) any_gem = 1;
      else all_gem = 0;
   }

   ChIP_t *ChIP = NULL;
   if ((any_sam && !all_sam) || (any_gem && !all_gem)) {
      fprintf(stderr, "%s\n", "different file formats.");
      return 1;

   } else if (any_sam && all_sam) {
      char * samfiles[argc-1];
      const unsigned int nfiles = argc-1;
      for (int i = 1; i < argc; i++) samfiles[i-1] = argv[i];
      ChIP = read_sam(samfiles, nfiles);
      if (ChIP == NULL) {
         fprintf(stderr, "error wile reading input\n");
         return 1;
      }

   } else if (any_gem && all_gem) {
      char * gemfiles[argc-1];
      const unsigned int nfiles = argc-1;
      for (int i = 1; i < argc; i++) gemfiles[i-1] = argv[i];
      ChIP = read_gem((const char **) gemfiles, nfiles);
      if (ChIP == NULL) {
         fprintf(stderr, "error wile reading input\n");
         return 1;
      }

   } else if (argc == 2 && !any_sam) {
      FILE *inputf = fopen(argv[1], "r");
      if (inputf == NULL) {
         fprintf(stderr, "file not found: %s\n", argv[1]);
         return 1;
      }
      ChIP = read_file(inputf);
      fclose(inputf);
   }

   // Do zerone.
   const unsigned int m = 3; // number of states.
   zerone_t *zerone = do_zerone(ChIP);
   if (zerone == NULL) return 1;

   // XXX Broken on any computer of the planet except one. XXX //
   // XXX Better put the data directly in the file. XXX
   
//   char * centerfn = "/home/pcusco/Zerone/classifier/SVM_18x1_center.csv";
//   char * scalefn  = "/home/pcusco/Zerone/classifier/SVM_18x1_scale.csv";
//   char * svfn     = "/home/pcusco/Zerone/classifier/SVM_200x18_sv.csv";
//   char * coefsfn  = "/home/pcusco/Zerone/classifier/SVM_200x1_coefs.csv";

//   double * coefs  = readmatrix(coefsfn, NSV, 1);
//   double * center = readmatrix(centerfn, DIM, 1);
//   double * scale  = readmatrix(scalefn, DIM, 1);
//   double * sv     = readmatrix(svfn, NSV, DIM);
//   if (coefs == NULL || center == NULL || scale == NULL || sv == NULL) {
//      fprintf(stderr, "could not read SVM file\n");
//      return 1;
//   }
//
   // Quality control.

//   double * feat = extractfeats(ChIP, zerone);
//   double * sfeat = zscale(feat, center, scale);

//   double QCscore = predict(sfeat, sv, coefs);
//   fprintf(stdout, "# QC score: %.3f\n", QCscore);
//   fprintf(stdout, "# advice: %s discretization.\n",
//         QCscore >= 0 ? "accept" : "reject");
  
   // Print results (Viretbi path and phi matrix).
   for (size_t i = 0 ; i < nobs(ChIP) ; i++) {
// The commented lines below print the output together with
// the data used for discretization.
//      fprintf(stdout, "%d\t%f\t%f\t%f", zerone->path[i],
//            zerone->phi[0+i*m], zerone->phi[1+i*m], zerone->phi[2+i*m]);
//      for (int j = 0 ; j < ChIP->r ; j++) {
//         fprintf(stdout, "\t%d", ChIP->y[j+i*ChIP->r]);
//      }
//      fprintf(stdout, "\n");
      fprintf(stdout, "%d\t%f\t%f\t%f\n", zerone->path[i],
            zerone->phi[0+i*m], zerone->phi[1+i*m], zerone->phi[2+i*m]);
   }


//char * featsfn = "/home/pcusco/zerone/classifier/SVM_946x18_features.csv";
//char * labelsfn = "/home/pcusco/zerone/classifier/SVM_946x1_labels.csv";
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
//   }
//      fprintf(stdout, "%d:%d\n", (int)labels[i], p);
//}
//fprintf(stderr, "sum: %d\n", sum);

//   free(center);
//   free(scale);
//   free(sv);
//   free(coefs);
//   free(feat);
//   free(sfeat);

   destroy_zerone_all(zerone); // Also frees ChIP.

   return 0;

}
