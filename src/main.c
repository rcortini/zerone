/*
** Copyright 2014 Guillaume Filion, Eduard Valera Zorita and Pol Cusco.
**
** File authors:
**  Guillaume Filion     (guillaume.filion@gmail.com)
**  Eduard Valera Zorita (polcusco@gmail.com)
**
** License: 
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
*/

#include <stdio.h>
#include <getopt.h>
#include "predict.h"
#include "parsam.h"
#include "parse.h"
#include "zerone.h"

#define VERSION "zerone-v1.0"

#define MAXNARGS 255

#define has_map(a) strcmp(".map", (a) + strlen(a) - 4) == 0 || \
   strcmp(".map.gz", (a) + strlen(a) - 7) == 0
#define has_sam(a) strcmp(".sam", (a) + strlen(a) - 4) == 0
#define has_bam(a) strcmp(".bam", (a) + strlen(a) - 4) == 0


const char *USAGE =
"\n"
"Usage:"
"  zerone [options] <input file 1> ... <input file n>\n"
"\n"
"    -0 --baseline <input file> given file is a baseline control\n"
"    -v --version: display version and exit\n";

void say_usage(void) { fprintf(stderr, "%s\n", USAGE); }
void say_version(void) { fprintf(stderr, VERSION "\n"); }


//  -------  Definitions of local (helper) functions  ------- //

void
parse_fname
(
   char ** names,
   char *  value,
   int  *  offset
)
// Several file names may be passed as comma-separated values.
{
   

   char *token;
   while ((token = strsep(&value, ",")) != NULL) {
      if (*offset >= MAXNARGS) {
         fprintf(stderr, "too many arguments\n");
         exit(EXIT_FAILURE);
      }
      names[*offset] = strndup(token, 256);
      if (names[*offset] == NULL) {
         fprintf(stderr, "memory error\n");
         exit(EXIT_FAILURE);
      }
      (*offset)++;
   }

   return;

}


int main(int argc, char **argv) {

   // Input file names (mock and ChIP).
   char *mock_fnames[MAXNARGS] = {0};
   char *ChIP_fnames[MAXNARGS] = {0};

   int n_mock_files = 0;
   int n_ChIP_files = 0;

   while(1) {
      int option_index = 0;
      static struct option long_options[] = {
         {"mock",     required_argument,  0, '0'},
         {"help",     no_argument,        0, 'h'},
         {"version",  no_argument,        0, 'v'},
         {0, 0, 0, 0}
      };

      int c = getopt_long(argc, argv, "0:hv",
            long_options, &option_index);

      // Done parsing //
      if (c == -1) break;

      switch (c) {
      case '0':
         parse_fname(mock_fnames, optarg, &n_mock_files);
         break;

      case 'h':
         say_usage();
         return EXIT_SUCCESS;

      case 'v':
         say_version();
         return EXIT_SUCCESS;

      default:
         // Cannot parse. //
         say_usage();
         return EXIT_FAILURE;

      }

   }

   // Now parse positional arguments (file names).
   while (optind < argc) {
      parse_fname(ChIP_fnames, argv[optind++], &n_ChIP_files);
   }


   fprintf(stderr, "ChIP files:\n");
   for (int i = 0 ; i < n_ChIP_files ; i++) {
      fprintf(stderr, "%s\n", ChIP_fnames[i]);
   }
   fprintf(stderr, "\nmock files:\n");
   for (int i = 0 ; i < n_mock_files ; i++) {
      fprintf(stderr, "%s\n", mock_fnames[i]);
   }
   fprintf(stderr, "\n");

   // Read files and build ChIP structure.
   int any_sam = 0;
   int all_sam = 1;
   // Are input files sam/bam?
   for (int i = 1 ; i < argc ; i++) {
      if (has_sam(argv[i]) || has_bam(argv[i])) any_sam = 1;
      else all_sam = 0;
   }

   int any_gem = 0;
   int all_gem = 1;
   // Are input files gem?
   for (int i = 1 ; i < argc ; i++) {
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

   for (int i = 0 ; i < MAXNARGS ; i++) {
      free(mock_fnames[i]);
      free(ChIP_fnames[i]);
   }

   return 0;

}
