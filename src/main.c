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
#include "parse.h"
#include "zerone.h"


const char *USAGE =
"\n"
"Usage:"
"  zerone [options] <input file 1> ... <input file n>\n"
"\n"
"    -0 --baseline <input file> given file is a baseline control\n"
"    -v --version: display version and exit\n";



#define VERSION "zerone-v1.0"
#define MAXNARGS 255


// Useful macros.
#define has_map(a) strcmp(".map", (a) + strlen(a) - 4) == 0 || \
   strcmp(".map.gz", (a) + strlen(a) - 7) == 0
#define has_sam(a) strcmp(".sam", (a) + strlen(a) - 4) == 0
#define has_bam(a) strcmp(".bam", (a) + strlen(a) - 4) == 0


typedef ChIP_t * (*parser_t) (char **, char **);



//  -----------  Definitions of local one-liners  ----------- //
void say_usage(void) { fprintf(stderr, "%s\n", USAGE); }
void say_version(void) { fprintf(stderr, VERSION "\n"); }


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
   char *mock_fnames[MAXNARGS+1] = {0};
   char *ChIP_fnames[MAXNARGS+1] = {0};

   int n_mock_files = 0;
   int n_ChIP_files = 0;

   // Parse options.
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

      // Done parsing named options. //
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

   // Process input files.
   ChIP_t *ChIP = parse_input_files(mock_fnames, ChIP_fnames);

   if (ChIP == NULL) {
      fprintf(stderr, "error wile reading input\n");
      exit(EXIT_FAILURE);
   }

   // Do zerone.
   //const unsigned int m = 3; // number of states.
   zerone_t *zerone = do_zerone(ChIP);

   if (zerone == NULL) {
      fprintf(stderr, "run time error (sorry)\n");
      exit(EXIT_FAILURE);
   }

   // Quality control.
   double QCscore = predict(ChIP, zerone);
   fprintf(stdout, "# QC score: %.3f\n", QCscore);
   fprintf(stdout, "# advice: %s discretization.\n",
         QCscore >= 0 ? "accept" : "reject");

   // Print results (Viterbi path and phi matrix).
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
            zerone->phi[0+i*3], zerone->phi[1+i*3], zerone->phi[2+i*3]);
   }

   destroy_zerone_all(zerone); // Also frees ChIP.

   for (int i = 0 ; i < MAXNARGS ; i++) {
      free(mock_fnames[i]);
      free(ChIP_fnames[i]);
   }

   return 0;

}
