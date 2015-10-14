/*
** Copyright 2014 Guillaume Filion, Eduard Valera Zorita and Pol Cusco.
**
** File authors:
**  Pol Cuscó Pons    (polcusco@gmail.com)
**  Guillaume Filion  (guillaume.filion@gmail.com)
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
#include "debug.h"
#include "parse.h"
#include "predict.h"
#include "zerone.h"


const char *USAGE =
"\n"
"Usage:"
"  zerone [options] <input file 1> ... <input file n>\n"
"\n"
"    -0 --mock: given file is a mock control\n"
"    -1 --chip: given file is a ChIP-seq experiment\n"
"    -l --list-output: output list of targets (default table)\n"
"\n"
"    -h --help: display this message and exit\n"
"    -v --version: display version and exit\n";



#define VERSION "zerone-v1.0"
#define MAXNARGS 255



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

   if (argc == 1) {
      say_usage();
      exit(EXIT_SUCCESS);
   }

   // Input file names (mock and ChIP).
   char *mock_fnames[MAXNARGS+1] = {0};
   char *ChIP_fnames[MAXNARGS+1] = {0};

   int n_mock_files = 0;
   int n_ChIP_files = 0;

   static int list_flag = 0;

   // Parse options.
   while(1) {
      int option_index = 0;
      static struct option long_options[] = {
         {"list-output", no_argument,       &list_flag, 1 },
         {"mock",        required_argument,         0, '0'},
         {"chip",        required_argument,         0, '1'},
         {"help",        no_argument,               0, 'h'},
         {"version",     no_argument,               0, 'v'},
         {0, 0, 0, 0}
      };

      int c = getopt_long(argc, argv, "0:1:hv",
            long_options, &option_index);

      // Done parsing named options. //
      if (c == -1) break;

      switch (c) {
      case 0:
         // A flag was set //
         break;

      case '0':
         parse_fname(mock_fnames, optarg, &n_mock_files);
         break;

      case '1':
         parse_fname(ChIP_fnames, optarg, &n_ChIP_files);
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

   debug_print("%s", "done parsing arguments\n");

   // Process input files.
   ChIP_t *ChIP = parse_input_files(mock_fnames, ChIP_fnames);

   if (ChIP == NULL) {
      fprintf(stderr, "error while reading input\n");
      exit(EXIT_FAILURE);
   }

   debug_print("%s", "done reading input files\n");

   // Do zerone.
   //const unsigned int m = 3; // number of states.
   zerone_t *Z = do_zerone(ChIP);

   if (Z == NULL) {
      fprintf(stderr, "run time error (sorry)\n");
      exit(EXIT_FAILURE);
   }

   // Quality control.
   double QCscore = zerone_predict(Z);
   fprintf(stdout, "# QC score: %.3f\n", QCscore);
   fprintf(stdout, "# advice: %s discretization.\n",
         QCscore >= 0 ? "accept" : "reject");

   // Print results (Viretbi path and phi matrix).
//   const int target = Z->map[2];
//   const int inter  = ChIP->map[1];
//   const int ground = ChIP->map[0];

   if (list_flag) {
      int wid = 0;
      int target = 0;
      for (int i = 0 ; i < ChIP->nb ; i++) {
         char *name = ChIP->nm + 32*i;
         for (int j = 0 ; j < ChIP->sz[i] ; j++) {
            if (!target && Z->path[wid] == target) {
               fprintf(stdout, "%s\t%d\t", name, 300*j + 1);
               target = 1;
            }
            else if (target && Z->path[wid] != target) {
               fprintf(stdout, "%d\n", 300*(j+1));
               target = 0;
            }
         }
         if (target) {
            fprintf(stdout, "%d\n", 300 * ChIP->sz[i]);
            target = 0;
         }
      }
   }

   else {
      int wid = 0;
      for (int i = 0 ; i < ChIP->nb ; i++) {
         char *name = ChIP->nm + 32*i;
         for (int j = 0 ; j < ChIP->sz[i] ; j++) {
            fprintf(stdout, "%s\t%d\t%d\t%d\n", name, 300*j + 1,
                  300*(j+1), Z->path[wid++]);
         }
      }
   }

   destroy_zerone_all(Z); // Also frees ChIP.

   for (int i = 0 ; i < MAXNARGS ; i++) {
      free(mock_fnames[i]);
      free(ChIP_fnames[i]);
   }

   return 0;

}
