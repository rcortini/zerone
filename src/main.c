/* Copyright 2015, 2016 Pol Cusco and Guillaume Filion

   This file is part of Zerone.

   Zerone is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Zerone is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Zerone. If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <strings.h>
#include <getopt.h>
#include "debug.h"
#include "parse.h"
#include "predict.h"
#include "zerone.h"


const char *USAGE =
"\n"
"USAGE:"
"  zerone [options] <input file 1> ... <input file n>\n"
"\n"
"  Input options\n"
"    -0 --mock: given file is a mock control\n"
"    -1 --chip: given file is a ChIP-seq experiment\n"
"    -w --window: window size in bp (default 300)\n"
"    -q --quality: minimum mapping quality (default 20)\n"
"\n"
"  Output options\n"
"    -l --list-output: output list of targets (default table)\n"
"    -c --confidence: print targets only with higher confidence\n"
"                     restricts intervals accordingly in list output\n"
"\n"
"  Other options\n"
"    -h --help: display this message and exit\n"
"    -v --version: display version and exit\n"
"\n"
"EXAMPLES:\n"
" zerone --mock file1.bam,file2.bam --chip file3.bam,file4.bam\n"
" zerone -l -0 file1.map -1 file2.map -1 file4.map\n"
" zerone -l -c.99 -w200 -0 file1.sam -1 file2.sam,file4.sam\n";


#define VERSION "zerone-v1.0"
#define MAXNARGS 255


// Snippets.
int check_strtoX (char *, char *);


//  -----------  Definitions of local one-liners  ----------- //
void say_usage(void) { fprintf(stderr, "%s\n", USAGE); }
void say_version(void) { fprintf(stderr, VERSION "\n"); }


//  -----------  Globals  ----------- //
int errno = 0;


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


void
print_format_warnings_if_needed
(
   char * mock_fnames[],
   char * ChIP_fnames[]
)
// Print warning about .bed and .wig file formats in the header
// of the output file.
{

   int has_bed = 0;
   int has_wig = 0;

   // Check if any input was in .bed format.
   for (int i = 0 ; i < MAXNARGS ; i++) {
      char *fname = mock_fnames[i];
      if (fname == NULL) break;
      int sz = strlen(fname);
      if ( strcmp(fname + sz - 4, ".bed") == 0 ||
           strcmp(fname + sz - 4, ".BED") == 0) {
         has_bed = 1;
      }
      if ( strcmp(fname + sz - 4, ".wig") == 0 ||
           strcmp(fname + sz - 4, ".WIG") == 0) {
         has_wig = 1;
      }
   }

   for (int i = 0 ; i < MAXNARGS ; i++) {
      char *fname = ChIP_fnames[i];
      if (fname == NULL) break;
      int sz = strlen(fname);
      if ( strcmp(fname + sz - 4, ".bed") == 0 ||
           strcmp(fname + sz - 4, ".BED") == 0) {
         has_bed = 1;
      }
      if ( strcmp(fname + sz - 4, ".wig") == 0 ||
           strcmp(fname + sz - 4, ".WIG") == 0) {
         has_wig = 1;
      }
   }

   if (has_bed || has_wig) {
      fprintf(stdout,
         "# -----------    WARNING    ------------\n"
         "# Input file(s) in .bed/.wig format.\n"
         "# Make sure that scores are read counts\n"
         "# and that PCR duplicates were removed.\n"
         "# Make sure that window size is constant\n"
         "# and matches Zerone window option -w.\n"
         "# --------------------------------------\n"
      );
   }

}



int main(int argc, char **argv) {

   debug_print("%s (DEBUG)\n", VERSION);

   if (argc == 1) {
      say_usage();
      exit(EXIT_SUCCESS);
   }

   // Input file names (mock and ChIP).
   char *mock_fnames[MAXNARGS+1] = {0};
   char *ChIP_fnames[MAXNARGS+1] = {0};

   int n_mock_files = 0;
   int n_ChIP_files = 0;
   int no_mock_specified = 1;
   int no_ChIP_specified = 1;

   static int list_flag = 0;
   static int minmapq = 20;
   static int window = 300;
   static int mock_flag = 1;
   static double minconf = 0.0;

   // Needed to check 'strtoul()'.
   char *endptr;

   // Parse options.
   debug_print("%s", "arguments:'\n");
   while(1) {
      int option_index = 0;
      static struct option long_options[] = {
         {"chip",        required_argument,          0, '1'},
         {"confidence",  required_argument,          0, 'c'},
         {"help",        no_argument,                0, 'h'},
         {"list-output", no_argument,       &list_flag,  1 },
         {"mock",        required_argument,          0, '0'},
         {"no-mock",     no_argument,       &mock_flag,  0 },
         {"quality",     required_argument,          0, 'q'},
         {"version",     no_argument,                0, 'v'},
         {"window",      required_argument,          0, 'w'},
         {0, 0, 0, 0}
      };

      int c = getopt_long(argc, argv, "0:1:c:hlq:vw:",
            long_options, &option_index);

      // Done parsing named options. //
      if (c == -1) break;

      switch (c) {
      case 0:
         break;

      case '0':
         debug_print("| mock files(s): %s\n", optarg);
         parse_fname(mock_fnames, optarg, &n_mock_files);
         no_mock_specified = 0;
         break;

      case '1':
         debug_print("| ChIP files(s): %s\n", optarg);
         parse_fname(ChIP_fnames, optarg, &n_ChIP_files);
         no_ChIP_specified = 0;
         break;

      case 'h':
         say_usage();
         return EXIT_SUCCESS;

      case 'l':
         list_flag = 1;
         break;

      case 'c':
         // Decode argument with 'strtod()'
         errno = 0;
         endptr = NULL;
         minconf = strtod(optarg, &endptr);
         if (!check_strtoX(optarg, endptr) || minconf < 0 || minconf > 1) {
            fprintf(stderr,
                  "zerone error: confidence must be "
                  "a float between 0 and 1\n");
            say_usage();
            return EXIT_FAILURE;
         }
         debug_print("| minconf: %f\n", minconf);
         break;

      case 'q':
         // Decode argument with 'strtoul()'
         errno = 0;
         endptr = NULL;
         minmapq = strtoul(optarg, &endptr, 10);
         if (!check_strtoX(optarg, endptr) ||
               minmapq < 0 || minmapq > 254) {
            fprintf(stderr,
                  "zerone error: minimum mapping quality must be "
                  "an integer between 0 and 254\n");
            say_usage();
            return EXIT_FAILURE;
         }
         debug_print("| minmapq: %d\n", minmapq);
         break;

      case 'v':
         say_version();
         return EXIT_SUCCESS;

      case 'w':
         window = atoi(optarg);
         if (window <= 0) {
            fprintf(stderr, "zerone error: window must be a "
                  "positive integer\n");
            say_usage();
            return EXIT_FAILURE;
         }
         debug_print("| window: %d\n", window);
         break;

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

   // Check options.
   if (no_mock_specified && mock_flag) {
      fprintf(stderr,
         "zerone error: specify a file for mock control experiment\n");
      say_usage();
      return EXIT_FAILURE;
   }
   if (no_ChIP_specified) {
      fprintf(stderr,
         "zerone error: specify a file for ChIP-seq experiment\n");
      say_usage();
      return EXIT_FAILURE;
   }

   // Process input files.
   zerone_parser_args_t args;
   args.window = window;
   args.minmapq = minmapq;

   ChIP_t *ChIP = parse_input_files(mock_fnames, ChIP_fnames, args);

   if (ChIP == NULL) {
      fprintf(stderr, "error while reading input\n");
      exit(EXIT_FAILURE);
   }

   // debug info //
   {
      debug_print("%s", "done reading input files\n");
      debug_print("%s", "ChIP:\n");
      debug_print("| r = %ld (dimension)\n", ChIP->r);
      debug_print("| nb = %d (block number)\n", ChIP->nb);
      for (int j = 0 ; j < ChIP->nb ; j++) {
         debug_print("| block %s (size: %d)\n",
               ChIP->nm + 32*j, ChIP->sz[j]);
      }
      // Sum reads of all blocks.
      size_t *nreads = calloc(ChIP->r, sizeof(size_t));
      if (nreads == NULL) {
         fprintf(stderr, "memory error\n");
         exit(EXIT_FAILURE);
      }
      for (int i = 0 ; i < nobs(ChIP) ; i++) {
         for (int j = 0 ; j < ChIP->r ; j++) {
            nreads[j] += ChIP->y[j + i*ChIP->r];
         }
      }
      debug_print("| aggregated mock: %ld reads\n", nreads[0]);
      for (int j = 0 ; j < ChIP->r-1 ; j++) {
         debug_print("| %s: %ld reads\n", ChIP_fnames[j], nreads[j+1]);
      }
      free(nreads);
   }

   // Do zerone.
   debug_print("%s", "starting zerone\n");
   zerone_t *Z = do_zerone(ChIP);

   if (Z == NULL) {
      fprintf(stderr, "run time error (sorry)\n");
      exit(EXIT_FAILURE);
   }

   // debug info //
   {
      debug_print("%s", "Q:\n");
      debug_print("%.3f %.3f %.3f\n", Z->Q[0], Z->Q[3], Z->Q[6]);
      debug_print("%.3f %.3f %.3f\n", Z->Q[1], Z->Q[4], Z->Q[7]);
      debug_print("%.3f %.3f %.3f\n", Z->Q[2], Z->Q[5], Z->Q[8]);

      debug_print("%s", "p:\n");
      for (int j = 0 ; j < 3 ; j++) {
         int off = 0;
         char debuf[512];
         for (int i = 0 ; i < Z->r+1 ; i++) {
            off += sprintf(debuf + off, "%.3f ", Z->p[i+j*(Z->r+1)]);
            if (off > 499) break;
         }
         debug_print("%s\n", debuf);
      }
   }

   // First things first.
   print_format_warnings_if_needed(mock_fnames, ChIP_fnames);

   // Quality control.
   double feat[5];
   double QC = zerone_qc(Z, feat);
   fprintf(stdout, "# QC score: %.3f\n", QC);
   fprintf(stdout, "# features: %.3f, %.3f, %.3f, %.3f, %.3f\n",
                           feat[0], feat[1], feat[2], feat[3], feat[4]);
   fprintf(stdout, "# advice: %s discretization.\n",
         QC >= 0 ? "accept" : "reject");

   // List output.
   if (list_flag) {

      int offset = 0;
      int target = 0;
      double best = 0.0;

      for (int i = 0 ; i < ChIP->nb ; i++) {
         char *name = ChIP->nm + 32*i;

         // Do not print the last bin because it may extend
         // beyond the limit of the chromosome.
         for (int j = 0 ; j < ChIP->sz[i]-1 ; j++) {
            // Toggle on target state.
            double conf = Z->phi[2+(offset+j)*3];
            if (!target && Z->path[offset+j] == 2 && conf > minconf) {
               fprintf(stdout, "%s\t%d\t", name, window*j + 1);
               best = conf;
               target = 1;
            }
            // Toggle off target state.
            else if (target) {
               // Update best score.
               if (conf > best) best = conf;
               if (Z->path[offset+j] != 2 || conf < minconf) {
                  fprintf(stdout, "%d\t%.5f\n", window*j, best);
                  best = 0.0;
                  target = 0;
               }
            }
         }
         // In case the end of the block is a target.
         if (target) {
            fprintf(stdout, "%d\t%.5f\n", window * ChIP->sz[i], best);
            best = 0.0;
            target = 0;
         }

         // Update offset.
         offset += ChIP->sz[i];

      }
   }

   // Table output.
   else {
      // Use 'offset' to navigate in the ChIP blocks.
      uint64_t offset = 0;
      // In case no mock was provided, skip the column.
      const int skipmock = mock_flag ? 0 : 1;

      for (int i = 0 ; i < ChIP->nb ; i++) {
         char *name = ChIP->nm + 32*i;

         // Do not print the last bin because it may extend
         // beyond the limit of the chromosome.
         for (int j = 0 ; j < ChIP->sz[i]-1 ; j++) {
            // Skip if 'confidence' too low.
            if (Z->phi[2+(offset+j)*3] < minconf) continue;
            fprintf(stdout, "%s\t%d\t%d\t%d", name, window*j + 1,
                    // Block name, window start, end, state.
                    window*(j+1), Z->path[offset+j] == 2 ? 1 : 0);
            for (int k = skipmock ; k < Z->ChIP->r ; k++) {
               fprintf(stdout, "\t%d",
                    // Read numbers of each file.
                    Z->ChIP->y[(offset+j)*Z->ChIP->r+k]);
            }
            fprintf(stdout, "\t%.5f\n",
                    // Confidence score.
                    Z->phi[2+(offset+j)*3]);
         }
         // End of the block. Update 'offset' before
         // local window number is reset to 0.
         offset += Z->ChIP->sz[i];
      }
   }


   destroy_zerone_all(Z); // Also frees ChIP.

   for (int i = 0 ; i < MAXNARGS ; i++) {
      free(mock_fnames[i]);
      free(ChIP_fnames[i]);
   }

   return 0;

}
