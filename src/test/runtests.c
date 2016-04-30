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

#include "unittest.h"

int
main(
   int     argc,
   char ** argv
)
{

   // Import test cases.
   extern test_case_t test_cases_hmm[];
   extern test_case_t test_cases_utils[];
   extern test_case_t test_cases_zinm[];
   extern test_case_t test_cases_parse[];
   extern test_case_t test_cases_predict[];
   extern test_case_t test_cases_zerone[];

   // Register test cases.
   const test_case_t *list_of_test_cases[] = {
      test_cases_hmm,
      test_cases_utils,
      test_cases_zinm,
      test_cases_parse,
      test_cases_predict,
      test_cases_zerone,
      NULL,
   };

   // Run the tests.
   return run_unittest(argc, argv, list_of_test_cases);

}
