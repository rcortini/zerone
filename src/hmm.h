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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <err.h>

#ifndef _HMM_HEADER_
#define _HMM_HEADER_

#define D double
#define I int
#define U unsigned int
#define V void
#define cD const double
#define cU const unsigned int

// function      ( 1    2   3    4    5    6   7   8  )
D  block_fwdb    (  U,  U, cU*,  D*,  D*,  D*, D*, D* );
I  block_viterbi ( cU, cU, cU*, cD*, cD*,  D*, I*     );
V  bwd           (  U,  U, cD*,  D*,  D*,  D*         );
D  fwd           (  U,  U, cD*, cD*,  D*              );
D  fwdb          (  U,  U, cD*, cD*,  D*,  D*, D*     );
V  viterbi       (  U,  U, cD*, cD*,  cD*, I*         );

#undef D
#undef I
#undef U
#undef V
#undef cD
#undef cU

#endif
