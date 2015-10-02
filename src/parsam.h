#define _GNU_SOURCE

#ifndef _PARSAM_H
#define _PARSAM_H

#include "zerone.h"

ChIP_t    * read_sam(char *fn[], unsigned int nfiles);

#endif
