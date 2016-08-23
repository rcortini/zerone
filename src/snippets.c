#include <ctype.h>
#include <stdlib.h>

int
check_strtoX
(
  char *nptr,
  char *endptr
)
// Convenience function to check whether 'strtoul()' and
// 'strtod()' could not parse their respective input
// arguments. Returns 0 in case of failure.
//
// For 'strtoul()' accept integer numbers separated by
// spaces or at the end of the string, like "345 " or "2366".
// Reject numbers next to puncutation or letters like "3.4",
// "3,4", "3e6" and "4a".
// 
// For 'strtod()' accept float numbers separated by spaces
// or at the end of the string, like "345.0 ", "3e4" or "2366".
// Reject numbers next to puncutation or letters like "3,4"
// or "4a".
{
   // If 'endptr' is NULL or if no character was
   // read, parsing the number has failed.
   if (endptr == NULL)    return 0;
   if (endptr == nptr)    return 0;
   // The first unread character is either
   // the terminator of the string or a space.
   if (*endptr == '\0')   return 1;
   if (isspace(*endptr))  return 1;
   // In all other cases, another character
   // was appended to the number: reject.
   return 0;
}
