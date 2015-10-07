#ifndef _DEBUG_H
#define _DEBUG_H

#ifdef DEBUG
#define XDBGX 1
#else
#define XDBGX 0
#endif

#define debug_print(fmt, ...) \
   do { if (XDBGX) fprintf(stderr, "%s:%d:%s(): " #fmt, \
         __FILE__, __LINE__, __func__, ##__VA_ARGS__); } while (0)

#endif
