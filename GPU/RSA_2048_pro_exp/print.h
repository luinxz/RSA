#ifndef _PRINT_H_
#define _PRINT_H_  

#include "types.h"

void print_int64(unsigned long *in, char chain[], uint n);
void print_ws64(unsigned long *in, char chain[], uint n);
void fprint_vecint(elt_v in[], char chain[], uint n, uint m);

#endif /* PRINT */ 
