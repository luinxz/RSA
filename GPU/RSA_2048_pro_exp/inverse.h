#ifndef _INVERSE_H_
#define _INVERSE_H_  

#include "types.h"


#define types_comp(a, b)	\
	(memcmp(a, b, sizeof(unsigned long)*NR)==0)
	
#define types_iseven(a) !(a[0]&0x1)

#define types_singleshiftright(ax) {\
	int k, sshift = 64-1;\
	for (k=0; k<NR-1; k++) {\
			ax[k] >>= 1;\
			ax[k]  ^= (ax[k+1]&1)<<sshift;\
	}\
	ax[k] >>= 1;\
}

void inverse(unsigned long *out, eltr a, unsigned long *p);
// void inverse(unsigned long *out, eltr a, unsigned long p, uint k);

#endif /* INVERSE */  
