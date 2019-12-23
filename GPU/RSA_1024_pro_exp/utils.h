#ifndef _UTILS_H_
#define _UTILS_H_ 

#include "types.h"

#define types_cpy(d, a, n){\
	uint i;\
	for(i=n;i--;) d[i]=a[i];\
}

#define set_zero(a, n){\
	uint i;\
	for(i=n;i--;) a[i]=0;\
}

// n = NR
inline unsigned int iszero(unsigned long *a, unsigned long n)
{
	uint i, f=1;
	
	for(i=n; i--;)
	{ 
		if(a[i])
		{
			f=0; 
			break; 
		}
	}
	
	if(f) return 1;
	else return 0;
}

// n = NR-1
inline unsigned int isone(unsigned long *a, unsigned long n)
{
	uint i, f=1;
	
	for(i=n; i>0; i--)
	{ if(a[i] != 0) f=0; }
	
	if(f && a[i]==1)
		return 1;
	else
		return 0;
}

// n = NR-1
inline int types_less(unsigned long *a, unsigned long *b, unsigned long n)
{
	if(n >= NI)
	{
		if(b[n] != 0)
		{ return 1; }
		else
		{ n = n - 1; }
	}
	
	while(a[n] == b[n] && n>0)
	{
		n--;
	}
	return a[n] < b[n];
}

void init_vectors(elt_v a[]);
void init_primes(eltr p);
void types_randgen(elt ax, uint deg);
void get_cs(uint c[]);
void get_prime(ulong *l);
void get_mu(ulong *mu);
void get_exp(rns_v *rns);
void recoding(rns_v *rns);

void get_modp64(rns_v *rns);

#endif /* UTILS */ 
