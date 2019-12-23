#ifndef _ARITH_H_
#define _ARITH_H_  

#include "types.h"

#define add256(c0,c1,c2,c3,c4,a0,a1,a2,a3)\
	asm("add %5, %0 \n\t    \
	     adc %6, %1 \n\t    \
	     adc %7, %2 \n\t    \
	     adc %8, %3 \n\t    \
	     adc $0, %4 \n\t"	\
	     :"+r"(c0),"+r"(c1),"+r"(c2),"+r"(c3),"+r"(c4)    \
	     : "r"(a0),"r"(a1),"r"(a2),"r"(a3)   \
	)
	
#define add256_a(c0,c1,c2,c3,a0,a1,a2,a3)\
	asm("add %5, %0 \n\t    \
	     adc %6, %1 \n\t    \
	     adc %7, %2 \n\t    \
	     adc %8, %3 \n\t    \
	     adc $0, %4 \n\t"	\
	     :"+r"(c0),"+r"(c1),"+r"(c2),"+r"(c3)    \
	     : "r"(a0),"r"(a1),"r"(a2),"r"(a3)   \
	)
	
#define add128(c0,c1,c2,a0,a1)\
	asm("add %3, %0 \n\t    \
	     adc %4, %1 \n\t    \
	     adc $0, %2 \n\t"	\
	     :"+r"(c0),"+r"(c1),"+r"(c2)   \
	     : "r"(a0),"r"(a1)   \
	)
	
#define sub256(c0,c1,c2,c3,c4,a0,a1,a2,a3)\
	asm("sub %5, %0 \n\t    \
	     sbb %6, %1 \n\t    \
	     sbb %7, %2 \n\t    \
	     sbb %8, %3 \n\t    \
	     sbb $0, %4 \n\t"    \
	     :"+r"(c0),"+r"(c1),"+r"(c2),"+r"(c3),"+r"(c4)    \
	     : "r"(a0),"r"(a1),"r"(a2),"r"(a3)   \
	)

#define sub128(c0,c1,a0,a1)\
	asm("sub %2, %0 \n\t    \
	     sbb %3, %1 \n\t"    \
	     :"+r"(c0),"+r"(c1)    \
	     : "r"(a0),"r"(a1)   \
	)
	
#define sub64(c0,c1,a0)\
	asm("sub %2, %0 \n\t    \
	     sbb $0, %1 \n\t"    \
	     :"+r"(c0),"+r"(c1)    \
	     : "r"(a0)   \
	)

#define add64(c0,c1,a0)\
	asm("add %2, %0 \n\t	\
	     adc $0, %1 \n\t"	\
	     :"+r"(c0), "+r"(c1)\
	     :"r"(a0)		\
	)
	
#define mul64x64(c0,c1,a0,b0) \
	asm("mov %3, %%rdx 	\t\n	\
	     mulx %2, %0, %1 	\t\n"	\
	     : "+r"(c0), "+r"(c1)	\
	     : "r"(a0), "r"(b0) \
	     : "%rdx")

void in_mult(unsigned long out[], unsigned long a[], unsigned long b[], uint s, uint t);
void in_mult_ex(unsigned long out[], unsigned long a[], unsigned long b[], uint s, uint t);
void in_add256(unsigned long a[], unsigned long b[], uint n);
void in_sub256(unsigned long a[], unsigned long b[], uint n);
void Fp_sub256(unsigned long a[], unsigned long b[], unsigned long p[], uint n);

void in_add256_b(unsigned long a[], unsigned long b[], ulong *f, uint n);

#endif /* ARITH */ 
