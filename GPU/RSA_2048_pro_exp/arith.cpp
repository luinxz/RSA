#include "arith.h"
#include "utils.h"
#include "print.h"

/*Arithmetic functions */

// Function to multiply two integers
void in_mult(unsigned long out[], unsigned long a[], unsigned long b[], uint s, uint t) 
{
	uint i, j, f;
	__uint128_t hl;
	unsigned long U, V;
	unsigned long l, h;
	unsigned long c[102]; //__attribute__((aligned(64))); 
	
	if(s+t > NR)
	{ set_zero(c, s+t); f=1;}
	else
	{ set_zero(c, NR); f=0;}
	
	for(i=0; i<s; i++) 
	{
		U = 0;
		for(j=0; j<t; j++) 
		{
			//(C, S)  = t[i+j] + a[j]*b[i] + C;
			hl = ((__uint128_t) a[i]*b[j]);
			l  = hl;
			h  = hl >> 64;
			l += U; 
			h += U > l ? 1 : 0;
			V   = l + c[i+j];
			U   = l > V ? h+1 : h;
			c[i+j]  = V;
		}
		c[i+j] = U;
	}
	
	if(f)
	{ types_cpy(out, c, s+t); }
	else
	{ types_cpy(out, c, NR); }
}

void in_mult_ex(unsigned long out[], unsigned long a[], unsigned long b[], uint s, uint t) 
{
	uint i, j;
	__uint128_t hl;
	unsigned long U, V;
	unsigned long l, h;
	unsigned long c[NR]; //__attribute__((aligned(64))); 
	
	for(i = 0; i < NR; i++) { c[i] = 0; }
	
	for(i=0; i<s; i++) 
	{
		U = 0;
		for(j=0; j<t-i; j++) 
		{
			//(C, S)  = t[i+j] + a[j]*b[i] + C;
			hl = ((__uint128_t) a[i]*b[j]);
			l  = hl;
			h  = hl >> 64;
			l += U; 
			h += U > l ? 1 : 0;
			V   = l + c[i+j];
			U   = l > V ? h+1 : h;
			c[i+j]  = V;
		}
		c[i+j] = U;
	}
	
	for(i = 0; i < NR; i++) { out[i] = c[i]; }
}

// Function that add at least 4 words
void in_add256(unsigned long a[], unsigned long b[], uint n)
{
	uint i; 
	
	for(i=0; i<n; i+=4)
	{
		add256(a[i], a[i+1], a[i+2], a[i+3], a[i+4], b[i], b[i+1], b[i+2], b[i+3]);
	}
}

void in_add256_b(unsigned long a[], unsigned long b[], ulong *f, uint n)
{
	uint i; 
	
	n = n - 4;
	for(i=0; i<n; i+=4)
	{
		add256(a[i], a[i+1], a[i+2], a[i+3], a[i+4], b[i], b[i+1], b[i+2], b[i+3]);
	}
	
	add256(a[i], a[i+1], a[i+2], a[i+3], *f, b[i], b[i+1], b[i+2], b[i+3]);
}

// Function that sub at least 4 words
void in_sub256(unsigned long a[], unsigned long b[], uint n)
{
	uint i;
	
	for(i=0; i<n; i+=4)
	{
		sub256(a[i],a[i+1],a[i+2],a[i+3],a[i+4],b[i],b[i+1],b[i+2],b[i+3]);
	}
}

void in_sub64x64(unsigned long a[], unsigned long b[], uint n)
{
	uint i;

	for (i = 0; i < n; ++i)
	{
		sub64(a[i],a[i+1],b[i]);
	}
}

void in_add64x64(unsigned long a[], unsigned long b[], uint n)
{
	uint i;

	for (i = 0; i < n; ++i)
	{
		add64(a[i],a[i+1],b[i]);
	}
}

void Fp_sub256(unsigned long a[], unsigned long b[], unsigned long p[], unsigned int n)
{
	unsigned long x[NR]; 
	
	// Is a gt b?
	if (types_less(a, b, n-1))  //out = x + p - y;
	{
		types_cpy(x, p, n);
		//in_sub256(x, b, n);
		//in_add256(a, x, n);
		in_sub64x64(x, b, n-1);
		in_add64x64(a, x, n-1);
	} 
	else  // out = x-y
	{
		//in_sub256(a, b, n);
		in_sub64x64(a, b, n-1);
	}
}
