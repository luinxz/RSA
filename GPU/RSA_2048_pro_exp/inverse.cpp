#include "types.h"
#include "utils.h"
#include "inverse.h"
#include "arith.h"


//Field inversion
// assumes that a is different than zero or p
void inverse(unsigned long *out, eltr a, unsigned long *p)
{
	uint i;
	uint n = NR-1;
	ulong x1[NR]; 
	ulong x2[NR]; 
	ulong u[NR]; 
	ulong v[NR]; 
	
	for(i=0; i < NR; i++)
	{
		
		x2[i] = 0;
		u[i] = a[i];
		v[i] = 0;
		x1[i] = 0;
	}
	x1[0] = 1;
	v[0] = p[0];
	
	if (iszero(a, NR) || types_comp(a, v)) 
	{
		printf("Error!!!!!!: trying to compute the inverse of zero\n");
		*out = 0;
		return;
	}
	
	while(!iszero(u, NR))
	{
		while (types_iseven(u)) 
		{ 
			types_singleshiftright(u);
			
			if (types_iseven(x1))
			{ 
				types_singleshiftright(x1);
				
			} else {
				add128(x1[0], x1[1], x1[2], p[0], p[1]); 
				types_singleshiftright(x1);
			}
		}
		
		while (types_iseven(v)) 
		{
			types_singleshiftright(v);
			
			if (types_iseven(x2))
			{ 
				types_singleshiftright(x2);
			} else {
				add128(x2[0], x2[1], x2[2], p[0], p[1]); 
				types_singleshiftright(x2);
			}
		}
		
		if (!types_less(u, v, n)) 
		{
			Fp_sub256(u, v, p, NR); 
			Fp_sub256(x1, x2, p, NR); 
		} 
		else 
		{
			Fp_sub256(v, u, p, NR); 
			Fp_sub256(x2, x1, p, NR); 
		}
	}
	
	if (isone(u, n)) 
	{
		if (types_less(p, x1, n)) 
		{ printf("Barrett*\n"); } 
		else
		{ *out = x1[0]; }
	} else { 
		if (types_less(p, x2, n)) 
		{ printf("Barrett+\n"); } 
		else 
		{ *out = x2[0]; }
	}
}
