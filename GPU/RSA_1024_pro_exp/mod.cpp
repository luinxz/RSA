#include "mod.h"
#include "arith.h"
#include "utils.h"
#include "print.h"

void get_mods(eltr out, eltr p64)
{
	unsigned long b[2] = {0x0000000000000000, 0x1};
	unsigned long tb[2];
	unsigned long tp[2];
	uint i;
	
	for(i = 0; i < NR; i++)
	{
		// We obtain mods
		tb[0] = b[0];
		tb[1] = b[1];
		tp[1] = 0;
		tp[0] = p64[i];
		sub128(tb[0], tb[1], tp[0], tp[1]);
		mod64(tb[0], p64[i]);
		
		out[i] = tb[0];
	}
}

void new_reduction(unsigned long out[], unsigned long in[], unsigned long p[], unsigned long mods[], uint wsize, uint index)
{
	unsigned long T;
	unsigned long twidth;
	int i;
	uint N;
	uint bytes = BASE-1;
	unsigned long G[NI];
	
	types_cpy(G, in, NI);
	N = wsize - 1;
	
	while(N > 0)
	{
		T = G[N];
		
		for(i = BASE; i--; )
		{
			twidth = T & MASK;
			twidth = twidth >> bytes;
			T = T << 1;
			
			while(twidth)
			{
				twidth = 0;
// 				T = T + mods; // 2^(width(Y)) mod Y
				add64(T, twidth, mods[index]);
			}
		}
		
// 		G[N-1] = G[N-1] + T
		add64(G[N-1], twidth, T);
		
		while(twidth)
		{
			twidth = 0;
			add64(G[N-1], twidth, mods[index]);
		}
		
		N = N-1;
	}
	
	mod64(G[0], p[index]);
	
	out[index] = G[0];
}

void new_red(ulong out[], ulong in[], ulong p[], ulong mods[], uint sw)
{
	ulong T[NI];
	ulong twidth;
	int i, j;
	uint k;
	uint N;
	uint bytes = BASE-1;
	ulong G[40]; // 4 * 8
	
	memset(G, 0, 40*sizeof(ulong));
	types_cpy(G, in, sw);
	
	k = sw / NI; 
	N = k * NI;
	
	while(k > 0)
	{
		for(i = 0; i < NI; i++) { T[i] = *(G+N+i); }
		
		for(i = DEGREE; i--; )
		{
			twidth = T[NI-1] & MASK;
			twidth = twidth >> 63;
			
			j = NI - 1;
			while(j > 0)
			{
				T[j] = T[j] << 1 | ((T[j-1] & MASK) >> 63);
				j--;
			}
			T[j] = T[j] << 1;
			
			while(twidth)
			{
				twidth = 0;
// 				T = T + mods; // 2^(width(Y)) mod Y
				in_add256_b(T, mods, &twidth, NI);
			}
		}
		
// 		G[N-1] = G[N-1] + T
		in_add256_b(G+(k-1)*NI, T, &twidth, NI);
		
		while(twidth)
		{
			twidth = 0;
			in_add256_b(G+(k-1)*NI, mods, &twidth, NI);
		}
		
		k--;
		N = k * NI;
	}
	
	while(types_less(p, G, 7))
	{
		in_sub256(G, p, NI);
	}
	
	memcpy(out, G, NI * sizeof(ulong));
}

void Barrett(ulong out[], ulong in[], ulong prime[], ulong mu[])
{
	ulong qh[97];
	ulong rs[NR];
	
	in_mult(qh, mu, in+47, 49, 49);
	
	in_mult_ex(rs, prime, qh+49, 48, 48);
	
	types_cpy(out, in, 48);
	
	in_sub256(out, rs, 48);
	
	while(types_less(prime, out, 47))
	{
		in_sub256(out, prime, 48);
	}
}

void Barrett_Ex(ulong out[], ulong in[], ulong prime[], ulong mu[])
{
	ulong qh[100];
	ulong rs[100];
	
	in_mult(qh, in+47, mu, 50, 49);
	
	in_mult_ex(rs, prime, qh+49, 48, 49);
	
	types_cpy(out, in, 49);
	
	in_sub256(out, rs, 49);
	
	while(types_less(prime, out, 48))
	{
		in_sub256(out, prime, 48);
	}
}
