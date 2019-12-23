#include <omp.h>
#include <math.h>
#include "rns.h"
#include "utils.h"
#include "arith.h"
#include "inverse.h"
#include "mod.h"
#include "print.h"

void get_prod_primes(unsigned long P[], unsigned long p64[])
{
	uint i;
	
	P[0] = p64[0];
	for(i = 1; i < NR; i++)
	{
		in_mult(P, P, &p64[i], i, 1);
	}
	
// 	print_int64(P, "M", R48);
}

void get_prods_primes(eltr_v pps[], eltr invps, eltr modinvspp, eltr p64)
{
	uint i, j;
	uint c;
	unsigned long t[NR]; //__attribute__((aligned(64)));
	
	set_zero(t, NR);
	
	for(i = 0; i < NR; i++)
	{
		pps[i].x[0] = 1;
		c=0;
		
		for(j = 0; j < NR; j++)
		{
			if(i != j) 
			{
				in_mult(pps[i].x, pps[i].x, p64+j, (c==0)?1:c, 1);
				c++;
			}
		}
		
		t[0] = p64[i];
		inverse(&invps[i], pps[i].x, t);
		
		modinvspp[i] = invps[i];
		while(modinvspp[i]>p64[i])
		{
			modinvspp[i] = modinvspp[i] - p64[i];
		}
	}
// 	print_int64(pps[0].x, "pps", NR);
// 	print_ws64(invps, "invps", NR);
}

void ppmodlmodp(eltr_v pplp[], eltr_v pps[], rns_v *rns)
{
	uint i,j;
	ulong t[NI];
	eltr_v pplp_t[NR]; 
	
// 	print_int64(rns->l, "l", NI);
// 	print_int64(pps[0].x, "pps", NR);
	
	for(i = 0; i < NR; i++)
	{
		// t = Pi mod l
		new_red(t, pps[i].x, rns->l, rns->modp64, NR-1);
		
// 		if(i == 0)
// 		{
// 			print_int64(t, "t", NI);
// 		}
		
		for(j = 0; j < NR; j++)
		{
			// t mod p64
			new_reduction(pplp[i].x, t, rns->p64, rns->modsp, NI, j);
		}
	}
	
// 	print_ws64(pplp[0].x, "pplp", NR);
	
	// Transponer
	for(i = 0; i < NR; i++)
	{
		for(j = 0; j < NR; j++)
			pplp_t[i].x[j] = pplp[j].x[i];
	}
	
	memcpy(pplp, pplp_t, NR * sizeof(eltr_v));
	
// 	printf("pplp_t:= [ ");
// 	for(i = 0; i < NR-1; i++)
// 		printf("0x%.16lX, ", pplp[0].x[i]);
// 	
// 	printf("0x%.16lX];\n", pplp[0].x[i]);
}

void alphaPlp(eltr_v alphaPlp[], rns_v *rns)
{
	uint alpha = 1;
	uint i, j;
	ulong v[NR+1];
	ulong valp[NR];
	
	set_zero(valp, NR);
	
// 	print_int64(rns->l, "l", NI);
// 	print_ws64(rns->p64, "p", NR);
	
	for(i = 0; i < NR-1; i++, alpha++)
	{
		set_zero(v, NR+1);
		valp[0] = alpha;
		in_mult(v, rns->P, valp, NR, 1);
		
// 		if(i == 93)
// 			print_int64(v, "v", 97);
		
		new_red(v, v, rns->l, rns->modp64, NR+1);
		
// 		if(i == 62)
// 			print_int64(v, "v", NI);
		
		for(j = 0; j < NR; j++)
		{
			new_reduction(alphaPlp[i].x, v, rns->p64, rns->modsp, NI, j);
		}
		
// 		if(i == 93)
// 			print_ws64(alphaPlp[i].x, "alphaPlp", NR);
	}
	
	for(j = 0; j < NR; j++)
	{
		alphaPlp[i].x[j] = 0;
	}
}

void convert3072_64(eltr_v out[], elt_v in[], rns_v *rns, eltr p64)
{
	uint i, j;
	
// 	print_int64(in[499].x, "a", N48);
	for(i=0; i<TESTS; i++)
	{
		for(j=0; j<NR; j++)
		{
			new_reduction(out[i].x, in[i].x, p64, rns->modsp, NI, j);
		}
	}
// 	print_ws64(out[499].x, "out", R48);
}

void calc_delta(rns_v *rns)
{
	uint i;
	float epsilon = 0;
	uint q = 14;
	
	memset(rns->delta, 0, NDELTA * sizeof(float));
	
	for(i = 0; i < NR; i++)
	{
		epsilon += rns->c[i];
	}
	
	epsilon = (float)(epsilon / pow(2, BASE));
	rns->delta[0] = NR * (float)((float) (pow(2, BASE - q)) / pow(2, BASE));
	rns->delta[0] = rns->delta[0] + epsilon;
	rns->kmin = BASE - q;
	rns->twotos = 2 << (q - 1);
}

// Pre-computing
void pre_comp(rns_v *rns, eltr p64, elt_v a[], elt_v b[], eltr_v ra[], eltr_v rb[])
{
	double ti, tf;
	struct timeval tm;
	
	gettimeofday(&tm, NULL);
	ti = tm.tv_sec+(tm.tv_usec/1000000.0);
	
	// Product p1, p2, ..., pn
	get_prod_primes(rns->P, p64);
	
	gettimeofday(&tm, NULL);
	tf = tm.tv_sec+(tm.tv_usec/1000000.0);
// 	printf("Primes-mult time: %f seg\n", tf - ti);
	
	// -----------------------------------------------------------------------
	
	gettimeofday(&tm, NULL);
	ti = tm.tv_sec+(tm.tv_usec/1000000.0);
	
	// Products of the prime N32-1 and inverses
	get_prods_primes(rns->pps, rns->invspp, rns->modinvspp, p64);
	
	// We obtain mods
	get_mods(rns->modsp, p64);
	
	get_modp64(rns);
	
	ppmodlmodp(rns->pplp, rns->pps, rns);
	
	alphaPlp(rns->alphaPlp, rns);
	
	calc_delta(rns);
	
	get_exp(rns);
	recoding(rns);
	
	gettimeofday(&tm, NULL);
	tf = tm.tv_sec+(tm.tv_usec/1000000.0);
// 	printf("Parameters time: %f seg\n", tf - ti);
	//---------------------------------------------------------------------------------
	
	gettimeofday(&tm, NULL);
	ti = tm.tv_sec+(tm.tv_usec/1000000.0);
	
	convert3072_64(ra, a, rns, p64);
	convert3072_64(rb, b, rns, p64);
	
	gettimeofday(&tm, NULL);
	tf = tm.tv_sec+(tm.tv_usec/1000000.0);
// 	printf("Pre-computing time: %f seg\n", tf - ti);
	
// 	print_ws64(ra[0].x, "a", R48);
// 	print_ws64(rb[0].x, "b", R48);
}
