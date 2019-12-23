#include "types.h"
#include "utils.h"
#include "print.h"

void types_randgen(elt ax, uint deg) 
{
	
	uint wsize = 64;
	int  i, j, E = (deg+wsize-1)/wsize, aux1 = deg % wsize;
	uint64_t aux;
	
	set_zero(ax, NI);

	E = !E ? 1: E;
	for(i = 0; i < E; i++) {
		for(j = 32; j > 0; j >>=2) {
			aux = rand();
			if (rand()>>15) aux |= 0x8000;
			ax[i] <<= 16;
			ax[i] |= aux;
		}
		
		ax[i] <<= 32;
		for(j = 32; j > 0; j >>=2) {
			aux = rand();
			if (rand()>>15) aux |= 0x8000;
			ax[i] <<= 16;
			ax[i] |= aux;
		}
	}
	if (aux1) {ax[E-1] &= ~(-(long int)((ulong)1 << aux1));}
}

void init_primes(eltr p)
{
	uint i;
	ulong primes[] = {0xFFFFFFFFFFFFFFFD, 0xFFFFFFFFFFFFFFFB, 0xFFFFFFFFFFFFFFF9, 0xFFFFFFFFFFFFFFF7, 0xFFFFFFFFFFFFFFF5, 0xFFFFFFFFFFFFFFF1, 0xFFFFFFFFFFFFFFEF, 0xFFFFFFFFFFFFFFDF, 
		0xFFFFFFFFFFFFFFDD, 0xFFFFFFFFFFFFFFD9, 0xFFFFFFFFFFFFFFD3, 0xFFFFFFFFFFFFFFD1, 0xFFFFFFFFFFFFFFCB, 0xFFFFFFFFFFFFFFC7, 0xFFFFFFFFFFFFFFC5, 0xFFFFFFFFFFFFFFC1, 
		0xFFFFFFFFFFFFFFB5, 0xFFFFFFFFFFFFFFB3, 0xFFFFFFFFFFFFFFAD, 0xFFFFFFFFFFFFFFA9, 0xFFFFFFFFFFFFFFA7, 0xFFFFFFFFFFFFFFA1, 0xFFFFFFFFFFFFFF9D, 0xFFFFFFFFFFFFFF97, 
		0xFFFFFFFFFFFFFF8F, 0xFFFFFFFFFFFFFF8B, 0xFFFFFFFFFFFFFF89, 0xFFFFFFFFFFFFFF83, 0xFFFFFFFFFFFFFF7F, 0xFFFFFFFFFFFFFF71, 0xFFFFFFFFFFFFFF6D, 0xFFFFFFFFFFFFFF67};
	
	for(i = 0; i < NR; i++)
	{
		p[i] = primes[i];
	}
}

void init_vectors(elt_v a[])
{
	uint i;
	
	for(i = TESTS; i--;)
	{
		types_randgen(a[i].x, DEGREE);
	}
}

void get_cs(uint c[])
{
	uint i;
	
	uint t[] = {0x3, 0x5, 0x7, 0x9, 0xB, 0xF, 0x11, 0x21, 0x23, 0x27, 0x2D, 0x2F, 0x35, 0x39, 0x3B, 0x3F, 0x4B, 0x4D, 0x53, 0x57, 0x59, 0x5F, 0x63, 0x69, 0x71, 0x75, 0x77, 0x7D, 0x81, 0x8F, 0x93, 0x99};
	
	for(i=0; i < NR; i++)
	{
		c[i] = t[i];
	}
}

void get_prime(ulong *l)
{
	uint i;
	ulong t[] = {0x346BBB2C439C1D6B, 0x9141930074847820, 0x82BAA847AEC58692, 0xE42EEA45AAF9FF7A, 0xDDD30E553D9DCF92, 0x2F31D96CA4EB3ABB, 0xD1D5C92C50BF302C, 0xDB98440F8196DDAB};
	
	for(i=NI; i--; )
	{
		l[i] = t[i];
	}
}

void get_mu(ulong *mu)
{
	uint i;
	
	ulong t[] = {0x5218EED06980AF6A, 0x79C7977CE2F50D6B, 0x39A2AE7746167113, 0x864C84A65E36CD58, 0xE38327B105C5EC71, 0x5873347B021D69D5, 
		0xE09E5E27A178BFA, 0xF275456ADC16289F, 0x45A4AC4245CDDF59, 0xF4137BB500BC54C6, 0x7870C81537EAFDDA, 0x97E0B8D6089EDAF0, 
		0x77E42D38BDC54009, 0xCFD02FF32A3254A7, 0xC0092DD3725C667B, 0x36709A8961BC72F0, 0x1705B8F8BE7B776C, 0x5F99BF1D1F3D6D61, 
		0x8B41D300E3ADB744, 0xEF620E64C0C0C287, 0x4B1E70EE98A67AD1, 0x9DC33B6E28B26810, 0x3488056394870FF6, 0x4474174183D285FC, 
		0x926D6F4DF89072E1, 0x7E2BEC7DF5402A07, 0x1558A192F8B4C5FD, 0xF9C0A7F3C3B3C381, 0x569E5F6FA6D9BBC7, 0xB441ACA60EC2748B, 
		0xE2028B4B84ED2838, 0xFF18BE6F88C28A0E, 0xD12858DE21A76211, 0x4923C06920AC00A5, 0x40A1D21EA99DA25E, 0xB548FB986C3BB578, 
		0xF5C341B03321A0E5, 0x289330F66B950510, 0x17BCCF1556E9EE75, 0xA42F47D8042AF346, 0x625CCB800E60992, 0xDB81DB0730B97709, 
		0x19AA961BA7173C29, 0xE7E4610D115A4D6A, 0x8A3927C767AEE82A, 0x516D7F0896550F57, 0x866BF84A14FB69F1, 0xB3EAFAC6259B8C29, 0x1};
	
	for(i=49; i--; )
	{
		mu[i] = t[i];
	}
}

void get_exp(rns_v *rns)
{
	uint i;
	
	ulong t[] = {0xCCC64C25214DC835, 0x416E20A5A82DFEBB, 0x7F8E974A0D898447, 0x2A9D3DB70F34E313, 0x18CE2FD1D9B9D3C9, 0xEA3907E86A050DA0, 0xA19973CF225EE94B, 0x8EE9B8551AECFD6B};
	
	
	for(i=NI; i--; )
	{
		rns->exp[i] = t[i];
	}
}

void recoding(rns_v *rns)
{
	uint i = 0, j;
	uint g = 1;
	uint m;
	uint d, dp;
	uint k = NB; 
	uint s;
	
	m = 1 << k;
	s = 64 - k;
	
	while(rns->exp[0] >= m + 1)
	{
		d = rns->exp[0] & 0xF;
		
		dp = d + g + m - 2;
		
		rns->w[i] = dp & 0xF;
		
		g = dp >> k;
		
		for(j = 0; j < NI-1; j++)
		{
			rns->exp[j] = (rns->exp[j] >> k) | (rns->exp[j+1] << s);
		}
		rns->exp[j] = rns->exp[j] >> k;
		
		i++;
	}
	
	rns->w[i] = rns->exp[0] + g - 2;
	rns->nw = i; // [0 ... n]
	
// 	printf("i: %d\n", i);
// 	
// 	printf("w:= [\n");
// 	for(j = 0; j < i; j++)
// 	{
// 		printf("%d, ", rns->w[j]);
// 	}
// 	printf("%d ];\n", rns->w[j]);
// 	exit(0);
}

void get_modp64(rns_v *rns)
{
	uint i;
	
	ulong t[] = {0xCB9444D3BC63E295, 0x6EBE6CFF8B7B87DF, 0x7D4557B8513A796D, 0x1BD115BA55060085, 0x222CF1AAC262306D, 0xD0CE26935B14C544, 0x2E2A36D3AF40CFD3, 0x2467BBF07E692254};
		
	for(i=NI; i--; )
	{
		rns->modp64[i] = t[i];
	}
}
