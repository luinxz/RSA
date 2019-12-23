#ifndef _MOD_H_
#define _MOD_H_  

#include "types.h"

#define BASE 64 
#define MASK 0x8000000000000000

#define mod64(m, p) {\
	while(m > p){ m=m-p; }\
}

void get_mods_multi_core(eltr out, eltr p64);
void get_mods(eltr out, eltr p64);
void new_reduction(unsigned long out[], unsigned long in[], unsigned long p[], unsigned long mods[], uint wsize, uint index);
void Barrett(ulong out[], ulong in[], ulong prime[], ulong mu[]);
void Barrett_Ex(ulong out[], ulong in[], ulong prime[], ulong mu[]);

void new_red(ulong out[], ulong in[], ulong p[], ulong mods[], uint sw);

#endif /* MOD */  
