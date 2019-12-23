#include "types.h"
#include "utils.h"
#include "rns.h"
#include "mul.h"
#include "print.h"

elt_v a[TESTS];
elt_v b[TESTS];
eltr_v ra[TESTS];
eltr_v rb[TESTS];
eltr_v rc[TESTS];

void mod_exp()
{
	int x;
	rns_v rns;
	
	srand(time(NULL));
	
	x = posix_memalign((void**)&rns, 64, sizeof(rns_v));
	
	init_primes(rns.p64);
	
	x = posix_memalign((void**)&a, 64, sizeof(elt_v)*TESTS);
	x = posix_memalign((void**)&b, 64, sizeof(elt_v)*TESTS);
	
	init_vectors(a);
	init_vectors(b);
	
	x = posix_memalign((void**)&ra, 64, sizeof(eltr_v)*TESTS);
	x = posix_memalign((void**)&rb, 64, sizeof(eltr_v)*TESTS);
	x = posix_memalign((void**)&rc, 64, sizeof(eltr_v)*TESTS);
	
	get_prime(rns.l);
	get_cs(rns.c);
	
	// pre-computing
	pre_comp(&rns, rns.p64, a, b, ra, rb);
	
	mul_gpu(rc, ra, rb, &rns);
}

int main(int argc, char *argv[])
{
	mod_exp();
	return 0;
}
