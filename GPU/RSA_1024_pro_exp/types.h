#ifndef _TYPES_H_
#define _TYPES_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <string.h>
#include <time.h>

#define DEGREE 512
#define TESTS 3000
#define NI 8
#define NR 19
#define NTHREADS 4
#define NDELTA 19

#define NB	4
#define SW	16 // 4
#define TW	128 // 512 / 4

typedef ulong elt[8] __attribute__((aligned(64)));
typedef ulong eltr[19] __attribute__((aligned(64)));

typedef struct elt_v
{
	elt x;
} elt_v;

typedef struct eltr_v
{
	eltr x;
} eltr_v;

typedef struct rns_v
{
	elt l;
	elt modp64;
	eltr p64;
	ulong mu[9];
	eltr P;
	eltr_v pps[NR];
	eltr invspp;
	eltr modinvspp;
	eltr modsp;
	uint c[NR];
	eltr_v pplp[NR];
	eltr_v alphaPlp[NR];
	float delta[NDELTA];
	uint kmin;
	uint twotos;
	// Binary
	ulong exp[NR];
	
	// m-ary
	ushort w[TW];
	ushort nw;
} rns_v;

extern elt_v a[TESTS];
extern elt_v b[TESTS];
extern eltr_v ra[TESTS];
extern eltr_v rb[TESTS];
extern eltr_v rc[TESTS];

#endif /* TYPES */ 
