#ifndef _TYPES_H_
#define _TYPES_H_

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <string.h>
#include <time.h>

#define DEGREE 1024
#define TESTS 3000
#define NI 16
#define NR 35
#define NTHREADS 4
#define NDELTA 35

#define NB	4
#define SW	16 // 4
#define TW	256 // 1024 / 4

typedef ulong elt[16] __attribute__((aligned(64)));
typedef ulong eltr[35] __attribute__((aligned(64)));

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
	ulong mu[25];
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
