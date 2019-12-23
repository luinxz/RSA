#ifndef _TYPESGPU_H_
#define _TYPESGPU_H_

#define align __attribute__ ((aligned (32)))

typedef align uint2 gelt[16];
typedef align uint2 geltx[35];

typedef struct gelt_v
{
	gelt x;
} gelt_v;

typedef struct geltx_v
{
	geltx x;
} geltx_v;

typedef struct dparams_v
{
	uint kmin;
	uint twotos;
} dparams_v;

#endif /* TYPESGPU */  
