#include "types.h"
#include "typesgpu.h"
#include "texture.h"
#include "utils.h"
#include "print.h"

#define NTASKS 2

#define __mulmod( c0, c1, a0, a1, b0, b1, pc ) \
	asm( "{\n\t"					\
	".reg .u32 dhi;" ".reg .u32 dlo;"		\
	".reg .u32 tx;" ".reg .u32 ty;" ".reg .u32 tz;"	\
	"mul.lo.u32     %0, %2, %4;"			\
	"mul.hi.u32     %1, %2, %4;"			\
	"mad.lo.cc.u32  %1, %3, %4, %1;"		\
	"madc.hi.u32    dlo, %3, %4, 0;"		\
	"mad.lo.cc.u32  %1, %2, %5, %1;"		\
	"madc.hi.cc.u32 dlo, %2, %5, dlo;"		\
	"madc.hi.u32    dhi, %3, %5, 0;"		\
	"mad.lo.cc.u32  dlo, %3, %5, dlo;"		\
	"addc.u32       dhi, dhi, 0;"			\
	"mul.lo.u32     tx, %6, dlo;"			\
	"mul.hi.u32     ty, %6, dlo;"			\
	"mad.lo.cc.u32  ty, %6, dhi, ty;"		\
	"madc.hi.u32    tz, %6, dhi, 0;"		\
	"add.cc.u32     %0, %0, tx;"			\
	"addc.cc.u32    %1, %1, ty;"			\
	"addc.u32       dlo, tz, 0;"			\
	"mad.lo.cc.u32  %0, dlo, %6, %0;"		\
	"addc.u32       %1, %1, 0;"			\
	"}"						\
	: "=r" (c0), "=r" (c1)			\
	: "r" (a0), "r" (a1),			\
	"r" (b0), "r" (b1), "r" (pc) )
	
#define __submod( c_hi, c_lo, a_hi, a_lo, b_hi, b_lo, pc ) \
	asm( "{\n\t"				\
	".reg .s32 t;"				\
	"sub.cc.u32 %1, %5, %3;" 		\
	"subc.cc.u32 %0, %4, %2;"		\
	"subc.u32 t, 0, 0;"			\
	"slct.u32.s32 t, 0, %6, t;"		\
	"sub.cc.u32 %1, %1, t;" 		\
	"subc.u32 %0, %0, 0;"			\
	"}"					\
	: "=r" (c_hi), "=r" (c_lo)		\
	: "r" (a_hi), "r" (a_lo), "r" (b_hi), "r" (b_lo), "r" (pc) )
	
#define __addmod( c0, c1, a0, a1, b0, b1, pc ) \
	asm( "{\n\t"				\
	".reg .u32     t;"			\
	"add.cc.u32    %0, %2, %4;" 		\
	"addc.cc.u32   %1, %3, %5;"		\
	"addc.u32 t,   0, 0;"			\
	"mad.lo.cc.u32 %0, %6, t, %0;"		\
	"addc.u32      %1, %1, 0;" 		\
	"}"					\
	: "=r" (c0), "=r" (c1)		\
	: "r" (a0), "r" (a1), "r" (b0), "r" (b1), "r" (pc) )

#define __mask( c0, a0, b0) \
	asm( "{\n\t"				\
	".reg .pred  p, q;"			\
	"setp.eq.u32 p|q, %1, %2;"		\
	"@p mov.u32  %0, -1;"			\
	"@q mov.u32  %0, 0;"			\
	"}"					\
	: "=r" (c0)				\
	: "r" (a0), "r" (b0) )
	
#define __and(c0, c1, a0, a1, b0) \
	asm( "{\n\t"				\
	"and.b32 %0, %2, %4;"			\
	"and.b32 %1, %3, %4;"			\
	"}"					\
	: "=r" (c0), "=r" (c1)			\
	: "r" (a0), "r" (a1), "r" (b0) )
	
#define __xor(c0, c1, a0, a1, b0, b1) \
	asm( "{\n\t"				\
	"xor.b32 %0, %2, %4;"			\
	"xor.b32 %1, %3, %5;"			\
	"}"					\
	: "=r" (c0), "=r" (c1)			\
	: "r" (a0), "r" (a1), "r" (b0), "r" (b1) )

__device__ void static inline mm(uint2 *out, const uint2 a, const uint2 b, const uint pc)
{
	__mulmod(out->x, out->y, a.x, a.y, b.x, b.y, pc);
}

__device__ void static inline submod(uint2 *out, const uint2 a, const uint2 b, const uint pc)
{
	__submod(out->y, out->x, a.y, a.x, b.y, b.x, pc);
}

__device__ void static inline gand(uint2 *out, const uint2 a, const uint b)
{
	__and(out->y, out->x, a.y, a.x, b);
}

__device__ void static inline gxor(uint2 *out, const uint2 a, const uint2 b)
{
	__xor(out->y, out->x, a.y, a.x, b.y, b.x);
}

__device__ __forceinline__ void warpsReduceInt(volatile uint2 *m, uint tid, uint pc)
{
	__addmod(m[tid].x, m[tid].y, m[tid].x, m[tid].y, m[tid+16].x, m[tid+16].y, pc);
	
	__addmod(m[tid].x, m[tid].y, m[tid].x, m[tid].y, m[tid+8].x, m[tid+8].y, pc);
	
	__addmod(m[tid].x, m[tid].y, m[tid].x, m[tid].y, m[tid+4].x, m[tid+4].y, pc);
	
	__addmod(m[tid].x, m[tid].y, m[tid].x, m[tid].y, m[tid+2].x, m[tid+2].y, pc);
	
	__addmod(m[tid].x, m[tid].y, m[tid].x, m[tid].y, m[tid+1].x, m[tid+1].y, pc);
}

__device__ __forceinline__ void warpsReduceFlo(volatile float *v, uint tid)
{
	v[tid] += v[tid + 16];
	v[tid] += v[tid + 8];
	v[tid] += v[tid + 4];
	v[tid] += v[tid + 2];
	v[tid] += v[tid + 1];
}

__device__ void mulmod(uint2* __restrict__ c, uint2 *d, uint2 a, uint2 b, uint2* __restrict__ spinv, uint2* __restrict__ spplp, 
		       uint* __restrict__ spc, float* __restrict__ sdelta, uint bid) 
{
	__shared__ float alphas[38]; 
	uint tid = threadIdx.x;
	uint2 z, g;
	uint ci;
	
	ci = spc[blockIdx.x];
	
	mm(&z, a, b, spc[tid]);
	
	mm(&g, z, spinv[tid], spc[tid]);
	
	mm(&d[tid], g, spplp[tid], ci);
	
	warpsReduceInt(d, tid, ci);
	
	alphas[tid] = sdelta[tid];
	alphas[tid] += ((float)(g.y >> 18) / 0x4000);
	
	warpsReduceFlo(alphas, tid);
	
	if(tid < 1) 
	{
		if(alphas[0] >= 1)
			submod(&d[tid], fetch_tex((int)(alphas[0]-1) * NR + blockIdx.x), d[tid], ci);
	}
	
	if(tid < 1) { c[bid] = d[0]; }
	
	//__syncthreads();
}

__device__ void inline gpusync(int flag, volatile int* __restrict__ in, volatile int* __restrict__ out)
{
	uint tid = threadIdx.x;
	uint bid = blockIdx.x;
	int old = -9999;
	
	if(tid < 1){ in[bid] = flag; }
	
	if(bid < 1)
	{
		while(old != flag) { old = in[tid]; }
		
		__syncthreads();
		
		out[tid] = flag;
	}
	
	if(bid > 0)
	{
		if(tid < 1) { while(old != flag) { old = out[bid]; } }
	}
	
	__syncthreads();
}

__device__ __forceinline__ void linearpassing(uint2 *aux, uint *mask, volatile uint2 *gpcomp, const uint index, const uint tid)
{
	uint i;
	uint2 t;
	
	__mask(mask[tid], index, tid);
	
	aux->x = 0;
	aux->y = 0;
	
	#pragma unroll 16
	for(i = 0; i < SW; i++)
	{
		gand(&t, *(gpcomp+i*NR+tid), mask[i]);
		
		gxor(aux, t, *aux);
	}
}

__global__ void expmod(uint2* __restrict__ c, uint2* __restrict__ a, uint2* __restrict__ b, uint2* __restrict__ pcomp, 
		       volatile ushort* __restrict__ w, uint2* __restrict__ gpinv, uint2* __restrict__ pplp, 
		       uint* __restrict__ pc, float* __restrict__ delta, uint nw, int *vals1, int *vals2) 
{
	__shared__ uint2 d[38];
	__shared__ uint2 spinv[NR];
	__shared__ uint2 spplp[NR];
	__shared__ uint spc[NR];
	__shared__ float sdelta[NR];
	__shared__ uint mask[NR];
	
	uint2 aux;
	uint tid = threadIdx.x;
	uint bid = blockIdx.x;
	uint i, k, s;
	int j;
	uint2 *t1, *t2;
	uint f;
	
	spinv[tid] = gpinv[tid];
	spplp[tid] = pplp[NR * bid + tid];
	spc[tid] = pc[tid];
	sdelta[tid] = delta[tid];
	
	for(i = 0; i < TESTS; i++)
	{
		k = i * NR; // a
		s = SW * k; // pcomp
		
		if(bid < 1)
		{ pcomp[s+tid] = a[k+tid]; }
		gpusync(-1, vals1, vals2);
		
		#pragma unroll 16
		for(j = 1; j < SW; j++)
		{
			t1 = pcomp+s+j*NR;
			t2 = pcomp+s+(j-1)*NR;
			
			mulmod(t1, d, t2[tid], a[k+tid], spinv, spplp, spc, sdelta, bid);
			gpusync(j, vals1, vals2);
		}
		
		j = nw;
		
		if(bid < 1)
		{
			s = SW * k + w[j] * NR;
			b[tid] = pcomp[s+tid];
		}
		
		gpusync(j, vals1, vals2);
		j--;
		
		t1 = c+k;
		t2 = b;
		f = 1;
		
		while(j >= 0)
		{
			#pragma unroll 4
			for(s = 0; s < NB; s++)
			{
				mulmod(t1, d, t2[tid], t2[tid], spinv, spplp, spc, sdelta, bid);
				gpusync(s, vals1, vals2);
				
				if(f)
				{
					t1 = b;
					t2 = c+k;
					f = 0;
				}
				else
				{
					t1 = c+k;
					t2 = b;
					f = 1;
				}
			}
			
			s = SW * k;
			linearpassing(&aux, mask, pcomp+s, w[j], tid);
			
			mulmod(t1, d, t2[tid], aux, spinv, spplp, spc, sdelta, bid);
			gpusync(TW, vals1, vals2);
			
			if(f)
			{
				t1 = b;
				t2 = c+k;
				f = 0;
			}
			else
			{
				t1 = c+k;
				t2 = b;
				f = 1;
			}
			
			j--;
		}
		
		if(f)
		{
			if(bid < 1) { c[k+tid] = b[tid]; }
		}
	}
}

inline void convert_l2i2_a(uint2 *out, eltr_v *in, uint n, uint m)
{
	uint i,j,k = 0;
	
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < m; j++)
		{
			out[k].x = in[i].x[j]; // Low
			out[k].y = in[i].x[j] >> 32; // High
			k++;
		}
	}
}

inline void convert_l2i2_b(uint2 *out, eltr in, uint n)
{
	uint i, k = 0;
	
	for(i = 0; i < n; i++)
	{
		out[k].x = in[i]; // Low
		out[k].y = in[i] >> 32; // High
		k++;
	}
}

inline void copy_pc(uint *h_pc, uint *pc, uint n)
{
	int i;
	
	for(i = n; i--; )
	{
		h_pc[i] = pc[i];
	}
}

inline void copy_delta(float *h_delta, float *delta, uint n)
{
	int i;
	
	for(i = n; i--; )
	{
		h_delta[i] = delta[i];
	}
}

inline void copy_w(ushort *h_w, ushort *w, uint n)
{
	int i;
	
	for(i = n; i--; )
	{
		h_w[i] = w[i];
	}
}

void fprint_vecintr(geltx_v in[], char chain[], uint n, uint m)
{
	uint i, j;
	FILE *pf;
	char tmp[20];
	
	sprintf(tmp, "%s.txt", chain);
	
	pf = fopen(tmp, "w");
	
	fprintf(pf, "%s:= [\n", chain);
	for(i=0; i < n-1; i++)
	{
		for(j = 0; j < m-1; j++)
		{
			fprintf(pf, "0x%.8X%.8X,", in[i].x[j].y, in[i].x[j].x);
		}
		
		fprintf(pf, "0x%.8X%.8X", in[i].x[j].y, in[i].x[j].x);
		fprintf(pf, ",\n");
	}
	
	for(j = 0; j < m-1; j++)
	{
		fprintf(pf, "0x%.8X%.8X,", in[i].x[j].y, in[i].x[j].x);
	}
	
	fprintf(pf, "0x%.8X%.8X", in[i].x[j].y, in[i].x[j].x);
	fprintf(pf, "\n];\n\n");
	
	fclose(pf);
}

void mul_gpu(eltr_v g_c[], eltr_v g_a[], eltr_v g_b[], rns_v *rns)
{
	// Variables in CPU
	uint i;
	
	// Time on GPU
	float msec = 0;
	cudaEvent_t start;
	cudaEvent_t stop;
	
	// Variables in CPU
	uint2 *h_a;
	uint2 *h_b;
	uint2 *h_c;
	uint2 *h_pinv;
	uint2 *h_pplp;
	uint2 *h_alphaPlp;
	float *h_delta;
	ushort *h_w;
	uint *h_pc;
	short *h_nw;
	
	// Variables in GPU
	uint2 *d_a, *d_b;
	uint2 *d_aux;
	uint2 *d_c;
	uint2 *d_pinv;
	uint *d_pc;
	uint2 *d_pplp;
	uint2 *d_alphaPlp;
	float *d_delta;
	uint2 *d_pcomp;
	ushort *d_w;
	int *d_vals1;
	int *d_vals2;
	
	cudaSetDevice(1);

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaStream_t stream[NTASKS];
	
	// We allocate memory in GPU
	cudaMalloc(&d_a, (TESTS * NR) * sizeof(uint2));
	cudaMalloc(&d_b, (TESTS * NR) * sizeof(uint2));
	cudaMalloc(&d_aux, (NR * NTASKS) * sizeof(uint2));
	cudaMalloc(&d_c, (TESTS * NR * NTASKS) * sizeof(uint2));
	cudaMalloc(&d_pc, (NR * NTASKS) * sizeof(uint));
	cudaMalloc(&d_pinv, (NR * NTASKS) * sizeof(uint2));
	cudaMalloc(&d_pplp, (NR * NR * NTASKS) * sizeof(uint2));
	cudaMalloc(&d_alphaPlp, (NR * NR * NTASKS) * sizeof(uint2));
	cudaMalloc(&d_delta, (NDELTA * NTASKS) * sizeof(float));
	cudaMalloc(&d_pcomp, (TESTS * NR * SW * NTASKS) * sizeof(uint2));
	cudaMalloc(&d_w, (TW * NTASKS) * sizeof(ushort));
	cudaMalloc(&d_vals1, (NR * NTASKS) * sizeof(int));
	cudaMalloc(&d_vals2, (NR * NTASKS) * sizeof(int));
	
	// The vector c is set to zero
	cudaMemset(d_c, 0, (TESTS * NR * NTASKS) * sizeof(uint2));
	
	cudaMallocHost(&h_a, (TESTS * NR) * sizeof(uint2));
	cudaMallocHost(&h_b, (TESTS * NR) * sizeof(uint2));
	cudaMallocHost(&h_c, (TESTS * NR * NTASKS) * sizeof(uint2)); 
	cudaMallocHost(&h_pc, NR * sizeof(uint));
	cudaMallocHost(&h_pinv, NR * sizeof(uint2));
	cudaMallocHost(&h_pplp, (NR * NR) * sizeof(uint2));
	cudaMallocHost(&h_alphaPlp, (NR * NR) * sizeof(uint2));
	cudaMallocHost(&h_delta, NDELTA * sizeof(float));
	cudaMallocHost(&h_w, TW * sizeof(ushort));
	cudaMallocHost(&h_nw, 2 * sizeof(ushort));
	
	convert_l2i2_a(h_a, g_a, TESTS, NR);
	convert_l2i2_a(h_b, g_b, TESTS, NR);
	convert_l2i2_b(h_pinv, rns->modinvspp, NR);
	convert_l2i2_a(h_pplp, rns->pplp, NR, NR);
	convert_l2i2_a(h_alphaPlp, rns->alphaPlp, NR, NR);
	copy_pc(h_pc, rns->c, NR);
	copy_delta(h_delta, rns->delta, NR);
	copy_w(h_w, rns->w, rns->nw+1);
	h_nw[0] = h_nw[1] = rns->nw;
	
	//Create 2 streams
	for (i = 0; i < NTASKS; i++) 
	{
		cudaStreamCreate(&stream[i]);
	}
	
	// Changing L1 cache configuration to 48KB - L1 and 16KB - Shared Memory 
	cudaFuncSetCacheConfig(expmod, cudaFuncCachePreferNone);
// 	cudaFuncSetCacheConfig(expmod, cudaFuncCachePreferShared);
// 	cudaFuncSetCacheConfig(expmod, cudaFuncCachePreferL1);
	
	// Stream 0
	cudaMemcpyAsync(d_a, h_a, (TESTS * NR) * sizeof(uint2), cudaMemcpyHostToDevice, stream[0]);
	cudaMemcpyAsync(d_w, h_w, (rns->nw + 1) * sizeof(ushort), cudaMemcpyHostToDevice, stream[0]);
	cudaMemcpyAsync(d_pc, h_pc, NR * sizeof(uint), cudaMemcpyHostToDevice, stream[0]);
	cudaMemcpyAsync(d_pinv, h_pinv, NR * sizeof(uint2), cudaMemcpyHostToDevice, stream[0]);
	cudaMemcpyAsync(d_pplp, h_pplp, (NR * NR) * sizeof(uint2), cudaMemcpyHostToDevice, stream[0]);
	cudaMemcpyAsync(d_delta, h_delta, NDELTA * sizeof(float), cudaMemcpyHostToDevice, stream[0]);
	cudaMemcpyAsync(d_alphaPlp, h_alphaPlp, (NR * NR) * sizeof(uint2), cudaMemcpyHostToDevice, stream[0]);
	
	// Stream 1
	cudaMemcpyAsync(d_b, h_b, (TESTS * NR) * sizeof(uint2), cudaMemcpyHostToDevice, stream[1]);
	cudaMemcpyAsync(d_w + TW, h_w, (rns->nw + 1) * sizeof(ushort), cudaMemcpyHostToDevice, stream[1]);
	cudaMemcpyAsync(d_pc + NR, h_pc, NR * sizeof(uint), cudaMemcpyHostToDevice, stream[1]);
	cudaMemcpyAsync(d_pinv + NR, h_pinv, NR * sizeof(uint2), cudaMemcpyHostToDevice, stream[1]);
	cudaMemcpyAsync(d_pplp + NR * NR, h_pplp, (NR * NR) * sizeof(uint2), cudaMemcpyHostToDevice, stream[1]);
	cudaMemcpyAsync(d_delta + NDELTA, h_delta, NDELTA * sizeof(float), cudaMemcpyHostToDevice, stream[1]);
	cudaMemcpyAsync(d_alphaPlp + NR * NR, h_alphaPlp, (NR * NR) * sizeof(uint2), cudaMemcpyHostToDevice, stream[1]);
	
	bind_tex(d_alphaPlp);
	
	// beginning time
	cudaEventRecord(start, 0);
	
	for(i = 0; i < NTASKS; i++) 
	{
		cudaStreamSynchronize(stream[i]);
	}
	
	expmod<<<19, 19, 0, stream[0]>>>(d_c, d_a, d_aux, d_pcomp, d_w, d_pinv, d_pplp, d_pc, d_delta, h_nw[0], d_vals1, d_vals2); 
	expmod<<<19, 19, 0, stream[1]>>>(d_c + TESTS * NR, d_b, d_aux + NR, d_pcomp + TESTS * NR * SW, d_w + TW, d_pinv + NR, d_pplp + NR * NR, d_pc + NR, d_delta + NDELTA, h_nw[1], d_vals1 + NR, d_vals2 + NR);
	
	for(i = 0; i < NTASKS; i++) 
	{
		cudaStreamSynchronize(stream[i]);
	}
	
	// End time
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&msec, start, stop);
	
	unbind_tex(d_alphaPlp);
	
// 	// Average time		Total time		All in milliseconds
	printf("%.16f\t\t%.16f\n", msec / TESTS, msec);
	
	cudaMemcpyAsync(h_c, d_c, (TESTS * NR) * sizeof(uint2), cudaMemcpyDeviceToHost, stream[0]);
	cudaMemcpyAsync(h_c + TESTS * NR, d_c + TESTS * NR, (TESTS * NR) * sizeof(uint2), cudaMemcpyDeviceToHost, stream[1]);
	
	for(i = 0; i < NTASKS; i++) 
	{
		cudaStreamSynchronize(stream[i]);
	}
	
	printf("w0:= [\n");
	for(i=0; i < NR-1; i++)
	{
		printf("0x%.8X%.8X, ", h_c[i].y, h_c[i].x);
	}
	printf("0x%.8X%.8X];\n", h_c[i].y, h_c[i].x);
	
	printf("w1:= [\n");
	for(i=TESTS*NR; i < TESTS*NR+NR-1; i++)
	{
		printf("0x%.8X%.8X, ", h_c[i].y, h_c[i].x);
	}
	printf("0x%.8X%.8X];\n", h_c[i].y, h_c[i].x);
	
	for (i = 0; i < NTASKS; i++) 
	{
		cudaStreamDestroy(stream[i]);
	}
	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_c);
	cudaFree(d_pc);
	cudaFree(d_pinv);
	cudaFree(d_pplp);
	cudaFree(d_alphaPlp);
	cudaFree(d_delta);
	cudaFree(d_pcomp);
	cudaFree(d_w);
	cudaFree(d_vals1);
	cudaFree(d_vals2);
	
	cudaFree(h_a);
	cudaFree(h_b);
	cudaFree(h_c);
	cudaFree(h_pc);
	cudaFree(h_pinv);
	cudaFree(h_pplp);
	cudaFree(h_alphaPlp);
	cudaFree(h_delta);
	cudaFree(h_w);
	cudaFree(h_nw);
}
