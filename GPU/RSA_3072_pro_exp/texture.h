#ifndef _TEXTURE_H_
#define _TEXTURE_H_ 

#include <cuda.h>

texture < uint2, 1 > tex_uint2;

inline void bind_tex(const uint2 *x)
{
	size_t offset = size_t(-1);
	cudaBindTexture(&offset, tex_uint2, x);
	if (offset != 0)
		printf("memory is not aligned, refusing to use texture cache\n");
}


inline void unbind_tex(const uint2 *x)
{
	cudaUnbindTexture(tex_uint2);
}


__inline__ __device__ uint2 fetch_tex(const uint &i)
{
	return tex1Dfetch(tex_uint2, i);
}

#endif /* TEXTURE */ 
