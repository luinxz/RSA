#include "types.h"

void print_int64(unsigned long *in, char chain[], uint n)
{
	uint i; 
	
	printf("%s:= 0x", chain);
	for(i = n-1; i > 0; i--)
	{
		printf("%.16lX", in[i]);
	}
	printf("%.16lX;\n", in[i]);
}

void print_ws64(unsigned long *in, char chain[], uint n)
{
	uint i; 
	
	printf("%s:= [", chain);
	for(i = 0; i < n-1; i++)
	{
		printf("0x%.16lX,", in[i]);
	}
	printf("0x%.16lX];\n", in[i]);
} 

void fprint_vecint(elt_v in[], char chain[], uint n, uint m)
{
	uint i, j;
	FILE *pf;
	char tmp[20];
	
	sprintf(tmp, "%s.txt", chain);
	
	pf = fopen(tmp, "w");
	
	fprintf(pf, "%s:= [\n", chain);
	for(i=0; i < n-1; i++)
	{
		fprintf(pf, "0x");
		
		for(j = m; j--; )
		{
			fprintf(pf, "%.16lX", in[i].x[j]);
		}
		
		fprintf(pf, ",\n");
	}
	
	fprintf(pf, "0x");
	
	for(j = m; j--; )
	{
		fprintf(pf, "%.16lX", in[i].x[j]);
	}
	
	fprintf(pf, "\n];\n\n");
	
	fclose(pf);
}

