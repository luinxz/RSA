PATHLIB=.
CFLAGS=-O3 -m64
GENCODE_SM= -arch=compute_61 -code=sm_61

all: main clean

main: print utils arith inverse mod rns gpu test
	nvcc $(GENCODE_SM) -Xptxas -O3 -Xcompiler -fopenmp print.o utils.o arith.o inverse.o mod.o rns.o mul.o main.o $(CFLAGS) -o main

print:
	g++ -c -I$(PATHLIB) $(CFLAGS) -o print.o print.cpp
	
utils:
	g++ -c -I$(PATHLIB) $(CFLAGS) -o utils.o utils.cpp
	
arith:
	g++ -c -I$(PATHLIB) $(CFLAGS) -o arith.o arith.cpp
	
inverse:
	g++ -c -I$(PATHLIB) $(CFLAGS) -o inverse.o inverse.cpp
	
mod:
	g++ -c -I$(PATHLIB) $(CFLAGS) -fopenmp -o mod.o mod.cpp
	
rns:
	g++ -c -I$(PATHLIB) $(CFLAGS) -fopenmp -o rns.o rns.cpp
	
gpu:
	nvcc -c -Xptxas -O3 -I$(PATHLIB) $(GENCODE_SM) mul.o mul.cu
	
test:
	g++ -c -I$(PATHLIB) $(CFLAGS) -o main.o main.cpp

clean:
	rm *.o
