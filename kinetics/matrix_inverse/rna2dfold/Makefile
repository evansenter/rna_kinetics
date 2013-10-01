# Makefile for FFTbor2D

CCFLAGS           = -c -std=c99 -pedantic -fopenmp -funroll-loops
LDFLAGS           = -lRNA -lgomp -llapack -o
BINDIR            = ~/bin # Change this to the BINDIR
CC                = gcc
GCC_VERSION      := $(shell expr `$(CC) -dumpversion`)
CC_MAJ_VERSION   := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 1` \* 10000)
CC_MIN_VERSION   := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 2` \* 100)
CC_PATCH_VERSION := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 3`)
GCC_NUM_VERSION  := $(shell expr $(CC_MAJ_VERSION) \+ $(CC_MIN_VERSION) \+ $(CC_PATCH_VERSION))
GCC_GTEQ_4.6.0   := $(shell expr $(GCC_NUM_VERSION) \>= 40600)

ifeq "$(GCC_GTEQ_4.6.0)" "1"
	CCFLAGS += -Ofast -march=native
else
	CCFLAGS += -O3
endif

RNA2DFoldKinetics : rna2dfold.o
	$(CC) rna2dfold.o $(LDFLAGS) RNA2DFoldKinetics
	
rna2dfold.o : rna2dfold.c data_structures.h
	$(CC) $(CCFLAGS) rna2dfold.c

clean:
	rm -f *.o RNA2DFoldKinetics

install:
	cp RNA2DFoldKinetics $(BINDIR)
	