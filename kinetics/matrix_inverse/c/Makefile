# Makefile for FFTbor2D

CCFLAGS           = -c -std=c99 -pedantic -fopenmp -funroll-loops
LDFLAGS           = -llapack -o
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

RNAmfpt : rna_mfpt.o energy_grid_mfpt.o
	$(CC) rna_mfpt.o energy_grid_mfpt.o $(LDFLAGS) RNAmfpt
	
rna_mfpt.o : rna_mfpt.c rna_mfpt.h energy_grid_mfpt.h
	$(CC) $(CCFLAGS) rna_mfpt.c
	
energy_grid_mfpt.o : energy_grid_mfpt.c energy_grid_mfpt.h
	$(CC) $(CCFLAGS) energy_grid_mfpt.c

clean:
	rm -f *.o RNAmfpt

install:
	cp RNAmfpt $(BINDIR)
	