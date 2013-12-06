CCFLAGS           = -I/usr/local/opt/lapack/include -c -std=c99 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q
LDFLAGS           = -lm -lgomp -L/usr/local/opt/lapack/lib -llapacke -lgslcblas -lgsl -o
BINDIR            = ~/bin # Change this to the BINDIR
LIBDIR            = ~/lib # Change this to the LIBDIR
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

RNAmfpt : rna_mfpt.o params.o energy_grid_mfpt.o
	$(CC) rna_mfpt.o params.o energy_grid_mfpt.o $(LDFLAGS) RNAmfpt
	ar cr shared/libMFPT.a energy_grid_mfpt.o params.o
	
rna_mfpt.o : rna_mfpt.c rna_mfpt.h energy_grid_mfpt.h constants.h params.h
	$(CC) $(CCFLAGS) rna_mfpt.c
	
params.o : params.c params.h
	$(CC) $(CCFLAGS) params.c
	
energy_grid_mfpt.o : energy_grid_mfpt.c energy_grid_mfpt.h constants.h params.h
	$(CC) $(CCFLAGS) energy_grid_mfpt.c
	
RNAspectral : spectral_grid.o
	$(CC) spectral_grid.o -lRNA $(LDFLAGS) RNAspectral
	
spectral_grid.o : spectral_grid.c spectral_grid.h constants.h
	$(CC) $(CCFLAGS) spectral_grid.c

clean:
	rm -f *.o shared/*.a RNAmfpt RNAspectral

install:
	cp RNAmfpt $(BINDIR)
	cp shared/libMFPT.a $(LIBDIR)
	