CCFLAGS           = -c -std=c99 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q
LDFLAGS           = -lm -lgomp -llapack -L/usr/local/include -llapacke -lgslcblas -lgsl -o
BINDIR            = ~/bin # Change this to the BINDIR
LIBDIR            = ~/lib # Change this to the LIBDIR
CC                = gcc-4.9
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

all: RNAmfpt

RNAmfpt: rna_mfpt.o mfpt_params.o energy_grid_mfpt.o mfpt_parser.o
	$(CC) rna_mfpt.o mfpt_params.o energy_grid_mfpt.o $(LDFLAGS) RNAmfpt
	ar cr ../shared/libMFPT.a mfpt_parser.o energy_grid_mfpt.o mfpt_params.o
	
rna_mfpt.o: rna_mfpt.c mfpt_parser.h energy_grid_mfpt.h ../shared/constants.h mfpt_params.h
	$(CC) $(CCFLAGS) rna_mfpt.c

mfpt_parser.o: mfpt_parser.c mfpt_parser.h
	$(CC) $(CCFLAGS) mfpt_parser.c
	
mfpt_params.o: mfpt_params.c mfpt_params.h
	$(CC) $(CCFLAGS) mfpt_params.c
	
energy_grid_mfpt.o: energy_grid_mfpt.c energy_grid_mfpt.h ../shared/constants.h mfpt_params.h
	$(CC) $(CCFLAGS) energy_grid_mfpt.c

clean:
	rm -f *.o ../shared/*.a RNAmfpt

install: RNAmfpt
	cp RNAmfpt $(BINDIR)
	cp ../shared/libMFPT.a $(LIBDIR)
	