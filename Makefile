# Makefile for computing FFTbor

CFLAGS = -c -O0 # -O3
LDFLAGS = -lfftw3 -lm -L.
PG = #-pg
BINDIR = /usr/local/bin
CC = g++

RNAbor: delta.o misc.o main.o
	 $(CC) $(PG) $(LDFLAGS) -g delta.o misc.o main.o -lRNA -o RNAbor
	 
main.o: main.c delta.h
	   $(CC) -Wall -W $(PG) $(CFLAGS) $(H) main.c

misc.o: misc.c misc.h
	 $(CC) -Wall -W $(PG) -g $(CFLAGS) $(H) misc.c

delta.o: delta.c delta.h
	$(CC) -Wall -W $(PG) -g $(CFLAGS) $(H) delta.c

clean:
	rm -f *.o RNAbor

install:
	cp RNAbor $(BINDIR)
	cp energy.par $(BINDIR)
