CFLAGS = -Wall -fPIC -c -O3 -fopenmp
LDFLAGS = -lm -lgomp -lRNA
CC = gcc

a.out : rnaeval_r_wrapper.o
	$(CC) *.o -shared -o rnaeval_r.so $(LDFLAGS)

rnaeval_r_wrapper.o : rnaeval_r_wrapper.c
	$(CC) $(CFLAGS) rnaeval_r_wrapper.c

clean :
	 rm -f *.o *.so
