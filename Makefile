CC=g++

# Configure number of threads
THREADS=8

# path to FLENS folder
FLENS=~/FLENS

CFLAGS=-I $(FLENS) -std=c++11 -lgsl -lgslcblas -O3 -fopenmp -DOMPTHREADS=$(THREADS) -Wall # -ggdb # -pg

%.o: %.cc *.h
	$(CC) -c -o $@ $< $(CFLAGS) 
	
pidefem: pidefem.o chi.o rhs_put.o rhs_call.o jump_mat.o diff_op.o
	$(CC) $^ -o $@ $(CFLAGS) 
	
clean:
	rm -f *.o pidefem
