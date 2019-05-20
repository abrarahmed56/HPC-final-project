all: fw fw-mpi

fw: floyd-warshall.c
	gcc -fopenmp -O3 -o floyd-warshall floyd-warshall.c

fw-mpi:  floyd-warshall-mpi.c
	mpicxx -fopenmp -O3 -o floyd-warshall-mpi floyd-warshall-mpi.c
