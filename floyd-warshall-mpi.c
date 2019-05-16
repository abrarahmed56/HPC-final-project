#include <stdio.h>
#include <mpi.h>

void create_matrix() {
  
}

int main(int argc, char* argv[]) {
  int N = 1000;
  int mpirank, num_processes;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  int num_rows_per_process = N / num_processes;

  if ( (N % p != 0) && mpirank == 0 ) {
    printf("N must be a factor of number of processes\n");
    MPI_Abort(MPI_COMM_WORLD, 0);
  }
  else {
    printf("works for now lol\n");
  }
  //int* distance = (int*) malloc(N*N*sizeof(int));

  // create_matrix();
  // parallel_floyd_warshall();
}
