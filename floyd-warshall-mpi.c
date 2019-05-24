#include <stdio.h>
#include <mpi.h>
#include <omp.h>
#include <math.h>
#include <limits.h>
#include <stdlib.h>

void create_matrix(int* distance, double PERCENT_INF, int N, int num_rows_per_process, int mpirank) {

  // hard-coded version
  /*
  int distance_vals[16] = {0, 500, -2, 500, 4, 0, 3, 500, 500, 500, 0, 2, 500, -1, 500, 0};
  //int* distance = (int*) malloc(N*N*sizeof(int));

  
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 4; j++) {
      int array_index = mpirank*num_rows_per_process*N + i*4 + j;
      int alloc_index = i*4 + j;
      *(distance+alloc_index) = distance_vals[array_index];
    }
  }
  */

  for (int i = 0; i < num_rows_per_process; i++) {
    for (int j = 0; j < N; j++) {
      int index = i*N + j;
      int i_in_entire_matrix = mpirank*num_rows_per_process + i;
      if (i_in_entire_matrix == j) {
	*(distance+index) = 0;
      }
      else {
	int prob = rand();
	if (prob < PERCENT_INF * RAND_MAX) {
	  *(distance+index) = INT_MAX;
	}
	else {
	  *(distance+index) = rand();
	}
      }
    }
  }
}

void floyd_warshall_mpi(int* distance, int N, int num_rows_per_process, int mpirank) {
  int* row_k = (int*) malloc(N*sizeof(int));
  int process_with_row_k, index;
  for (int k = 0; k < N; k++) {
    process_with_row_k = k / num_rows_per_process;
    if (mpirank == process_with_row_k) {
      index = k % num_rows_per_process;
      for (int i = 0; i < N; i++) {
	*(row_k+i) = *(distance+(index*N)+i);
      }
    }
    //MPI_Bcast(row_k, N, MPI_INT, process_with_row_k, MPI_COMM_WORLD);

    #pragma omp parallel for
    for (int i = 0; i < num_rows_per_process; i++) {
      for (int j = 0; j < N; j++) {
	int new_val = *(distance+(i*N)+k) + *(row_k+j);
	if (new_val < *(distance+(i*N)+j) && new_val > 0) {
	  *(distance+(i*N)+j) = new_val;
	}
      }
    }
  }
  free(row_k);
}

void print_matrix(int* distance) {
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 4; j++) {
      printf("%d, ", *(distance+i*4+j));
    }
    printf("\n");
  }
  printf("\n");
}

int main(int argc, char* argv[]) {
  //printf("threads: %d\n", omp_get_num_threads());
  int STARTING_N = 10000;
  int ENDING_N = 20000;
  int N_INCR = 5000;
  int STARTING_NUM_THREADS = 1;
  int ENDING_NUM_THREADS = 1;
  int THREAD_INCR = 1;
  int num_threads = STARTING_NUM_THREADS;
  int N = STARTING_N;
  double PERCENT_INF = 0.5;
  int mpirank, num_processes;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  
  while (N <= ENDING_N && num_threads <= ENDING_NUM_THREADS) {
    omp_set_num_threads(num_threads);

    int num_rows_per_process = N / num_processes;
    
    if ( (N % num_processes != 0) && mpirank == 0 ) {
      printf("N must be a factor of number of processes\n");
      MPI_Abort(MPI_COMM_WORLD, 0);
    }
    //MPI_Barrier(MPI_COMM_WORLD);
    
    int* distance = (int*) malloc(N*num_rows_per_process*sizeof(int));
    
    create_matrix(distance, PERCENT_INF, N, num_rows_per_process, mpirank);
    MPI_Barrier(MPI_COMM_WORLD);
    double tt = MPI_Wtime();
    floyd_warshall_mpi(distance, N, num_rows_per_process, mpirank);
    /*
      MPI_Barrier(MPI_COMM_WORLD);
      if (mpirank == 0) {
      printf("creation of rows 1 and 2:\n");
      print_matrix(distance);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      if (mpirank == 1) {
      printf("creation of rows 3 and 4:\n");
      print_matrix(distance);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      
      if (mpirank == 0) {
      printf("floyd-warshall on rows 1 and 2:\n");
      print_matrix(distance);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      if (mpirank == 1) {
      printf("floyd-warshall on rows 3 and 4:\n");
      print_matrix(distance);
      }
    */
    tt = MPI_Wtime() - tt;
    free(distance);
    
    MPI_Barrier(MPI_COMM_WORLD);
    if (mpirank == 0) {
      //printf("Completed Floyd-Warshall MPI in %f seconds\n", tt);
      printf("N = %d, num threads = %d, num processes = %d: %f seconds\n", N, num_threads, num_processes, tt);
    }
    
    if (N < ENDING_N) {
      N += N_INCR;
    }
    else {
      printf("\n");
      num_threads += THREAD_INCR;
      N = STARTING_N;
    }
    
  }

  MPI_Finalize();
  return 0;

  // parallel_floyd_warshall();
}
