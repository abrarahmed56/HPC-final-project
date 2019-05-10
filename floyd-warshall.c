#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int* makeCopy(int* original_array, int N) {
  int* array_copy = (int*) malloc(N*N*sizeof(int));
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      *(array_copy+(i*N)+j) = *(original_array+(i*N)+j);
    }
  }
  return array_copy;
}

int* floyd_warshall_serial(int* distance, int N) {
  int* distanceCopy = makeCopy(distance, N);
  for (int k = 0; k < N; k++) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
	if (*(distanceCopy+(i*N)+k) + *(distanceCopy+(k*N)+j) < *(distanceCopy+(i*N)+j)) {
	  *(distanceCopy+(i*N)+j) = *(distanceCopy+(i*N)+k) + *(distanceCopy+(k*N)+j);
	}
      }
    }
  }
  return distanceCopy;
}

int* floyd_warshall_omp(int* distance, int N) {
  int* distanceCopy = makeCopy(distance, N);
  for (int k = 0; k < N; k++) {
#pragma omp parallel
    {
    #pragma omp for
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
	if (*(distanceCopy+(i*N)+k) + *(distanceCopy+(k*N)+j) < *(distanceCopy+(i*N)+j)) {
	  *(distanceCopy+(i*N)+j) = *(distanceCopy+(i*N)+k) + *(distanceCopy+(k*N)+j);
	}
      }
    }
  }
  }
  return distanceCopy;
}

void initialize_matrix(int* distance, int PERCENT_INF, int N) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int index = i*N + j;
      if (i == j) {
	*(distance+index) = 0;
      }
      else {
	int prob = rand();
	if (prob < PERCENT_INF * RAND_MAX) {
	  *(distance+index) = -1;
	}
	else {
	  *(distance+index) = rand();
	}
      }
    }
  }  
}

void time_functions(int* distance, int N) {
  clock_t start, end;
  double cpu_time_used_serial, cpu_time_used_omp;
  start = clock();
  int* serial_matrix = floyd_warshall_serial(distance, N);
  end = clock();
  cpu_time_used_serial = ((double) (end - start));
  printf("Time for serial code: %f\n", (cpu_time_used_serial/CLOCKS_PER_SEC));
  start = clock();
  int* parallel_matrix = floyd_warshall_omp(distance, N);
  end = clock();
  cpu_time_used_omp = ((double) (end - start));
  printf("Time for omp code: %f\n", cpu_time_used_omp/CLOCKS_PER_SEC);
  double error = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int index = i*N + j;
      error += fabs(*(serial_matrix+index) - *(parallel_matrix+index));
    }
  }
  printf("Error: %f\n", error);
  if (cpu_time_used_omp < cpu_time_used_serial) {
    printf("omp is faster\n");
  }
  else if (cpu_time_used_omp > cpu_time_used_serial) {
    printf("serial is faster\n");
  }
  else {
    printf("miraculously, same time for serial and omp\n");
  }
  free(serial_matrix);
  free(parallel_matrix);

}

int main (int argc, char *argv[]) 
{
  int N = 1000;
  int* distance = (int*) malloc(N*N*sizeof(int));
  double PERCENT_INF = 0;
  while (PERCENT_INF <= 1) {
    printf("About %.0f%% of the vertices are connected\n", (PERCENT_INF*100));
    initialize_matrix(distance, PERCENT_INF, N);
    time_functions(distance, N);
    PERCENT_INF += 0.25;
  }
  free(distance);
}
