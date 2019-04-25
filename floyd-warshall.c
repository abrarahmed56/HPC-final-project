#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double* makeCopy(double* original_array, int N) {
  double* array_copy = (double*) malloc(N*N*sizeof(double));
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      *(array_copy+(i*N)+j) = *(original_array+(i*N)+j);
    }
  }
  return array_copy;
}

double* floyd_warshall_serial(double* distance, int N) {
  double* distanceCopy = makeCopy(distance, N);
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

double* floyd_warshall_omp(double* distance, int N) {
  double* distanceCopy = makeCopy(distance, N);
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

int main (int argc, char *argv[]) 
{
  int N = 500;
  double* distance = (double*) malloc(N*N*sizeof(double));
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int index = i*N + j;
      if (i == j) {
	*(distance+index) = 0;
      }
      else {
	*(distance+index) = rand();
      }
    }
  }

  clock_t start, end;
  double cpu_time_used;
  start = clock();
  double* serial_matrix = floyd_warshall_serial(distance, N);
  end = clock();
  cpu_time_used = ((double) (end - start));
  printf("Time for serial code: %f\n", cpu_time_used);
  start = clock();
  double* parallel_matrix = floyd_warshall_omp(distance, N);
  end = clock();
  cpu_time_used = ((double) (end - start));
  printf("Time for omp code: %f\n", cpu_time_used);
  double error = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int index = i*N + j;
      error += fabs(*(serial_matrix+index) - *(parallel_matrix+index));
    }
  }
  printf("Error: %f\n", error);
  free(distance);
  free(serial_matrix);
  free(parallel_matrix);
}
