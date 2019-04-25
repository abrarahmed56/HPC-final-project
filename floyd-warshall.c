#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void floyd_warshall_serial(double* distance, int N) {
  for (int k = 0; k < N; k++) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
	if (*(distance+(i*N)+k) + *(distance+(k*N)+j) < *(distance+(i*N)+j)) {
	  *(distance+(i*N)+j) = *(distance+(i*N)+k) + *(distance+(k*N)+j);
	}
      }
    }
  }
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%f, ", *(distance+(i*N)+j));
    }
    printf("\n");
  }
}

void floyd_warshall_omp(double* distance, int N) {
  for (int k = 0; k < N; k++) {
#pragma omp parallel
    {
    #pragma omp for
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
	if (*(distance+(i*N)+k) + *(distance+(k*N)+j) < *(distance+(i*N)+j)) {
	  *(distance+(i*N)+j) = *(distance+(i*N)+k) + *(distance+(k*N)+j);
	}
      }
    }
  }
  }
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%f, ", *(distance+(i*N)+j));
    }
    printf("\n");
  }
}

int main (int argc, char *argv[]) 
{
  int N = 4;
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
  floyd_warshall_serial(distance, N);
  end = clock();
  cpu_time_used = ((double) (end - start));
  printf("Time for serial code: %f\n", cpu_time_used);
  start = clock();
  floyd_warshall_omp(distance, N);
  end = clock();
  cpu_time_used = ((double) (end - start));
  printf("Time for omp code: %f\n", cpu_time_used);
}
