#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void floyd_warshall_serial() {
  int N = 4;
  //double* adjacency_matrix = (double*) malloc(N*N*sizeof(double));
  //hardcode for now
  double distance[4][4] = {{0, INFINITY, -2, INFINITY}, {4, 0, 3, INFINITY}, {INFINITY, INFINITY, 0, 2}, {INFINITY, -1, INFINITY, 0}};
  for (int k = 0; k < N; k++) {
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
	if (distance[i][k] + distance[k][j] < distance[i][j]) {
	  distance[i][j] = distance[i][k] + distance[k][j];
	}
      }
    }
  }
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%f, ", distance[i][j] );
    }
    printf("\n");
  }
}

void floyd_warshall_omp() {
  int N = 4;
  //double* adjacency_matrix = (double*) malloc(N*N*sizeof(double));
  //hardcode for now
  double distance[4][4] = {{0, INFINITY, -2, INFINITY}, {4, 0, 3, INFINITY}, {INFINITY, INFINITY, 0, 2}, {INFINITY, -1, INFINITY, 0}};
  for (int k = 0; k < N; k++) {
#pragma omp parallel
    {
    #pragma omp for
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
	if (distance[i][k] + distance[k][j] < distance[i][j]) {
	  distance[i][j] = distance[i][k] + distance[k][j];
	}
      }
    }
  }
  }
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%f, ", distance[i][j] );
    }
    printf("\n");
  }
}

int main (int argc, char *argv[]) 
{
  clock_t start, end;
  double cpu_time_used;
  start = clock();
  floyd_warshall_serial();
  end = clock();
  cpu_time_used = ((double) (end - start));
  printf("Time for serial code: %f\n", cpu_time_used);
  start = clock();
  floyd_warshall_omp();
  end = clock();
  cpu_time_used = ((double) (end - start));
  printf("Time for omp code: %f\n", cpu_time_used);
}
