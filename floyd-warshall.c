#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

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
	int new_val = *(distanceCopy+(i*N)+k) + *(distanceCopy+(k*N)+j);
	if (new_val < *(distanceCopy+(i*N)+j) && new_val > 0) {
	  *(distanceCopy+(i*N)+j) = new_val;
	}
      }
    }
  }
  return distanceCopy;
}

int* floyd_warshall_omp(int* distance, int N) {
  int* distanceCopy = makeCopy(distance, N);
  for (int k = 0; k < N; k++) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
      //printf("thread %d of %d\n", omp_get_thread_num(), omp_get_num_threads());
      for (int j = 0; j < N; j++) {
	int new_val = *(distanceCopy+(i*N)+k) + *(distanceCopy+(k*N)+j);
	if (new_val < *(distanceCopy+(i*N)+j) && new_val > 0) {
	  *(distanceCopy+(i*N)+j) = new_val;
	}
      }
    }
  }
  return distanceCopy;
}

void initialize_matrix(int* distance, double PERCENT_INF, int N) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int index = i*N + j;
      if (i == j) {
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

void print_matrix(int* distance, int N) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      printf("%d, ", *(distance+i*N+j));
    }
    printf("\n");
  }
  printf("\n");
}

void write_matrix(int* distance, int N, FILE* fd) {
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      fprintf(fd, "%d, ", *(distance+i*N+j));
    }
    fprintf(fd, "\n");
  }
  fprintf(fd, "\n");
}

void time_functions(int* distance, int N, double PERCENT_INF) {
  clock_t start_serial, end_serial;
  double start_omp, end_omp;
  double cpu_time_used_serial, cpu_time_used_omp;

  start_serial = clock();
  int* serial_matrix = floyd_warshall_serial(distance, N);
  end_serial = clock();
  cpu_time_used_serial = ((double) (end_serial - start_serial));
  cpu_time_used_serial = cpu_time_used_serial/CLOCKS_PER_SEC;
  printf("Time for serial code: %f\n", (cpu_time_used_serial));

  FILE* serial_matrix_file = NULL;
  char filename_serial[256];
  snprintf(filename_serial, 256, "serial_matrix_%f_percent.txt", PERCENT_INF);
  serial_matrix_file = fopen(filename_serial, "w");
  if (serial_matrix_file == NULL) {
    printf("Error opening file\n");
  }
  write_matrix(serial_matrix, N, serial_matrix_file);
  fclose(serial_matrix_file);

  start_omp = omp_get_wtime();
  int* parallel_matrix = floyd_warshall_omp(distance, N);
  end_omp = omp_get_wtime();
  cpu_time_used_omp = end_omp - start_omp;

  printf("Time for omp code: %f\n", cpu_time_used_omp);

  FILE* omp_matrix_file = NULL;
  char filename_omp[256];
  snprintf(filename_omp, 256, "omp_matrix_%f_percent.txt", PERCENT_INF);
  omp_matrix_file = fopen(filename_omp, "w");
  if (omp_matrix_file == NULL) {
    printf("Error opening file\n");
  }
  write_matrix(parallel_matrix, N, omp_matrix_file);

  fclose(omp_matrix_file);
  int error = 0;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int index = i*N + j;
      error += abs(*(serial_matrix+index) - *(parallel_matrix+index));
    }
  }

  printf("Error: %d\n", error);

  if (cpu_time_used_omp < cpu_time_used_serial) {
    double factor = cpu_time_used_serial / cpu_time_used_omp;
    printf("omp is faster by a factor of %f\n", factor);
  }
  else if (cpu_time_used_omp > cpu_time_used_serial) {
    double factor = cpu_time_used_omp / cpu_time_used_serial;
    printf("serial is faster by a factor of %f\n", factor);
  }
  else {
    printf("miraculously, same time for serial and omp\n");
  }
  free(serial_matrix);
  free(parallel_matrix);

}

int main (int argc, char *argv[]) 
{
  omp_set_num_threads(4);
  int N = 2000;//1000000;

  // In case we need to re-evaluate with a hard-coded graph:
  /*
  int N = 4;

  //int distance_vals[16] = {0, INF, -2, INF, 4, 0, 3, INF, INF, INF, 0, 2, INF, -1, INF, 0};
  int distance_vals[16] = {0, 1, 4, 1, 3, 0, 2, 500, 4, 1, 0, 3, 2, 500, 500, 0};
  int* distance = (int*) malloc(N*N*sizeof(int));

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      int index = i*4 + j;
      *(distance+index) = distance_vals[index];
    }
  }
  */

  int* distance = (int*) malloc(N*N*sizeof(int));
  double PERCENT_INF = 0;
  while (PERCENT_INF <= 1) {
    printf("About %.0f%% of the vertices are connected\n", (PERCENT_INF*100));
    initialize_matrix(distance, PERCENT_INF, N);
    FILE* initial_matrix_file = NULL;
    char filename[256];
    snprintf(filename, 256, "initial_matrix_%f_percent.txt", PERCENT_INF);
    initial_matrix_file = fopen(filename, "w");
    if (initial_matrix_file == NULL) {
      printf("Error opening file\n");
      return 1;
    }
    write_matrix(distance, N, initial_matrix_file);
    fclose(initial_matrix_file);

    time_functions(distance, N, PERCENT_INF);
    PERCENT_INF += 0.25;
  }
  free(distance);
}
