#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include "rna_mfpt.h"
#include "energy_grid_mfpt.h"

#define DEBUG 0

int main(int argc, char* argv[]) {
  unsigned long line_count;
  int i, j;
  int* k;
  int* l;
  double mfpt;
  double* p;
  double** transition_matrix;
  
  if (argc != 2) {
    fprintf(stderr, "Call with a csv file.\n");
    return 0;
  }
  
  line_count = count_lines(argv[1]);
  
  if (!line_count) {
    fprintf(stderr, "%s appears to have no data.\n", argv[1]);
    return 0;
  }
  
  k = (int*)malloc(line_count * sizeof(int));
  l = (int*)malloc(line_count * sizeof(int));
  p = (double*)malloc(line_count * sizeof(double));
  
  populate_arrays(argv[1], k, l, p);
  
  #if DEBUG
    printf("Input data:\n");
    for (i = 0; i < line_count; ++i) {
      printf("%d\t%d\t%.8f\n", k[i], l[i], p[i]);
    }
  #endif
    
  transition_matrix = convertEnergyGridToTransitionMatrix(p, line_count);
  
  #if DEBUG
    printf("Transition matrix:\n");
    printf("i\tj\t(x, y)\t(a, b)\tp((x, y) -> (a, b))\n");
  
    for (i = 0; i < line_count; ++i) {
      for (j = 0; j < line_count; ++j) {
        printf("%d\t%d\t(%d, %d)\t=>\t(%d, %d)\t%.8f\n", i, j, k[i], l[i], k[j], l[j], transition_matrix[i][j]);
      }
    
      printf("\n");
    }
  #endif
  
  #if DEBUG
    mfpt = computeMFPT(k, l, transition_matrix, line_count, 1);
  #else
    mfpt = computeMFPT(k, l, transition_matrix, line_count, 0);
  #endif

  printf("%.8f\n", mfpt);
  
  return 0;
}

unsigned long count_lines(char* file_path) {
  FILE *file = fopen(file_path, "r");
  int c;
  unsigned long line_count = 0;
  
  if (file == NULL) {
    fprintf(stderr, "File not found.\n");
    fclose(file);
    return 0;
  }
  
  while ((c = fgetc(file)) != EOF) {
    if (c == '\n') {
      line_count++;
    }
  }
  
  fclose(file);    
  
  return line_count;
}

void populate_arrays(char* file_path, int* k, int* l, double* p) {
  int i = 0;
  FILE *file = fopen(file_path, "r");
  char *token;
  char line[1024];
  
  while (fgets(line, 1024, file)) {
    token = strtok(line, ",");
    k[i] = atoi(token);
    token = strtok(NULL, ",");
    l[i] = atoi(token);
    token = strtok(NULL, ",");
    p[i] = atof(token);
    
    i++;
  }
  
  fclose(file);
}
