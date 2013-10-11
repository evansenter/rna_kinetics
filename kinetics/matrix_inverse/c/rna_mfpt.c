#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include "rna_mfpt.h"
#include "energy_grid_mfpt.h"

#define DEBUG 0

short ENERGY_BASED, TRANSITION_MATRIX_INPUT;

int main(int argc, char* argv[]) {
  unsigned long line_count;
  int i, j;
  int* k;
  int* l;
  double mfpt;
  double* p;
  double** transition_matrix;
  
  parse_args(argc, argv);
  
  line_count = count_lines(argv[argc - 1]);
  
  if (!line_count) {
    fprintf(stderr, "%s appears to have no data.\n", argv[argc - 1]);
    return 0;
  }
  
  k = (int*)malloc(line_count * sizeof(int));
  l = (int*)malloc(line_count * sizeof(int));
  p = (double*)malloc(line_count * sizeof(double));
  
  populate_arrays(argv[argc - 1], k, l, p);
  
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

void parse_args(int argc, char* argv[]) {
  int i;
  
  ENERGY_BASED            = 0;
  TRANSITION_MATRIX_INPUT = 0;
  
  if (argc < 2) {
    usage();
  }
  
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      if (strcmp(argv[i], "-E") == 0) {
        if (i == argc - 1) {
          usage();
        }
        ENERGY_BASED = 1;
      } else if (strcmp(argv[i], "-T") == 0) {
        if (i == argc - 1) {
          usage();
        }
        TRANSITION_MATRIX_INPUT = 1;
      } else {
        usage();
      }
    }
  }
}

void usage() {
  fprintf(stderr, "RNAmfpt [options] input_csv\n\n");
    
  fprintf(stderr, "where input_csv is a CSV file (with *no* header) of the format:\n");
  fprintf(stderr, "k_0,l_0,p_0\n");
  fprintf(stderr, "k_1,l_1,p_1\n");
  fprintf(stderr, "...\n");
  fprintf(stderr, "k_n,l_n,p_n\n\n");

  fprintf(stderr, "Options include the following:\n");
  fprintf(stderr, "-E\tenergy-based transitions, the default is disabled. If this flag is provided, the transition from state a to b will be calculated as (min(1, p_b - p_a) / n) rather than (min(1, p_b / p_a) / n)\n");
  fprintf(stderr, "-T\ttransition matrix input, the default is disabled. If this flag is provided, the input is expected to be a transition probability matrix, rather than a 2D energy grid. In this case, the first two columns in the CSV file are row-order indices into the transition probability matrix, and the third (final) column is the transition probability of that cell.\n");
  
  exit(0);
}
