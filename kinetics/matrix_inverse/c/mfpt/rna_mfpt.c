#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../shared/constants.h"
#include "params.h"
#include "mfpt_parser.h"
#include "energy_grid_mfpt.h"

int main(int argc, char* argv[]) {
  int i, j, line_count, row_length;
  int* k;
  int* l;
  double mfpt;
  double* p;
  double** transition_matrix;
  GlobalParameters parameters;
  
  parameters = parse_args(argc, argv);
  line_count = count_lines(argv[argc - 1]);
  
  if (!line_count) {
    fprintf(stderr, "%s appears to have no data.\n", argv[argc - 1]);
    return 0;
  }
  
  k = malloc(line_count * sizeof(int));
  l = malloc(line_count * sizeof(int));
  p = malloc(line_count * sizeof(double));
  
  populate_arrays(argv[argc - 1], k, l, p, parameters);
  
  #ifdef SUPER_HEAVY_DEBUG
    printf("\nInput data:\n");
    for (i = 0; i < line_count; ++i) {
      printf("%d\t%d\t%.8f\n", k[i], l[i], p[i]);
    }
  #endif
  
  // We already have a transition matrix, this is the easy case. Just need to find MFPT.
  if (parameters.transition_matrix_input) {
    // We need to infer the dimensions of the transition matrix.
    row_length = 0;
    for (i = 0; i < line_count; ++i) {
      row_length = k[i] > row_length ? k[i] : row_length;
      row_length = l[i] > row_length ? l[i] : row_length;
    }
    row_length++;
    
    transition_matrix = malloc(row_length * sizeof(double*));
    for (i = 0; i < row_length + 1; ++i) {
      transition_matrix[i] = calloc(row_length, sizeof(double));
    }
    
    for (i = 0; i < line_count; ++i) {
      transition_matrix[k[i]][l[i]] = p[i];
    }
    
    #if SUPER_HEAVY_DEBUG
      printf("Transition matrix:\n");
      printf("(x)\t(y)\tp(x . y)\n");
      
      for (i = 0; i < line_count; ++i) {
        printf("(%d)\t=>\t(%d)\t%.8f\n", k[i], l[i], transition_matrix[k[i]][l[i]]);
      }
    #endif
  // We have an energy grid, this requires converting the energy grid into a transition matrix data structure before finding MFPT.
  } else {
    row_length        = line_count;
    transition_matrix = convert_energy_grid_to_transition_matrix(&k, &l, &p, &row_length, parameters);
    
    #if SUPER_HEAVY_DEBUG
      printf("Transition matrix:\n");
      printf("i\tj\t(x, y)\t(a, b)\tp((x, y) . (a, b))\n");

      for (i = 0; i < row_length; ++i) {
        for (j = 0; j < row_length; ++j) {
          printf("%d\t%d\t(%d, %d)\t=>\t(%d, %d)\t%.8f\n", i, j, k[i], l[i], k[j], l[j], transition_matrix[i][j]);
        }
  
        printf("\n");
      }
    #endif
  }
  
  mfpt = compute_mfpt(k, l, transition_matrix, row_length, parameters);
  printf("%+.8f\n", mfpt);
  
  free(transition_matrix);
  free(k);
  free(l);
  free(p);
  
  return 0;
}

int count_lines(char* file_path) {
  FILE *file = fopen(file_path, "r");
  int c;
  int line_count = 0;
  
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

void populate_arrays(char* file_path, int* k, int* l, double* p, GlobalParameters parameters) {
  int i = 0;
  FILE *file = fopen(file_path, "r");
  char *token;
  char line[1024];
  
  while (fgets(line, 1024, file)) {
    token = strtok(line, ",");
    k[i]  = atoi(token);
    token = strtok(NULL, ",");
    l[i]  = atoi(token);
    token = strtok(NULL, ",");
    p[i]  = atof(token);
    
    if (!parameters.energy_based && (p[i] < 0 || p[i] > 1)) {
      fprintf(stderr, "Error: line number %d (0-indexed) in the input doesn't satisfy 0 <= probability (%+1.2f) <= 1. Did you forget the -E flag?\n\n", i, p[i]);
      usage();
    }
    
    i++;
  }
  
  fclose(file);
}
