#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include "../shared/constants.h"
#include "spectral_params.h"
#include "spectral_grid.h"
#include "spectral_functions.h"

extern SPECTRAL_PARAMS parameters;

double* convert_structures_to_transition_matrix(SOLUTION* all_structures, int num_structures) {  
  int i, j;
  double col_sum;
  double* transition_matrix = malloc(num_structures * num_structures * sizeof(double));
  
  for (i = 0; i < num_structures; ++i) {
    col_sum = 0;
    
    for (j = 0; j < num_structures; ++j) {
      if (i != j) {
        if (parameters.use_min) {
          transition_matrix[i + num_structures * j] = MIN(1, exp(-((double)all_structures[j].energy - (double)all_structures[i].energy) / RT));
        } else {
          transition_matrix[i + num_structures * j] = exp(-((double)all_structures[j].energy - (double)all_structures[i].energy) / RT);
        }
        
        col_sum += transition_matrix[i + num_structures * j];
        
        #ifdef INSANE_DEBUG
          printf("%d\t%d\t%.4e\n", i, j, transition_matrix[i + num_structures * j]);
        #endif
      }
    }
    
    #ifdef INSANE_DEBUG
      printf("%d col_sum:\t%.4e\n\n", i, col_sum);
    #endif
    
    transition_matrix[i * num_structures + i] = -col_sum;
  }
  
  return transition_matrix;
}

EIGENSYSTEM convert_transition_matrix_to_eigenvectors(double* transition_matrix, int num_structures) {
  int i, j;
  
  EIGENSYSTEM eigensystem = {
    .values          = malloc(num_structures * sizeof(double)),
    .vectors         = malloc(num_structures * num_structures * sizeof(double)),
    .inverse_vectors = malloc(num_structures * num_structures * sizeof(double))
  };

  gsl_matrix_view matrix_view      = gsl_matrix_view_array(transition_matrix, num_structures, num_structures);  
  gsl_vector_complex *eigenvalues  = gsl_vector_complex_alloc(num_structures);
  gsl_matrix_complex *eigenvectors = gsl_matrix_complex_alloc(num_structures, num_structures);
  
  gsl_eigen_nonsymmv_workspace *workspace = gsl_eigen_nonsymmv_alloc(num_structures);
  gsl_eigen_nonsymmv_params(1, workspace);
  gsl_eigen_nonsymmv(&matrix_view.matrix, eigenvalues, eigenvectors, workspace);
  gsl_eigen_nonsymmv_free(workspace);
  
  
  for (i = 0; i < num_structures; ++i) {
    eigensystem.values[i]               = GSL_REAL(gsl_vector_complex_get(eigenvalues, i));
    gsl_vector_complex_view eigenvector = gsl_matrix_complex_column(eigenvectors, i);
    
    for (j = 0; j < num_structures; ++j) {
      eigensystem.vectors[i + num_structures * j] = GSL_REAL(gsl_vector_complex_get(&eigenvector.vector, j));
    }
  }
  
  free(transition_matrix);
  gsl_vector_complex_free(eigenvalues);
  gsl_matrix_complex_free(eigenvectors);
  
  return eigensystem;
}

void invert_matrix(EIGENSYSTEM eigensystem, int num_structures) {
  int i, j, signum;
  
  gsl_matrix *matrix_to_invert = gsl_matrix_alloc(num_structures, num_structures);
  gsl_matrix *inversion_matrix = gsl_matrix_alloc(num_structures, num_structures);
  gsl_permutation *permutation = gsl_permutation_alloc(num_structures);
  
  for (i = 0; i < num_structures; ++i) {
    for (j = 0; j < num_structures; ++j) {
      gsl_matrix_set(matrix_to_invert, i, j, eigensystem.vectors[i * num_structures + j]);
    }
  }
  
  gsl_linalg_LU_decomp(matrix_to_invert, permutation, &signum);
  gsl_linalg_LU_invert(matrix_to_invert, permutation, inversion_matrix);
  
  for (i = 0; i < num_structures; ++i) {
    for (j = 0; j < num_structures; ++j) {
      eigensystem.inverse_vectors[i * num_structures + j] = gsl_matrix_get(inversion_matrix, i, j);
    }
  }
  
  gsl_matrix_free(matrix_to_invert);
  gsl_matrix_free(inversion_matrix);
  gsl_permutation_free(permutation);
}

double probability_at_time(EIGENSYSTEM eigensystem, double timepoint, int start_index, int target_index, int num_structures) {
  // This function is hard-wired to only consider the kinetics for folding from a distribution where p_{0}(start_index) == 1.
  
  int i;
  double cumulative_probability = 0;
  
  for (i = 0; i < num_structures; ++i) {
    cumulative_probability += 
      eigensystem.vectors[target_index * num_structures + i] * 
      eigensystem.inverse_vectors[i * num_structures + start_index] * 
      exp(eigensystem.values[i] * timepoint);
  }
  
  return cumulative_probability;
}

void print_array(char* title, double* matrix, int length) {
  int i;

  printf("%s\n", title);

  for (i = 0; i < length; ++i) {
    printf("%+.4f\n", matrix[i]);
  }

  printf("\n");
}


void print_matrix(char* title, double* matrix, int length) {
  int i, j;

  printf("%s\n", title);

  for (i = 0; i < length; ++i) {
    for (j = 0; j < length; ++j) {
      printf("%+.4f\t", matrix[i * length + j]);
    }
    printf("\n");
  }

  printf("\n");
}

int populate_arrays(char* file_path, int** k, int** l, double** p) {
  FILE *file = fopen(file_path, "r");
  int i, c, line_count = 0;
  char *token;
  char line[1024];
  
  if (file == NULL) {
    fprintf(stderr, "File not found.\n");
    fclose(file);
    return -1;
  }
  
  while ((c = fgetc(file)) != EOF) {
    if (c == '\n') {
      line_count++;
    }
  }
  
  *k = malloc(line_count * sizeof(int));
  *l = malloc(line_count * sizeof(int));
  *p = malloc(line_count * sizeof(double));
  
  rewind(file);
  
  while (fgets(line, 1024, file)) {
    token = strtok(line, ",");
    (*k)[i]  = atoi(token);
    token = strtok(NULL, ",");
    (*l)[i]  = atoi(token);
    token = strtok(NULL, ",");
    (*p)[i]  = atof(token);
    
    i++;
  }
  
  fclose(file);
  return line_count;
}
