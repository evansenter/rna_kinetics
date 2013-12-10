#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include "constants.h"
#include "spectral_params.h"
#include "spectral_grid.h"

#define TIMING(start, stop, task) printf("Time in ms for %s: %.2f\n", task, (double)(((stop.tv_sec * 1000000 + stop.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec)) / 1000.0));

extern double temperature;
  
int main(int argc, char* argv[]) {
  #ifdef DEBUG
    struct timeval fullStart, fullStop, start, stop;
    gettimeofday(&fullStart, NULL);
    gettimeofday(&start, NULL);
  #endif
    
  SPECTRAL_PARAMS parameters;
  parameters = parse_args(argc, argv);
  
  int i, seq_length, from_index = -1, to_index = -1, empty_index = -1, num_structures = 0;
  double step_counter;
  double* transition_matrix;
  EIGENSYSTEM eigensystem;
  
  char* sequence  = parameters.sequence;
  seq_length      = strlen(sequence);
  char* empty_str = malloc((seq_length + 1) * sizeof(char));
  char* mfe_str   = malloc((seq_length + 1) * sizeof(char));
  
  for (i = 0; i < seq_length; ++i) {
    empty_str[i] = '.';
  }
  empty_str[i] = '\0';
  
  (double)fold(sequence, mfe_str);
  
  SOLUTION* all_structures = subopt(sequence, empty_str, 1000000, NULL);
  while (all_structures[num_structures].structure != NULL) {
    num_structures++;
  }
  
  for (i = 0; i < num_structures; ++i) {
    if (parameters.start_structure != NULL) {
      if (!strcmp(parameters.start_structure, all_structures[i].structure)) {
        from_index = i;
      }
    } else {
      if (!strcmp(empty_str, all_structures[i].structure)) {
        parameters.start_structure = all_structures[i].structure;
        from_index                 = i;
      }
    }
    
    if (parameters.end_structure != NULL) {
      if (!strcmp(parameters.end_structure, all_structures[i].structure)) {
        to_index = i;
      }
    } else {
      if (!strcmp(mfe_str, all_structures[i].structure)) {
        parameters.end_structure = all_structures[i].structure;
        to_index                 = i;
      }
    }
    
    if (!strcmp(empty_str, all_structures[i].structure)) {
      empty_index = i;
    }
  }
  
  if (parameters.verbose) {
    printf("sequence:\t%s\n", sequence);
    printf("start:\t\t%s\t(%d)\n", parameters.start_structure, from_index);
    printf("stop:\t\t%s\t(%d)\n", parameters.end_structure, to_index);
    printf("num str:\t%d\n", num_structures);
  }
  
  #ifdef SUPER_HEAVY_DEBUG
    for (i = 0; i < num_structures; ++i) {
      printf("%s\t%+.2f\t%d\n", all_structures[i].structure, all_structures[i].energy, i);
    }
  
    printf("\n");
  #endif
    
  #ifdef DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "initialization")
    gettimeofday(&start, NULL);
  #endif
  
  transition_matrix  = convert_structures_to_transition_matrix(all_structures, num_structures);
  
  #ifdef INSANE_DEBUG
    print_matrix("transition_matrix", transition_matrix, num_structures);
  #endif
  
  #ifdef DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "convert_structures_to_transition_matrix")
    gettimeofday(&start, NULL);
  #endif
  
  eigensystem = convert_transition_matrix_to_eigenvectors(transition_matrix, num_structures);
  
  #ifdef DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "convert_transition_matrix_to_eigenvectors")
    gettimeofday(&start, NULL);
  #endif
  
  invert_matrix(eigensystem, num_structures);
  
  #ifdef DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "invert_matrix")
    gettimeofday(&start, NULL);
  #endif
  
  #ifdef INSANE_DEBUG
    print_array("eigensystem.values", eigensystem.values, num_structures);
    print_matrix("eigensystem.vectors", eigensystem.vectors, num_structures);
    print_matrix("eigensystem.inverse_vectors", eigensystem.inverse_vectors, num_structures);
  #endif
    
  for (step_counter = 0; step_counter <= parameters.end_time; step_counter += parameters.step_size) {
    printf("%f\t%+.8f\n", step_counter, probability_at_time(eigensystem, step_counter, from_index, to_index, num_structures));
  }
  
  #ifdef DEBUG
    gettimeofday(&stop, NULL);
    TIMING(start, stop, "generate_population_points")
    gettimeofday(&fullStop, NULL);
    TIMING(fullStart, fullStop, "total")
  #endif
  
  return 0;
}
  
double* convert_structures_to_transition_matrix(SOLUTION* all_structures, int num_structures) {  
  int i, j;
  double col_sum;
  double* transition_matrix = malloc(num_structures * num_structures * sizeof(double));
  
  for (i = 0; i < num_structures; ++i) {
    col_sum = 0;
    
    for (j = 0; j < num_structures; ++j) {
      if (i != j) {
        transition_matrix[i + num_structures * j] = exp(-((double)all_structures[j].energy - (double)all_structures[i].energy) / RT);
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
