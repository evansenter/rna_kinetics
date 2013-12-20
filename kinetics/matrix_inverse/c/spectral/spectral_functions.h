#ifndef SPECTRAL_FUNCTIONS_H
#define SPECTRAL_FUNCTIONS_H

typedef struct {
  double *values;       
  double *vectors;
  double *inverse_vectors;
} EIGENSYSTEM;

double* convert_structures_to_transition_matrix(SOLUTION*, int, int);
EIGENSYSTEM convert_transition_matrix_to_eigenvectors(double*, int);
void invert_matrix(EIGENSYSTEM, int);
double probability_at_time(EIGENSYSTEM, double, int, int, int);
void print_array(char*, double*, int);
void print_matrix(char*, double*, int);
void free_eigensystem(EIGENSYSTEM);
int populate_arrays(char*, int**, int**, double**);

#endif