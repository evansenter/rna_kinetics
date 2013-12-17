typedef struct {
  float energy;       
  char *structure;    
} SOLUTION;

typedef struct {
  double *values;       
  double *vectors;
  double *inverse_vectors;
} EIGENSYSTEM;

extern float fold_par(char*, char*, paramT*, int, int);
extern SOLUTION* subopt_par(char*, char*, paramT*, int, int, int, FILE*);

double* convert_structures_to_transition_matrix(SOLUTION*, int);
EIGENSYSTEM convert_transition_matrix_to_eigenvectors(double*, int);
void invert_matrix(EIGENSYSTEM, int);
double probability_at_time(EIGENSYSTEM, double, int, int, int);
void print_array(char*, double*, int);
void print_matrix(char*, double*, int);
int populate_arrays(char*, int**, int**, double**);
