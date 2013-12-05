char** get_all_structures_for_sequence(char* sequence);
double* get_all_energies_for_structures(char* sequence, char** structures);
double** convert_sequence_to_transition_matrix(char* sequence);

typedef struct {
  float energy;       
  char *structure;    
} SOLUTION;

extern SOLUTION* subopt(char*, char*, int, FILE*);