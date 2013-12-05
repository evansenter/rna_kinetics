typedef struct {
  float energy;       
  char *structure;    
} SOLUTION;

extern SOLUTION* subopt(char*, char*, int, FILE*);

double** convert_structures_to_transition_matrix(SOLUTION*, int);