#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "spectral_grid.h"
  
int main(int argc, char* argv[]) {
  int i, seq_length, num_structures = 0;
  char* sequence  = argv[1];
  seq_length      = strlen(sequence);
  char* empty_str = (char*)malloc((seq_length + 1) * sizeof(char));
  
  for (i = 0; i < seq_length; ++i) {
    empty_str[i] = '.';
  }
  empty_str[seq_length] = '\0';
  
  SOLUTION* subopt_structures = subopt(sequence, empty_str, 1000000, NULL);
  while (subopt_structures[num_structures].structure != NULL) {
    num_structures++;
  }
  
  printf("%s\n", argv[1]);
  printf("%s\n", empty_str);
  printf("%d\n", num_structures);
  
  for (i = 0; i < num_structures; ++i) {
    printf("%s\t%+.2f\n", subopt_structures[i].structure, subopt_structures[i].energy);
  }
  
  return 0;
}
  
char** get_all_structures_for_sequence(char* sequence) {
  
}

double* get_all_energies_for_structures(char* sequence, char** structures) {
  
}
  
double** convert_sequence_to_transition_matrix(char* sequence) {
  
}
