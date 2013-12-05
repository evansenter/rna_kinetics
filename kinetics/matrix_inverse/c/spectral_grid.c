#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "constants.h"
#include "spectral_grid.h"
  
int main(int argc, char* argv[]) {
  int i, seq_length, num_structures = 0;
  double** transition_matrix;
  
  char* sequence  = argv[1];
  seq_length      = strlen(sequence);
  char* empty_str = malloc((seq_length + 1) * sizeof(char));
  
  for (i = 0; i < seq_length; ++i) {
    empty_str[i] = '.';
  }
  empty_str[seq_length] = '\0';
  
  SOLUTION* all_structures = subopt(sequence, empty_str, 1000000, NULL);
  while (all_structures[num_structures].structure != NULL) {
    num_structures++;
  }
  
  printf("%s\n", argv[1]);
  printf("%s\n", empty_str);
  printf("%d\n", num_structures);
  
  for (i = 0; i < num_structures; ++i) {
    printf("%s\t%+.2f\n", all_structures[i].structure, all_structures[i].energy);
  }
  
  transition_matrix = convert_structures_to_transition_matrix(all_structures, num_structures);
  
  return 0;
}
  
double** convert_structures_to_transition_matrix(SOLUTION* all_structures, int num_structures) {  
  int i, j;
  double row_sum;
  double** transition_matrix = malloc(num_structures * sizeof(double*));
  
  for (i = 0; i < num_structures; ++i) {
    transition_matrix[i] = malloc(num_structures * sizeof(double));
    row_sum              = 0;
    
    for (j = 0; j < num_structures; ++j) {
      if (i != j) {
        transition_matrix[i][j] = MIN(1., exp(-((double)all_structures[j].energy - (double)all_structures[i].energy) / RT));
        row_sum                += transition_matrix[i][j];
      }
    }
    
    transition_matrix[i][i] = -row_sum;
  }
  
  for (i = 0; i < num_structures; ++i) {
    for (j = 0; j < num_structures; ++j) {
      printf("%+.2f\t", transition_matrix[i][j]);
    }
    printf("\n");
  }
  
  return transition_matrix;
}
