#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include "constants.h"
#include "rna_mfpt.h"
#include "energy_grid_mfpt.h"

short ENERGY_BASED, TRANSITION_MATRIX_INPUT, PSEUDOINVERSE, SINGLE_BP_MOVES_ONLY, HASTINGS;
int START_STATE, END_STATE, SEQ_LENGTH;
double EPSILON, RT = 1e-3 * 1.9872041 * (273.15 + 37);

int main(int argc, char* argv[]) {
  unsigned long line_count, row_length;
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
  
  if (TRANSITION_MATRIX_INPUT) {
    row_length = 0;
    for (i = 0; i < line_count; ++i) {
      row_length = k[i] > row_length ? k[i] : row_length;
      row_length = l[i] > row_length ? l[i] : row_length;
    }
    row_length++;
    
    transition_matrix = (double**)malloc(row_length * sizeof(double*));
    for (i = 0; i < row_length + 1; ++i) {
      transition_matrix[i] = (double*)calloc(row_length, sizeof(double));
    }
    
    for (i = 0; i < line_count; ++i) {
      transition_matrix[k[i]][l[i]] = p[i];
    }
    
    #if HEAVY_DEBUG
      printf("Transition matrix:\n");
      printf("(x)\t(y)\tp(x -> y)\n");
      
      for (i = 0; i < line_count; ++i) {
        printf("(%d)\t=>\t(%d)\t%.8f\n", k[i], l[i], transition_matrix[k[i]][l[i]]);
      }
    #endif
  } else {
    row_length = line_count;
    transition_matrix = convertEnergyGridToTransitionMatrix(k, l, p, row_length);
    
    #if HEAVY_DEBUG
      printf("Transition matrix:\n");
      printf("i\tj\t(x, y)\t(a, b)\tp((x, y) -> (a, b))\n");

      for (i = 0; i < line_count; ++i) {
        for (j = 0; j < line_count; ++j) {
          printf("%d\t%d\t(%d, %d)\t=>\t(%d, %d)\t%.8f\n", i, j, k[i], l[i], k[j], l[j], transition_matrix[i][j]);
        }
  
        printf("\n");
      }
    #endif
  }
  
  mfpt = computeMFPT(k, l, transition_matrix, row_length, PSEUDOINVERSE ? &pseudoinverse : &inverse);

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
    
    if (!ENERGY_BASED && (p[i] < 0 || p[i] > 1)) {
      fprintf(stderr, "Error: line number %d (0-indexed) in the input doesn't satisfy 0 <= probability (%+1.2f) <= 1. Did you forget the -E flag?\n\n", i, p[i]);
      usage();
    }
    
    i++;
  }
  
  fclose(file);
}

void parse_args(int argc, char* argv[]) {
  int i, error = 0;
  
  ENERGY_BASED            = 0;
  TRANSITION_MATRIX_INPUT = 0;
  PSEUDOINVERSE           = 0;
  SINGLE_BP_MOVES_ONLY    = 0;
  HASTINGS                = 0;
  SEQ_LENGTH              = 0;
  EPSILON                 = 0;
  START_STATE             = -1;
  END_STATE               = -1;
  
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
      } else if (strcmp(argv[i], "-P") == 0) {
        if (i == argc - 1) {
          usage();
        }
        PSEUDOINVERSE = 1;
      } else if (strcmp(argv[i], "-X") == 0) {
        if (i == argc - 1) {
          usage();
        }
        SINGLE_BP_MOVES_ONLY = 1;
      } else if (strcmp(argv[i], "-H") == 0) {
        if (i == argc - 1) {
          usage();
        }
        HASTINGS = 1;
      } else if (strcmp(argv[i], "-A") == 0) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%d", &START_STATE)) {
          usage();
        } else if (START_STATE < 0 || (END_STATE >= 0 && START_STATE == END_STATE)) {
          usage();
        }
      } else if (strcmp(argv[i], "-Z") == 0) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%d", &END_STATE)) {
          usage();
        } else if (END_STATE < 0 || (START_STATE >= 0 && START_STATE == END_STATE)) {
          usage();
        }
      } else if (strcmp(argv[i], "-N") == 0) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%d", &SEQ_LENGTH)) {
          usage();
        } else if (SEQ_LENGTH <= 0) {
          usage();
        }
      } else if (strcmp(argv[i], "-O") == 0) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%lf", &EPSILON)) {
          usage();
        } else if (SEQ_LENGTH <= 0) {
          usage();
        }
      } else {
        usage();
      }
    }
  }
  
  if (SEQ_LENGTH && !EPSILON) {
    EPSILON = 1e-8;
  }
  
  #ifdef DEBUG
    printf("ENERGY_BASED\t\t%hd\n", ENERGY_BASED);
    printf("TRANSITION_MATRIX_INPUT\t%hd\n", TRANSITION_MATRIX_INPUT);
    printf("PSEUDOINVERSE\t\t%hd\n", PSEUDOINVERSE);
    printf("SINGLE_BP_MOVES_ONLY\t%hd\n", SINGLE_BP_MOVES_ONLY);
    printf("HASTINGS\t\t%hd\n", HASTINGS);
    printf("SEQ_LENGTH\t\t%d\n", SEQ_LENGTH);
    printf("EPSILON\t\t\t%.15f\n", EPSILON);
    printf("START_STATE\t\t%d\n", START_STATE);
    printf("END_STATE\t\t%d\n", END_STATE);
  #endif
  
  if (TRANSITION_MATRIX_INPUT && (HASTINGS || SEQ_LENGTH || EPSILON || SINGLE_BP_MOVES_ONLY)) {
    fprintf(stderr, "Error: If the -T flag is provided, -H, -N, -O and -X are not permitted!\n");
    error++;
  }
  
  if (TRANSITION_MATRIX_INPUT && !(START_STATE >= 0 && END_STATE >= 0)) {
    fprintf(stderr, "Error: If the -T flag is provided, -A and -Z must be explicitly set!\n");
    error++;
  }
  
  if (HASTINGS && !SINGLE_BP_MOVES_ONLY) {
    fprintf(stderr, "Error: If the -H flag is provided, -X must be explicitly set!\n");
    error++;
  }
  
  if (error) {
    fprintf(stderr, "\n");
    usage();
  }
}

void usage() {
  fprintf(stderr, "RNAmfpt [options] input_csv\n\n");
    
  fprintf(stderr, "where input_csv is a CSV file (with *no* header) of the format:\n");
  fprintf(stderr, "k_0,l_0,p_0\n");
  fprintf(stderr, "...,...,...\n");
  fprintf(stderr, "k_n,l_n,p_n\n\n");

  fprintf(stderr, "Options include the following:\n");
  
  fprintf(stderr, "-A\tstart state, the default is -1 (inferred from input data as the first row in the CSV whose entry in the first column is 0). If provided, should indicate the 0-indexed line in the input CSV file representing the start state.\n");
    
  fprintf(stderr, "-E\tenergy-based transitions, the default is disabled. If this flag is provided, the transition from state a to b will be calculated as (min(1, p_b - p_a) / n) rather than (min(1, p_b / p_a) / n).\n");
  
  fprintf(stderr, "-H\tHastings adjustment, the default is disabled. If this flag is provided, the input must be in the form of an energy grid, and only diagonally adjacent moves are permitted (in the all-to-all transition case, N(X) / N(Y) == 1). Calculating N(X) and N(Y) will respect grid boundaries and the triangle equality, and the basepair distance between the two structures for kinetics is inferred from the energy grid.\n");
  
  fprintf(stderr, "-N\tsequence length, the default is disabled. This flag represents the sequence length of the sequence on which kinetics is being performed. It is used in conjunction with the -O flag to ensure that the graph is fully connected. If -O is not explicitly set, it is taken to be 1e-8.\n");
  
  fprintf(stderr, "-O\tepsilon, the default is disabled. This flag should be a %%f-parseable epsilon value added to *all* positions in the energy grid (including valid positions not present in the input), which is then renormalized and used to ensure that the graph is fully connected.\n");
  
  fprintf(stderr, "-P\tpseudoinverse, the default is disabled. If this flag is provided, the Moore-Penrose pseudoinverse is computed for the transition probability matrix, rather than the true inverse.\n");
  
  fprintf(stderr, "-T\ttransition matrix input, the default is disabled. If this flag is provided, the input is expected to be a transition probability matrix, rather than a 2D energy grid. In this case, the first two columns in the CSV file are row-order indices into the transition probability matrix, and the third (final) column is the transition probability of that cell.\n");
  
  fprintf(stderr, "-X\tsingle basepair moves, the default is disabled. If this flag is provided, the input must be in the form of an energy grid, and only diagonally adjacent moves are permitted. This option makes the assumption that the input is *not* a transition probability matrix already, and the input energy grid already satisfies the triangle inequality / parity condition.\n");

  fprintf(stderr, "-Z\tend state, the default is -1 (inferred from input data as the first row in the CSV whose entry in the second column is 0). If provided, should indicate the 0-indexed line in the input CSV file representing the end state.\n");
  
  fprintf(stderr, "\nProgram returns -1 (resp. -2) if the start state (resp. end state) probability is 0. Otherwise returns the MFPT as predicted by matrix inversion.\n");
  
  exit(0);
}
