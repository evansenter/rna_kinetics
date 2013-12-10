#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spectral_params.h"

extern double temperature;

SPECTRAL_PARAMS init_params() {
  SPECTRAL_PARAMS parameters = {
    .verbose         = 0,
    .sequence        = NULL,
    .start_structure = NULL,
    .end_structure   = NULL,
    .end_time        = 1e-1,
    .step_size       = 1e-4
  };
  
  return parameters;
}

SPECTRAL_PARAMS parse_args(int argc, char* argv[]) {
  int i;
  SPECTRAL_PARAMS parameters;
  
  if (argc < 1) {
    usage();
  }
  
  parameters = init_params();
  
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      if (!strcmp(argv[i], "--seq")) {
        if (i == argc - 1) {
          usage();
        } else {
          parameters.sequence = argv[++i];
        }
      } else if (!strcmp(argv[i], "--from")) {
        if (i == argc - 1) {
          usage();
        } else {
          parameters.start_structure = argv[++i];
        }
      } else if (!strcmp(argv[i], "--to")) {
        if (i == argc - 1) {
          usage();
        } else {
          parameters.end_structure = argv[++i];
        }
      } else if (!strcmp(argv[i], "--end-time")) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%lf", &(parameters.end_time))) {
          usage();
        }
      } else if (!strcmp(argv[i], "--step-size")) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%lf", &(parameters.step_size))) {
          usage();
        }
      } else if (!strcmp(argv[i], "--temperature")) {
        if (i == argc - 1) {
          usage();
        } else if (!sscanf(argv[++i], "%lf", &temperature)) {
          usage();
        }
      } else if (!strcmp(argv[i], "-v")) {
        parameters.verbose = 1;
      } else {
        usage();
      }
    }
  }
  
  if (parameters.verbose) {
    debug_parameters(parameters);
  }
    
  if (error_handling(parameters)) {
    usage();
  }
  
  return parameters;
}

int error_handling(SPECTRAL_PARAMS parameters) {
  int error = 0;
  
  if (parameters.sequence == NULL) {
    fprintf(stderr, "Error: No sequence provided.\n");
    error++;
  }
  
  // Also need error handling for the size of the structures, if provided.
  
  if (error) {
    fprintf(stderr, "\n");
  }
  
  return error;
}

void debug_parameters(SPECTRAL_PARAMS parameters) {
  printf("parameters.sequence\t\t%s\n",      parameters.sequence);
  printf("parameters.start_structure\t%s\n", parameters.start_structure == NULL ? "empty" : parameters.start_structure);
  printf("parameters.end_structure\t%s\n",   parameters.end_structure == NULL ? "mfe" : parameters.end_structure);
  printf("parameters.end_time\t\t%.2e\n",    parameters.end_time);
  printf("parameters.step_size\t\t%.2e\n",   parameters.step_size);
  printf("temperature\t\t\t%.1f\n",          temperature);
}

void usage() {
  fprintf(stderr, "RNAspectral [options]\n\n");
    
  exit(0);
}
