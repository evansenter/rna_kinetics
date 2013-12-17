#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "spectral_params.h"

extern double temperature;

SPECTRAL_PARAMS init_spectral_params() {
  SPECTRAL_PARAMS parameters = {
    .verbose          = 0,
    .sequence         = NULL,
    .energy_grid_file = NULL,
    .start_structure  = NULL,
    .end_structure    = NULL,
    .temperature      = 37.,
    .start_time       = -11,
    .end_time         = 0,
    .step_size        = 1e-1,
    .lonely_bp        = 0,
    .energy_cap       = 1,
    .use_min          = 1,
    .eigen_only       = 0
  };
  
  return parameters;
}

SPECTRAL_PARAMS parse_spectral_args(int argc, char* argv[]) {
  int i;
  SPECTRAL_PARAMS parameters;
  
  if (argc < 1) {
    spectral_usage();
  }
  
  parameters = init_spectral_params();
  
  for (i = 1; i < argc; i++) {
    if (argv[i][0] == '-') {
      if (!strcmp(argv[i], "-seq")) {
        if (i == argc - 1) {
          spectral_usage();
        } else {
          parameters.sequence = argv[++i];
        }
      } else if (!strcmp(argv[i], "-from")) {
        if (i == argc - 1) {
          spectral_usage();
        } else {
          parameters.start_structure = argv[++i];
        }
      } else if (!strcmp(argv[i], "-to")) {
        if (i == argc - 1) {
          spectral_usage();
        } else {
          parameters.end_structure = argv[++i];
        }
      } else if (!strcmp(argv[i], "-energy_grid_file")) {
        if (i == argc - 1) {
          spectral_usage();
        } else {
          parameters.energy_grid_file = argv[++i];
        }
      } else if (!strcmp(argv[i], "-start_time")) {
        if (i == argc - 1) {
          spectral_usage();
        } else if (!sscanf(argv[++i], "%lf", &(parameters.start_time))) {
          spectral_usage();
        }
      } else if (!strcmp(argv[i], "-end_time")) {
        if (i == argc - 1) {
          spectral_usage();
        } else if (!sscanf(argv[++i], "%lf", &(parameters.end_time))) {
          spectral_usage();
        }
      } else if (!strcmp(argv[i], "-step_size")) {
        if (i == argc - 1) {
          spectral_usage();
        } else if (!sscanf(argv[++i], "%lf", &(parameters.step_size))) {
          spectral_usage();
        }
      } else if (!strcmp(argv[i], "-temperature")) {
        if (i == argc - 1) {
          spectral_usage();
        } else if (!sscanf(argv[++i], "%lf", &(parameters.temperature))) {
          spectral_usage();
        }
        temperature = parameters.temperature;
      } else if (!strcmp(argv[i], "-lonely_bp")) {
        parameters.lonely_bp = 1;
      } else if (!strcmp(argv[i], "-no_energy_cap")) {
        parameters.energy_cap = 0;
      } else if (!strcmp(argv[i], "-no_min_transition")) {
        parameters.use_min = 0;
      } else if (!strcmp(argv[i], "-eigen_only")) {
        parameters.eigen_only = 1;
      } else if (!strcmp(argv[i], "-verbose")) {
        parameters.verbose = 1;
      } else {
        spectral_usage();
      }
    }
  }
  
  if (parameters.verbose) {
    debug_spectral_parameters(parameters);
  }
    
  if (spectral_error_handling(parameters)) {
    spectral_usage();
  }
  
  return parameters;
}

int spectral_error_handling(SPECTRAL_PARAMS parameters) {
  int error = 0;
  
  // if (parameters.sequence == NULL) {
  //   fprintf(stderr, "Error: No sequence provided.\n");
  //   error++;
  // }
  // 
  // Also need error handling for the size of the structures, if provided.
  
  if (error) {
    fprintf(stderr, "\n");
  }
  
  return error;
}

void debug_spectral_parameters(SPECTRAL_PARAMS parameters) {
  printf("parameters.sequence\t\t%s\n",       parameters.sequence);
  printf("parameters.start_structure\t%s\n",  parameters.start_structure == NULL ? "empty" : parameters.start_structure);
  printf("parameters.end_structure\t%s\n",    parameters.end_structure == NULL ? "mfe" : parameters.end_structure);
  printf("parameters.energy_grid_file\t%s\n", parameters.energy_grid_file == NULL ? "none" : parameters.energy_grid_file);
  printf("parameters.temperature\t\t%.1f\n",  parameters.temperature);
  printf("parameters.start_time\t\t%.2e\n",   parameters.start_time);
  printf("parameters.end_time\t\t%.2e\n",     parameters.end_time);
  printf("parameters.step_size\t\t%.2e\n",    parameters.step_size);
  printf("parameters.lonely_bp\t\t%s\n",      parameters.lonely_bp ? "No" : "Yes");
  printf("parameters.energy_cap\t\t%s\n",     parameters.energy_cap ? "Yes" : "No");
  printf("parameters.eigen_only\t\t%s\n",     parameters.eigen_only ? "Yes" : "No");
  printf("parameters.use_min\t\t%s\n",        parameters.use_min ? "MIN(1, exp(-(E(j) - E(i)) / RT))" : "exp(-(E(j) - E(i)) / RT)");
  printf("temperature\t\t\t%.1f\n",           temperature);
}

void spectral_usage() {
  fprintf(stderr, "RNAspectral [options]\n\n");
    
  exit(0);
}
