typedef struct {
  int verbose;
  char* sequence;
  char* start_structure;
  char* end_structure;
  double start_time;
  double end_time;
  double step_size;
} SPECTRAL_PARAMS;

SPECTRAL_PARAMS init_params();
SPECTRAL_PARAMS parse_args(int, char*[]);
int error_handling(SPECTRAL_PARAMS);
void debug_parameters(SPECTRAL_PARAMS);
void usage();
