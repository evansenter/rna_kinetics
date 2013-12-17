#ifndef PARAMS_H
#define PARAMS_H

typedef struct {
  int start_state;
  int end_state;
  int sequence_length;
  int bp_dist;
  short energy_based;
  short transition_matrix_input;
  short pseudoinverse;
  short single_bp_moves_only;
  short hastings;
  short verbose;
  double additive_epsilon;
  double distributed_epsilon;
} MFPT_PARAMETERS;

MFPT_PARAMETERS init_mfpt_params();
MFPT_PARAMETERS parse_mfpt_args(int, char*[]);
int mfpt_error_handling(MFPT_PARAMETERS);
void debug_mfpt_parameters(MFPT_PARAMETERS);
void mfpt_usage();

#endif
