#ifndef PARAMS_H
#define PARAMS_H

typedef struct GlobalParameters {
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
} GlobalParameters;

GlobalParameters init_params();
GlobalParameters parse_args(int, char*[]);
int error_handling(GlobalParameters);
void debug_parameters(GlobalParameters);
void usage();

#endif
