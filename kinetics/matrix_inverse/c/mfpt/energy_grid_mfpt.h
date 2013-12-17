#ifndef ENERGY_GRID_MFPT_H
#define ENERGY_GRID_MFPT_H

double** convert_energy_grid_to_transition_matrix(int**, int**, double**, int*, MFPT_PARAMETERS);
double compute_mfpt(int*, int*, double**, int, MFPT_PARAMETERS);
double* inverse(double*, int);
double* pseudoinverse(double*, int);
int number_of_permissible_single_bp_moves(int, int, int*, int*, int);
double transition_rate_from_probabilities(double, double, double);
double transition_rate_from_energies(double, double, double);
double transition_rate_from_probabilities_with_hastings(double, double, double, double);
double transition_rate_from_energies_with_hastings(double, double, double, double);

#endif
