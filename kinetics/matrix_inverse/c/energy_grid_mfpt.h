#ifndef ENERGY_GRID_MFPT_H
#define ENERGY_GRID_MFPT_H

double** convertEnergyGridToTransitionMatrix(int*, int*, double*, int);
double computeMFPT(int*, int*, double**, int, double* (*)(double*, int));
double* inverse(double*, int);
double* pseudoinverse(double*, int);
int numSingleBpMoves(int, int, int*, int*, int, int);
double transitionRateFromProbabilities(double, double, double);
double transitionRateFromEnergies(double, double, double);
double transitionRateFromProbabilitiesWithHastings(double, double, double, double);
double transitionRateFromEnergiesWithHastings(double, double, double, double);

#endif
