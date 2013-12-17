#ifndef SPECTRAL_GRID_H
#define SPECTRAL_GRID_H

typedef struct {
  float energy;       
  char *structure;    
} SOLUTION;

extern float fold_par(char*, char*, paramT*, int, int);
extern SOLUTION* subopt_par(char*, char*, paramT*, int, int, int, FILE*);

#endif
