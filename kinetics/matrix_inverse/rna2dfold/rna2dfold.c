// Call with: ./RNA2DFoldKinetics GGGGGCCCCC ".........." "(((....)))"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data_structures.h"

#define DEBUG 1
#define FLT_OR_DBL double
#define VRNA_GQUAD_MAX_STACK_SIZE 7
#define VRNA_GQUAD_MAX_LINKER_LENGTH 15

extern void read_parameter_file(const char energyfile[]);
extern TwoDpfold_vars *get_TwoDpfold_variables(const char *seq, const char *structure1, char *structure2, int circ);
extern TwoDpfold_solution *TwoDpfoldList(TwoDpfold_vars *vars, int distance1, int distance2);
extern void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
extern void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);

void computeTransitionMatrix(int *, int, double *, double **);
double computeMFPT(int *, int, double **, int, int);
void inverse(double*, int);

int main(int argc, char *argv[]) {
  int i, length, maxDistance1 = -1, maxDistance2 = -1;
  double temperature = 37., Q = 0., kT = (temperature + K0) * GASCONST / 1000.0; /* in Kcal */
  char *ParamFile = NULL, *sequence = argv[1], *structure1 = argv[2], *structure2 = argv[3];
  TwoDpfold_vars *q_vars;
  TwoDpfold_solution *pf_s;
  
  if (argc != 4) {
    printf("./RNA2DFoldKinetics SEQ STR1 STR2");
  }
  
  read_parameter_file("rna_turner1999.par");
  
  q_vars          = get_TwoDpfold_variables(sequence, structure1, structure2, 0);
  q_vars->dangles = 0;
  pf_s            = TwoDpfoldList(q_vars, -1, -1);

  for (i = 0, length = 0; pf_s[i].k != INF; i++, length++) {
    Q += pf_s[i].q;
  }
  
  #ifdef DEBUG
    printf("%s\t%f\n", "init_temp", q_vars->init_temp);
    printf("%s\t%f\n", "temperature", q_vars->temperature);
    printf("%s\t%d\n", "circ", q_vars->circ);
    printf("%s\t%d\n", "dangles", q_vars->dangles);
    printf("non-zero entries: %d\n", length);
  #endif

  for (i = 0; pf_s[i].k != INF; i++) {
    printf("%d\t%d\t%+.8f\n", pf_s[i].k, pf_s[i].l, (float)(pf_s[i].q) / (float)Q);
  }
  
  // for (i = 0; pf_s[i].k != INF; ++i) {
  //   nonZeroIndices[nonZeroCount++] = 2 * i + offset;
  //   
  //   if (x == 0 && y == inputStructureDist) {
  //     nonZeroA = 1;
  //   }
  //   
  //   if (x == inputStructureDist && y == 0) {
  //     nonZeroB = 1;
  //   }
  // }
  // 
  // inputStrsAreNonZero = nonZeroA && nonZeroB;
  // 
  // if (inputStrsAreNonZero) {
  //   transitionProbabilities = (double **)malloc(nonZeroCount * sizeof(double *));
  //   for (i = 0; i < nonZeroCount; ++i) {
  //     transitionProbabilities[i] = (double *)xcalloc(nonZeroCount, sizeof(double));
  //   }
  // 
  //   computeTransitionMatrix(nonZeroIndices, nonZeroCount, probabilities, transitionProbabilities);
  //   double mfpt = computeMFPT(nonZeroIndices, nonZeroCount, transitionProbabilities, inputStructureDist, rowLength);
  //   
  //   printf("%08f\n", mfpt);
  // 
  //   for (i = 0; i < nonZeroCount; ++i) {
  //     delete[] transitionProbabilities[i];
  //   }
  //   delete[] transitionProbabilities;
  // } else {
  //   printf("INFINITY\n");
  // }
  // 
  // delete[] nonZeroIndices;

  return 0;
}

void computeTransitionMatrix(int *nonZeroIndices, int nonZeroCount, double *probabilities, double **transitionProbabilities) {
  int i, j;
  double rowSum;
  
  for (i = 0; i < nonZeroCount; ++i) {
    rowSum = 0;
    
    for (j = 0; j < nonZeroCount; ++j) { 
      if (i != j) {
        // transitionProbabilities[i][j] = MIN2(
        //   1., 
        //   probabilities[nonZeroIndices[j]] / probabilities[nonZeroIndices[i]]
        // ) / (nonZeroCount - 1);
        
        rowSum += transitionProbabilities[i][j];
      }
    }
  
    transitionProbabilities[i][i] = abs(1 - rowSum);
  }
}

double computeMFPT(int *nonZeroIndices, int nonZeroCount, double **transitionProbabilities, int inputStructureDist, int rowLength) {
  int i, j, x, y, remappedIndexForStrA = -1, remappedIndexForStrB = -1, transitionMatrixRowLength = nonZeroCount - 1;
  double mpftFromAtoB;
  
  double *mfpt            = (double*)malloc(transitionMatrixRowLength * sizeof(double));
  double *inversionMatrix = (double*)malloc((int)pow((double)transitionMatrixRowLength, 2.) * sizeof(double));
  
  for (i = 0; i < nonZeroCount; ++i) {
    if (nonZeroIndices[i] == inputStructureDist * rowLength) {
      remappedIndexForStrB = i;
    } else if (nonZeroIndices[i] == inputStructureDist) {
      remappedIndexForStrA = i;
    }
  }
  
  if (remappedIndexForStrA < 0) {
    printf("Something has gone horribly wrong. We can not find which transition probabilities correspond to the starting state.");
    return -1;
  }
  
  if (remappedIndexForStrB < 0) {
    printf("Something has gone horribly wrong. We can not find which transition probabilities correspond to the finish state.");
    return -1;
  }
  
  // If remappedIndexForStrA > remappedIndexForStrB, we need to shift to the left by one because the remappedIndexForStrB
  // row / column is being removed.
  remappedIndexForStrA = (remappedIndexForStrA > remappedIndexForStrB ? remappedIndexForStrA - 1 : remappedIndexForStrA);
  
  #if MFPT_DEBUG
    printf("remappedIndexForStrA: %d\n", remappedIndexForStrA);
    printf("remappedIndexForStrB: %d\n", remappedIndexForStrB);
  
    printf("Mapping values:\nx\tRMOI(x)\n");
    for (i = 0; i < nonZeroCount; ++i) {
      printf("%d\t%d\n", i, nonZeroIndices[i]);
    }
    
    printf("Transition matrix values:\nx\ty\tRMOI(i)\tRMOI(j)\tinversionMatrix[x][y]\n");
  #endif
  
  for (i = 0; i < nonZeroCount; ++i) {
    for (j = 0; j < nonZeroCount; ++j) { 
      if (i != remappedIndexForStrB && j != remappedIndexForStrB) {
        x = (i > remappedIndexForStrB ? i - 1 : i);
        y = (j > remappedIndexForStrB ? j - 1 : j);
        
        // Be VERY careful changing anything here. nonZeroIndices has a mapping of all the positions in the energy grid (using 
        // row major ordering with a row length of rowLength) that have non-zero probabilities. Of these, we throw out anything
        // at base pair distance 0 from the second structure (the target of the MFPT calculation) and maximally distant from the
        // first structure. Because of this, there's a chunk of indices that need to get shifted to the left by one, to keep the
        // array tight (this is what x, y are doing). Hence, x and y are used for indexing into inversionMatrix and i, j are 
        // used for indexing into transitionProbabilities.
        inversionMatrix[x * transitionMatrixRowLength + y] = (i == j ? 1 - transitionProbabilities[i][j] : -transitionProbabilities[i][j]);
        
        #if MFPT_DEBUG
          printf("%d\t%d\t%d\t%d\t%+.08f\n", x, y, nonZeroIndices[i], nonZeroIndices[j], inversionMatrix[x * transitionMatrixRowLength + y]);
        #endif
      }
    }
  }
  
  inverse(inversionMatrix, transitionMatrixRowLength);
  
  for (i = 0; i < transitionMatrixRowLength; ++i) {
    for (j = 0; j < transitionMatrixRowLength; ++j) {
      mfpt[i] += inversionMatrix[i * transitionMatrixRowLength + j];
    }
  }
  
  #if MFPT_DEBUG
    for (i = 0; i < transitionMatrixRowLength; ++i) {
      printf("%d\t%+.08f\n", i, mfpt[i]);
    }
  #endif
    
  mpftFromAtoB = mfpt[remappedIndexForStrA];
  free(inversionMatrix);
  free(mfpt);
  
  return mpftFromAtoB;
}

void inverse(double* A, int N) {
  // http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
  int *IPIV = (int*)malloc((N + 1) * sizeof(int));
  int LWORK = N*N;
  double *WORK = (double*)malloc(LWORK * sizeof(double));
  int INFO;

  dgetrf_(&N,&N,A,&N,IPIV,&INFO);
  dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

  free(IPIV);
  free(WORK);
}