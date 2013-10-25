#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "constants.h"
#include "energy_grid_mfpt.h"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define SMLSIZ 25

#ifdef __cplusplus
  extern "C" {
    void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
    void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);  
    int dgelsd_(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *s, double *rcond, int *rank, double *work, int *lwork, int *iwork, int *info);
  }
#else
  extern void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  extern void dgetri_(int* N, double* A, int* lda, int* IPIV, double* WORK, int* lwork, int* INFO);
  extern int dgelsd_(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *s, double *rcond, int *rank, double *work, int *lwork, int *iwork, int *info);
#endif
  
extern double RT, EPSILON;
extern short ENERGY_BASED, SINGLE_BP_MOVES_ONLY, HASTINGS;
extern int START_STATE, END_STATE, SEQ_LENGTH;

double** convertEnergyGridToTransitionMatrix(int* k, int* l, double* p, int length) {
  int i, j, bpDist, numX, numY, distFromK = -1, distFromL = -1;
  double rowSum;
  double* numAdjacentMoves         = (double*)malloc(length * sizeof(double));
  double** transitionProbabilities = (double**)malloc(length * sizeof(double*));
  
  if (SINGLE_BP_MOVES_ONLY) {
    for (i = 0; i < length; ++i) {
      if (k[i] == 0) {
        distFromL = l[i];
      }
      
      if (l[i] == 0) {
        distFromK = k[i];
      }
    }
    
    if (distFromK == distFromL && distFromK >= 0) {
      bpDist = distFromK;
    } else {
      fprintf(stderr, "Can't infer the input structure distances for the energy grid. We found (0, %d) and (%d, 0).\n", distFromL, distFromK);
      exit(0);
    }
    
    for (i = 0; i < length; ++i) {
      numAdjacentMoves[i] = (double)numSingleBpMoves(k[i], l[i], k, l, bpDist, length);
      
      #ifdef DEBUG
        printf("numSingleBpMoves(%d, %d, bpDist = %d, length = %d):\t%f\n", k[i], l[i], bpDist, length, numAdjacentMoves[i]);
      #endif
    }
  }
  
  for (i = 0; i < length; ++i) {
    rowSum                     = 0.;
    transitionProbabilities[i] = (double*)calloc(length, sizeof(double));
      
    for (j = 0; j < length; ++j) {
      if (i != j) {
        if (SINGLE_BP_MOVES_ONLY) {
          if ((int)abs(k[i] - k[j]) == 1 && (int)abs(l[i] - l[j]) == 1) {
            #ifdef HEAVY_DEBUG
              printf("Single bp move: %d, %d => %d, %d\n", k[i], l[i], k[j], l[j]);
            #endif
            
            if (HASTINGS) {
              if (ENERGY_BASED) {
                transitionProbabilities[i][j] = transitionRateFromEnergiesWithHastings(p[i], p[j], numAdjacentMoves[i], numAdjacentMoves[j]);
              } else {
                transitionProbabilities[i][j] = transitionRateFromProbabilitiesWithHastings(p[i], p[j], numAdjacentMoves[i], numAdjacentMoves[j]);
              }
            } else {
              if (ENERGY_BASED) {
                transitionProbabilities[i][j] = transitionRateFromEnergies(p[i], p[j], numAdjacentMoves[i]);
              } else {
                transitionProbabilities[i][j] = transitionRateFromProbabilities(p[i], p[j], numAdjacentMoves[i]);
              }
            }
          }
        } else {
          if (ENERGY_BASED) {
            transitionProbabilities[i][j] = transitionRateFromEnergies(p[i], p[j], (double)(length - 1));
          } else {
            transitionProbabilities[i][j] = transitionRateFromProbabilities(p[i], p[j], (double)(length - 1));
          }
        }
        
        rowSum += transitionProbabilities[i][j];
      }
    }
    
    transitionProbabilities[i][i] = 1 - rowSum;
  }
  
  return transitionProbabilities;
}

double computeMFPT(int* k, int* l, double **transitionProbabilities, int length, double* (*invert)(double*, int)) {
  int i, j, x, y, startIndex, endIndex, inversionMatrixRowLength = length - 1;
  double mfptFromStart, rowSum;
  
  if (START_STATE == -1) {
    for (i = 0, startIndex = -1; i < length; ++i) {
      if (k[i] == 0) {
        startIndex = i;
      }
    }
  } else {
    startIndex = START_STATE;
  }
  
  if (END_STATE == -1) {
    for (i = 0, endIndex = -1; i < length; ++i) {
      if (l[i] == 0) {
        endIndex = i;
      }
    }
  } else {
    endIndex = END_STATE;
  }
  
  #ifdef DEBUG
    printf("startIndex:\t%d\n", startIndex);
    printf("endIndex:\t%d\n", endIndex);
  #endif
  
  if (startIndex < 0) {
    #ifdef DEBUG
      fprintf(stderr, "We can not find any position in the energy grid correspondent to the starting state.\n");
    #endif
    return -1;
  }
  
  if (endIndex < 0) {
    #ifdef DEBUG
      fprintf(stderr, "We can not find any position in the energy grid correspondent to the stopping state.\n");
    #endif
    return -2;
  }
  
  // If startIndex > endIndex, we need to shift to the left by one because the endIndex row / column is being removed.
  if (startIndex > endIndex) {
    startIndex--;
  }
  
  double *mfpt            = (double*)calloc(inversionMatrixRowLength, sizeof(double));
  double *inversionMatrix = (double*)malloc((int)pow((double)inversionMatrixRowLength, 2.) * sizeof(double));
  
  #ifdef HEAVY_DEBUG
    printf("Inversion matrix:\n");
    printf("i\tj\tx\ty\tinversionMatrix[x, y]\n");
  #endif
  
  for (i = 0; i < length; ++i) {
    for (j = 0; j < length; ++j) { 
      if (i != endIndex && j != endIndex) {
        x = (i > endIndex ? i - 1 : i);
        y = (j > endIndex ? j - 1 : j);
        
        // Be VERY careful changing anything here. We throw out anything at base pair distance 0 (endIndex) from the second structure (the target of the MFPT calculation) and maximally distant from the first structure. Because of this, there's a chunk of indices that need to get shifted to the left by one, to keep the array tight (this is what x, y are doing). Hence, x and y are used for indexing into inversionMatrix and i, j are used for indexing into transitionProbabilities.
        inversionMatrix[x * inversionMatrixRowLength + y] = (i == j ? 1 - transitionProbabilities[i][j] : -transitionProbabilities[i][j]);
        
        #ifdef HEAVY_DEBUG
          printf("%d\t%d\t%d\t%d\t%f\n", i, j, x, y, inversionMatrix[x * inversionMatrixRowLength + y]);
        #endif
      }
    }
    
    #ifdef HEAVY_DEBUG
      printf("\n");
    #endif
  }
  
  #ifdef DEBUG
    printf("Transition matrix (diagonal):\n");
    printf("i\tj\ttransition(i => j)\trow sum\n");
    for (i = 0; i < length; ++i) {
      rowSum = 0;
      for (j = 0; j < length; ++j) {
        rowSum += transitionProbabilities[i][j];
      }
    
      printf("%d, %d\t%f\t%f\n", i, i, transitionProbabilities[i][i], rowSum);
    }
  
    printf("Matrix *to be* inverted (diagonal):\n");
    printf("i\tj\tinversion matrix(i => j)\trow sum\n");
    for (i = 0; i < inversionMatrixRowLength; ++i) {
      rowSum = 0;
      for (j = 0; j < inversionMatrixRowLength; ++j) {
        rowSum += inversionMatrix[i * inversionMatrixRowLength + j];
      }
    
      printf("%d, %d\t%f\t%f\n", i, i, inversionMatrix[i * inversionMatrixRowLength + i], rowSum * inversionMatrixRowLength);
    }
  #endif
  
  inversionMatrix = (*invert)(inversionMatrix, inversionMatrixRowLength);
  
  for (i = 0; i < inversionMatrixRowLength; ++i) {
    for (j = 0; j < inversionMatrixRowLength; ++j) {
      mfpt[i] += inversionMatrix[i * inversionMatrixRowLength + j];
    }
    
    #ifdef DEBUG
      // The business with this i < endIndex stuff is inorder to ensure that the output MFPT debug indices are representative of the input data.
      printf("MFPT %-8d%f\n", i < endIndex ? i : i + 1, mfpt[i]);
    #endif
  }
    
  mfptFromStart = mfpt[startIndex];
  free(mfpt);
  free(inversionMatrix);
  
  return mfptFromStart;
}

double* inverse(double* a, int size) {
  // http://stackoverflow.com/questions/3519959/computing-the-inverse-of-a-matrix-using-lapack-in-c
  int* ipiv    = (int*)malloc((size + 1) * sizeof(int));
  int lwork    = size * size;
  double* work = (double*)malloc(lwork * sizeof(double));
  int info;

  dgetrf_(&size, &size, a, &size, ipiv, &info);
  
  #ifdef DEBUG
    printf("info (dgetrf):\t%d\n", info);
  #endif
  
  dgetri_(&size, a, &size, ipiv, work, &lwork, &info);
  
  #ifdef DEBUG
    printf("info (dgetri):\t%d\n", info);
  #endif

  free(ipiv);
  free(work);
  
  return a;
}

double* pseudoinverse(double* a, int size) {
  // Least-squares fit solution to B - Ax, where (in this case) A is square and B is the identity matrix.
  int i, j, m, n, nrhs, lda, ldb, rank, nlvl, lwork, liwork, info;
  double rcond;

  m     = size;
  n     = m;
  nrhs  = m;
  lda   = m;
  ldb   = m;
  rcond = -1.;
  
  // NLVL = MAX(0, INT(LOG_2(MIN(M, N) / (SMLSIZ + 1))) + 1)
  nlvl = MAX(0, (int)(log2(MIN(m, n) / (SMLSIZ + 1)) + 1));
  
  // LWORK [when M >= N] = MAX(1, 12 * M + 2 * M * SMLSIZ + 8 * M * NLVL + M * NRHS + (SMLSIZ + 1) ** 2)
  lwork = MAX(1, 12 * m + 2 * m * SMLSIZ + 8 * m * nlvl + m * nrhs + (int)pow((double)(SMLSIZ + 1), 2.));
  
  // LIWORK = 3 * MIN(M, N) * NLVL + 11 * MIN(M, N)
  liwork = 3 * MIN(m, n) * nlvl + 11 * MIN(m, n);
  
  double* b = (double*)calloc(ldb * nrhs, sizeof(double));
  for (i = 0; i < ldb; ++i) {
    b[i * nrhs + i] = 1.;
  }
  double* s    = (double*)malloc(MIN(m, n) * sizeof(double));
  double* work = (double*)malloc(MAX(1, lwork) * sizeof(double));
  int* iwork   = (int*)malloc(MAX(1, liwork) * sizeof(int));
  
  #ifdef DEBUG
    printf("nlvl:\t%d\n", nlvl);
    printf("lwork:\t%d\n", lwork);
    printf("liwork:\t%d\n", liwork);
  #endif
  
  dgelsd_(&n, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
  
  #ifdef DEBUG
    printf("info:\t%d\n", info);
  #endif
  
  #ifdef HEAVY_DEBUG
    printf("least-squares fit:\t");
    for (i = 0; i < ldb; ++i) {
      printf("%+.8f\t", b[i * nrhs + i]);
    }
    printf("\n");
  
    printf("s:\t");
    for (i = 0; i < n; ++i) {
      printf("%+.8f\t", s[i]);
    }
    printf("\n");
  
    printf("rank:\t%d\n", rank);
  #endif
  
  free(s);
  free(work);
  free(iwork);
  
  return(b);
}

int numSingleBpMoves(int x, int y, int* k, int* l, int bpDist, int length) {
  int j, a, b, numMoves = 0;
  
  for (j = 0; j < length; ++j) {
    a = k[j];
    b = l[j];
    
    if (
      // Because N(x, y) is restricted to non-zero probabilities, we *assume* the input data satisfies the triangle inequality and bounds.
      // a + b >= bpDist &&
      // a + bpDist >= b &&
      // b + bpDist >= a &&
      (int)abs(x - a) == 1 && (int)abs(y - b) == 1
    ) {
      numMoves++;
    }
  }
  
  return numMoves;
}

double transitionRateFromProbabilities(double from, double to, double numFrom) {
  return MIN(1., to / from) / numFrom;
}

double transitionRateFromEnergies(double from, double to, double numFrom) {
  return MIN(1., exp(-(to - from) / RT)) / numFrom;
}

double transitionRateFromProbabilitiesWithHastings(double from, double to, double numFrom, double numTo) {
  return MIN(1., (numFrom / numTo) * (to / from)) / numFrom;
}

double transitionRateFromEnergiesWithHastings(double from, double to, double numFrom, double numTo) {
  return MIN(1., (numFrom / numTo) * exp(-(to - from) / RT)) / numFrom;
}