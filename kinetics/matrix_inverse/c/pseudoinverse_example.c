// M       (input) INTEGER   
//         The number of rows of A. M >= 0.   
// 
// N       (input) INTEGER   
//         The number of columns of A. N >= 0.   
// 
// NRHS    (input) INTEGER   
//         The number of right hand sides, i.e., the number of columns   
//         of the matrices B and X. NRHS >= 0.   
// 
// A       (input) DOUBLE PRECISION array, dimension (LDA,N)   
//         On entry, the M-by-N matrix A.   
//         On exit, A has been destroyed.   
// 
// LDA     (input) INTEGER   
//         The leading dimension of the array A.  LDA >= max(1,M).   
// 
// B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)   
//         On entry, the M-by-NRHS right hand side matrix B.   
//         On exit, B is overwritten by the N-by-NRHS solution   
//         matrix X.  If m >= n and RANK = n, the residual   
//         sum-of-squares for the solution in the i-th column is given   
//         by the sum of squares of elements n+1:m in that column.   
// 
// LDB     (input) INTEGER   
//         The leading dimension of the array B. LDB >= max(1,max(M,N)).   
// 
// S       (output) DOUBLE PRECISION array, dimension (min(M,N))   
//         The singular values of A in decreasing order.   
//         The condition number of A in the 2-norm = S(1)/S(min(m,n)).   
// 
// RCOND   (input) DOUBLE PRECISION   
//         RCOND is used to determine the effective rank of A.   
//         Singular values S(i) <= RCOND*S(1) are treated as zero.   
//         If RCOND < 0, machine precision is used instead.   
// 
// RANK    (output) INTEGER   
//         The effective rank of A, i.e., the number of singular values   
//         which are greater than RCOND*S(1).   
// 
// WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))   
//         On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   
// 
// LWORK   (input) INTEGER   
//         The dimension of the array WORK. LWORK must be at least 1.   
//         The exact minimum amount of workspace needed depends on M,   
//         N and NRHS. As long as LWORK is at least   
//             12*N + 2*N*SMLSIZ + 8*N*NLVL + N*NRHS + (SMLSIZ+1)**2,   
//         if M is greater than or equal to N or   
//             12*M + 2*M*SMLSIZ + 8*M*NLVL + M*NRHS + (SMLSIZ+1)**2,   
//         if M is less than N, the code will execute correctly.   
//         SMLSIZ is returned by ILAENV and is equal to the maximum   
//         size of the subproblems at the bottom of the computation   
//         tree (usually about 25), and   
//            NLVL = MAX( 0, INT( LOG_2( MIN( M,N )/(SMLSIZ+1) ) ) + 1 )   
//         For good performance, LWORK should generally be larger.   
// 
//         If LWORK = -1, then a workspace query is assumed; the routine   
//         only calculates the optimal size of the WORK array, returns   
//         this value as the first entry of the WORK array, and no error   
//         message related to LWORK is issued by XERBLA.   
// 
// IWORK   (workspace) INTEGER array, dimension (MAX(1,LIWORK))   
//         IWORK dimension should be at least 3*MIN(M,N)*NLVL + 11*MIN(M,N),
//         where NLVL = MAX( 0, INT( LOG_2( MIN(M,N)/(SMLSIZ+1) ) )+1 )
//         and SMLSIZ = 25  
// 
// INFO    (output) INTEGER   
//         = 0:  successful exit   
//         < 0:  if INFO = -i, the i-th argument had an illegal value.   
//         > 0:  the algorithm for computing the SVD failed to converge;   
//               if INFO = i, i off-diagonal elements of an intermediate   
//               bidiagonal form did not converge to zero.

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define SMLSIZ 25

extern int dgelsd_(int *m, int *n, int *nrhs, double *a, int *lda, double *b, int *ldb, double *s, double *rcond, int *rank, double *work, int *lwork, int *iwork, int *info);

int main (int argc, const char * argv[]) {
  int i, j, m, n, nrhs, lda, ldb, rank, nlvl, lwork, liwork, info;
  double rcond;

  m     = 5000;
  n     = m;
  nrhs  = m;
  lda   = m;
  ldb   = m;
  rcond = -1.;
  
  // NLVL = MAX(0, INT(LOG_2(MIN(M, N) / (SMLSIZ + 1))) + 1)
  nlvl = MAX(0, (int)(log2(MIN(m, n) / (SMLSIZ + 1)) + 1));
  
  // LWORK [when M >= N] = MAX(1, 12 * M + 2 * M * SMLSIZ + 8 * M * NLVL + M * NRHS + (SMLSIZ + 1) ** 2)
  lwork = MAX(1, 12 * m + 2 * m * SMLSIZ + 8 * m * nlvl + m * nrhs + (int)pow((double)(SMLSIZ + 1), 2.));
  
  // LIWORK = 3*MIN(M,N)*NLVL + 11*MIN(M,N)
  liwork = 3 * MIN(m, n) * nlvl + 11 * MIN(m, n);
  
  double* a = (double*)malloc(lda * n * sizeof(double));
  for (i = 0; i < lda * n; ++i) {
    a[i] = (double)rand();
  }
  
  double* b = (double*)calloc(ldb * nrhs, sizeof(double));
  for (i = 0; i < ldb; ++i) {
    b[i * nrhs + i] = 1.;
  }
  double* s    = (double*)malloc(MIN(m, n) * sizeof(double));
  double* work = (double*)malloc(MAX(1, lwork) * sizeof(double));
  int* iwork   = (int*)malloc(MAX(1, liwork) * sizeof(int));
  
  printf("nlvl:\t%d\n", nlvl);
  printf("lwork:\t%d\n", lwork);
  printf("liwork:\t%d\n", liwork);
  
  // for (i = 0; i < n; ++i) {
  //   for (j = 0; j < nrhs; ++j) {
  //     printf("%+.8f\t", b[i * nrhs + j]);
  //   }
  //   printf("\n");
  // }
  
  dgelsd_(&n, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork, iwork, &info);
  
  printf("info:\t%d\n", info);
  
  printf("least-squares fit:\t");
  for (i = 0; i < ldb; ++i) {
    printf("%+.8f\t", b[i * nrhs + i]);
  }
  printf("\n");
  
  // for (i = 0; i < n; ++i) {
  //   for (j = 0; j < nrhs; ++j) {
  //     printf("%+.8f\t", b[i * nrhs + j]);
  //   }
  //   printf("\n");
  // }
  // 
  // printf("s:\t");
  // for (i = 0; i < n; ++i) {
  //   printf("%+.8f\t", s[i]);
  // }
  // printf("\n");
  
  printf("rank:\t%d\n", rank);
  
  free(s);
  free(work);
  free(iwork);
  
  return(info);
}