
///@author Milan Zaloudek

#include "smallCodeFunction.h"


double * allocateVector(int n) {
  double *A;
  A = (double *) malloc(n*sizeof(double));
  if (A == NULL) {printf("Allocation Error\n"); exit(1);};
  return A;
};


double ** allocateMatrix(int n, int m) {
  double **A;
  A = (double **) malloc(n*sizeof(double*));
  if (A == NULL) {printf("Allocation Error\n"); exit(1);};
  A[0] = (double *) malloc(n*m*sizeof(double));
  if (A[0] == NULL) {printf("Allocation Error\n"); exit(1);};
  for (int i = 1; i<n; i++)
    A[i] = A[0] + i*m;
  return A;
}


double *** allocateTensor(int n, int m, int k) {
  double ***A;
  A = (double ***) malloc(n*sizeof(double**));
  A[0] = (double **) malloc(n*m*sizeof(double*));
  for (int i=0; i<n; i++)
    A[i] = allocateMatrix(m, k);
  return A;
}


ONE_NODE *allocateCell(int nb) {
  ONE_NODE *A;
  A = (ONE_NODE *) malloc(nb*sizeof(ONE_NODE));
  if (A == NULL) {printf("Allocation Error\n"); exit(1);};
  return A;
}


ONE_TETRA_ELEMENT *allocateTetra(int nb) {
  ONE_TETRA_ELEMENT *A;
  A = (ONE_TETRA_ELEMENT *) malloc(nb*sizeof(ONE_TETRA_ELEMENT));
  if (A == NULL) {printf("Allocation Error\n"); exit(1);};
  return A;
}


ONE_TRS *allocateTRS(int nb) {
  ONE_TRS *A;
  A = (ONE_TRS *) malloc(nb*sizeof(ONE_TRS));
  if (A == NULL) {printf("Allocation Error\n"); exit(1);};
  return A;
}
