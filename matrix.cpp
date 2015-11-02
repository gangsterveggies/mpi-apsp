#include "matrix.h"

int Matrix::min(int a, int b)
{
  return (a < b) ? a : b;
}

int Matrix::valid(int vl, int a, int b)
{
  return vl || (a == b);
}

void Matrix::local_Multiply(int** mat_A, int** mat_B, int** mat_C, int sz)
{
  int i, j, k;
  for (i = 0; i < sz; i++)
    for (j = 0; j < sz; j++)
      for (k = 0; k < sz; k++)
        if (mat_A[i][k] != -1 && mat_B[k][j] != -1 && mat_C[i][j] != -1)
          mat_C[i][j] = min(mat_C[i][j], mat_A[i][k] + mat_B[k][j]);
        else if (mat_A[i][k] != -1 && mat_B[k][j] != -1)
          mat_C[i][j] = mat_A[i][k] + mat_B[k][j];
}
