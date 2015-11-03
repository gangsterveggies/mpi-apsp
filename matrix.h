#ifndef _MATRIX_
#define _MATRIX_

#include "common.h"

class Matrix
{
 private:
  static int valid(int vl, int a, int b);

 public:
  static int min(int a, int b);
  static void local_Multiply(int** mat_A, int** mat_B, int** mat_C, int sz);
};

#endif _MATRIX_
