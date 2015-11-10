#ifndef _SEQFOX_
#define _SEQFOX_

#include "common.h"
#include "matrix.h"

class SeqFox
{
 private:
  static double start, finish;
  static int nodes;
  static int **global_Matrix;

  static void read_Input();
  static void setup_Grid();
  static void calculate_APSP();
  static void print_Result();

 public:
  static void APSP(int print_result, int print_time);
};

#endif _FOX_
