#ifndef _SEQFLOYD_
#define _SEQFLOYD_

#include "common.h"
#include "matrix.h"

class SeqFloyd
{
 private:
  static double start, finish;
  static int nodes;
  static int **global_Matrix;

  static void read_Input();
  static void relax(int &a, int &b, int &c);
  static void calculate_APSP();
  static void print_Result();

 public:
  static void APSP(int print_result, int print_time);
};

#endif _SEQFLOYD_
