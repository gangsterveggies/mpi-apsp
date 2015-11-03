#ifndef _FLOYD_
#define _FLOYD_

#include "common.h"
#include "matrix.h"

class Floyd
{
 private:
  static double start, finish;
  static int num_Proc, global_Rank;
  static int nodes, q;
  static int **global_Matrix;
  static int **local_Matrix;
  static int *local_Matrix_K;

  static void read_Input();
  static void setup_Grid();
  static void distribute_Input();
  static void multiply_Matrices(MPI_Datatype matrix_Type);
  static void relax(int &a, int &b, int &c);
  static void calculate_APSP();
  static void collect_Results();
  static void sequential_calculate_APSP();
  static void print_Result();

 public:
  static void APSP(int print_result, int print_time);
};

#endif _FLOYD_
