#ifndef _FOX_
#define _FOX_

#include "common.h"
#include "matrix.h"

class Fox
{
 private:
  static double start, finish;
  static int num_Proc, sq_Proc, global_Rank;
  static int nodes, q, org_Nodes;
  static int **global_Matrix;
  static int **local_Matrix;
  static int **local_Matrix_A;
  static int **local_Matrix_B;
  static int **local_Matrix_C;

  static MPI_Comm whole_Grid, collumn_Grid, row_Grid;
  static int collumn, row;
  static int grid_Rank, collumn_Rank, row_Rank;

  static void read_Input();
  static void setup_Grid();
  static void distribute_Input();
  static void multiply_Matrices(MPI_Datatype matrix_Type);
  static void calculate_APSP();
  static void collect_Results();
  static void sequential_calculate_APSP();
  static void print_Result();

 public:
  static void APSP(int print_result, int print_time);
};

#endif _FOX_
