#include "common.h"
#include "matrix.h"
#include "fox.h"
#include "floyd.h"

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  if (argc > 1 && strcmp(argv[1], "-fw") == 0)
    Floyd::APSP(PRINT_MODE, TIME_MODE);
  else
    Fox::APSP(PRINT_MODE, TIME_MODE);

  MPI_Finalize();

  return 0;
}
