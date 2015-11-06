#include "common.h"
#include "matrix.h"
#include "fox.h"
#include "floyd.h"

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  int i;
  int to_Floyd = 0, to_Print = 0, to_Time = 0;  

  for (i = 1; i < argc; i++)
    if (strcmp(argv[i], "-fw") == 0)
      to_Floyd = 1;
    else if (strcmp(argv[i], "-tm") == 0)
      to_Time = 1;
    else if (strcmp(argv[i], "-pr") == 0)
      to_Print = 1;

  if (to_Floyd)
    Floyd::APSP(to_Print, to_Time);
  else
    Fox::APSP(to_Print, to_Time);

  MPI_Finalize();

  return 0;
}
