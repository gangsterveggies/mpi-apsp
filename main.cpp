#include "common.h"
#include "matrix.h"
#include "fox.h"
#include "floyd.h"
#include "seqfox.h"
#include "seqfloyd.h"

void print_Help()
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == MASTER)
  {
    printf("\t\t---- Parallel APSP ----\n");
    printf("\tFilipe Figueiredo, Pedro Paredes\n\n");
    printf("usage:\n");
    printf("\tmpirun -np <np> [--hostfile <hostfile>] fox [arguments]\n\n");
    printf("Available arguments:\n");
    printf("\t -h\t\tdisplay this help file\n"
           "\t -fw\t\tuse Floyd-Warshall algorithm (default is Fox)\n"
           "\t -sq\t\trun sequential version of selected algorithm (default is parallel)\n"
           "\t -pr\t\tdo not print output of computation (default is to print)\n"
           "\t -tm\t\tprint time spent in computation (default is no print)\n"
           "\t -i <infile>\tread input from file <infile>\n");
  }
}

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  int i;
  int to_Floyd = 0, to_Seq = 0, to_Print = 1, to_Time = 0;  

  for (i = 1; i < argc; i++)
    if (strcmp(argv[i], "-fw") == 0)
      to_Floyd = 1;
    else if (strcmp(argv[i], "-sq") == 0)
      to_Seq = 1;
    else if (strcmp(argv[i], "-tm") == 0)
      to_Time = 1;
    else if (strcmp(argv[i], "-pr") == 0)
      to_Print = 0;
    else if (strcmp(argv[i], "-i") == 0)
    {
      freopen(argv[i + 1], "r", stdin);
      i++;
    }
    else if (strcmp(argv[i], "-h") == 0)
    {
      print_Help();

      MPI_Finalize();
      exit(0);
    }

  if (to_Seq && !to_Floyd)
    SeqFox::APSP(to_Print, to_Time);
  else if (to_Seq && to_Floyd)
    SeqFloyd::APSP(to_Print, to_Time);
  else if (to_Floyd)
    Floyd::APSP(to_Print, to_Time);
  else
    Fox::APSP(to_Print, to_Time);

  MPI_Finalize();

  return 0;
}
