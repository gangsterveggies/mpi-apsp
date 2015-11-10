#include "seqfox.h"

double SeqFox::start, SeqFox::finish;
int SeqFox::nodes;
int **SeqFox::global_Matrix;

void SeqFox::read_Input()
{
  int i, j;
  scanf("%d", &nodes);

  int* tmp_Matrix = (int*) malloc(nodes * nodes * sizeof(int));
  global_Matrix = (int**) malloc(nodes * sizeof(int*));
  for (i = 0; i < nodes; i++)
    global_Matrix[i] = &(tmp_Matrix[i * nodes]);

  for (i = 0; i < nodes; i++)
    for (j = 0; j < nodes; j++)
      scanf("%d", &global_Matrix[i][j]);

  for (i = 0; i < nodes; i++)
    for (j = 0; j < nodes; j++)
      if (i != j && global_Matrix[i][j] == 0)
        global_Matrix[i][j] = -1;

  for (i = 0; i < nodes; i++)
    global_Matrix[i][i] = 0;
}

void SeqFox::calculate_APSP()
{
  int run;
  for (run = 1; run <= nodes; run <<= 1)
    Matrix::local_Multiply(global_Matrix, global_Matrix, global_Matrix, nodes);
}

void SeqFox::print_Result()
{
  int i, j;
  for (i = 0; i < nodes; i++)
    for (j = 0; j < nodes; j++)
      if (global_Matrix[i][j] == -1)
        global_Matrix[i][j] = 0;

  for (i = 0; i < nodes; i++)
  {
    for (j = 0; j < nodes; j++)
      printf("%d%c", global_Matrix[i][j], j == nodes - 1 ? '\n' : ' ');
  }
}

void SeqFox::APSP(int print_result, int print_time)
{
  int global_Rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &global_Rank);

  if (global_Rank == 0)
    printf("Running using sequential Fox\n");

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  if (global_Rank == 0)
  {
    read_Input();
    calculate_APSP();
  }

  MPI_Barrier(MPI_COMM_WORLD);
  finish = MPI_Wtime();

  if (print_result && global_Rank == 0)
    print_Result();

  if (print_time && global_Rank == 0)
    printf("Execution time: %0.3lf seconds\n", finish - start);
}
