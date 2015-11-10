#include "seqfloyd.h"

double SeqFloyd::start, SeqFloyd::finish;
int SeqFloyd::nodes;
int **SeqFloyd::global_Matrix;

void SeqFloyd::read_Input()
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


void SeqFloyd::relax(int &a, int &b, int &c)
{
  if (a != -1 && b != -1 && c != -1)
    a = Matrix::min(a, b + c);
  else if (b != -1 && c != -1)
    a = b + c;
}

void SeqFloyd::calculate_APSP()
{
  int i, j, k;
  int dest_Rank = 0, curr_Done = 0;

  for (k = 0; k < nodes; k++)
    for (i = 0; i < nodes; i++)
      for (j = 0; j < nodes; j++)
        relax(global_Matrix[i][j], global_Matrix[i][k], global_Matrix[k][j]);
}

void SeqFloyd::print_Result()
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

void SeqFloyd::APSP(int print_result, int print_time)
{
  int global_Rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &global_Rank);

  if (global_Rank == 0)
    printf("Running using sequential Floyd-Warshall\n");

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

