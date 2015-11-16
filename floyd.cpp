#include "floyd.h"

double Floyd::start, Floyd::finish;
int Floyd::num_Proc, Floyd::global_Rank, Floyd::q;
int Floyd::nodes;
int **Floyd::global_Matrix;
int **Floyd::local_Matrix;
int *Floyd::local_Matrix_K;

void Floyd::read_Input()
{
  if (global_Rank != MASTER)
    return;

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

void Floyd::setup_Grid()
{
  MPI_Bcast(&nodes, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

  q = nodes / num_Proc + ((nodes - num_Proc * (nodes / num_Proc)) > global_Rank);

  int i;
  int* tmp_Matrix = (int*) malloc(nodes * q * sizeof(int));
  local_Matrix = (int**) malloc(q * sizeof(int*));
  for (i = 0; i < q; i++)
    local_Matrix[i] = &(tmp_Matrix[i * nodes]);

  local_Matrix_K = (int*) malloc(nodes * sizeof(int));
}

void Floyd::distribute_Input()
{
  if (global_Rank == MASTER)
  {
    int i, j, acc = q;
    for (i = 0; i < num_Proc; i++)
    {
      int dest_Rank = i;

      if (dest_Rank == MASTER)
      {
        int x, y;
        for (y = 0; y < q; y++)
          for (x = 0; x < nodes; x++)
            local_Matrix[y][x] = global_Matrix[y][x];
        continue;
      }

      int tq = nodes / num_Proc + ((nodes - num_Proc * (nodes / num_Proc)) > dest_Rank);
      MPI_Send(&(global_Matrix[acc][0]), nodes * tq, MPI_INT, dest_Rank, SEND_TAG, MPI_COMM_WORLD);
      acc += tq;
    }
  }
  else
    MPI_Recv(&(local_Matrix[0][0]), nodes * q, MPI_INT, MASTER, SEND_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void Floyd::relax(int &a, int &b, int &c)
{
  if (a != -1 && b != -1 && c != -1)
    a = Matrix::min(a, b + c);
  else if (b != -1 && c != -1)
    a = b + c;
}

void Floyd::calculate_APSP()
{
  int i, j, k;
  int dest_Rank = 0, curr_Done = 0;

  for (k = 0; k < nodes; k++, curr_Done++)
  {
    int tq = nodes / num_Proc + ((nodes - num_Proc * (nodes / num_Proc)) > dest_Rank);
    if (tq == curr_Done)
    {
      curr_Done = 0;
      dest_Rank++;
    }

    if (dest_Rank == global_Rank)
      MPI_Bcast(&(local_Matrix[k % q][0]), nodes, MPI_INT, dest_Rank, MPI_COMM_WORLD);
    else
      MPI_Bcast(&(local_Matrix_K[0]), nodes, MPI_INT, dest_Rank, MPI_COMM_WORLD);

    for (i = 0; i < q; i++)
      for (j = 0; j < nodes; j++)
        if (dest_Rank == global_Rank)
          relax(local_Matrix[i][j], local_Matrix[i][k], local_Matrix[k % q][j]);
        else
          relax(local_Matrix[i][j], local_Matrix[i][k], local_Matrix_K[j]);
  }
}

void Floyd::collect_Results()
{
  if (global_Rank == MASTER)
  {
    int i, j, acc = q;
    for (i = 0; i < num_Proc; i++)
    {
      int dest_Rank = i;

      if (dest_Rank == MASTER)
      {
        int x, y;
        for (y = 0; y < q; y++)
          for (x = 0; x < nodes; x++)
            global_Matrix[y][x] = local_Matrix[y][x];
        continue;
      }

      int tq = nodes / num_Proc + ((nodes - num_Proc * (nodes / num_Proc)) > dest_Rank);
      MPI_Recv(&(global_Matrix[acc][0]), nodes * tq, MPI_INT, dest_Rank, SEND_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      acc += tq;
    }
  }
  else
    MPI_Send(&(local_Matrix[0][0]), nodes * q, MPI_INT, MASTER, SEND_TAG, MPI_COMM_WORLD);
}

void Floyd::print_Result()
{
  if (global_Rank != MASTER)
    return;

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

void Floyd::APSP(int print_result, int print_time)
{
  MPI_Comm_size(MPI_COMM_WORLD, &num_Proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &global_Rank);

  if (global_Rank == MASTER)
    printf("Running using parallel Floyd with %d processes\n", num_Proc);

  read_Input();

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  setup_Grid();
  distribute_Input();
  calculate_APSP();
  collect_Results();

  MPI_Barrier(MPI_COMM_WORLD);
  finish = MPI_Wtime();

  if (print_result)
    print_Result();

  if (print_time && global_Rank == MASTER)
    printf("Execution time: %0.3lf seconds\n", finish - start);
}
