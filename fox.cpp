#include "fox.h"

double Fox::start, Fox::finish;
int Fox::num_Proc, Fox::sq_Proc, Fox::global_Rank;
int Fox::nodes, Fox::q, Fox::org_Nodes;
int **Fox::global_Matrix;
int **Fox::local_Matrix;
int **Fox::local_Matrix_A;
int **Fox::local_Matrix_B;
int **Fox::local_Matrix_C;

MPI_Comm Fox::whole_Grid, Fox::collumn_Grid, Fox::row_Grid;
int Fox::collumn, Fox::row;
int Fox::grid_Rank, Fox::collumn_Rank, Fox::row_Rank;

void Fox::read_Input()
{
  if (global_Rank != MASTER)
    return;

  int i, j;
  scanf("%d", &nodes);

  int tmp = (int)sqrt((double) num_Proc);
  org_Nodes = nodes;

  if (nodes % tmp != 0)
  {
    nodes = (nodes + tmp - 1) / tmp * tmp;
    fprintf(stderr, "Warning: matrix size not divisible by number of nodes\n");
    fprintf(stderr, "         filled extra space to fullfil contrain (%d)\n", nodes);
  }

  int* tmp_Matrix = (int*) malloc(nodes * nodes * sizeof(int));
  global_Matrix = (int**) malloc(nodes * sizeof(int*));
  for (i = 0; i < nodes; i++)
    global_Matrix[i] = &(tmp_Matrix[i * nodes]);

  for (i = 0; i < org_Nodes; i++)
    for (j = 0; j < org_Nodes; j++)
      scanf("%d", &global_Matrix[i][j]);

  for (i = 0; i < org_Nodes; i++)
    for (j = 0; j < org_Nodes; j++)
      if (i != j && global_Matrix[i][j] == 0)
        global_Matrix[i][j] = -1;

  for (i = 0; i < org_Nodes; i++)
    global_Matrix[i][i] = 0;

  for (i = 0; i < nodes; i++)
    for (j = 0; j < nodes; j++)
      if (i >= org_Nodes || j >= org_Nodes)
        global_Matrix[i][j] = -1;
}

void Fox::setup_Grid()
{
  MPI_Bcast(&nodes, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  sq_Proc = (int)sqrt((double) num_Proc);
  q = nodes / sq_Proc;

  if (sq_Proc * sq_Proc != num_Proc)
  {
    if (global_Rank == MASTER)
      fprintf(stderr, "Non perfect square number of processes chosen\n");
    exit(1);
  }

  if (nodes % sq_Proc != 0)
  {
    if (global_Rank == MASTER)
      fprintf(stderr, "Matrix size not divisible by number of processes\n");
    exit(1);
  }

  int i;
  int* tmp_Matrix = (int*) malloc(q * q * sizeof(int));
  local_Matrix = (int**) malloc(q * sizeof(int*));
  for (i = 0; i < q; i++)
    local_Matrix[i] = &(tmp_Matrix[i * q]);

  tmp_Matrix = (int*) malloc(q * q * sizeof(int));
  local_Matrix_A = (int**) malloc(q * sizeof(int*));
  for (i = 0; i < q; i++)
    local_Matrix_A[i] = &(tmp_Matrix[i * q]);

  tmp_Matrix = (int*) malloc(q * q * sizeof(int));
  local_Matrix_B = (int**) malloc(q * sizeof(int*));
  for (i = 0; i < q; i++)
    local_Matrix_B[i] = &(tmp_Matrix[i * q]);

  tmp_Matrix = (int*) malloc(q * q * sizeof(int));
  local_Matrix_C = (int**) malloc(q * sizeof(int*));
  for (i = 0; i < q; i++)
    local_Matrix_C[i] = &(tmp_Matrix[i * q]);

  int dimensions[2] = {sq_Proc, sq_Proc};
  int period[2] = {1, 1};
  MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, period, 1, &whole_Grid);
  MPI_Comm_rank(whole_Grid, &grid_Rank);

  int tcoordinates[2];
  MPI_Cart_coords(whole_Grid, grid_Rank, 2, tcoordinates);
  row = tcoordinates[0];
  collumn = tcoordinates[1];

  int free_Dimensions[2] = {0, 1};
  MPI_Cart_sub(whole_Grid, free_Dimensions, &row_Grid);

  free_Dimensions[0] = 1, free_Dimensions[1] = 0;
  MPI_Cart_sub(whole_Grid, free_Dimensions, &collumn_Grid);
}

void Fox::distribute_Input()
{
  MPI_Datatype matrix_Type;
  MPI_Type_vector(q, q, nodes, MPI_INT, &matrix_Type);
  MPI_Type_commit(&matrix_Type);

  if (global_Rank == MASTER)
  {
    int i, j;
    int coords[2];
    for (i = 0; i < nodes; i += q)
      for (j = 0; j < nodes; j += q)
      {
        int dest_Rank;
        coords[0] = i / q, coords[1] = j / q;
        MPI_Cart_rank(whole_Grid, coords, &dest_Rank);

        if (dest_Rank == MASTER)
        {
          int x, y;
          for (y = 0; y < q; y++)
            for (x = 0; x < q; x++)
              local_Matrix[y][x] = global_Matrix[y][x];
          continue;
        }

        MPI_Send(&(global_Matrix[i][j]), 1, matrix_Type, dest_Rank, SEND_TAG, whole_Grid);
      }
  }
  else
    MPI_Recv(local_Matrix[0], q * q, MPI_INT, MASTER, SEND_TAG, whole_Grid, MPI_STATUS_IGNORE);

  MPI_Type_free(&matrix_Type);
}

void Fox::multiply_Matrices(MPI_Datatype matrix_Type)
{
  int i, j;
  for (i = 0; i < q; i++)
    for (j = 0; j < q; j++)
    {
      local_Matrix_B[i][j] = local_Matrix[i][j];
      local_Matrix_C[i][j] = local_Matrix[i][j];
    }

  int top_Rank, bot_Rank;
  int tmp[1] = {(row + 1) % sq_Proc};
  MPI_Cart_rank(collumn_Grid, tmp, &bot_Rank);
  tmp[0] = (row - 1 + sq_Proc) % sq_Proc;
  MPI_Cart_rank(collumn_Grid, tmp, &top_Rank);

  int run;
  for (run = 0; run < sq_Proc; run++)
  {
    int u = (run + row) % sq_Proc;
    int u_Rank;
    tmp[0] = u;
    MPI_Cart_rank(row_Grid, tmp, &u_Rank);

    if (u == collumn)
    {
      MPI_Bcast(local_Matrix[0], 1, matrix_Type, u_Rank, row_Grid);
      Matrix::local_Multiply(local_Matrix, local_Matrix_B, local_Matrix_C, q);
    }
    else
    {
      MPI_Bcast(local_Matrix_A[0], 1, matrix_Type, u_Rank, row_Grid);
      Matrix::local_Multiply(local_Matrix_A, local_Matrix_B, local_Matrix_C, q);
    }

    if (run != sq_Proc - 1)
      MPI_Sendrecv_replace(local_Matrix_B[0], 1, matrix_Type, top_Rank, SEND_TAG, bot_Rank, SEND_TAG, collumn_Grid, MPI_STATUS_IGNORE);
  }

  for (i = 0; i < q; i++)
    for (j = 0; j < q; j++)
      local_Matrix[i][j] = local_Matrix_C[i][j];
}

void Fox::calculate_APSP()
{
  MPI_Datatype matrix_Type;
  MPI_Type_vector(q * q, 1, 1, MPI_INT, &matrix_Type);
  MPI_Type_commit(&matrix_Type);

  int run;
  for (run = 1; run <= nodes; run <<= 1)
    multiply_Matrices(matrix_Type);

  MPI_Type_free(&matrix_Type);
}

void Fox::collect_Results()
{
  MPI_Datatype matrix_Type;
  MPI_Type_vector(q, q, nodes, MPI_INT, &matrix_Type);
  MPI_Type_commit(&matrix_Type);

  if (global_Rank == MASTER)
  {
    int i, j;
    int coords[2];
    for (i = 0; i < nodes; i += q)
      for (j = 0; j < nodes; j += q)
      {
        int dest_Rank;
        coords[0] = i / q, coords[1] = j / q;
        MPI_Cart_rank(whole_Grid, coords, &dest_Rank);

        if (dest_Rank == MASTER)
        {
          int x, y;
          for (y = 0; y < q; y++)
            for (x = 0; x < q; x++)
              global_Matrix[y][x] = local_Matrix[y][x];
          continue;
        }

        MPI_Recv(&(global_Matrix[i][j]), 1, matrix_Type, dest_Rank, SEND_TAG, whole_Grid, MPI_STATUS_IGNORE);
      }
  }
  else
    MPI_Send(local_Matrix[0], q * q, MPI_INT, MASTER, SEND_TAG, whole_Grid);

  MPI_Type_free(&matrix_Type); 
}

void Fox::sequential_calculate_APSP()
{
  int* tmp_Matrix = (int*) malloc(nodes * nodes * sizeof(int));
  int** cp = (int**) malloc(nodes * sizeof(int*));
  int i, j, run;
  for (i = 0; i < nodes; i++)
    cp[i] = &(tmp_Matrix[i * nodes]);

  for (i = 0; i < nodes; i++)
    for (j = 0; j < nodes; j++)
      cp[i][j] = global_Matrix[i][j];

  for (run = 1; run < 2; run++)
  {
    Matrix::local_Multiply(global_Matrix, global_Matrix, cp, nodes);
              
    for (i = 0; i < nodes; i++)
      for (j = 0; j < nodes; j++)
        global_Matrix[i][j] = cp[i][j];

  }
}

void Fox::print_Result()
{
  if (global_Rank != MASTER)
    return;

  nodes = org_Nodes;

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

void Fox::APSP(int print_result, int print_time)
{
  MPI_Comm_size(MPI_COMM_WORLD, &num_Proc);
  MPI_Comm_rank(MPI_COMM_WORLD, &global_Rank);

  MPI_Barrier(MPI_COMM_WORLD);
  start = MPI_Wtime();

  read_Input();
  setup_Grid();
  distribute_Input();
  calculate_APSP();
  collect_Results();

  MPI_Barrier(MPI_COMM_WORLD);
  finish = MPI_Wtime();

  if (print_result)
    print_Result();

  if (print_time && global_Rank == 0)
    printf("Execution time: %0.3lf seconds\n", finish - start);
}
