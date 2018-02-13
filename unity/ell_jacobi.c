#include "ell.h"
#include <stdio.h>

#define NRM  "\x1B[0m"
#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"

#define N 100

int main(void)
{
  printf("ell.c/ell.h\n");

  FILE *fl;
  double *x = malloc(N*sizeof(double));
  double *b = malloc(N*sizeof(double));
  double m_e[4] = { 1, -1, -1, 1 };

  ell_matrix m;
  ell_solver solver;
  ell_init(&m, N, N, 4);
  int ix[2], iy[2];
  for (int i = 0 ; i < N-1 ; i++) {
    b[i] = 1.0;
    ix[0] = i; ix[1] = i+1;
    iy[0] = i; iy[1] = i+1;
    ell_add_vals(&m, ix, 2, iy, 2, m_e);
  }
  ell_set_zero_row(&m, 0  , 1);
  ell_set_zero_row(&m, N-1, 1);
  ell_set_zero_col(&m, 0  , 1);
  ell_set_zero_col(&m, N-1, 1);

//  printf("\nm=\n");
//  ell_print_full(&m);

  for (int i = 0 ; i < N ; i++)
    x[i] = 0.0;
  solver.max_its = 20;
  solver.min_tol = 1.0e-5;
  ell_solve_cg(&solver, &m, b, x);

  fl = fopen("sol_cg.dat","w");
  for (int i = 0 ; i < N ; i++)
    fprintf(fl,"%lf\n", x[i]);
  fclose(fl);
  printf(GRN"cg\nerr = %lf\nits = %d\n"NRM, solver.err, solver.its);

  for (int i = 0 ; i < N ; i++)
    x[i] = 0.0;
  solver.max_its = 100000;
  solver.min_tol = 1.0e-5;
  ell_solve_jacobi(&solver, &m, b, x);

  fl = fopen("sol_jacobi.dat","w");
  for (int i = 0 ; i < N ; i++)
    fprintf(fl,"%lf\n", x[i]);
  fclose(fl);
  printf(GRN"jacobi\nerr = %lf\nits = %d\n"NRM, solver.err, solver.its);

  return 0;
}
