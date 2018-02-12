#ifndef ELL_H_
#define ELL_H_

#include <stdlib.h> 
#include <stdio.h>
#include <math.h>

#define NRM  "\x1B[0m"
#define RED  "\x1B[31m"
#define GRN  "\x1B[32m"

#define MAX_ITS 100000;
#define MIN_ERR 1.0e-10;

typedef struct ell_matrix_ {
  int nrow;
  int ncol;
  int nnz;
  int *cols;
  double *vals;
} ell_matrix;

typedef struct ell_solver_ {
  int max_its;
  int its;
  double min_err;
  double err;
} ell_solver;

int ell_init (ell_matrix *m, int nrow, int ncol, int nnz);
int ell_set_val (ell_matrix *m, int row, int col, double val);
int ell_add_val (ell_matrix *m, int row, int col, double val);
int ell_add_vals (ell_matrix *m, int *ix, int nx, int *iy, int ny, double *vals);
int ell_mvp (ell_matrix *m, double *x, double *y);
int ell_get_val (ell_matrix *m, int row, int col, double *val);
int ell_solve_jacobi (ell_solver *solver, ell_matrix * m, double *b, double *x);
int ell_set_zero_row (ell_matrix *m, int row, double diag_val);
int ell_set_zero_col (ell_matrix *m, int col, double diag_val);
int ell_set_zero_mat (ell_matrix * m);
int ell_print_full (ell_matrix * m);

#endif
