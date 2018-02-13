#include "ell.h"

int ell_init (ell_matrix * m, int nrow, int ncol, int nnz)
{
  if (m == NULL) return 1;
  m->nnz = nnz;
  m->nrow = nrow;
  m->ncol = ncol;
  m->cols = malloc((nrow*nnz) * sizeof(int));
  m->vals = malloc((nrow*nnz) * sizeof(double));
  if (m->vals == NULL || m->cols == NULL) return 2;
  for (int i = 0 ; i < (nrow*nnz) ; i++) m->cols[i] = -1;
  for (int i = 0 ; i < (nrow*nnz) ; i++) m->vals[i] = +0;
  return 0;
}

int ell_set_val (ell_matrix * m, int row, int col, double val)
{
  if (row >= m->nrow || col >= m->ncol) {
    printf(RED "ell error: row %d or col %d greater than the dimension of the matrix\n" NRM, row, col);
    return 1;
  }
  if (row < 0 || col < 0) {
    printf(RED "ell error: negative values in row %d or col %d\n" NRM, row, col);
    return 2;
  }
  int j = 0;
  while (j < m->nnz) {
    if (m->cols[(row*m->nnz) + j] == -1) {
      m->cols[(row*m->nnz) + j] = col;
      m->vals[(row*m->nnz) + j] = val;
      return 0;
    } else if (m->cols[(row*m->nnz) + j] == col) {
      m->vals[(row*m->nnz) + j] = val;
      return 0;
    }
    j++;
  }
  if (j == m->nnz) {
    printf(RED "ell error: not enought space to store value in row %d and col %d\n" NRM, row, col);
    return 3;
  }
  return 4;
}

int ell_add_val (ell_matrix * m, int row, int col, double val)
{
  if (row >= m->nrow || col >= m->ncol) {
    printf(RED "ell error: row %d or col %d greater than the dimension of the matrix\n" NRM, row, col);
    return 1;
  }
  if (row < 0 || col < 0) {
    printf(RED "ell error: negative values in row %d or col %d\n" NRM, row, col);
    return 2;
  }
  int j = 0;
  while (j < m->nnz) {
    if (m->cols[(row*m->nnz) + j] == -1) {
      m->cols[(row*m->nnz) + j] = col;
      m->vals[(row*m->nnz) + j] = val;
      return 0;
    } else if (m->cols[(row*m->nnz) + j] == col) {
      m->vals[(row*m->nnz) + j] += val;
      return 0;
    }
    j++;
  }
  if (j == m->nnz) {
    printf(RED "ell error: not enought space to add value in row %d and col %d\n" NRM, row, col);
    return 3;
  }
  return 4;
}

int ell_add_vals (ell_matrix *m, int *ix, int nx, int *iy, int ny, double *vals)
{
  if (m == NULL || ix == NULL || iy == NULL || vals == NULL) return 1;
  for (int i = 0 ; i < nx ; i++) {
    for (int j = 0 ; j < ny ; j++) {
      ell_add_val (m, ix[i], iy[j], vals[i*nx + j]);
    }
  }
  return 0;
}

int ell_set_zero_row (ell_matrix *m, int row, double diag_val)
{
  if (m == NULL) return 1;
  int j = 0;
  while (j < m->nnz) {
    if (m->cols[(row*m->nnz) + j] == -1) {
      return 0;
    } else if (m->cols[(row*m->nnz) + j] == row) {
      m->vals[(row*m->nnz) + j] = diag_val;
    } else {
      m->vals[(row*m->nnz) + j] = 0.0;
    }
    j++;
  }
  return 0;
}

int ell_set_zero_col (ell_matrix *m, int col, double diag_val)
{
  if (m == NULL) return 1;
  for (int i = 0 ; i < m->nrow ; i++) {
    int j = 0;
    while (j < m->nnz) {
      if (m->cols[(i*m->nnz) + j] == -1) {
	return 0;
      } else if (m->cols[(i*m->nnz) + j] == col) {
	m->vals[(i*m->nnz) + j] = (i == col) ? diag_val : 0.0;
      }
      j++;
    }
  }
  return 0;
}

int ell_set_zero_mat (ell_matrix * m)
{
  for (int i = 0 ; i < m->nrow ; i++) ell_set_zero_row (m, i, 0.0);
  return 0;
}


int ell_mvp (ell_matrix * m, double *x, double *y)
{
  //  y = m * x
  if (m == NULL || x == NULL || y == NULL) return 1;
  for (int i = 0 ; i < m->nrow ; i++) {
    y[i] = 0;
    int j = 0;
    while (j < m->nnz) {
      if (m->cols[(i*m->nnz) + j] == -1) break;
      y[i] += m->vals[(i*m->nnz) + j] * x[m->cols[(i*m->nnz) + j]];
      j++;
    }
  }
  return 0;
}

int ell_solve_jacobi (ell_solver *solver, ell_matrix * m, double *b, double *x)
{
  // A = K - N  
  // K = diag(A)
  // N_ij = -a_ij for i!=j  and =0 if i=j
  // x_(i) = K^-1 * ( N * x_(i-1) + b )

  if (m == NULL || b == NULL || x == NULL) return 1;
  
  double *k = malloc (m->nrow * sizeof(double)); // K = diag(A)
  double *e_i = malloc(m->nrow * sizeof(double));

  for (int i = 0 ; i < m->nrow ; i++) {
    ell_get_val (m, i, i, &k[i]);
    k[i] = 1 / k[i];
  }

  int its = 0;
  int max_its = solver->max_its;
  double err;
  double min_tol = solver->min_tol;

  while (its < max_its) {

    err = 0;
    int i = 0;
    while (i < m->nrow) {
      double aux = 0.0; // sum_(j!=i) a_ij * x_j
      int j = 0;
      while (j < m->nnz) {
        if (m->cols[i*m->nnz + j] == -1) break;
        if (m->cols[i*m->nnz + j] != i)
	  aux += m->vals[i*m->nnz + j] * x[m->cols[i*m->nnz + j]];
	j++;
      }
      x[i] = k[i] * (-1*aux + b[i]);
      i++;
    }

    err = 0;
    ell_mvp (m, x, e_i);
    for (int i = 0 ; i < m->nrow ; i++){
      e_i[i] -= b[i];
      err += e_i[i] * e_i[i];
    }
    err = sqrt(err); if (err < min_tol) break;
    its ++;
  }
  solver->err = err;
  solver->its = its;
  return 0;
}

int ell_get_val (ell_matrix * m, int row, int col, double *val)
{
  if (row >= m->nrow || col >= m->ncol) {
    printf(RED "ell_get_val: row %d or col %d greater than the dimension of the matrix\n" NRM, row, col);
    return 1;
  }
  if (row < 0 || col < 0) {
    printf(RED "ell_get_val: negative values in row %d or col %d\n" NRM, row, col);
    return 2;
  }
  int j = 0;
  while (j < m->nnz) {
    if (m->cols[(row * m->nnz) + j] == -1) {
      *val = 0.0;
      return 0;
    } else if (m->cols[(row * m->nnz) + j] == col) {
      *val = m->vals[(row * m->nnz) + j];
      return 0;
    }
    j++;
  }
  return 3;
}

int ell_print_full (ell_matrix * m)
{
  if (m == NULL) return 1;
  if (m->vals == NULL || m->cols == NULL) return 2;
  double val;
  for (int i = 0 ; i < m->nrow ; i++) {
    for (int j = 0 ; j < m->ncol ; j++) {
      printf("%lf%s",(ell_get_val(m, i, j, &val) == 0)?val:0.0, (j == m->ncol - 1) ? "\n" : " ");
    }
  }
  return 0;
}
