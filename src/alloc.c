#include "micro.h"

int alloc_memory(void) {

  int ierr = 0;
  int nnz = (mesh_struct.dim == 2)? 18:81;
  int nn = mesh_struct.nn;
  int nx = mesh_struct.nx;
  int ny = mesh_struct.ny;
  int dim = mesh_struct.dim;

  int nrows = 0;
  switch (params.fe2_bc) {
    case BC_USTRAIN:
      nrows = nn * dim;
      break;
    case BC_PER_LM:
      nrows = (nn + (nx - 2) + (ny - 2)) * dim;
      break;
    case BC_PER_MS:
      nrows = (nn - nx - ny + 2) * dim;
      break;
    case BC_USTRESS:
      nrows = nn * dim + nvoi;
      break;
    default:
      return 1;
  }

  x_ell = malloc(nrows*sizeof(double));
  dx_ell = malloc(nrows*sizeof(double));
  res_ell = malloc(nrows*sizeof(double));
  ell_init(&jac_ell, nrows, nrows, nnz);

  elem_index = malloc(dim*mesh_struct.npe*sizeof(int));
  elem_nods = malloc(mesh_struct.npe*sizeof(int));
  elem_disp = malloc(dim*mesh_struct.npe*sizeof(double));
  stress_gp = malloc(nvoi*sizeof(double));
  strain_gp = malloc(nvoi*sizeof(double));
  elem_strain = malloc(mesh_struct.nelm*nvoi*sizeof(double));
  elem_stress = malloc(mesh_struct.nelm*nvoi*sizeof(double));
  elem_energy = malloc(mesh_struct.nelm*sizeof(double));
  elem_type = malloc(mesh_struct.nelm*sizeof(int));

  struct_bmat = malloc(nvoi*sizeof(double**));
  for (int i = 0 ; i < nvoi ; i++) {
    struct_bmat[i] = malloc(mesh_struct.npe*dim*sizeof(double*));
    for (int j = 0 ; j < mesh_struct.npe*dim ; j++)
      struct_bmat[i][j] = malloc(ngp*sizeof(double));
  }

  flags.allocated = true;

  return ierr;
}
