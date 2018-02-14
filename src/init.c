#include "micro.h"

int init_variables_1(void)
{
  params.multis_method = MULTIS_NULL;
  params.fe2_bc = BC_NULL;
  params.nl_max_its = 2;
  params.nl_min_norm = 1.0e-4;
  params.rho = 1.0e7;
  params.solver = SOL_PETSC;

  solver.type = SOLVER_PETSC;

  flags.coupled = false;
  flags.linear_materials = false;
  flags.allocated = false;
  flags.print_pvtu = false;
  flags.print_vectors = false;
  flags.print_matrices = false;

  return 0;
}

int init_variables_2(void)
{

  for (int gp = 0; gp < ngp ; gp++) {
    for (int is = 0; is < mesh_struct.npe ; is++) {
      if (dim == 2) {
	struct_bmat[0][is*dim + 0][gp] = struct_dsh[is][0][gp];
	struct_bmat[0][is*dim + 1][gp] = 0;
	struct_bmat[1][is*dim + 0][gp] = 0;
	struct_bmat[1][is*dim + 1][gp] = struct_dsh[is][1][gp];
	struct_bmat[2][is*dim + 0][gp] = struct_dsh[is][1][gp];
	struct_bmat[2][is*dim + 1][gp] = struct_dsh[is][0][gp];
      }
    }
  }

  return 0;
}
