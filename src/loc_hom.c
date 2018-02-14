#include "micro.h"

int localize_strain (double *strain)
{
  clock_t start, end;
  double time_ass_b = 0.0, time_ass_A = 0.0, time_sol = 0.0;

  set_disp_0(strain);

  int nl_its = 0;
  double res_norm = params.nl_min_norm * 10;
  while (nl_its < params.nl_max_its && res_norm > params.nl_min_norm) {

    start = clock();
    assembly_res_ell(&res_norm, strain);
    end = clock();
    time_ass_b = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf(GREEN "|b| = %lf" NORMAL "\n", res_norm);

    if (res_norm < params.nl_min_norm) break;

    negative_res();

    start = clock();
    assembly_jac_ell();
    end = clock();
    time_ass_A = ((double) (end - start)) / CLOCKS_PER_SEC;

    start = clock();
    solve();
    end = clock();
    time_sol = ((double) (end - start)) / CLOCKS_PER_SEC;
 
    add_x_dx();

    nl_its ++;
  }

  printf(BLUE "time ass_b = %lf" NORMAL "\n", time_ass_b);
  printf(BLUE "time ass_A = %lf" NORMAL "\n", time_ass_A);
  printf(BLUE "time sol   = %lf" NORMAL "\n", time_sol);
  return 0;
}

int homogenize_stress (double *strain, double *stress)
{
  for (int i = 0; i < nvoi; i++) stress[i] = strain[i] = 0.0;

  for (int e = 0 ; e < mesh_struct.nelm ; e++) {
    for (int gp = 0; gp < ngp; gp++) {
      get_strain(e, gp, strain_gp);
      get_stress(e, gp, strain_gp, stress_gp);
      for (int i = 0; i < nvoi; i++) {
	stress[i] += stress_gp[i] * struct_wp[gp];
	strain[i] += strain_gp[i] * struct_wp[gp];
      }
    }
  }

  for (int i = 0; i < nvoi; i++) {
    stress[i] /= mesh_struct.vol;
    strain[i] /= mesh_struct.vol;
  }
  return 0;
}

int get_ctang (double *strain, double *ctang)
{
  int ierr;
  double *eps_2 = malloc(nvoi * sizeof(double));
  double *sig_1 = malloc(nvoi * sizeof(double));
  double *sig_2 = malloc(nvoi * sizeof(double));
  double *eps_dummy = malloc(nvoi * sizeof(double));

  for (int i = 0 ; i < nvoi ; i++)
    printf("%s%lf%s", (i == 0)?"calc stress in for strain: ":" ", strain[i], (i == (nvoi-1))?"\n":"");

  for (int i = 0 ; i < nvoi ; i++)
    eps_2[i] = strain[i];

  ierr = localize_strain(strain);
  ierr = homogenize_stress(eps_dummy, sig_1); // we can check with eps_dummy if the localization have been done well

  for (int i = 0 ; i < nvoi ; i++) {

    printf("exp %d\n", i);
    for (int j = 0 ; j < nvoi ; j++)
      eps_2[j] = strain[j];
    eps_2[i] = eps_2[i] + DELTA_EPS;

    ierr = localize_strain(eps_2);
    ierr = homogenize_stress(eps_dummy, sig_2);

    for (int j = 0 ; j < nvoi ; j++)
      ctang[j*nvoi + i] = (sig_2[j] - sig_1[j]) / (eps_2[i] - strain[i]);

    if (flags.print_pvtu == true) {
      get_elem_properties();
      char filename[64]; sprintf(filename,"micro_exp%d",i);
      ierr = micro_pvtu(filename);
    }
  }
  return ierr;
}

int set_disp_0(double *strain_mac)
{
  int nn = mesh_struct.nn;
  int dim = mesh_struct.dim;
  double disp[2];
  double coor[2];

  if (params.fe2_bc == BC_USTRAIN) {
    for (int n = 0 ; n < (nn*dim) ; n++) x_ell[n] = 0.0;
    for (int n = 0 ; n < (mesh_struct.ny - 2) ; n++) { 
      // X0
      coor[0] = 0.0;
      coor[1] = (n + 1)*mesh_struct.hy;
      strain_x_coord(strain_mac, coor, disp);
      for (int d = 0; d < dim ; d++) 
	x_ell[mesh_struct.nods_x0[n]*dim + d] = disp[d];
      // X1
      coor[0] = mesh_struct.lx;
      coor[1] = (n + 1)*mesh_struct.hy;
      strain_x_coord(strain_mac, coor, disp);
      for (int d = 0; d < dim ; d++) 
	x_ell[mesh_struct.nods_x1[n]*dim + d] = disp[d];
    }
    for (int n = 0 ; n < (mesh_struct.nx - 2) ; n++) { 
      // Y0
      coor[0] = (n + 1)*mesh_struct.hx;
      coor[1] = 0.0;
      strain_x_coord(strain_mac, coor, disp);
      for (int d = 0; d < dim ; d++) 
	x_ell[mesh_struct.nods_y0[n]*dim + d] = disp[d];
      // Y1
      coor[0] = (n + 1)*mesh_struct.hx;
      coor[1] = mesh_struct.ly;
      strain_x_coord(strain_mac, coor, disp);
      for (int d = 0; d < dim ; d++) 
	x_ell[mesh_struct.nods_y1[n]*dim + d] = disp[d];
    }

    coor[0] = 0.0;
    coor[1] = 0.0;
    strain_x_coord(strain_mac, coor, disp);
    for (int d = 0; d < dim ; d++) 
      x_ell[mesh_struct.nod_x0y0*dim + d] = disp[d];

    coor[0] = mesh_struct.lx;
    coor[1] = 0.0;
    strain_x_coord(strain_mac, coor, disp);
    for (int d = 0; d < dim ; d++) 
      x_ell[mesh_struct.nod_x1y0*dim + d] = disp[d];

    coor[0] = mesh_struct.lx;
    coor[1] = mesh_struct.ly;
    strain_x_coord(strain_mac, coor, disp);
    for (int d = 0; d < dim ; d++) 
      x_ell[mesh_struct.nod_x1y1*dim + d] = disp[d];

    coor[0] = 0.0;
    coor[1] = mesh_struct.ly;
    strain_x_coord(strain_mac, coor, disp);
    for (int d = 0; d < dim ; d++) 
      x_ell[mesh_struct.nod_x0y1*dim + d] = disp[d];
  }

  return 0;
}

int negative_res(void)
{
  int nn = mesh_struct.nn;
  int dim = mesh_struct.dim;
  for (int i = 0 ; i < (nn*dim) ; i++)
    res_ell[i] *= -1.0;
  return 0;
}

int add_x_dx (void)
{
  int nn = mesh_struct.nn;
  for (int i = 0 ; i < (nn*dim) ; i++)
    x_ell[i] += dx_ell[i];
  return 0;
}

int solve (void)
{
  ell_solver solver;
  solver.max_its = 400;
  solver.min_tol = 1.0e-5;
  ell_solve_cg(&solver, &jac_ell, res_ell, dx_ell);
  return 0;
}

int strain_x_coord (double *strain, double *coord, double *u)
{
  if (dim == 2) {
    u[0] = strain[0]   * coord[0] + strain[2]/2 * coord[1];
    u[1] = strain[2]/2 * coord[0] + strain[1]   * coord[1];
  }
  return 0;
}

int get_averages(double *strain_ave, double *stress_ave)
{
  for (int i = 0; i < nvoi; i++) {
    strain_ave[i] = 0.0;
    stress_ave[i] = 0.0;
  }

  for (int e = 0 ; e < mesh_struct.nelm ; e++) {

    for (int gp = 0; gp < ngp; gp++) {

      get_strain(e, gp, strain_gp);
      get_stress(e, gp, strain_gp, stress_gp);

      for (int i = 0; i < nvoi; i++) {
	stress_ave[i] += stress_gp[i] * struct_wp[gp];
	strain_ave[i] += strain_gp[i] * struct_wp[gp];
      }
    }
  }

  for (int i = 0; i < nvoi; i++) {
    stress_ave[i] /= mesh_struct.vol;
    strain_ave[i] /= mesh_struct.vol;
  }
  return 0;
}

int get_elem_properties(void)
{
  for (int e = 0; e < mesh_struct.nelm; e++) {

    double strain_aux[MAX_NVOIGT];
    double stress_aux[MAX_NVOIGT];
    for (int i = 0; i < nvoi; i++) {
      strain_aux[i] = 0.0;
      stress_aux[i] = 0.0;
    }

    for (int gp = 0 ; gp < ngp ; gp++) {

      get_strain( e , gp, strain_gp );
      get_stress( e , gp, strain_gp, stress_gp );
      for (int v = 0 ; v < nvoi ; v++ ) {
	strain_aux[v] += strain_gp[v] * struct_wp[gp];
	stress_aux[v] += stress_gp[v] * struct_wp[gp];
      }

    }
    for (int v = 0 ; v < nvoi ; v++) {
      elem_strain[ e*nvoi + v ] = strain_aux[v] / mesh_struct.vol_elm;
      elem_stress[ e*nvoi + v ] = stress_aux[v] / mesh_struct.vol_elm;
    }
  }
  return 0;
}

int get_elem_disp(int e, double *disp_e)
{
  int npe = mesh_struct.npe;
  int dim = mesh_struct.dim;

  mesh_struct_get_elem_indeces(&mesh_struct, e, elem_index);

  for (int i = 0; i < (npe*dim); i++) disp_e[i] = x_ell[elem_index[i]];
}

int get_strain(int e, int gp, double *strain_gp)
{
  int npe = mesh_struct.npe;
  int dim = mesh_struct.dim;
  get_elem_disp(e, elem_disp);

  for (int v = 0; v < nvoi; v++) {
    strain_gp[v] = 0.0;
    for (int i = 0 ; i < (npe*dim) ; i++)
      strain_gp[v] += struct_bmat[v][i][gp] * elem_disp[i];
  }
  return 0;
}


int get_stress(int e, int gp, double *strain_gp, double *stress_gp)
{
  char *word_to_search;

  switch (elem_type[e]) {

    case ID_FIBER:
      word_to_search = strdup("FIBER");
      break;

    case ID_MATRIX:
      word_to_search = strdup("MATRIX");
      break;

    default:
      return 1;
  }

  material_t  *mat_p;
  node_list_t *pm = material_list.head;

  while (pm != NULL) {
    mat_p = (material_t *)pm->data;
    if (strcmp(mat_p->name, word_to_search) == 0) break;
    pm = pm->next;
  }

  if (pm == NULL) return 1;

  return material_get_stress(mat_p, dim, strain_gp, stress_gp);
}

int get_c_tan(const char *name, int e, int gp, double *strain_gp, double *c_tan)
{
  char *word_to_search;

  if (name != NULL)
    word_to_search = strdup(name);
  else if ( elem_type[e] == ID_FIBER )
    word_to_search = strdup("FIBER");
  else if ( elem_type[e] == ID_MATRIX )
    word_to_search = strdup("MATRIX");

  material_t  *mat_p;
  node_list_t *pm = material_list.head;

  while (pm != NULL){
    mat_p = (material_t *)pm->data;
    if (strcmp(mat_p->name, word_to_search) == 0) break;
    pm = pm->next;
  }
  if (pm == NULL) return 1;

  return material_get_c_tang(mat_p, dim, strain_gp, c_tan);
}

int get_elem_centroid(int e, int dim, double *centroid)
{
  if (dim == 2) {
    centroid[0] = (e % mesh_struct.nex + 0.5) * mesh_struct.hx;
    centroid[1] = (e / mesh_struct.nex + 0.5) * mesh_struct.hy;
  }
  return 0;
}
