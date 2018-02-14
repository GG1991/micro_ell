#ifndef MICRO_H
#define MICRO_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_linalg.h>
#include "list.h"
#include "material.h"
#include "myio.h"
#include "vtk.h"
#include "micro_struct.h"
#include "solvers.h"
#include "mesh.h"
#include "fem.h"
#include "ell.h"

enum {MULTIS_NULL, MULTIS_MIXS, MULTIS_MIXP, MULTIS_FE2};
enum {BC_NULL, BC_USTRAIN, BC_USTRESS, BC_PER_MS, BC_PER_LM};
enum {SOL_PETSC, SOL_ELL};

#define DELTA_EPS 0.005

#define GREEN  "\x1B[32m"
#define RED    "\x1B[31m"
#define BLUE   "\x1B[34m"
#define NORMAL "\x1B[0m"

#define MAX_NVOIGT   6

#define NBUF         256

int rank_mic;
int nproc_mic;
int dim;
int nvoi;

int ngp;
int *elem_index;
int *elem_nods;
int *elem_type;
double **struct_sh;
double ***struct_dsh;
double *struct_wp;
double ***struct_bmat;
double *elem_disp;
double *strain_gp;
double *stress_gp;
double *elem_strain;
double *elem_stress;
double *elem_energy;

ell_matrix jac_ell;
double *x_ell;
double *res_ell;
double *dx_ell;

typedef struct {

  int multis_method;
  int fe2_bc;
  int nl_max_its;
  int solver;

  double nl_min_norm;
  double rho;

} params_t;

typedef struct {

  bool coupled;
  bool allocated;
  bool linear_materials;
  bool print_pvtu;
  bool print_vectors;
  bool print_matrices;

} flags_t;

extern params_t params;
extern flags_t flags;

double center_domain[3];

int micro_print_info( void );
int micro_pvtu( char *name );

int micro_check_material_and_elem_type(list_t *material_list, int *elem_type, int nelm);

int get_local_elem_index(int e, int *loc_index);
int get_global_elem_index(int e, int *glo_elem_index);
int get_stress(int e , int gp, double *strain_gp, double *stress_gp);
int get_strain(int e , int gp, double *strain_gp);
int get_c_tan(const char * name, int e, int gp, double *strain_gp, double *c_tan);
int get_elem_centroid( int e, int dim, double *centroid);
int strain_x_coord( double * strain , double * coord , double *u);
int get_node_local_coor(int n, double *coord);
int get_node_ghost_coor(int n, double *coord);
int get_local_elem_node(int e, int *n_loc);
int local_to_global_index(int local);
int get_averages(double * strain_ave, double *stress_ave);
int get_elem_type(int e, int *type);
int get_elem_properties(void);

int assembly_jac_ell(void);
int assembly_res_ell(double *norm, double *strain_mac);

int init_variables_1(void);
int init_variables_2(void);
int finalize(void);

int alloc_memory(void);

int comm_line_set_flags(void);

int localize_strain (double *strain);
int homogenize_stress (double *strain, double *stress);
int get_ctang (double *strain, double *ctang);

int set_disp_0(double *strain_mac);
int negative_res(void);
int add_x_dx(void);
int solve(void);

#endif
