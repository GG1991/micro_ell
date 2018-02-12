#include "mesh.h"

int mesh_struct_init(int dim, int *sizes, double *length, mesh_struct_t *mesh_struct)
{
  mesh_struct->dim = dim;

  mesh_struct->lx = length[0];
  mesh_struct->ly = length[1];
  mesh_struct->lz = (dim == 3) ? length[2] : 1;
  mesh_struct->nx = sizes[0];
  mesh_struct->ny = sizes[1];
  mesh_struct->nz = (mesh_struct->dim == 3) ? sizes[2] : 1;
  mesh_struct->nn = mesh_struct->nx * mesh_struct->ny * mesh_struct->nz;

  mesh_struct->nex = (mesh_struct->nx - 1);
  mesh_struct->ney = (mesh_struct->ny - 1);
  mesh_struct->nez = (mesh_struct->dim == 3) ? (mesh_struct->nz - 1) : 1;
  mesh_struct->nelm = mesh_struct->nex * mesh_struct->ney * mesh_struct->nez;

  mesh_struct->hx = mesh_struct->lx / mesh_struct->nex;
  mesh_struct->hy = mesh_struct->ly / mesh_struct->ney;
  mesh_struct->hz = (mesh_struct->dim == 2) ? 1 : (mesh_struct->lz / mesh_struct->nez);

  mesh_struct->vol = mesh_struct->lx * mesh_struct->ly * mesh_struct->lz;
  mesh_struct->vol_elm = mesh_struct->hx * mesh_struct->hy * mesh_struct->hz;

  mesh_struct->npe = (mesh_struct->dim == 2) ? 4 : 8;

  if (dim == 2) {
    mesh_struct->nnods_boundary = 2*mesh_struct->nx + 2*(mesh_struct->ny - 2);
    mesh_struct->nods_x0 = malloc((mesh_struct->ny - 2) * sizeof(int));
    mesh_struct->nods_x1 = malloc((mesh_struct->ny - 2) * sizeof(int));
    mesh_struct->nods_y0 = malloc((mesh_struct->nx - 2) * sizeof(int));
    mesh_struct->nods_y1 = malloc((mesh_struct->nx - 2) * sizeof(int));
    mesh_struct->coor_x0 = malloc((mesh_struct->ny - 2)*mesh_struct->dim*sizeof(double));
    mesh_struct->coor_x1 = malloc((mesh_struct->ny - 2)*mesh_struct->dim*sizeof(double));
    mesh_struct->coor_y0 = malloc((mesh_struct->nx - 2)*mesh_struct->dim*sizeof(double));
    mesh_struct->coor_y1 = malloc((mesh_struct->nx - 2)*mesh_struct->dim*sizeof(double));
  }
  mesh_struct->boundary_nods = malloc(mesh_struct->nnods_boundary * sizeof(int));
  mesh_struct->boundary_coord = malloc(mesh_struct->nnods_boundary * mesh_struct->dim * sizeof(double));
  mesh_struct->boundary_indeces = malloc(mesh_struct->nnods_boundary * mesh_struct->dim * sizeof(int));

  if (dim == 2) {

    /*  ________
     *  |   2   |
     *  |3     4|
     *  |___1___|
     */
    mesh_struct->nod_x0y0 = 0;
    mesh_struct->nod_x1y0 = mesh_struct->nx - 1;
    mesh_struct->nod_x1y1 = mesh_struct->nx*mesh_struct->ny - 1;
    mesh_struct->nod_x0y1 = mesh_struct->nx*(mesh_struct->ny - 1);

    int node_count = 0;
    for (int n = 0; n < mesh_struct->nx; n++) { // y = 0
      mesh_struct->boundary_nods[node_count] = n;
      mesh_struct->boundary_coord[node_count*dim + 0] = n * mesh_struct->hx;
      mesh_struct->boundary_coord[node_count*dim + 1] = 0.0;
      node_count++;
    }

    for (int n = 0; n < mesh_struct->nx; n++) { // y = ly
      mesh_struct->boundary_nods[node_count] = (mesh_struct->ny - 1) * mesh_struct->nx + n;
      mesh_struct->boundary_coord[node_count*dim + 0] = n * mesh_struct->hx;
      mesh_struct->boundary_coord[node_count*dim + 1] = mesh_struct->ly;
      node_count++;
    }

    for (int n = 0; n < (mesh_struct->ny - 2); n++) { // x = 0
      mesh_struct->boundary_nods[node_count] = (n + 1) * mesh_struct->nx;
      mesh_struct->boundary_coord[node_count*dim + 0] = 0.0;
      mesh_struct->boundary_coord[node_count*dim + 1] = (n + 1)*mesh_struct->hy;
      node_count++;
    }

    for (int n = 0; n < (mesh_struct->ny - 2); n++) { // x = lx
      mesh_struct->boundary_nods[node_count] = (n + 2) * mesh_struct->nx - 1;
      mesh_struct->boundary_coord[node_count*dim + 0] = mesh_struct->lx;
      mesh_struct->boundary_coord[node_count*dim + 1] = (n + 1)*mesh_struct->hy;
      node_count++;
    }

    for (int n = 0; n < mesh_struct->nx - 2; n++) { // y = 0 & y = ly
      mesh_struct->nods_y0[n] = n + 1;
      mesh_struct->nods_y1[n] = (mesh_struct->ny - 1) * mesh_struct->nx + n + 1;
      mesh_struct->coor_y0[n*dim + 0] = (n + 1)*mesh_struct->hx;
      mesh_struct->coor_y0[n*dim + 1] = 0.0;
      mesh_struct->coor_y1[n*dim + 0] = (n + 1)*mesh_struct->hx;
      mesh_struct->coor_y1[n*dim + 1] = mesh_struct->ly;
    }
    for (int n = 0; n < mesh_struct->ny - 2; n++) { // x = 0 & x = lx
      mesh_struct->nods_x0[n] = (n + 1) * mesh_struct->nx;
      mesh_struct->nods_x1[n] = (n + 2) * mesh_struct->nx - 1;
      mesh_struct->coor_x0[n*dim + 0] = 0.0;
      mesh_struct->coor_x0[n*dim + 1] = (n + 1)*mesh_struct->hy;
      mesh_struct->coor_x1[n*dim + 0] = mesh_struct->lx;
      mesh_struct->coor_x1[n*dim + 1] = (n + 1)*mesh_struct->hy;
    }
  }

  for (int i = 0; i < mesh_struct->nnods_boundary ; i++)
    for (int d = 0; d < mesh_struct->dim; d++)
      mesh_struct->boundary_indeces[i*mesh_struct->dim + d] = mesh_struct->boundary_nods[i] * mesh_struct->dim + d;

  return 0;
}

int mesh_struct_get_node_coord(mesh_struct_t *mesh_struct, int node, double *coord)
{
  if (mesh_struct->dim == 2) {
    coord[0] = (node % mesh_struct->nx) * mesh_struct->hx;
    coord[1] = (node / mesh_struct->nx) * mesh_struct->hy;
  }
  return 0;
}

int mesh_struct_get_elem_nods(mesh_struct_t *mesh_struct, int e, int *elem_nods)
{
  if (mesh_struct->dim == 2) {
    int n0 = (e % mesh_struct->nex) + (e / mesh_struct->nex) * mesh_struct->nx;
    elem_nods[0] = n0;
    elem_nods[1] = n0 + 1;
    elem_nods[2] = n0 + mesh_struct->nx + 1;
    elem_nods[3] = n0 + mesh_struct->nx;
  }
  return 0;
}

int mesh_struct_get_elem_indeces(mesh_struct_t *mesh_struct, int e, int *elem_indeces)
{
  if (mesh_struct->dim == 2) {
    for (int d = 0; d < mesh_struct->dim; d++) {
      int n0 = (e % mesh_struct->nex) + (e / mesh_struct->nex) * mesh_struct->nx;
      elem_indeces[0*mesh_struct->dim + d] = (n0) * mesh_struct->dim + d;
      elem_indeces[1*mesh_struct->dim + d] = (n0 + 1) * mesh_struct->dim + d;
      elem_indeces[2*mesh_struct->dim + d] = (n0 + mesh_struct->nx + 1) * mesh_struct->dim + d;
      elem_indeces[3*mesh_struct->dim + d] = (n0 + mesh_struct->nx) * mesh_struct->dim + d;
    }
  }
  return 0;
}
