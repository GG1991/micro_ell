#ifndef MESH_H
#define MESH_H

#include "list.h"
#include "util.h"
#include "myio.h"
#include "geometry.h"

#define NBUF 256

#define MAX_ADJ_NODES 30
#define MAX_NUM_OF_BOUNDARIES 4
#define MAX_NPE 8
#define MAX_DIM 3

typedef struct{

  int dim;
  int nx, ny, nz;
  int nex, ney, nez;
  int nelm;
  int nn;
  int npe;
  int nnods_boundary;
  int *boundary_nods;
  int *boundary_indeces;
  int *nods_x0;
  int *nods_x1;
  int *nods_y0;
  int *nods_y1;
  int nod_x0y0;
  int nod_x1y0;
  int nod_x1y1;
  int nod_x0y1;
  double *boundary_coord;
  double lx, ly, lz;
  double hx, hy, hz;
  double vol;
  double vol_elm;
  double *coor_x0;
  double *coor_x1;
  double *coor_y0;
  double *coor_y1;

}mesh_struct_t;

extern mesh_struct_t mesh_struct;

int mesh_struct_init(int dim, int *sizes, double *length, mesh_struct_t *mesh_struct);
int mesh_struct_get_node_coord(mesh_struct_t *mesh_struct, int node, double *coord);
int mesh_struct_get_elem_nods(mesh_struct_t *mesh_struct, int e, int *elem_nods);
int mesh_struct_get_elem_indeces(mesh_struct_t *mesh_struct, int e, int *elem_indeces);

#endif
