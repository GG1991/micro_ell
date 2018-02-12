#include "micro.h"

static char help[] =
"micro multiscale code \n"
"-homo_ts     : c =  vi ci + vm cm            (serial) \n"
"-homo_tp     : c = (vi ci^-1 + vm cm^-1)^-1  (parallel) \n"
"-homo_us     : homogenization using uniform strains approach \n"
"-struct_n [<nx,ny>] if dim = 2 \n"
"-print_matrices \n"
"-print_vectors \n"
"-print_pvtu \n";

params_t params;
flags_t flags;
solver_t solver;
mesh_struct_t mesh_struct;

#define CHECK(error, message) {\
  if (error != 0) {\
    printf( "%s\n", message);\
    goto end;}}

int main(int argc, char **argv)
{
  int ierr;
  bool found;

  myio_comm_line_init(argc, argv, &command_line);

  myio_comm_line_search_option(&command_line, "-help", &found);
  if (found == true) {
    printf("%s", help);
    goto end;
  }

  init_variables_1();

  ierr = comm_line_set_flags();
  if(ierr != 0) CHECK(ierr, "error in reading command line");

  nvoi = (dim == 2) ? 3 : 6;

  char string_buf[NBUF];
  ierr = myio_comm_line_get_string(&command_line, "-micro_struct", string_buf, &found);
  ierr = micro_struct_init(dim, string_buf, &micro_struct);

  int nval_found, nval_expect = (dim == 2) ? 2 : 3;
  int values_int[3];
  myio_comm_line_get_int_array(&command_line, "-struct_n", values_int, nval_expect, &nval_found, &found);

  if (found == true) {
    if (nval_found != nval_expect) {
      printf(RED "-struct_n should include %d arguments\n" NORMAL, nval_expect);
      return 1;
    }
  }else{
    printf(RED "-struct_n is request\n" NORMAL, nval_expect);
    return 1;
  }

  mesh_struct_init(dim, values_int, micro_struct.size, &mesh_struct);

  ngp = (dim == 2) ? 4 : 8;

  ierr = material_fill_list_from_command_line(&command_line, &material_list);
  CHECK(ierr, "error parsing material from command line");

  printf(GREEN
      "--------------------------------------------------\n"
      "  MICRO: START\n"
      "--------------------------------------------------" NORMAL "\n\n");

  ierr = alloc_memory();

  ierr = micro_struct_init_elem_type(&micro_struct, dim, mesh_struct.nelm, &get_elem_centroid, elem_type);
  CHECK(ierr, "error initializing elem_type array using micro_struct");

  ierr = micro_check_material_and_elem_type(&material_list, elem_type, mesh_struct.nelm);
  CHECK(ierr, "error checking elem_type and material_list");

  double h[3]; h[0] = mesh_struct.hx; h[1] = mesh_struct.hy; h[2] = mesh_struct.hz;
  fem_init();
  fem_init_struct(&struct_sh, &struct_dsh, &struct_wp, h, mesh_struct.dim);

  init_variables_2();

  printf("\nctang_0 = \n");
  for (int i = 0 ; i < nvoi ; i++) {
    for (int j = 0 ; j < nvoi ; j++)
      printf("%e ", (fabs(params.c_tangent_linear[i*nvoi+j])>1.0) ? params.c_tangent_linear[i*nvoi+j] : 0.0);
    printf("\n");
  }
  printf("\n");

  micro_print_info();

end:

  printf(GREEN
      "--------------------------------------------------\n"
      "  MICRO: FINISH\n"
      "--------------------------------------------------" NORMAL "\n");

end_no_message:

  ierr = finalize();

  return ierr;
}


int micro_check_material_and_elem_type(list_t *material_list, int *elem_type, int nelm) {

  char *word_to_search;

  for (int e ; e < nelm ; e++ ) {

    switch(elem_type[e]) {

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
    node_list_t *pm = material_list->head;

    while (pm != NULL) {
      mat_p = (material_t *)pm->data;
      if (strcmp(mat_p->name, word_to_search) == 0) break;
      pm = pm->next;
    }

    if (pm == NULL) return 1;
  }

  return 0;
}
