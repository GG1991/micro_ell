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
  if(ierr != 0) goto end;

  nvoi = (dim == 2) ? 3 : 6;

  char string_buf[NBUF];
  ierr = myio_comm_line_get_string(&command_line, "-micro_struct", string_buf, &found);
  ierr = micro_struct_init(dim, string_buf, &micro_struct);

  int nval_found, nval_expect = (dim == 2) ? 2 : 3;
  int values_int[3];
  myio_comm_line_get_int_array(&command_line, "-struct_n", values_int, nval_expect, &nval_found, &found);

  if (found == true) {
    if (nval_found != nval_expect) {
      printf(RED "-struct_n should include %d arguments\n" NOR, nval_expect);
      return 1;
    }
  }else{
    printf(RED "-struct_n is request\n" NOR, nval_expect);
    return 1;
  }

  mesh_struct_init(dim, values_int, micro_struct.size, &mesh_struct);

  ngp = (dim == 2) ? 4 : 8;

  ierr = material_fill_list_from_command_line(&command_line, &material_list);
  if (ierr != 0) goto end;

  printf(GRE
      "--------------------------------------------------\n"
      "  MICRO: START\n"
      "--------------------------------------------------" NOR"\n");

  ierr = alloc_memory();

  ierr = micro_struct_init_elem_type(&micro_struct, dim, mesh_struct.nelm, &get_elem_centroid, elem_type);
  if (ierr != 0) goto end;

  ierr = micro_check_material_and_elem_type(&material_list, elem_type, mesh_struct.nelm);
  if (ierr != 0) goto end;

  double h[3]; h[0] = mesh_struct.hx; h[1] = mesh_struct.hy; h[2] = mesh_struct.hz;
  fem_init();
  fem_init_struct(&struct_sh, &struct_dsh, &struct_wp, h, mesh_struct.dim);

  init_variables_2();

  double strain[3] = {0.0, 0.0, 0.0};
  double ctang[9];
  get_ctang (strain, ctang);

  printf("\nctang = \n");
  for (int i = 0 ; i < (nvoi*nvoi) ; i++)
    printf("%e%s", (fabs(ctang[i])>1.0) ? ctang[i] : 0.0, ((i+1) % nvoi == 0) ? "\n":" ");

end:

  printf(GRE
      "--------------------------------------------------\n"
      "  MICRO: FINISH\n"
      "--------------------------------------------------" NOR"\n");

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
