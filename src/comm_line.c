#include "micro.h"

int comm_line_set_flags(void)
{
  bool found;

  myio_comm_line_get_int(&command_line, "-dim", &dim, &found); if (found == false) { printf("-dim not given on command line.\n"); return 1; }

  myio_comm_line_search_option(&command_line, "-bc_ustrain", &found); if (found == true) params.fe2_bc = BC_USTRAIN;
  myio_comm_line_search_option(&command_line, "-bc_per_ms", &found); if (found == true) params.fe2_bc = BC_PER_MS;
  myio_comm_line_search_option(&command_line, "-bc_per_lm", &found); if (found == true) params.fe2_bc = BC_PER_LM;
  myio_comm_line_search_option(&command_line, "-bc_ustress", &found); if (found == true) params.fe2_bc = BC_USTRESS;
  if (params.fe2_bc == BC_NULL) { printf("comm_line.c : A BC should be specify : (-bc_ustrain -bc_ustress -bc_ustress)\n"); return 2; }

  myio_comm_line_search_option(&command_line, "-print_matrices", &found); if (found == true) flags.print_matrices = true;
  myio_comm_line_search_option(&command_line, "-print_vectors", &found); if (found == true) flags.print_vectors = true;
  myio_comm_line_search_option(&command_line, "-print_pvtu", &found); if (found == true) flags.print_pvtu = true;

  myio_comm_line_get_int(&command_line, "-nl_max_its", &params.nl_max_its, &found);
  myio_comm_line_get_double(&command_line, "-nl_min_norm", &params.nl_min_norm, &found);

  return 0;
}
