#include "gmsh.h"


int gmsh_get_node_index(const char *mesh_n, const char *bou_name, int nnods, int *nods_arr, int dim, int *n, int **ix) {

  FILE *fm = fopen(mesh_n,"r"); if (fm == NULL) return 1;

  list_t nod_list;
  list_init(&nod_list, sizeof(int), &gmsh_funcmp_int_a);

  int id_s = gmsh_which_id(mesh_n, bou_name); if (id_s < 0) return 1;
  int flag = 0;
  char buf[NBUF_GMSH];

  while (fgets(buf, NBUF_GMSH, fm) != NULL) {

    char *data = strtok(buf," \n");
    if (flag == 0) {

      if (strcmp(data,"$Elements") == 0) {
	fgets(buf, NBUF_GMSH, fm);
	flag = 1;
      }

    }else{

      data = strtok(NULL, " \n");

      if (gmsh_is_surf(atoi(data), dim) != 0) {

	int npe = gmsh_npe(atoi(data));
	data = strtok(NULL," \n");
	int ntag = atoi(data);
	data = strtok(NULL," \n");
	int id = atoi(data);

	if (id == id_s) {

	  int d = 1;
	  while (d < ntag) {
	    data = strtok(NULL, " \n");
	    d++;
	  }
	  int i = 0;
	  while (i < npe) {
	    data = strtok(NULL, " \n");
	    int ns = atoi(data) - 1;
	    int *nf = bsearch(&ns, nods_arr, nnods, sizeof(int), &gmsh_funcmp_int_b);
	    if (nf != NULL)
	      list_insert_se(&nod_list, nf);
	    i++;
	  }
	}

      }else{

	*n = nod_list.sizelist;
	*ix = malloc( *n * sizeof(int));

	node_list_t *pn = nod_list.head;
	int i = 0;
	while (pn != NULL) {
	  (*ix)[i++] = *(int *)pn->data;
	  pn = pn->next;
	}
	list_clear(&nod_list);
	fclose(fm);
	return 0;
      }
    }
  }

  return 1;
}


int gmsh_which_id(const char *mesh_n, const char *name) {

  FILE         *fm = fopen(mesh_n,"r"); if ( fm == NULL ) return 1;
  int          id;
  int          i = 0, total;
  int          flag = 0;
  char         buf[NBUF_GMSH], *data;

  while (fgets(buf, NBUF_GMSH, fm) != NULL) {

    data = strtok(buf," \n");

    if (flag == 0) {
      if (strcmp(data,"$PhysicalNames") == 0) {
	fgets( buf, NBUF_GMSH, fm );
	data  = strtok(buf," \n");
	total = atoi(data);
	flag = 1;
      }
    }else{
      data = strtok( NULL, " \n" );
      id   = atoi(data);
      data = strtok( NULL, " \"\n" );
      if ( strcmp( name , data ) == 0 )
	return id;
      i++;
    }
    if (i == total) break;

  }

  return -1;
}


int gmsh_get_physical_list(char *mesh_n, list_t *physical_list) {

  FILE  *fm;
  int   i = 0, total;
  int   flag = 0;
  char  buf[NBUF_GMSH], *data;

  physical_t physical;

  fm = fopen(mesh_n,"r"); if (fm == NULL) return 1;

  while (fgets(buf,NBUF_GMSH,fm) != NULL) {

    data=strtok(buf," \n");
    if (strcmp(data,"$PhysicalNames") == 0) {

      fgets( buf, NBUF_GMSH, fm );
      data  = strtok(buf," \n");
      total = atoi(data);
      flag  = 1;

    }else if (flag == 1) {

      physical.dim = atoi(data);
      data=strtok(NULL," \n");
      physical.id = atoi(data);
      data=strtok(NULL," \"\n");
      physical.name = strdup( data );

      list_insertlast(physical_list, &physical);
      i ++;

      if (i == total) break;
    }
  }
  return 0;
}


int gmsh_read_mesh(MPI_Comm COMM, const char *gmsh_file, gmsh_mesh_t *gmsh_mesh) {

  int rank, nproc;
  MPI_Comm_size(COMM, &nproc);
  MPI_Comm_rank(COMM, &rank);

  FILE *fm = fopen(gmsh_file,"r"); if (fm == NULL) return 1;

  char buf[NBUF_GMSH];
  int offset;

  myio_file_get_offset_line_start_word(gmsh_file, "$Nodes", &offset);
  fseek(fm, offset, SEEK_SET);
  fgets(buf, NBUF_GMSH, fm);
  fgets(buf, NBUF_GMSH, fm);

  char *str_token = strtok(buf, " \n");
  gmsh_mesh->nnods = atoi(str_token);
  gmsh_mesh->coord = malloc(gmsh_mesh->nnods*gmsh_mesh->dim*sizeof(double));
  for (int i = 0 ; i < gmsh_mesh->nnods ; i++) {
    fgets(buf, NBUF_GMSH, fm);
    char *str_token = strtok(buf, " \n");
    if (strcmp(str_token, "$EndNodes") == 0) break;
    for (int d = 0 ; d < 3 ; d++) {
      str_token = strtok(NULL, " \n");
      gmsh_mesh->coord[i*gmsh_mesh->dim + d] = atof(str_token);
    }
  }

  myio_file_get_offset_line_start_word(gmsh_file, "$Elements", &offset);
  fseek(fm, offset, SEEK_SET);
  fgets(buf, NBUF_GMSH, fm);
  fgets(buf, NBUF_GMSH, fm);

  int nelm_total = 0, nelm_surf = 0;

  while (fgets(buf, NBUF_GMSH, fm) != NULL) {
    char *str_token = strtok(buf, " \n");
    if (strcmp(str_token, "$EndElements") == 0) break;
    str_token = strtok(NULL, " \n");
    if (gmsh_is_vol_elm(gmsh_mesh->dim, atoi(str_token)) != 0)
      nelm_total ++;
    else
      nelm_surf ++;
  }

  gmsh_mesh->nelm_total = nelm_total;
  gmsh_mesh->nelm_surf = nelm_surf;
  gmsh_mesh->elem_per_proc = malloc(nproc*sizeof(int));
  gmsh_mesh->nelm_dist = malloc((nproc + 1)*sizeof(int));

  for (int i = 0 ; i < nproc ; i++)
    gmsh_mesh->elem_per_proc[i] = gmsh_mesh->nelm_total/nproc + ((i < gmsh_mesh->nelm_total % nproc) ? 1 : 0);

  gmsh_mesh->nelm_local = gmsh_mesh->elem_per_proc[rank];

  gmsh_mesh->nelm_dist[0] = 0;
  for (int i = 1 ; i < nproc + 1 ; i++)
    gmsh_mesh->nelm_dist[i] = gmsh_mesh->nelm_dist[i-1] + gmsh_mesh->elem_per_proc[i-1];

  gmsh_mesh->eptr = malloc((gmsh_mesh->nelm_local + 1)*sizeof(int));
  gmsh_mesh->elem_centroid = malloc(gmsh_mesh->nelm_local * gmsh_mesh->dim * sizeof(double));
  gmsh_mesh->elm_id = malloc(gmsh_mesh->nelm_local*sizeof(int));

  fseek(fm, offset, SEEK_SET);
  fgets(buf, NBUF_GMSH, fm);
  fgets(buf, NBUF_GMSH, fm);

  for (int i = 0 ; i < gmsh_mesh->nelm_surf ; i++)
    fgets(buf, NBUF_GMSH, fm);

  for (int i = 0 ; i < gmsh_mesh->nelm_dist[rank] ; i++)
    fgets(buf, NBUF_GMSH, fm);

  gmsh_mesh->eptr[0] = 0;
  for (int i = 0 ; i < gmsh_mesh->nelm_local ; i++) {
    fgets(buf, NBUF_GMSH, fm);
    char *str_token = strtok(buf, " \n");
    str_token = strtok(NULL, " \n");
    int npe = gmsh_npe(atoi(str_token));
    gmsh_mesh->eptr[i+1] = gmsh_mesh->eptr[i] + npe;
  }

  gmsh_mesh->eind = malloc(gmsh_mesh->eptr[gmsh_mesh->nelm_local] * sizeof(int));

  fseek(fm, offset, SEEK_SET);
  fgets(buf, NBUF_GMSH, fm);
  fgets(buf, NBUF_GMSH, fm);

  for (int i = 0 ; i < gmsh_mesh->nelm_surf ; i++)
    fgets(buf, NBUF_GMSH, fm);

  for (int i = 0 ; i < gmsh_mesh->nelm_dist[rank] ; i++)
    fgets(buf, NBUF_GMSH, fm);

  int start_eind = 0;
  for (int i = 0 ; i < gmsh_mesh->nelm_local ; i++) {
    fgets(buf, NBUF_GMSH, fm);
    char *str_token = strtok(buf, " \n");

    str_token = strtok(NULL, " \n");
    int npe = gmsh_npe(atoi(str_token));

    str_token = strtok(NULL, " \n");
    int ntag = atoi(str_token);

    str_token = strtok(NULL, " \n");
    gmsh_mesh->elm_id[i] = atoi(str_token);

    for (int j = 0 ; j < ntag - 1 ; j++)
      str_token = strtok(NULL, " \n");

    for (int j = 0 ; j < npe ; j++) {
      str_token = strtok(NULL, " \n");
      gmsh_mesh->eind[start_eind + j] = atoi(str_token);
    }
    start_eind += npe;
  }

  fclose(fm);

  return 0;
}


int gmsh_funcmp_int_a(void *a, void *b) {
     if (*(int*)a > *(int*)b) return  1;
     if (*(int*)a < *(int*)b) return -1;
     return 0;
}


int gmsh_funcmp_int_b(const void *a, const void *b) {
     if (*(int*)a > *(int*)b) return  1;
     if (*(int*)a < *(int*)b) return -1;
     return 0;
}


int gmsh_is_surf( int code , int dim) {
  if (dim == 2)
    return (code == 1 || code == 15) ? 1 : 0;
  else if (dim == 3)
    return (code == 1 || code == 2 || code == 3 || code == 15) ? 1 : 0;
  return 1;
}


int gmsh_is_vol_elm(int dim, int code) {

  if (dim == 2)
    return (code == 2 || code == 3) ? 1 : 0;

  else if (dim == 3)
    return (code == 4 || code == 5 || code == 6) ? 1 : 0;

  return -1;
}


int gmsh_npe(int code) {

    switch(code) {
	case 1:
	    return 2;
	case 2:
	    return 3;
	case 3:
	    return 4;
	case 4:
	    return 4;
	case 5:
	    return 8;
	case 6:
	    return 6;
	case 15:
	    return 1;
	default:
	    return -1;
    }
}
