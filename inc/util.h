#ifndef UTIL_H
#define UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int strbin2dec(char *str);
int util_is_in_vector(int val, int *vector, int size);
int util_clean_and_sort_vector(int *in_vector, int n_in, int **out_vector, int *n_out);
int util_cmpfunc(const void * a, const void * b);
int util_sort_vector_intersec(int *array_1, int n1, int *array_2, int n2, int **inter_arr, int *n_inter);

#endif
