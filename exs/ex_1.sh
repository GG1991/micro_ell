#!/bin/bash

#../micro \
gdb -x file.gdb --args ../micro \
    -struct_n 4,4 \
    -dim 2 \
    -material "MATRIX MAT_ELASTIC 1.0e7 1.0e6 0.3","FIBER MAT_ELASTIC 1.0e7 1.0e7 0.3" \
    -micro_struct "fiber_cilin 3.0 3.0 1 1 0.75 0.0 0.0" \
    -nl_max_its 2 \
    -bc_ustrain \
    -print_pvtu

#-bc_periodic \
#-bc_ustrain \
#-bc_ustress \

#-micro_struct  "size"[dim]     \
#               "nx_fib"        \
#               "ny_fib"        \
#               "radio"         \
#               "desv"[2]
#
#-micro_struct  "fiber_line     \
#               "size"[dim]     \
#               "ntype"         \
#               "nfib[ntype]"   \
#               "tetha[ntype]"  \
#               "seps[ntype]"   \
#               "width[ntype]"  \
#               "desv[ntype]"

#-micro_struct  "fiber_cilin 3.0 3.0 1 1 0.75 0.0 0.0" \
#-micro_struct  "fiber_line 3.0 3.0 1 1 0.0 0.4 1.0 0.0" \
#-micro_struct  "fiber_line 3.0 3.0 2 9 9 0.785398 2.35619 0.4 0.4 0.2 0.2 0.0 0.0" \
