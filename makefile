
SRC_DIR= ./src
OBJ_DIR= ./obj
INC_DIR= ./inc

OPT = 0
ifeq ($(OPT),1)
CFLAGS=-O3
else
CFLAGS=-g -O0 -lm -std=gnu11 -fopenmp
endif

MIC = ${INC_DIR}/micro.h \
      ${INC_DIR}/micro_struct.h \
      ${INC_DIR}/comm.h \
      ${INC_DIR}/material.h \
      ${INC_DIR}/myio.h \
      ${INC_DIR}/vtk.h \
      ${INC_DIR}/geometry.h \
      ${INC_DIR}/trace.h

OBJ  = ${OBJ_DIR}/main.o \
       ${OBJ_DIR}/homogenize.o \
       ${OBJ_DIR}/micro_struct.o \
       ${OBJ_DIR}/init.o \
       ${OBJ_DIR}/finalize.o \
       ${OBJ_DIR}/alloc.o \
       ${OBJ_DIR}/assembly.o \
       ${OBJ_DIR}/output.o \
       ${OBJ_DIR}/comm_line.o \
       ${OBJ_DIR}/util.o \
       ${OBJ_DIR}/ell.o \
       ${OBJ_DIR}/mesh.o \
       $(OBJ_DIR)/material.o \
       $(OBJ_DIR)/vtk.o \
       $(OBJ_DIR)/solvers.o \
       $(OBJ_DIR)/myio.o \
       $(OBJ_DIR)/geometry.o \
       ${OBJ_DIR}/fem.o \
       $(OBJ_DIR)/list.o

.PHONY: clean all

all: micro

LDFLAG = -lgsl -lgslcblas -lm

micro: ${OBJ}
	gcc -o micro $^ ${LDFLAG} ${CFLAGS}

${OBJ_DIR}/%.o: ${SRC_DIR}/%.c ${INC_DIR}
	gcc -c -o $@ $< ${CFLAGS} -I${INC_DIR}

clean:	    
	rm -f obj/* micro
