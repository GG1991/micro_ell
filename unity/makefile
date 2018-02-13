SRC_DIR= ../src
OBJ_DIR= ../obj
HEA_DIR= ../inc

OPT = 0
ifeq ($(OPT),1)
CFLAGS=-O3
else
CFLAGS=-g -O0 -lm -std=gnu11
endif

all: ell_jacobi

ell_jacobi: ell_jacobi.o $(OBJ_DIR)/ell.o
	gcc -o ell_jacobi $^ -lm ${CFLAGS}

ell_jacobi.o: ell_jacobi.c $(HEA_DIR)/ell.h
	gcc -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11 ${CFLAGS}

$(OBJ_DIR)/ell.o: ../src/ell.c $(HEA_DIR)/ell.h
	gcc -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11 ${CFLAGS}

clean:	    
	rm -f $(OBJ_DIR)/* *.o ell_jacobi

run:
	./ell_jacobi
