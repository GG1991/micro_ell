SRC_DIR= ../src
OBJ_DIR= ../obj
HEA_DIR= ../inc

OPT = 0
ifeq ($(OPT),1)
CFLAGS=-O3
else
CFLAGS=-g -O0 -lm -std=gnu11
endif

all: unity_ell

unity_ell: unity_ell.o $(OBJ_DIR)/ell.o
	gcc -o unity_ell $^ -lm ${CFLAGS}

unity_ell.o: unity_ell.c $(HEA_DIR)/ell.h
	gcc -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11 ${CFLAGS}

$(OBJ_DIR)/ell.o: ../src/ell.c $(HEA_DIR)/ell.h
	gcc -c -g -o $@ $< -I$(HEA_DIR) -std=gnu11 ${CFLAGS}

clean:	    
	rm -f $(OBJ_DIR)/* *.o unity_ell

run:
	./unity_ell
