
SRC_C   = $(wildcard *.c)
DIR_C   = $(notdir $(SRC_C))
OBJ_C  = $(patsubst %.c,%.o,$(DIR_C))
BIN_C  = $(patsubst %.c,%,$(DIR_C))

default: lib $(BIN_C)

lib:
	@(cd ../; make)

include ../Makefile.inc

arms: arms.o ../src/libitsol.a
iluc: iluc.o ../src/libitsol.a
iluk: iluk.o ../src/libitsol.a
ilut: ilut.o ../src/libitsol.a
vbiluk: vbiluk.o ../src/libitsol.a
vbilut: vbilut.o ../src/libitsol.a
solver: solver.o ../src/libitsol.a
simplest: simplest.o ../src/libitsol.a

clean:
	rm -f *.o $(BIN_C)
