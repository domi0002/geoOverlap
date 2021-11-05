
.SUFFIXES: .cxx

default: lib
all: lib

include ../Makefile.inc
include Makefile.dep

SRC_C   = $(wildcard *.c)
DIR_C   = $(notdir $(SRC_C))
OBJS_C  = $(patsubst %.c,%.o,$(DIR_C))

SRC_F   = $(wildcard *.f)
DIR_F   = $(notdir $(SRC_F))
OBJS_F  = $(patsubst %.f,%.o,$(DIR_F))

OBJS     := $(OBJS_C) $(OBJS_F)

lib: lib$(NLIB).a

lib$(NLIB).a: ${OBJS}
	${AR} -ruv $@ $(OBJS)
	${RANLIB} $@

.PHONY: clean

clean:
	@rm -fv lib*.a
	@rm -fv *.o 
