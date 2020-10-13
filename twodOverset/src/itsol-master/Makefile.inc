NLIB            = itsol
LIB_SUFFIX      = .a

CPP	            = /usr/bin/gcc -E
CPPFLAGS        = -I/home/chandard/OpenFOAM/chandard-v1812/run/Overset1D/1DMatrix/laplacianFoam/twoDOverset/itsol-master/include  -I/usr/local/include

CC              = /usr/bin/gcc 
CFLAGS          = -g -O2 -Wall -Wextra -Wmissing-declarations -Wmissing-prototypes -Wimplicit-function-declaration -Wunused -Wuninitialized -Wvla -Wformat -Wformat-security -Wmain -Wsequence-point -Wshadow -Wbad-function-cast -Wcomments -Wdeclaration-after-statement -Waddress -Wpointer-sign -Wempty-body -Wnested-externs -Wold-style-declaration -Wold-style-definition -Wmissing-field-initializers -Wmaybe-uninitialized -Wno-return-local-addr -Wreturn-type -Wsequence-point -Wundef -Wunknown-pragmas -Wunused-label

FC              = /usr/bin/gfortran 
FCFLAGS         = -g -O2
FCLIBS          =  -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/7/../../.. -lgfortran -lm -lquadmath

LINKER          = /usr/bin/gcc 
LDFLAGS         = -L/home/chandard/OpenFOAM/chandard-v1812/run/Overset1D/1DMatrix/laplacianFoam/twoDOverset/itsol-master/src  -L/usr/local/lib
LIBS            = -litsol -L/usr/local/lib  -lm  -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/7/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/7/../../.. -lgfortran -lm -lquadmath -lm

AR              = ar
RANLIB          = ranlib

prefix		= /home/chandard/OpenFOAM/chandard-v1812/run/Overset1D/1DMatrix/laplacianFoam/twoDOverset/itsol-master
libdir      = /home/chandard/OpenFOAM/chandard-v1812/run/Overset1D/1DMatrix/laplacianFoam/twoDOverset/itsol-master/lib
includedir  = /home/chandard/OpenFOAM/chandard-v1812/run/Overset1D/1DMatrix/laplacianFoam/twoDOverset/itsol-master/include
RPATH_FLAG	= -Wl,-rpath,

%.o : %.c
	${CC} ${CFLAGS} ${CPPFLAGS} -c $*.c
	@echo

%.o : %.f
	${FC} ${FCFLAGS} -c $*.f
	@echo

% : %.o
	${LINKER} ${LDFLAGS} -o $@ $*.o ${LIBS}
	@echo

% : %.c
	${LINKER} ${LDFLAGS} ${CFLAGS} ${CPPFLAGS} -c $*.c ${LIBS}
	@echo

