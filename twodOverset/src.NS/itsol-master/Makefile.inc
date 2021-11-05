NLIB            = itsol
LIB_SUFFIX      = .a

CPP	            = /usr/bin/gcc -E
CPPFLAGS        = -I/users/3054576/OpenFOAM/3054576-v1812/geoOverlap/twodOverset/src/itsol-master/include  -I/usr/local/include

CC              = /usr/bin/gcc 
CFLAGS          = -g -O2 -Wall -Wextra -Wmissing-declarations -Wmissing-prototypes -Wimplicit-function-declaration -Wunused -Wuninitialized -Wvla -Wformat -Wformat-security -Wmain -Wsequence-point -Wshadow -Wbad-function-cast -Wcomments -Wdeclaration-after-statement -Waddress -Wpointer-sign -Wempty-body -Wnested-externs -Wold-style-declaration -Wold-style-definition -Wmissing-field-initializers -Wmaybe-uninitialized -Wno-return-local-addr -Wreturn-type -Wsequence-point -Wundef -Wunknown-pragmas -Wunused-label

FC              = /usr/bin/gfortran 
FCFLAGS         = -g -O2
FCLIBS          =  -L/opt/gridware/depots/54e7fb3c/el7/pkg/libs/lapack/3.5.0/gcc-4.8.5/lib -L/opt/gridware/depots/54e7fb3c/el7/pkg/libs/blas/3.6.0/gcc-4.8.5 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../.. -lgfortran -lm -lquadmath

LINKER          = /usr/bin/gcc 
LDFLAGS         = -L/users/3054576/OpenFOAM/3054576-v1812/geoOverlap/twodOverset/src/itsol-master/src  -L/opt/gridware/depots/54e7fb3c/el7/pkg/libs/lapack/3.5.0/gcc-4.8.5/lib -L/opt/gridware/depots/54e7fb3c/el7/pkg/libs/blas/3.6.0/gcc-4.8.5 -L/usr/local/lib
LIBS            = -litsol -L/usr/local/lib  -lm  -L/opt/gridware/depots/54e7fb3c/el7/pkg/libs/lapack/3.5.0/gcc-4.8.5/lib -L/opt/gridware/depots/54e7fb3c/el7/pkg/libs/blas/3.6.0/gcc-4.8.5 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.8.5/../../.. -lgfortran -lm -lquadmath -lm

AR              = ar
RANLIB          = ranlib

prefix		= /users/3054576/OpenFOAM/3054576-v1812/geoOverlap/twodOverset/src/itsol-master
libdir      = /users/3054576/OpenFOAM/3054576-v1812/geoOverlap/twodOverset/src/itsol-master/lib
includedir  = /users/3054576/OpenFOAM/3054576-v1812/geoOverlap/twodOverset/src/itsol-master/include
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

