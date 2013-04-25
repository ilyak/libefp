# C compiler
CC= clang

# Fortran compiler
FC= gfortran46

# lapack libraries
LAPACK= -lalapack -lcblas -lf77blas -latlas -lgfortran -lmpi

# install prefix
PREFIX= /home/ilya/libefp

# fragment library path
FRAGLIB= $(PREFIX)/libefp/fraglib

# additional C flags
MYCFLAGS= -DWITH_MPI -std=c99 -g -Weverything -Wno-vla -Wno-sign-conversion -Wno-padded -Wno-format-nonliteral -Wno-disabled-macro-expansion -Wno-cast-align -I/usr/local/include -I/usr/local/mpi/openmpi/include

# additional Fortran flags
MYFFLAGS= -g

# additional linker flags
MYLDFLAGS= -L/usr/local/lib -L/usr/local/lib/gcc46 -L/usr/local/mpi/openmpi/lib

# ar command
AR= ar rcu

# ranlib command
RANLIB= ranlib

# rm command
RM= rm -f

# tags program
CTAGS= ctags

# tags file
TAGS= tags
