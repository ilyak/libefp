# C compiler
CC= mpicc

# Fortran compiler
FC= gfortran

# install prefix
PREFIX= /usr/local

# fragment library path
FRAGLIB= $(PREFIX)/share/libefp/fraglib

# additional link libraries
MYLIBS= -lalapack -lf77blas -latlas -lgfortran

# additional linker flags
MYLDFLAGS=

# additional C flags
MYCFLAGS= -DEFP_USE_MPI -std=c99 -O3 -g -fopenmp

# additional Fortran flags
MYFFLAGS= -g

# ar command
AR= ar rc

# ranlib command
RANLIB= ranlib

# rm command
RM= rm -f
