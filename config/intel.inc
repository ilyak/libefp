# C compiler
CC= icc

# Fortran compiler
FC= ifort

# install prefix
PREFIX= /usr/local

# fragment library path
FRAGLIB= $(PREFIX)/share/libefp/fraglib

# additional link libraries
MYLIBS= -lifcore

# additional linker flags
MYLDFLAGS=

# additional C flags
MYCFLAGS= -std=c99 -O3 -g -mkl=sequential -fopenmp

# additional Fortran flags
MYFFLAGS= -g
