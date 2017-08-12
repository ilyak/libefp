# C compiler
CC= cc

# Fortran compiler
FC= gfortran

# install prefix
PREFIX= /home/ilya/libefp

# fragment library path
FRAGLIB= $(PREFIX)/libefp/fraglib

# additional link libraries
MYLIBS= -llapack -lblas -lg2c -lgfortran

# additional linker flags
MYLDFLAGS= -L/usr/local/lib

# additional C flags
MYCFLAGS= -Wall -Wextra -g -I/usr/local/include

# additional Fortran flags
MYFFLAGS= -g
