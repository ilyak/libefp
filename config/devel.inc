# C compiler
CC= cc

# Fortran compiler
FC= gfortran

# install prefix
PREFIX= /home/ilya/libefp

# fragment library path
FRAGLIB= $(PREFIX)/fraglib

# additional link libraries
MYLIBS= -llapack -lblas -lgfortran

# additional linker flags
MYLDFLAGS= -L/usr/local/lib

# additional C flags
MYCFLAGS= -Wall -Wextra -g -I/usr/local/include

# additional Fortran flags
MYFFLAGS= -g
