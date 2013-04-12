# C compiler
CC= clang

# Fortran compiler
FC= gfortran47

# lapack libraries
LAPACK= -lalapack -lcblas -lf77blas -latlas -lgfortran

# install prefix
PREFIX= /home/ilya/libefp

# fragment library path
FRAGLIB= $(PREFIX)/libefp/fraglib

# additional C flags
MYCFLAGS= -std=c99 -g -Weverything -Wno-sign-conversion -Wno-padded -Wno-format-nonliteral -Wno-disabled-macro-expansion -Wno-gnu -I/usr/local/include

# additional Fortran flags
MYFFLAGS= -g

# additional linker flags
MYLDFLAGS= -L/usr/local/lib -L/usr/local/lib/gcc47

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
