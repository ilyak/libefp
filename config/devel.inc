# C compiler
CC= clang

# Fortran compiler
FC= egfortran

# install prefix
PREFIX= /home/ilya/libefp

# fragment library path
FRAGLIB= $(PREFIX)/libefp/fraglib

# additional link libraries
MYLIBS= -llapack -lblas -lg2c -lgfortran -lmpi

# additional linker flags
MYLDFLAGS= -L/usr/local/lib -L/usr/local/lib/openmpi

# additional C flags
MYCFLAGS= -DWITH_MPI -std=c99 -g -fcolor-diagnostics -Weverything -Wno-unreachable-code -Wno-c11-extensions -Wno-vla -Wno-sign-conversion -Wno-padded -Wno-format-nonliteral -Wno-disabled-macro-expansion -Wno-cast-align -Wno-cast-qual -I/usr/local/include

# additional Fortran flags
MYFFLAGS= -g

# ar command
AR= ar rc

# ranlib command
RANLIB= ranlib

# rm command
RM= rm -f

# tags program
CTAGS= exctags

# tags file
TAGS= tags
