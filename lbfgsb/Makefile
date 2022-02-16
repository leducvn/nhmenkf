NAME = liblbfgsb.a
include ../Make.inc
#MODULE = -module ../../module
#INCLUDES = -I/home/leduc/qc/include -I/home/leduc/opt/mpich2/include
#LIBPATHS = -L/home/leduc/qc/lib -L/home/leduc/opt/mpich2/lib
#INCLUDES = -I. -I/opt/sgi/mpt/mpt-2.04/include -I/home/leduc/qc/include
#LIBPATHS = -L/opt/sgi/mpt/mpt-2.04/lib -L/home/leduc/qc/lib
#FC = /home/leduc/opt/mpich2/bin/mpif90
#FFLAGS = -fp-model precise -fimf-arch-consistency=true
#FC = ifort
#FFLAGS = -traceback -O3 -w -auto -fpic -convert big_endian -assume byterecl $(INCLUDES) $(MODULE)
#F90 = ifort
#F90 = /home/leduc/opt/mpich2/bin/mpif90
#F90FLAGS = -traceback -O3 -w -auto -fpic -convert big_endian -assume byterecl $(INCLUDES) $(MODULE)
#CC = /home/leduc/opt/mpich2/bin/mpicc
#CC = icc
#CFLAGS = -traceback -O3 -w -fPIC -DC_FORTRAN_LOWER_ $(INCLUDES)
#LIBPATH = ../../lib

OBJS = blas.o linpack.o timer.o lbfgsb.o 

.SUFFIXES: .o .c .f .f90

.c.o:
	$(CC) $(CFLAGS) -c $<

.f.o:
	$(FC) $(FFLAGS) -c $<

.f90.o:
	$(F90) $(F90FLAGS) -c $<

all: $(OBJS)
	rm -f $(LIBPATH)/$(NAME)
	ar q $(LIBPATH)/$(NAME) $(OBJS)

clean:
	rm -f *.o *.mod

