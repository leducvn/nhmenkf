NAME = libeigenexa.a
include ../Make3.inc
#MODULE = -module ../../module
#INCLUDES = -I/home/leduc/qc/include -I/home/leduc/opt/mpich2/include
#LIBPATHS = -L/home/leduc/qc/lib -L/home/leduc/opt/mpich2/lib
#INCLUDES = -I. -I/opt/sgi/mpt/mpt-2.04/include -I/home/leduc/qc/include
#LIBPATHS = -L/opt/sgi/mpt/mpt-2.04/lib -L/home/leduc/qc/lib
#FC = /home/leduc/opt/mpich2/bin/mpif90
#FC = ifort
#FFLAGS = -traceback -O3 -w -auto -fpic -convert big_endian -assume byterecl $(INCLUDES) $(MODULE)
#F90 = /home/leduc/opt/mpich2/bin/mpif90
#F90 = ifort
#F90FLAGS = -traceback -O3 -w -auto -fpic -convert big_endian -assume byterecl $(INCLUDES) $(MODULE)
#CC = /home/leduc/opt/mpich2/bin/mpicc
#CC = icc
#CFLAGS = -traceback -O3 -w -fPIC -DC_FORTRAN_LOWER_ $(INCLUDES)
#LIBPATH = ../../lib

OBJS = miscC.o eigen_devel.o comm.o eigen_blacs.o CSTAB.o eigen_libs.o bisect2.o bisect.o eigen_dc.o \
       mx_pdlaed3.o mx_pdlaed2.o mx_pdlaedz.o mx_pdlaed1.o mx_pdlaed0.o \
       my_pdlaed3.o my_pdlaed2.o my_pdlaedz.o my_pdlaed1.o lapack_eigen.o my_pdlaed0.o \
       dlacpy.o my_pdlasrt.o mx_pdstedc.o my_pdsxedc.o dc2.o \
       dc_redist1.o dc_redist1_r.o dc_redist2.o dc_redist2_r.o dlaed6_common.o dlaed6_init.o dcx.o \
       eigen_house.o eigen_prd_t8.o eigen_prd_t7.o eigen_prd_t6_3.o eigen_prd_t5.o eigen_prd_t4x.o eigen_prd_t2.o \
       eigen_t1.o eigen_prd.o eigen_scaling.o eigen_trbak.o trbakwy4_body.o trbakwy4.o \
       eigen_trd_t8.o eigen_trd_t7.o eigen_trd_t6_3.o eigen_trd_t5x.o eigen_trd_t5.o eigen_trd_t4.o eigen_trd_t2.o \
       eigen_trd.o eigen_s.o eigen_sx.o

.SUFFIXES: .o .c .F

.c.o:
	$(CC) $(CFLAGS) -c $<

.F.o:
	$(FC) $(FFLAGS) -c $<

all: $(OBJS)
	rm -f $(LIBPATH)/$(NAME)
	ar q $(LIBPATH)/$(NAME) $(OBJS)

clean:
	rm -f *.o *.mod

