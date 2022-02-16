NAME = ../../bin/grd2mfin
include ../Make2.inc
#INCLUDES =
#LIBPATHS =
#LIBS =
#CC = icc
#CFLAGS = -O3 -w
#F90 = ifort
#F90FLAGS = -O3 -w -auto -convert big_endian -assume byterecl $(INCLUDES)
#LINK = ifort

OBJS = variable.o grd2mfinlib.o main.o

.SUFFIXES: .o .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

all: $(OBJS)
	$(LD) $(LDFLAGS) -o $(NAME) $(OBJS) $(LIBPATHS) $(LIBS)

clean:
	rm -f *.o *.mod

