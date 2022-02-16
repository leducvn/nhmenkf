INCLUDES = -I. -I${HOME}/nhm/include 
MODULE = 
LIBPATHS = -L${HOME}/nhm/lib
LIBPATH = ../../lib
FC = mpifrtpx
#FFLAGS = -f s -Kvisimpact,ocl,ilfunc,preex,array_private,auto,nooptmsg $(INCLUDES) $(MODULE)
FFLAGS = -O3 -w -Kparallel $(INCLUDES) $(MODULE)
F90 = mpifrtpx
#F90FLAGS = -f s -Kvisimpact,ocl,ilfunc,preex,array_private,auto,nooptmsg -Cfpp -DR_KIND=8 $(INCLUDES) $(MODULE)
F90FLAGS = -O3 -w -Kparallel -Cfpp -DSWP_ON -DR_KIND=8 $(INCLUDES) $(MODULE)
CC = mpifccpx
#CFLAGS = -Kvisimpact,ocl,ilfunc,preex,array_private,nooptmsg
CFLAGS = -O3 -w -Kparallel -DC_FORTRAN_LOWER_
LD = mpifrtpx
#LDFLAGS = -Kvisimpact,ocl,ilfunc,preex,array_private,auto,nooptmsg
LDFLAGS = -O3 -w -Kparallel
