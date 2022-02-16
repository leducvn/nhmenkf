INCLUDES = -I.
MODULE = -M ../../module
LIBPATHS = -L${HOME}/enkf/lib
LIBPATH = ../../lib
FC = mpifrtpx
#FFLAGS = -O2 w -Kocl,fast,simd,swp,fp_precision,SVE -SCALAPACK -SSL2BLAMP $(INCLUDES) $(MODULE)
FFLAGS = -O3 w -Kparallel,ocl,ilfunc,preex,array_private,fp_precision,SVE -SCALAPACK -SSL2BLAMP $(INCLUDES) $(MODULE)
F90 = mpifrtpx
#F90FLAGS = -O2 -w -Kocl,fast,simd,swp,fp_precision -Cfpp -DSWP_ON -DR_KIND=8 $(INCLUDES) $(MODULE)
F90FLAGS = -O3 -w -Kparallel,ocl,ilfunc,preex,array_private,fp_precision -Cfpp -DSWP_ON -DR_KIND=8 $(INCLUDES) $(MODULE)
CC = mpifccpx
#CFLAGS = -O2 -w -Kocl,fast,simd,swp,fp_precision -DC_FORTRAN_LOWER_
CFLAGS = -O3 -w -Kparallel,ocl,ilfunc,preex,array_private,fp_precision -DC_FORTRAN_LOWER_
LD = mpifrtpx
#LDFLAGS = -O2 -w -Kocl,fast,simd,swp,fp_precision
LDFLAGS = -O3 -w -Kparallel,ocl,ilfunc,preex,array_private,fp_precision
