# SUNTANSHOME=/Users/fringer/suntans-periodic/main
SUNTANSHOME=../../main
include $(SUNTANSHOME)/Makefile.in

ifneq ($(MPIHOME),)
  CC = $(MPIHOME)/bin/mpicc
  MPIDEF = 
  MPIINC = -I$(MPIHOME)
else
  CC = mpicc
  MPIDEF =
  MPIINC = 
endif

ifneq ($(PARMETISHOME),)
  PARMETISINC = -I$(PARMETISHOME)/ParMETISLib
endif
ifneq ($(NETCDF4HOME),)
  NETCDFINCLUDE=-I$(NETCDF4HOME)/include
  NETCDFDEF = -DUSENETCDF
else
  NETCDFLIBDIR = 
  NETCDFINCLUDE = 
  NETCDFLD = 
  NETCDFSRC=
endif

LD = $(CC) 
CFLAGS = 
MATHLIB = -lm

EXEC = lee_wave
OBJS = 
SUN = $(SUNTANSHOME)/sun
INCLUDES = -I$(SUNTANSHOME) $(MPIINC) $(PARMETISINC) $(NETCDFINCLUDE) 
DEFS = $(MPIDEF) $(NETCDFDEF)
NUMPROCS = 40

all:	data

grid: clobber data
	sh $(EXEC).sh $(NUMPROCS)  2>&1 |  tee output.txt

test: clobber data	
	sh $(EXEC).sh $(NUMPROCS)  2>&1 |  tee output.txt
	rm -f ./data/*.dat.*

data:	$(SUN)

.c.o:	
	$(LD) -c $(INCLUDES) $(DEFS) $*.c

$(SUN):	initialization.o boundaries.o state.o sources.o
	cp initialization.o boundaries.o state.o sources.o $(SUNTANSHOME)
	make -C $(SUNTANSHOME) 2>&1 |  tee output_compile.txt

clean:
	rm -f *.o ../../main/*.o ../../main/sun

clobber: clean
	rm -rf ./*.txt *~ \#*\# PI* $(EXEC) gmon.out data rundata/*~
