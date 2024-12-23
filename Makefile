# Fortran compiler
FC = ifx

# Target executable
TARGET = main

# Source files
SRCS = global.f90 geometry.f90 dsa.f90 vtk.f90 mra.f90 main.f90

# Object files
OBJS = $(SRCS:.f90=.o)

# Module files
MOD_FILES = $(SRCS:.f90=.mod)

# Base flags
FFLAGS = -fpp -heap-arrays -g -traceback -check bounds -debug
#FFLAGS = -fpp -heap-arrays -Ofast -xhost

# Run clean before building to remove .o and .mod files
all: clean $(TARGET)

# Build the target executable
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -qmkl=parallel -qopenmp $(OBJS) -o $(TARGET)

# Compile the source files
# mkl_spblas.o : mkl_spblas.f90
# 	$(FC) $(FFLAGS) -c mkl_spblas.f90 -o mkl_spblas.o

global.o : global.f90
	$(FC) $(FFLAGS) -c global.f90 -o global.o

geometry.o : geometry.f90
	$(FC) $(FFLAGS) -qopenmp -c geometry.f90 -o geometry.o 

dsa.o : dsa.f90
	$(FC) $(FFLAGS) -qopenmp -c dsa.f90 -o dsa.o

vtk.o : vtk.f90
	$(FC) $(FFLAGS) -qopenmp -c vtk.f90 -o vtk.o

mra.o : mra.f90
	$(FC) $(FFLAGS) -qmkl=parallel -qopenmp -c mra.f90 -o mra.o

main.o : main.f90
	$(FC) $(FFLAGS) -c main.f90 -o main.o

clean:
	rm -f *.mod *.o $(TARGET)
.PHONY: all clean