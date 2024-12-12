FC = ifx
TARGET = test
SRCS = utility.F90 geometry.F90 data_structures.F90 algorithms.F90 vtk.F90 tetrahedral.F90 test.F90
OBJS = $(SRCS:.f90=.o)
OBJS := $(OBJS:.F90=.o)
MOD_FILES = $(SRCS:.f90=.mod)
MOD_FILES := $(MOD_FILES:.F90=.mod)

# Base flags
FFLAGS = -fpp -qopenmp -heap-arrays -g -traceback -check bounds -debug
# FFLAGS = -fpp -qopenmp -heap-arrays -Ofast -xHost

# Run clean before building to remove .o and .mod files
all: clean $(TARGET)

# Build the target from object files
$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) -o $@ $(OBJS)
	@echo $(OMP_STATUS)
	@rm -f $(OBJS) $(MOD_FILES)  # Clean up .o and .mod files after compilation

# Compile .f90 and .F90 files into .o
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.F90
	$(FC) $(FFLAGS) -c $< -o $@

# Clean rule to remove all .o, .mod, and the executable
clean:
	rm -f $(OBJS) $(MOD_FILES) $(TARGET)

.PHONY: all clean