# Define the compiler
FC = ifx

# Define the flags for the compiler for both release and debug builds
RELEASE_FLAGS = -Ofast
DEBUG_FLAGS = -g -O0 -check bounds -check uninit -traceback -real-size 128 -fimf-precision=high -fp-model=precise
# Define the source files
SRC = error_handling.f90 data_structures.f90 mesh_types.f90 geometry.f90 delaunay2.f90 io.f90 meshf.f90

# Define the object files for both release and debug builds
RELEASE_OBJ = $(SRC:.f90=.o)
DEBUG_OBJ = $(SRC:.f90=.debug.o)

# Define the output executable for both release and debug builds
RELEASE_OUT = meshf
DEBUG_OUT = meshf_debug

# Default target
all: $(RELEASE_OUT) post_build

# Debug target
debug: FFLAGS = $(DEBUG_FLAGS)
debug: OBJ = $(DEBUG_OBJ)
debug: OUT = $(DEBUG_OUT)
debug: $(DEBUG_OUT) post_build

# Rule to link everything and produce the release executable
$(RELEASE_OUT): $(RELEASE_OBJ)
	$(FC) $(RELEASE_FLAGS) -o $@ $^

# Rule to link everything and produce the debug executable
$(DEBUG_OUT): $(DEBUG_OBJ)
	$(FC) $(DEBUG_FLAGS) -o $@ $^

# Rule to compile each source file for release
%.o: %.f90
	$(FC) $(RELEASE_FLAGS) -c $< -o $@

# Rule to compile each source file for debug
%.debug.o: %.f90
	$(FC) $(DEBUG_FLAGS) -c $< -o $@

post_build: 
	make clean

clean:
	rm -f *.mod
	rm -f *.o
	rm -f *.i90
	rm -f *.debug.o
