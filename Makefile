# Makefile
# Automates the compilation and building of the MFT solver. See documentation on Implicit Rules.

# Compiler name
CXX=icpc# C++ compiler; used in implicit rules.
CC=$(CXX)# C compiler; is used as the linker in implicit rules.
# Compiler flags (add -fopenmp to compilation and linking for OpenMP). Used in implicit compiling rules.
CXXFLAGS=-std=c++14 -DMKL_ILP64 -I${MKLROOT}/include -I${NETCDF_LIB}
# Linker flags (add -fopenmp to compilation and linking for OpenMP). Used in implicit linking rules.
LDFLAGS=-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -L${NETCDF_LIB}
LDLIBS=-lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -lnetcdf

# Compiler name
#CXX=g++-9
# Flags for including in path search
#INC_FLAGS=-I${LAPACKE_INC} -I${NETCDF_INC}
# Compiler flags (add -fopenmp to compilation and linking for OpenMP)
#CXXFLAGS=-std=c++14 $(INC_FLAGS)
# Linker flags (add -fopenmp to compilation and linking for OpenMP)
#LDLIBS=-llapacke -lnetcdf
# Flags for linking with libraries (place after all object files)
#LDFLAGS=-L${LAPACKE_LIB} -L${NETCDF_LIB} $(LDLIBS)

# List of object files and header files belonging to modules
MODULES=ticktock alloc array_init IO kspace math diag ncio
HEADERS=$(patsubst %,%.h,$(MODULES))
OBJECTS=$(patsubst %,%.o,$(MODULES))
TESTS=$(patsubst %,%_test,$(MODULES))
TESTOBJECTS=$(patsubst %,%_test.o,$(MODULES))

# all: Default target; empty
.PHONY: all
all: help

# #######################################################################################
# ham1
# #######################################################################################
ham1.o: ham1.cc ham1.h $(HEADERS)
driver_ham1.o: driver_ham1.cc ham1.h $(HEADERS)
driver_ham1: driver_ham1.o ham1.o $(OBJECTS)

ham1_test.o: ham1_test.cc ham1.h IO.h kspace.h
ham1_test: ham1_test.o ham1.o $(OBJECTS)

.PHONY: clean_ham1
clean_ham1:
	rm -f ham1.o driver_ham1.o driver_ham1 ham1_test.o ham1_test

# #######################################################################################
# Modules
# #######################################################################################
ticktock.o: ticktock.cc ticktock.h
alloc.o: alloc.cc alloc.h
array_init.o: array_init.cc array_init.h
IO.o: IO.cc IO.h
kspace.o: kspace.cc kspace.h alloc.h array_init.h
math.o: math.cc math.h
diag.o: diag.cc diag.h
ncio.o: ncio.cc ncio.h

alloc_test.o: alloc_test.cc alloc.h
array_init_test.o: array_init_test.cc array_init.h
IO_test.o: IO_test.cc IO.h alloc.h
kspace_test.o: kspace_test.cc kspace.h
math_test.o: math_test.cc math.h
diag_test.o: diag_test.cc diag.h
ncio_test.o: ncio_test.cc ncio.h


# #######################################################################################
# Tests
# #######################################################################################
alloc_test: alloc_test.o alloc.o
array_init_test: array_init_test.o array_init.o
IO_test: IO_test.o IO.o alloc.o
kspace_test: kspace_test.o kspace.o alloc.o array_init.o
math_test: math_test.o math.o
diag_test: diag_test.o diag.o
ncio_test: ncio_test.o ncio.o


# #######################################################################################
# Test targets
# #######################################################################################
## ut_alloc
.PHONY: ut_alloc
ut_alloc: alloc_test
	./$<
## ut_array_init
.PHONY: ut_array_init
ut_array_init: array_init_test
	./$<
## ut_IO
.PHONY: ut_IO
ut_IO: IO_test
	./$<
## ut_kspace
.PHONY: ut_kspace
ut_kspace: kspace_test
	./$<
## ut_math
.PHONY: ut_math
ut_math: math_test
	./$<
## ut_diag
.PHONY: ut_diag
ut_diag: diag_test
	./$<
## ut_ncio
.PHONY: ut_ncio
ut_ncio: ncio_test
	./$<
## ut_ham1
.PHONY: ut_ham1
ut_ham1: ham1_test
	./$<


# #######################################################################################
# Clean rules
# #######################################################################################

.PHONY: clean_objectfiles
clean_objectfiles:
	rm -f $(OBJECTS)

.PHONY: clean_unittests
clean_unittests:
	rm -f $(TESTS) $(TESTOBJECTS)

## clean: Removes everything but source files
.PHONY: clean
clean: clean_objectfiles clean_unittests clean_ham1

## help: Shows targets and their descriptions
.PHONY: help
help: Makefile
	@sed -n 's/^##//p' $<


# #######################################################################################
# Notes on implicit rules
# #######################################################################################

# Best guess for the implicit rule compiling .cc into .o (object) files:
# %.o : %.cc (...)
#	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@
# Note that the first prerequisite is the source code.
# Best guess for the implicit rule linking .o files into executables:
# % : (...)
#	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)