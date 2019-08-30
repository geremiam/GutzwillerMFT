# Makefile
# Automates the compilation and building of the MFT solver

# Compiler name
CXX=icpc
# Flags for including in path search
INC_FLAGS=-I${LAPACKE_INC} -I${NETCDF_INC}
# Compiler flags (add -fopenmp to compilation and linking for OpenMP)
CXXFLAGS=-std=c++14
# Linker flags (add -fopenmp to compilation and linking for OpenMP)
LDFLAGS=-L${LAPACKE_LIB} -L${NETCDF_LIB}
# Flags for linking with libraries (place after all object files)
LDLIBS=-llapacke -lnetcdf
# List of object files and header files belonging to modules
OBJECTS=alloc.o array_init.o IO.o kspace.o math.o
HEADERS=alloc.h array_init.h IO.h kspace.h math.h


## all: Default target; empty
.PHONY: all
all: help

# #######################################################################################
# DRIVER_HAM1

# #######################################################################################
# MODULE HAM1

# #######################################################################################
# MODULE ALLOC

# Creation of the alloc.o object file, which depends on the module's header file 
# and its source file
alloc.o: alloc.cc alloc.h
	${CXX} $(CXXFLAGS) -c alloc.cc -o alloc.o

# Creation of the alloc_test.o object file, which depends on its source file and 
# on the alloc.h header file.
alloc_test.o: alloc_test.cc alloc.h
	${CXX} $(CXXFLAGS) -c alloc_test.cc -o alloc_test.o

# Linking of the alloc_test.o object and alloc.o object files to create 
# the executable file.
alloc_test: alloc_test.o alloc.o
	${CXX} alloc_test.o alloc.o -o alloc_test

## ut_alloc: Runs the testing suite for the module alloc
.PHONY: ut_alloc
ut_alloc: alloc_test # Runs the testing suite's executable
	./alloc_test

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: ut_alloc_clean
ut_alloc_clean:
	rm -f alloc_test.o alloc.o alloc_test
# #######################################################################################
# MODULE ARRAY_INIT

# Creation of the object file
array_init.o: array_init.cc array_init.h
	${CXX} $(CXXFLAGS) -c array_init.cc -o array_init.o

# Creation of the alloc_test.o object file
array_init_test.o: array_init_test.cc array_init.h
	${CXX} $(CXXFLAGS) -c array_init_test.cc -o array_init_test.o

# Linking of the alloc_test.o object and alloc.o object files
array_init_test: array_init_test.o array_init.o
	${CXX} array_init_test.o array_init.o -o array_init_test

## ut_array_init: Runs the testing suite for the module array_init
.PHONY: ut_array_init
ut_array_init: array_init_test # Runs the testing suite's executable
	./array_init_test

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: array_init_clean
array_init_clean:
	rm -f array_init_test.o array_init.o array_init_test
# #######################################################################################
# MODULE IO

# IO.o object file depends on header file, source file, and all included header files
IO.o: IO.cc IO.h
	${CXX} $(CXXFLAGS) -c IO.cc -o IO.o

# IO_test.o object file depends on source file and IO.h header
IO_test.o: IO_test.cc IO.h alloc.h
	${CXX} $(CXXFLAGS) -c IO_test.cc -o IO_test.o

# Testing suite executable depends on IO_test.o and IO.o
IO_test: IO_test.o IO.o alloc.o
	${CXX} IO_test.o IO.o alloc.o -o IO_test

## ut_IO: Runs the testing suite for the module IO
.PHONY: ut_IO
ut_IO: IO_test # Runs the testing suite's executable
	./IO_test

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: IO_clean
IO_clean:
	rm -f IO_test.o IO.o alloc.o IO_test
# #######################################################################################
# MODULE KSPACE

# kspace.o object file depends on header file, source file, and all included header files
kspace.o: kspace.cc kspace.h alloc.h array_init.h
	${CXX} $(CXXFLAGS) -c kspace.cc -o kspace.o

# kspace_test.o object file depends on source file and header
kspace_test.o: kspace_test.cc kspace.h
	${CXX} $(CXXFLAGS) -c kspace_test.cc -o kspace_test.o

# Testing suite executable depends on all linked object files
kspace_test: kspace_test.o kspace.o alloc.o array_init.o
	${CXX} kspace_test.o kspace.o alloc.o array_init.o -o kspace_test

## ut_kspace: Runs the testing suite for the module kspace
.PHONY: ut_kspace
ut_kspace: kspace_test # Runs the testing suite's executable
	./kspace_test

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: kspace_clean
kspace_clean:
	rm -f kspace_test.o kspace.o alloc.o array_init.o kspace_test
# #######################################################################################
# MODULE MATH

# Creation of the object file
math.o: math.cc math.h
	${CXX} $(CXXFLAGS) -c math.cc -o math.o

# Creation of the math_test.o object file
math_test.o: math_test.cc math.h
	${CXX} $(CXXFLAGS) -c math_test.cc -o math_test.o

# Linking of the math_test.o object and math.o object files
math_test: math_test.o math.o
	${CXX} math_test.o math.o -o math_test

## ut_math: Runs the testing suite for the module math
.PHONY: ut_math
ut_math: math_test # Runs the testing suite's executable
	./math_test

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: math_clean
math_clean:
	rm -f math_test.o math.o math_test
# #######################################################################################
# #######################################################################################

## clean: Removes module object files as well as driver object files and executables
.PHONY: clean
clean: #driver_ham4_clean
	rm -f $(OBJECTS)

## ut_clean: Runs clean rules for all unit tests
.PHONY: ut_clean
ut_clean: ut_alloc_clean \
          IO_clean \
          array_init_clean\
          kspace_clean\
          math_clean#\
          #ham4_clean

## help: Shows targets and their descriptions
.PHONY: help
help: Makefile
	@sed -n 's/^##//p' $<