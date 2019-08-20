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
OBJECTS=alloc.o init_routines.o chempot.o math_routines.o misc_routines.o diag_routines.o kspace.o nc_IO.o IO.o
HEADERS=alloc.h init_routines.h chempot.h math_routines.h misc_routines.h diag_routines.h kspace.h nc_IO.h IO.h


## all: Default target; empty
.PHONY: all
all: help

# #######################################################################################
# DRIVER_HAM4

## driver_ham4: Builds the final executable for driver_ham4
# Linking of the object files into the final executable. Depends on all .o files.
driver_ham4: driver_ham4.o ham4.o $(OBJECTS)
	${CXX} $(LDFLAGS) -o driver_ham4 driver_ham4.o ham4.o $(OBJECTS) $(LDLIBS)

# Creation of the driver_ham4.o object file, which depends on the header files
driver_ham4.o: driver_ham4.cc ham4.h $(HEADERS)
	${CXX} $(CXXFLAGS) $(INC_FLAGS) -c -o driver_ham4.o driver_ham4.cc

# Deletion of the object file and executable file for this driver
.PHONY: driver_ham4_clean
driver_ham4_clean:
	rm -f driver_ham4 driver_ham4.o ham4.o

# #######################################################################################
# MODULE ALLOC

# Creation of the alloc.o object file, which depends on the module's header file 
# and its source file
alloc.o: alloc.cc alloc.h
	${CXX} $(CXXFLAGS) -c alloc.cc -o alloc.o

## ut_alloc: Runs the testing suite for the module alloc
.PHONY: ut_alloc
ut_alloc: alloc_test # Runs the testing suite's executable
	./alloc_test

# Creation of the alloc_test.o object file, which depends on its source file and 
# on the alloc.h header file.
alloc_test.o: alloc_test.cc alloc.h
	${CXX} $(CXXFLAGS) -c alloc_test.cc -o alloc_test.o

# Linking of the alloc_test.o object and alloc.o object files to create 
# the executable file.
alloc_test: alloc_test.o alloc.o
	${CXX} alloc_test.o alloc.o -o alloc_test

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: ut_alloc_clean
ut_alloc_clean:
	rm -f alloc_test.o alloc.o alloc_test

# #######################################################################################
# MODULE INIT_ROUTINES

# Creation of the init_routines.o object file, which depends on the module's header file 
# and its source file
init_routines.o: init_routines.cc init_routines.h
	${CXX} $(CXXFLAGS) -c init_routines.cc -o init_routines.o

## ut_init_routines: Runs the testing suite for the module init_routines
.PHONY: ut_init_routines
ut_init_routines: init_routines_test # Runs the testing suite's executable
	./init_routines_test

# Linking of the init_routines_test.o object and init_routines.o object files to create 
# the executable file.
init_routines_test: init_routines_test.o init_routines.o
	${CXX} init_routines_test.o init_routines.o -o init_routines_test

# Creation of the init_routines_test.o object file, which depends on its source file and 
# on the init_routines.h header file.
init_routines_test.o: init_routines_test.cc init_routines.h
	${CXX} $(CXXFLAGS) -c init_routines_test.cc -o init_routines_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: ut_init_routines_clean
ut_init_routines_clean:
	rm -f init_routines_test.o init_routines.o init_routines_test

# #######################################################################################
# MODULE CHEMPOT

# Creation of the chempot.o object file, which depends on header file and source file as 
# well as math_routines.h
chempot.o: chempot.cc chempot.h math_routines.h
	${CXX} $(CXXFLAGS) -c chempot.cc -o chempot.o

## ut_math_routines: Runs the testing suite for the module math_routines
.PHONY: ut_chempot
ut_chempot: chempot_test # Runs the testing suite's executable
	./chempot_test

# Executable file from linking of chempot_test.o, chempot.o, and math_routines.o
chempot_test: chempot_test.o chempot.o math_routines.o
	${CXX} $(LDFLAGS) chempot_test.o chempot.o math_routines.o -o chempot_test

# Creation of the chempot_test.o object file, which depends on its source file and 
# on the chempot.h header file.
chempot_test.o: chempot_test.cc chempot.h
	${CXX} $(CXXFLAGS) -c chempot_test.cc -o chempot_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: chempot_clean
chempot_clean:
	rm -f chempot_test.o chempot.o math_routines.o chempot_test

# #######################################################################################
# MODULE MATH_ROUTINES

# Creation of the math_routines.o object file, which depends on the module's header file 
# and its source file
math_routines.o: math_routines.cc math_routines.h
	${CXX} $(CXXFLAGS) -c math_routines.cc -o math_routines.o

## ut_math_routines: Runs the testing suite for the module math_routines
.PHONY: ut_math_routines
ut_math_routines: math_routines_test # Runs the testing suite's executable
	./math_routines_test

# Linking of the math_routines_test.o object and math_routines.o object files to create 
# the executable file.
math_routines_test: math_routines_test.o math_routines.o
	${CXX} math_routines_test.o math_routines.o -o math_routines_test

# Creation of the math_routines_test.o object file, which depends on its source file and 
# on the math_routines.h header file.
math_routines_test.o: math_routines_test.cc math_routines.h
	${CXX} $(CXXFLAGS) -c math_routines_test.cc -o math_routines_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: ut_math_routines_clean
ut_math_routines_clean:
	rm -f math_routines_test.o math_routines.o math_routines_test

# #######################################################################################
# MODULE MISC_ROUTINES

# Creation of the misc_routines.o object file, which depends on the module's header file 
# and its source file
misc_routines.o: misc_routines.cc misc_routines.h
	${CXX} $(CXXFLAGS) -c misc_routines.cc -o misc_routines.o

## ut_misc_routines: Runs the testing suite for the module misc_routines
.PHONY: ut_misc_routines
ut_misc_routines: misc_routines_test # Runs the testing suite's executable
	./misc_routines_test

# Linking of the misc_routines_test.o object and misc_routines.o object files to create 
# the executable file.
misc_routines_test: misc_routines_test.o misc_routines.o
	${CXX} misc_routines_test.o misc_routines.o -o misc_routines_test

# Creation of the misc_routines_test.o object file, which depends on its source file and 
# on the misc_routines.h header file.
misc_routines_test.o: misc_routines_test.cc misc_routines.h
	${CXX} $(CXXFLAGS) -c misc_routines_test.cc -o misc_routines_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: misc_routines_clean
misc_routines_clean:
	rm -f misc_routines_test.o misc_routines.o misc_routines_test

# #######################################################################################
# MODULE DIAG_ROUTINES

# diag_routines.o object file depends on header file and source file
# Must include LAPACKE_INC in the path search
diag_routines.o: diag_routines.cc diag_routines.h
	${CXX} $(CXXFLAGS) $(INC_FLAGS) -c diag_routines.cc -o diag_routines.o

## ut_diag_routines: Runs the testing suite for the module diag_routines
.PHONY: ut_diag_routines
ut_diag_routines: diag_routines_test # Runs the testing suite's executable
	./diag_routines_test

# Testing suite executable depends on diag_routines_test.o diag_routines.o
# Must include LAPACKE_LIB in path search and lapacke library
diag_routines_test: diag_routines_test.o diag_routines.o
	${CXX} $(LDFLAGS) -o diag_routines_test diag_routines_test.o diag_routines.o $(LDLIBS)

# diag_routines_test.o object file depends on source file and diag_routines.h header file
diag_routines_test.o: diag_routines_test.cc diag_routines.h
	${CXX} $(CXXFLAGS) -c diag_routines_test.cc -o diag_routines_test.o

# Clean target for this unit test
.PHONY: ut_diag_routines_clean
ut_diag_routines_clean:
	rm -f diag_routines_test.o diag_routines.o diag_routines_test

# #######################################################################################
# MODULE KSPACE

# kspace.o object file depends on header file, source file, and all included header files
kspace.o: kspace.cc kspace.h init_routines.h alloc.h
	${CXX} $(CXXFLAGS) -c kspace.cc -o kspace.o

## ut_kspace: Runs the testing suite for the module kspace
.PHONY: ut_kspace
ut_kspace: kspace_test # Runs the testing suite's executable
	./kspace_test

# Testing suite executable depends on kspace_test.o, kspace.o, and the other modules used 
# in the source code.
kspace_test: kspace_test.o kspace.o init_routines.o alloc.o
	${CXX} kspace_test.o kspace.o init_routines.o alloc.o -o kspace_test

# kspace_test.o object file depends on source file and kspace.h header
kspace_test.o: kspace_test.cc kspace.h
	${CXX} $(CXXFLAGS) -c kspace_test.cc -o kspace_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: kspace_clean
kspace_clean:
	rm -f init_routines.o alloc.o kspace_test.o kspace.o kspace_test

# #######################################################################################
# MODULE NC_IO

# nc_IO.o object file depends on header file, source file, and all included header files
nc_IO.o: nc_IO.cc nc_IO.h
	${CXX} $(CXXFLAGS) $(INC_FLAGS) -c nc_IO.cc -o nc_IO.o

## ut_nc_IO: Runs the testing suite for the module nc_IO
.PHONY: ut_nc_IO
ut_nc_IO: nc_IO_test # Runs the testing suite's executable
	./nc_IO_test

# Testing suite executable depends on nc_IO_test.o, Z.o, and the other modules used 
# in the source code.
nc_IO_test: nc_IO_test.o nc_IO.o
	${CXX} $(LDFLAGS) -o nc_IO_test nc_IO_test.o nc_IO.o $(LDLIBS)

# nc_IO_test.o object file depends on source file and nc_IO.h header
nc_IO_test.o: nc_IO_test.cc nc_IO.h
	${CXX} $(CXXFLAGS) -c nc_IO_test.cc -o nc_IO_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: ut_nc_IO_clean
ut_nc_IO_clean:
	rm -f nc_IO_test.o nc_IO.o nc_IO_test

# #######################################################################################
# MODULE IO

# IO.o object file depends on header file, source file, and all included header files
IO.o: IO.cc IO.h
	${CXX} $(CXXFLAGS) -c IO.cc -o IO.o

## ut_IO: Runs the testing suite for the module IO
.PHONY: ut_IO
ut_IO: IO_test # Runs the testing suite's executable
	./IO_test

# Testing suite executable depends on IO_test.o and IO.o
IO_test: IO_test.o IO.o alloc.o
	${CXX} IO_test.o IO.o alloc.o -o IO_test

# IO_test.o object file depends on source file and IO.h header
IO_test.o: IO_test.cc IO.h alloc.h
	${CXX} $(CXXFLAGS) -c IO_test.cc -o IO_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: IO_clean
IO_clean:
	rm -f IO_test.o IO.o alloc.o IO_test

# #######################################################################################
# MODULE HAM4

# ham4.o object file depends on header file, source file, and all included header 
# files
ham4.o: ham4.cc ham4.h alloc.h init_routines.h math_routines.h misc_routines.h chempot.h kspace.h diag_routines.h
	${CXX} $(CXXFLAGS) $(INC_FLAGS) -c ham4.cc -o ham4.o

## ut_ham4: Runs the testing suite for the module ham4
.PHONY: ut_ham4
ut_ham4: ham4_test # Runs the testing suite's executable
	./ham4_test

# Testing suite executable depends on ham4_test.o, ham4.o, and the other 
# modules used in the source code.
ham4_test: ham4_test.o ham4.o IO.o alloc.o init_routines.o math_routines.o misc_routines.o chempot.o kspace.o diag_routines.o
	${CXX} $(LDFLAGS) -o ham4_test ham4_test.o IO.o ham4.o alloc.o init_routines.o math_routines.o misc_routines.o chempot.o kspace.o diag_routines.o $(LDLIBS)

# ham4_test.o object file depends on source file and IO.h header
ham4_test.o: ham4_test.cc ham4.h alloc.h IO.h
	${CXX} $(CXXFLAGS) -c ham4_test.cc -o ham4_test.o

# Deletion of the object files and executable files pertaining to this unit test.
.PHONY: ham4_clean
ham4_clean:
	rm -f ham4_test.o ham4.o IO.o alloc.o init_routines.o math_routines.o misc_routines.o chempot.o kspace.o diag_routines.o ham4_test

# #######################################################################################

## clean: Removes module object files as well as driver object files and executables
.PHONY: clean
clean: driver_ham4_clean
	rm -f $(OBJECTS)

## ut_clean: Runs clean rules for all unit tests
.PHONY: ut_clean
ut_clean: ut_alloc_clean \
          ut_init_routines_clean \
          chempot_clean \
          ut_math_routines_clean \
          misc_routines_clean \
          ut_diag_routines_clean \
          kspace_clean \
          ut_nc_IO_clean \
          IO_clean \
          ham4_clean

## help: Shows targets and their descriptions
.PHONY: help
help: Makefile
	@sed -n 's/^##//p' $<