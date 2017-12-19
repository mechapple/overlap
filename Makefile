# Load the common configuration file
include /users/pala2/lore/voro++-0.4.6/config.mk

# List of executables
EXECUTABLES= radical_hmx create_assembly_stl gccm_single gccm_generate get_intersection_mesh read_stl

# Makefile rules -L/users/pala2/lore/gmp/lib -lgmp
all: $(EXECUTABLES)

radical_hmx: radical_hmx.cc
	$(CXX) $(CFLAGS) $(E_INC) $(E_LIB) -o radical_hmx radical_hmx.cc -lvoro++

create_assembly_stl: create_assembly_stl.cpp
	g++ -std=c++11 -std=c++0x -L/users/pala2/lore/boost_1_59_0 -L/users/pala2/lore/boost_1_59_0/stage/lib -g -O3 -o create_assembly_stl create_assembly_stl.cpp -lboost_timer -lboost_system -lboost_filesystem -lm

read_stl: read_stl.cpp
	g++ -std=c++11 -std=c++0x -g -O3 -o read_stl read_stl.cpp -lm

gccm_single: GCCM_SingleCrystal1.f90
	gfortran -o gccm_single GCCM_SingleCrystal1.f90 
	
gccm_generate: GCCM_GenerateCrystal1.f90
	gfortran -o gccm_generate GCCM_GenerateCrystal1.f90

get_intersection_mesh: get_intersection_mesh.cpp
	g++ -o get_intersection_mesh -frounding-math -std=c++11 get_intersection_mesh.cpp \
		-I/users/pala2/lore/cgal/include -I/users/pala2/lore/libigl/include -I/users/pala2/lore/eigen \
		-L/users/pala2/lore/boost_1_59_0/stage/lib -lboost_thread -lboost_system \
		-L/users/pala2/lore/cgal/lib -lCGAL \
		-L/users/pala2/lore/mpfr/lib -lmpfr \
		-L/users/pala2/lore/gmp/lib -lgmp 

clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
