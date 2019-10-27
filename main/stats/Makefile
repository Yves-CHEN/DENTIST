# Directory of the target
OUTPUT = lib/libStat.a
#
# # Compiler
CXX = g++
AR = ar
LIBTOOL = libtool --mode=link cc -static

# EIGEN library
mklRoot = /gpfs1/scratch/90days/uqzzhen4/local/intel/compilers_and_libraries/linux/mkl
PKG_CPPFLAGS =  -DEIGEN_NO_DEBUG -DUSEDOUBLE  -DNDEBUG -fpic  -g -O3   -fopenmp  -std=gnu++11  -Wno-deprecated

SRC =  dcdflib.cpp
OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT) 

$(OUTPUT) :  ${OBJ}
	$(AR) -crf $(OUTPUT) $(OBJ) 
	ranlib  $(OUTPUT) 


.cpp.o : 
	$(CXX) $(PKG_CPPFLAGS)   -c $< -o $@ $(HDR)

.SUFFIXES : .cpp .c .o $(SUFFIXES)
FORCE:
clean: 
	rm -f *.o

.PHONY: all
