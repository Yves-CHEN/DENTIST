# Directory of the target
OUTPUT = libStat.a
#
# Compiler
CXX = g++
AR = ar
LIBTOOL = libtool --mode=link cc -static
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
	rm -f *.a

.PHONY: all
