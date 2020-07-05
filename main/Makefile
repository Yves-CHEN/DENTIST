# Directory of the target
OUTPUT = libmain.a
# # Compiler
CXX = g++
AR = ar
LIBTOOL = libtool --mode=link cc -static
# Require BOOST lib
# Require EIGEN library

PKG_CPPFLAGS =  -m64 -I./stats/  -I$(BOOST_PATH)  -I$(MKLROOT)/include -I$(EIGEN3_INCLUDE_DIR)  -I. -DEIGEN_NO_DEBUG -DUSEDOUBLE    -DNDEBUG   -fpic  -g -O3   -fopenmp  -std=gnu++11  -Wno-deprecated -DEIGEN_USE_MKL_ALL -DDEBUG_DENTIST
PKG_LIBS =    $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a  $(MKLROOT)/lib/intel64/libmkl_gnu_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a


# export EIGEN3_INCLUDE_DIR="$local/include/eigen3
# EIGEN_PATH

CXXFLAGS = -w -O3 -m64 -fopenmp -I $(EIGEN3_INCLUDE_DIR)  -I $(MKLROOT)/include
HDR += *.h 
SRC =  bfileOperations.cpp \
	storage.cpp \
	calInversion.tests.cpp \
	fastIO.cpp             \
	getGeno.cpp            \
	inversion.cpp          \
	LDinspect.cpp

SUBDIRS = stats

OBJ = $(SRC:.cpp=.o)



all : $(OUTPUT) 

${SUBDIRS} :
	$(MAKE) -C ${SUBDIRS}

$(OUTPUT) : ${SUBDIRS}
	$(AR) -crf $(OUTPUT) $(OBJ)  ${SUBDIRS}/*.o
	ranlib  $(OUTPUT) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(PKG_CPPFLAGS)  -c $*.cpp 
.SUFFIXES : .cpp .c .o $(SUFFIXES)


$(OUTPUT) : $(OBJ)

FORCE:

clean: 
	rm -f *.o
	rm -f *.a
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done



.PHONY: all $(SUBDIRS)
