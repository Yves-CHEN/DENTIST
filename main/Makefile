# Directory of the target
OUTPUT = ./lib/libLDinspect.a
#
# # Compiler
CXX = g++
AR = ar
LIBTOOL = libtool --mode=link cc -static

BOOST_PATH = /home/uqwche11/90days/tool-sources/boost_1_69_0

# EIGEN library
mklRoot = /gpfs1/scratch/90days/uqzzhen4/local/intel/compilers_and_libraries/linux/mkl

PKG_CPPFLAGS =  -m64 -I./stats/  -I${BOOST_PATH}  -I$(mklRoot)/include -I$(EIGEN3_INCLUDE_DIR)  -I. -DEIGEN_NO_DEBUG -DUSEDOUBLE    -DNDEBUG   -fpic  -g -O3   -fopenmp  -std=gnu++11  -Wno-deprecated -DEIGEN_USE_MKL_ALL -DDEBUG_DENTIST
PKG_LIBS =   lib/libHeader.a $(mklRoot)/lib/intel64/libmkl_intel_lp64.a  $(mklRoot)/lib/intel64/libmkl_gnu_thread.a $(mklRoot)/lib/intel64/libmkl_core.a



CXXFLAGS = -w -O3 -m64 -fopenmp -I $(EIGEN_PATH)  -I $(mklRoot)/include

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
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done
			

run : 
	$(CXX) -o tt  tt.cpp ./lib/libLDinspect.a   -DUSEDOUBLE -g3 -m64  -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl 
	./tt


install:
	#cp lib/libLDinspect.a /home/uqwche11/work/SMR_source/smr0.866/lib/
	cp lib/libLDinspect.a  ../lib/

test:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir; \
	done

.PHONY: all $(SUBDIRS)
