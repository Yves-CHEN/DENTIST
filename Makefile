# Directory of the target
dirs = builts


# Name of the executable.
OUTPUT = ./builts/DENTIST.tmp2

# Compiler
#CXX =  /gpfs1/scratch/90days/uqzzhen4/local/.local/bin/g++-7
#CXX = ~/90days/utils/bin/g++
CXX = g++

#mklRoot = /gpfs1/scratch/90days/uqzzhen4/local/intel/compilers_and_libraries/linux/mkl
#mklRoot = /gpfs1/scratch/90days/uqzzhen4/local/intel/compilers_and_libraries/linux/mkl

main = main

# EIGEN library
#EIGEN_PATH = $(EIGEN3_INCLUDE_DIR)

#BOOST_PATH = /home/uqwche11/90days/tool-sources/boost_1_69_0

# Intel MKL library
#MKL_PATH = /opt/intel/mkl

# Compiler flags
#CXXFLAGS = -w -O3 -m64 -fopenmp -I $(EIGEN_PATH) -DDEBUG -g 
#-I $(EIGEN_PATH)
CXXFLAGS = -w -O3 -fopenmp  -DEIGEN_NO_DEBUG 

LIB +=  -lz -Wl,-lm -ldl  
#LIB += -lz -Wl, -lm -ldl

# PKG_CPPFLAGS =  -m64 -I${MKLROOT}/include -I${EIGEN3_INCLUDE_DIR} -I/home/uqwche11/utils/R-3.2.2/lib64/R/include  -I. -DUSEDOUBLE -g3 -fopenmp  -std=gnu++11  -Wno-deprecated
# PKG_LIBS =  -Llib -lLDinspect -DUSEDOUBLE -g3 -m64  -Wl,--start-group /opt/intel/composer_xe_2017.4/mkl/lib/intel64/libmkl_intel_lp64.a /opt/intel/composer_xe_2017.4/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/composer_xe_2017.4/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl -L /home/uqwche11/utils/R-3.2.2/lib64/R/lib/ -lR -lRblas -lRlapack


PKG_CPPFLAGS =  -m64 -DEIGEN_NO_DEBUG -DNDEBUG   -fpic  -g -O2     -I${BOOST_PATH}   -I${MKLROOT}/include  -DMKL_ILP64  -I${EIGEN3_INCLUDE_DIR}  -I. -DUSEDOUBLE -g3 -fopenmp  -std=gnu++11  -Wno-deprecated -DEIGEN_USE_MKL_ALL
PKG_LIBS = -static -L$(main) -l$(main) -DUSEDOUBLE -g3 -m64 -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group  -lgomp -lpthread -lz -lm -ldl  -DNDEBUG -DEIGEN_USE_MKL_ALL









HDR += DENTIST.h \
	 invoker.h 

SRC = DENTIST.cpp


OBJ = $(SRC:.cpp=.o)

all : $(OUTPUT)

$(dirs): 
	@echo "Creating $@ dirs"
	mkdir -p $@




$(main)/$(main).a :
	$(MAKE) -C ${main}


$(OUTPUT) : | $(dirs)  $(main)/$(main).a
	$(CXX)  -o $(OUTPUT)   $(OBJ) $(PKG_LIBS) $(CXXFLAGS) $(LIB) 

$(OBJ) : $(HDR)

.cpp.o : 
	$(CXX) $(CXXFLAGS) $(PKG_CPPFLAGS) -c $*.cpp
.SUFFIXES : .cpp .c .o $(SUFFIXES)

$(OUTPUT) : $(OBJ)


run:
	export LD_LIBRARY_PATH=/home/uqwche11/utils/R-3.2.2/lib64/R/lib/:$LD_LIBRARY_PATH
	time /home/uqwche11/30days/simulation.UK10K/multipleSites.dup2/DENTIST.new2  --gwas-summary /home/uqwche11/30days/simulation.UK10K//anaHeight/summaryData/height_UKB_FULLSAMPLE_Julia.txt \
	--bfile test.chr22 \
	--out ~/30days/test --thread-num 8  --maf 0.0001 \
	--degree-of-QC 0.03 --lambda 0.1 \
	--wind-dist 2000000 \
	--target   rs3827330 \
	--extract <(echo rs3827330)
run2:
	export LD_LIBRARY_PATH=/home/uqwche11/utils/R-3.2.2/lib64/R/lib/:$LD_LIBRARY_PATH
	time ./DENTIST --gwas-summary /home/uqwche11/30days/simulation.UK10K//anaHeight/summaryData/height_UKB_FULLSAMPLE_Julia.txt \
	--bfile test.chr22 \
	--out ~/30days/test.2 --thread-num 8  --maf 0.0001 \
	--degree-of-QC 0.03 --lambda 0.1 \
	--wind-dist 2000000 \
	--target   rs3827330
run3:
	bash -c "time ./DENTIST  --gwas-summary /home/uqwche11/30days/simulation.UK10K//anaHeight/summaryData/height_UKB_FULLSAMPLE_Julia.txt --bfile test.chr22 --out ~/30days/test --thread-num 8  --maf 0.0001 --degree-of-QC 0.03 --lambda 0.1 --wind-dist 2000000 --target   rs3827330 --extract <(echo rs3827330)"



FORCE:

clean: 
	rm -f *.o
	for dir in $(main); do \
		$(MAKE) -C $$dir clean; \
	done



