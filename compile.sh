#!/bin/bash - 
#===============================================================================
#
#          FILE: compile.sh
# 
#         USAGE: ./compile.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: YOUR NAME (), 
#  ORGANIZATION: 
#       CREATED: 09/11/2018 14:27
#      REVISION:  ---
#===============================================================================

set -o nounset                              # Treat unset variables as an error

export local="/gpfs1/scratch/90days/uqzzhen4/local/.local2/"
export MKLROOT="/gpfs1/scratch/90days/uqzzhen4/local/intel/compilers_and_libraries/linux/mkl"# GCC settings
export PATH="$local/bin:$PATH"
export EIGEN3_INCLUDE_DIR="$local/include/eigen3"
export BOOST_LIB="$local/include"
export CC=gcc-7
export CXX=g++-7export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8


