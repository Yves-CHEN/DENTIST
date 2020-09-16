#ifndef __HEADERS__
#define __HEADERS__
#include <algorithm>
#include <assert.h>
#include <cstdio>
#include <cmath>
#include <cstring>
//#include <Eigen/Dense> 
#include <stdarg.h>  
#include <fstream>
#include <iostream>
#include <map>
#include <numeric> 
#include <omp.h>
#include <string> 
#include <vector> 
#include <cstdint>
#include "main/storage.h"


#include "utils.h"

using namespace std;

//typedef  float LDType;
//
typedef unsigned int  uint;
typedef unsigned char uchar;
typedef double        LDType;
typedef long long     int64;
typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper> LDType2;


const uint maxSummaryRowSize = 1000000;


#ifdef DEBUG_DENTIST
#  define D(x) x
#else
#  define D(x)
#endif // DEBUG


#endif
