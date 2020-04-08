#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <math.h>
#include <bitset>
#include <list>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>


#include <omp.h>
#include <mkl.h>
//#include <stdio.h>
//#include <string.h>
#include "sortedHeap.h"
#include <thread>
#include <time.h>
#include <assert.h>
//#include "utils.h"
#include <numeric>
#include <algorithm>
#include <Eigen/Dense> 
#include "storage.h"


// This defines the processing of bfile binnary in single bype / double byte.
typedef unsigned char dataType;
//typedef unsigned short dataType;

typedef unsigned int  uint;
typedef unsigned long long     uint64;
typedef          long long      int64;
typedef unsigned char uchar;

//  #ifndef CONST_VARS
//  #define CONST_VARS
//  const float lamda = 0.001; // To avoid inversion fails. 
//  #endif

#ifdef DEBUG_DENTIST
#  define D(x) x
#else
#  define D(x)
#endif // DEBUG

