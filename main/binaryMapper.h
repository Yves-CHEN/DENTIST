#ifndef BINARY_MAPPER_H
#define BINARY_MAPPER_H

// #include "dataTypeDef.h"
#include<cmath>
#include <bitset>
typedef unsigned char uchar;
typedef unsigned int uint;
typedef unsigned short dataType;

template <class T>
class BinaryMapper
{
    //***********************************************************************************
    /// Creating map for a bype.
    /// number of ones,  00 coded for 0 in  a additive model, 11 for 2, 10 or 01 for 1
    /// This step assumes no missingness. (10 for missing.)
    ///    mapper : 65536 = 4 ^ 1byte
    // invCountOnes is a "function (aa)  aa = ~aa; bitset<8> (aa).count"
    // mapper2   is a "function (aa) aa = ~aa; aa & (aa<<1) & 0xaa" 
    // mapper3   is a "function (aa)  aa & (aa<<1) & 0xaa" 
    // countOnes is a "function (aa)  bitset<8> (aa).count"
    // mapper5   is a "function (aa)  aa = ~aa ; ( ((aa ^ aa <<1) & 0xaa ) >>1 ) *3"
    //
    //***********************************************************************************
    
    size_t sizeOfMap = size_t (pow( 4 , (sizeof(T) * 8) ));
    uint *  invCountOnes ;
    uint *  mapper2      ; 
    uint *  mapper3      ; 
    uint *  countOnes    ; 
    uchar*  mapper5      ; 
    uchar*  markMissing  ; 

public:
    ~BinaryMapper <T>()
    {
        delete[] invCountOnes ;
        delete[] mapper2      ;
        delete[] mapper3      ;
        delete[] countOnes    ;
        delete[] mapper5      ;
        delete[] markMissing  ;
    };
    BinaryMapper <T>()
    {
        sizeOfMap = (size_t) (pow( 2 , (sizeof(T) * 8) ));

        std::cout << sizeOfMap << std::endl;
        invCountOnes = new uint   [sizeOfMap ] ;
        mapper2      = new uint   [sizeOfMap ]; 
        mapper3      = new uint   [sizeOfMap ]; 
        countOnes    = new uint   [sizeOfMap ]; 
        mapper5      = new uchar  [sizeOfMap ]; 
        markMissing  = new uchar  [sizeOfMap ]; 
#pragma omp parallel for 
        for (unsigned int i = 0; i < sizeOfMap ; i ++)
            countOnes[i] = (std::bitset<sizeof(T)*8> (i)).count();
        T screener = T (0xaaaaaaaaaaaaaaaa); /// the 64bit value automatically truncated to the T size.
#pragma omp parallel for 
        for (unsigned int i = 0; i < sizeOfMap ; i ++)
        {
            T aa = (T)(i);
            invCountOnes [i] = countOnes[aa];
            mapper2[i]       = countOnes[(T)(aa & (aa << 1 ) & screener )];
            mapper3[i]       = countOnes[(T)((T)(i) & ((T)(i) << 1 ) & screener) ];
            markMissing[i]   =  ((aa ^ aa <<1) & screener & ~aa ) ;
            markMissing[i]  = ~( markMissing[i] | (markMissing[i] >> 1) );
        }
        std::cout << "Running" << std::endl;
    };
public:
    inline void print(){};


};


//class TT
//{
//public:
//    static  BinnaryMapper<short>  mapper;;
//};
//

//BinnaryMapper<short>  TT::mapper;







#endif

