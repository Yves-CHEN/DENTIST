#include "dataTypeDef.h"

//#include <variant>

#include <bitset>
//#include "binaryMapper.h"

// ****************************************************************************
/// The "theMarkIdx" contains the indices for the SNPs in a sequential order.
//  
//      1. Performs LD calculation for given theMarkIdx. The function load
//          all the specified markers into memory from the bed file.
//          That is N marks between <the first SNP, last SNP>. Therefore 
//          this function assumes all the SNPs are from a same small regions 
//          (ideally winSize <20Mb).
//
//      2. The function can inverted the genotypes for specified SNPs (toAvert),
//         so as to change the coding allele. 
//
//
//      3. Deal with missingness. For a pair of SNPs, it takes the common 
//         individuals with no missing genotypes.
//
//  Todos: 
//     2. Sample selection. Currently it does not allow " --keep samplesIDs.txt "
///// Important: the size LD matrix can exceed uint (1e9)
//
//

template <class T>
int calcLDFromBfile_old (std::string bedFile, uint nSample, long nMarker, uint* theMarkIdx, long sizeOfMarkIdx, uint* toAvert, int cutoff, int* ncpus, T* result)
{
    // **************************************************************
    ///                    set timer
    // **************************************************************
    struct timespec start, finish;
    double elapsed;
    
    // **************************************************************
    ///                    set number of cpus
    // **************************************************************
    int nProcessors = omp_get_max_threads();
    if(*ncpus < nProcessors) nProcessors = *ncpus;
    omp_set_num_threads( nProcessors );
        //omp_set_num_threads( 1);
    D(printf("[info] Calc LD  based on %d cpus \n", nProcessors););
    // **************************************************************
    //                     Variables
    // **************************************************************
    uint perMakerSize = ceil ( (nSample) / 4.0  );
    uint nBlanks   = ( 4 - (nSample) % 4  ) % 4; //
    long lSize =0;
    /// headerCoding is 9 bytes in total for plink 1.9 bedfile format, 3 bytes in older version.
    int const nByteHeader = 9;
    uchar headerCoding[nByteHeader +1 ] = { 0x6c, 0x1b, 0x01, 0xdc, 0x0f, 0xe7, 0x0f, 0x6b, 0x01};
    int const nByteHeader_older = 3;
    std::string formatVersion = "";
    int nThrowAway = 0;
    uchar  headerBuf[nByteHeader ] ;
    uint nKeptSample = nSample;
    std::vector<std::vector<uchar> > GENO      (sizeOfMarkIdx, std::vector<uchar> ( nSample, 0) );
    std::vector<double>         GENO_VAR (sizeOfMarkIdx, -1 );
    std::vector<double>         GENO_Ex  (sizeOfMarkIdx, -1 );
    std::vector<double>         GENO_sum11  (sizeOfMarkIdx, -1 ); // sum of sample with genotype 11
    if(cutoff > sizeOfMarkIdx  ) cutoff = sizeOfMarkIdx;

    // **************************************************************
    //                1. validation of the bfile version
    // 2. validation of the bed file size, given the number of markers and samples
    // **************************************************************
    FILE* bedFileReader = fopen ( bedFile.c_str()  , "rb" );
    if (bedFileReader ==NULL) {fputs ("File not found error",stderr); }
    fseek (bedFileReader, 0 , SEEK_END);
    lSize = ftell (bedFileReader);
    rewind (bedFileReader);
    /// This checks the version of bed file. Examine if there is damage to
    /// bedfile, in which case, the estimated bed file size would be inconsistent with the
    /// acutal size.
    fread (&headerBuf,1, nByteHeader, bedFileReader);
    if(!memcmp(headerCoding, headerBuf, nByteHeader)) {printf("[info]This bed file is plink 1.9 bedfile format. (Newer) \n"); formatVersion="1.9"; nThrowAway = nByteHeader;};
    if(!memcmp(headerCoding, headerBuf, nByteHeader_older)) {printf("[info]This bed file is plink 1.0 bedfile format. (Older)\n"); formatVersion="1.0"; nThrowAway = nByteHeader_older;};
    if(lSize  != long(perMakerSize * nMarker + nThrowAway) )
    {
        printf("[error] The size of bedFile %ld is inconsistenty with the estimated %u basd on %u samples and %d markers. \n", lSize, perMakerSize * nMarker + nThrowAway, perMakerSize, nMarker);
    }
    ////////////////////////////////////////////////////
    /// Creating map for a bype.
    /// number of ones,  00 coded for 0 in  a additive model, 11 for 2, 10 or 01 for 1
    /// This step assumes no missingness. (10 for missing.)
    ///    mapper : 65536 = 4 ^ 1byte
    size_t const sizeOfMap = size_t (pow( 4 , (sizeof(char) * 8) ));
    // invCountOnes is a "function (aa)  aa = ~aa; bitset<8> (aa).count"
    // mapper2   is a "function (aa) aa = ~aa; aa & (aa<<1) & 0xaa" 
    // mapper3   is a "function (aa)  aa & (aa<<1) & 0xaa" 
    // countOnes is a "function (aa)  bitset<8> (aa).count"
    // mapper5   is a "function (aa)  aa = ~aa ; ( ((aa ^ aa <<1) & 0xaa ) >>1 ) *3"
    uint  invCountOnes  [sizeOfMap ] = {0};
    uint  mapper2[sizeOfMap ]        = {0}; 
    uint  mapper3[sizeOfMap ]        = {0}; 
    uint  countOnes[sizeOfMap ]      = {0}; 
    uchar mapper5[sizeOfMap ]        = {0}; 
    uchar markMissing[sizeOfMap ]        = {0}; 
#pragma omp parallel for 
    for (unsigned int i = 0; i < sizeOfMap ; i ++)
        countOnes[i] = std::bitset<8> ((unsigned char)(i)).count();
#pragma omp parallel for 
    for (unsigned int i = 0; i < sizeOfMap ; i ++)
    {
        //unsigned char aa = ~(unsigned char)(i);
        unsigned char aa = (unsigned char)(i);
        invCountOnes [i] = countOnes[aa];
        mapper2[i]       = countOnes[(uchar)(aa & (aa << 1 ) & 0xaa )];
        mapper3[i]       = countOnes[(uchar)((unsigned char)(i) & ((unsigned char)(i) << 1 ) & 0xaa) ];
        //mapper5[i]       = ( ((aa ^ aa <<1) & 0xaa ) >>1 ) *3;
        markMissing[i]   =  ((aa ^ aa <<1) & 0xaa & ~aa ) ;
        markMissing[i]  = ~( markMissing[i] | (markMissing[i] >> 1) );
    }

    // ******************************************************
    ///     1. Calculating the var(SNP_i), E(SNP_i), E(SNP_i^2)
    //      2. Avert the genotype if specified.
    // ******************************************************
    unsigned long loadSize = perMakerSize * sizeof(char) * (theMarkIdx[sizeOfMarkIdx-1] - theMarkIdx[0] +1) ;
    uchar* bufferAllMarkers = new unsigned char [loadSize ];

    D(printf("[info] Buffer size is %d Mb. \n", int(loadSize/1e6)););
    
    fseek (bedFileReader , perMakerSize * sizeof(char) * (theMarkIdx[0]) + nThrowAway, SEEK_SET );
    fread (bufferAllMarkers, 1, loadSize, bedFileReader);
    uchar* bufferMaker = NULL;  // This is pointer to the memory of a particular marker.
    for (unsigned int i = 0; i < sizeOfMarkIdx; i ++)
    {

////        uint nMissing = 0;
////        double sum = 0, sum2 = 0;
////        double  E_x=0, E_x_square=0;
        bufferMaker = perMakerSize * (theMarkIdx[i] - theMarkIdx[0]) + bufferAllMarkers;
////
////        // avert the genotypes. For 00 and 11, it would ~ geno. For 01 and 10, it is more difficult.
////        if(toAvert[i]) {
////            for (unsigned int j =0; j< perMakerSize; j ++) 
////            {
////                bufferMaker[j] = (bufferMaker[j]  | ~mapper5[bufferMaker[j]]) & (mapper5[bufferMaker[j]] | ~bufferMaker[j] ) ;
////            } 
////            uchar trunc = 0xff >> nBlanks *2;
////            bufferMaker[perMakerSize -1] = trunc & bufferMaker[perMakerSize -1];
////        }
////#pragma omp parallel for reduction(+:sum, sum2)
////        for (unsigned int j =0; j< perMakerSize; j ++)
////        {
////            nMissing += countOnes [markMissing[bufferMaker[j]] ]; 
////            sum += invCountOnes [(unsigned int )(bufferMaker[j]  | markMissing[bufferMaker[j]])  ];
////            //sum += invCountOnes [(unsigned int )(bufferMaker[j])  ];
////            sum2 += mapper2[(unsigned int )(bufferMaker[j] )];
////        }
////        sum2        = sum2 - nBlanks;
////        E_x         = double(sum - nBlanks *2  ) / (nKeptSample - nMissing);
////        E_x_square  = double(sum    - nBlanks*2  + sum2*2 )/ (nKeptSample -nMissing);
////        GENO_VAR[i] =  E_x_square - E_x * E_x ;
////        GENO_Ex[i]  =  E_x ;
////        GENO_sum11[i]  =  sum2 ;
////        //cout << "Var: " << GENO_VAR[i] << endl;
////        //cout << "MAF: " << E_x / 2 << endl;
////         printf("nMissing : %d sum: %f\n", nMissing, sum);
////        //copy(bufferMaker, bufferMaker + nSample, GENO[i].begin());
        for (unsigned int j =0; j<perMakerSize; j ++)
            GENO[i][j] =  bufferMaker[j];

    }

    
    clock_gettime(CLOCK_MONOTONIC, &start);

    // *****************************************************************
    //
    //        Calculating cov(SNP_i, SNP_j),  LD(SNP_i, SNP_j)
    //
    //  When there is no missingness no-missingness, 
    //     cov(SNP_i, SNP_j) =  SNP_i,k & SNP_j,k 
    //                            + (SNP_i,k==2) + (SNP_j,k ==2) 
    //                            - (SNP_i,k == 2 & SNP_i,j == 0)
    //    Notably, - padding bits (when nSample %mod% 4 != 0) are set to 0.
    // *****************************************************************
    //
    //





#pragma omp parallel for
    for (unsigned int i = 0; i < sizeOfMarkIdx; i ++)
    {
        for (unsigned int j =i+1; j< sizeOfMarkIdx; j ++)
        {

            int sign = 1;
            if(toAvert[i] != toAvert[j]) 
            {
                sign = -1;
            }
            uint nMissing = 0;
            double   E_i = 0;
            double   E_j = 0;
            double   sum11_i = 0;
            double   sum11_j = 0;
            double   sum_i = 0;
            double   sum_j = 0;

            double   E_i_sq = 0;
            double   E_j_sq = 0;
            double var_i = 0 ;
            double var_j = 0 ;
            if (j > i + cutoff)  break;
            double sum_XY = 0;
            double sum_0011 = 0; // 00 OR 11

            uchar marker = 0  ;
            unsigned int k =0;
            for (; k < perMakerSize; k ++)
            {

                marker = markMissing[GENO [i][k]]  & markMissing[GENO [j][k]]   ;
                nMissing += countOnes[(uchar)(~marker)];
                sum_i    += countOnes[GENO[i][k] & marker];
                sum_j    += countOnes[GENO[j][k] & marker];
                sum11_i  += mapper2[(unsigned int )(GENO[i][k] & marker)];
                sum11_j  += mapper2[(unsigned int )(GENO[j][k] & marker)];

                // GENO [i][k] &=  markMissing[GENO [i][k]];
                // GENO [j][k] &=  markMissing[GENO [j][k]];
                GENO [i][k] &=  marker;
                GENO [j][k] &=  marker;
                sum_XY   +=  countOnes  [ (uchar) (GENO [i][k] & GENO [j][k]) ];
                sum_0011 +=  mapper3[ (uchar) (GENO [i][k] ^ GENO [j][k]) ];
            }

            //marker = markMissing[GENO [i][k]]  | markMissing[GENO [j][k]]   ;
            nMissing = nMissing /2;

        // sum11_i     = sum11_i  + mapper2[(unsigned int )(GENO[i][k] | marker)] - nBlanks;
        // sum11_j     = sum11_j  + mapper2[(unsigned int )(GENO[j][k] | marker)] - nBlanks;
            E_i         = double(sum_i  ) / (nKeptSample - nMissing);
            E_j         = double(sum_j  ) / (nKeptSample - nMissing);
            E_i_sq     = double(sum_i   + sum11_i*2 )/ (nKeptSample -nMissing);
            E_j_sq     = double(sum_j   + sum11_j*2 )/ (nKeptSample -nMissing);
            var_i =  E_i_sq - E_i * E_i ;
            var_j =  E_j_sq - E_j * E_j ;



            //uchar trunc = 0xff >> nBlanks *2;
            //sum_XY   +=  invCountOnes  [ (uchar) (GENO [i][perMakerSize -1] | GENO [j][perMakerSize -1] | ~trunc) ];
            //sum_0011 +=  mapper3[ (uchar) ( (GENO [i][perMakerSize -1] ^ GENO [j][perMakerSize -1]) & trunc) ];
            sum_XY += sum11_i + sum11_j - sum_0011;
            //// printf("nMissing cov : %d sum_0011 : %f \n", nMissing, sum_0011);
            /// double cov_XY = sum_XY / (nSample - nMissing) - GENO_Ex[i] * GENO_Ex[j];
            /// double LD = cov_XY / sqrt ( (GENO_VAR[i] * GENO_VAR[j]) );
            double cov_XY = sum_XY / (nSample - nMissing) - E_i * E_j;
            double LD = cov_XY / sqrt ( (var_i * var_j) ) * sign;

            /// printf("sum_XY %f, sum_11=%f, sum_11=%f, cov: %f var1:  %f var2: %f E_i = %f, Ej= %f\n", sum_XY, sum11_i, sum11_j, cov_XY, var_i, var_j, E_i, E_j);
            if(sizeof(*result) == 1) // char
            {
                char aa =  char(LD*250);
                printf("(%f, %f)\n",LD, (0L |255 & aa)/250.0  );
                result[sizeOfMarkIdx * i + j] = char(LD * 250);
                result[sizeOfMarkIdx * j + i] = char(LD * 250);
            }
            else if(sizeof(*result) == 2) // short
            {
                result[sizeOfMarkIdx * i + j] = short(LD * 1000);
                result[sizeOfMarkIdx * j + i] = short(LD * 1000);
     
            }
            else
            {
                result[sizeOfMarkIdx * i + j] = LD ;
                result[sizeOfMarkIdx * j + i] = LD ;
            }
        }
    }


                                               


    delete[] bufferAllMarkers ;
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    D(printf("[info] Step2:  Time elapsed is %f. \n", elapsed ););
    fclose (bedFileReader);
    return 0;

}




template <class T>
int calcLDFromBfile (std::string bedFile, uint nSample, long nMarker, uint* theMarkIdx, long sizeOfMarkIdx, uint* toAvert, int cutoff, int* ncpus, T* result, int* jump)
{
    // **************************************************************
    ///                    set timer
    // **************************************************************
    struct timespec start, finish;
    double elapsed;
    
    // **************************************************************
    ///                    set number of cpus
    // **************************************************************
    int nProcessors = omp_get_max_threads();
    if(*ncpus < nProcessors) nProcessors = *ncpus;
    omp_set_num_threads( nProcessors );
        //omp_set_num_threads( 1);
    D(printf("[info] Calc LD  based on %d cpus \n", nProcessors););
    // **************************************************************
    //                     Variables
    // **************************************************************
    uint processSamplePerRound = sizeof(dataType)*8 /2;
    //uint perMakerSize = ceil ( (nSample) / 1.0 / processSamplePerRound );
    uint perMakerSizeOrig = ceil ( (nSample) / 1.0 / 4);
    uint perMakerSize = int( (nSample) / 1.0 / processSamplePerRound );
    uint nBlanks   = ( processSamplePerRound - (nSample) % processSamplePerRound  ) % processSamplePerRound; 
    long lSize =0;
    /// headerCoding is 9 bytes in total for plink 1.9 bedfile format, 3 bytes in older version.
    int const nByteHeader = 9;
    uchar headerCoding[nByteHeader +1 ] = { 0x6c, 0x1b, 0x01, 0xdc, 0x0f, 0xe7, 0x0f, 0x6b, 0x01};
    int const nByteHeader_older = 3;
    std::string formatVersion = "";
    int nThrowAway = 0;
    uchar  headerBuf[nByteHeader ] ;
    //uint nKeptSample = nSample;
    uint nKeptSample = perMakerSize * processSamplePerRound;
    //// std::vector<std::vector<dataType> > GENO_old      (sizeOfMarkIdx, std::vector<dataType> ( nSample, 0) );
    std::vector<double>         GENO_VAR (sizeOfMarkIdx, -1 );
    std::vector<double>         GENO_Ex  (sizeOfMarkIdx, -1 );
    std::vector<double>         GENO_sum11  (sizeOfMarkIdx, -1 ); // sum of sample with genotype 11
    if(cutoff > sizeOfMarkIdx  ) cutoff = sizeOfMarkIdx;

    // **************************************************************
    //                1. validation of the bfile version
    // 2. validation of the bed file size, given the number of markers and samples
    // **************************************************************
    FILE* bedFileReader = fopen ( bedFile.c_str()  , "rb" );
    if (bedFileReader ==NULL) {fputs ("File not found error",stderr); }
    fseek (bedFileReader, 0 , SEEK_END);
    lSize = ftell (bedFileReader);
    rewind (bedFileReader);
    /// This checks the version of bed file. Examine if there is damage to
    /// bedfile, in which case, the estimated bed file size would be inconsistent with the
    /// acutal size.
    fread (&headerBuf,1, nByteHeader, bedFileReader);
    if(!memcmp(headerCoding, headerBuf, nByteHeader)) {printf("[info]This bed file is plink 1.9 bedfile format. (Newer) \n"); formatVersion="1.9"; nThrowAway = nByteHeader;};
    if(!memcmp(headerCoding, headerBuf, nByteHeader_older)) {printf("[info]This bed file is plink 1.0 bedfile format. (Older)\n"); formatVersion="1.0"; nThrowAway = nByteHeader_older;};
    if(lSize  != long(perMakerSizeOrig * nMarker + nThrowAway) )
    {
        printf("[error] The size of bedFile %ld is inconsistenty with the estimated %u basd on %u samples and %d markers. \n", lSize, perMakerSizeOrig * nMarker + nThrowAway, perMakerSizeOrig, nMarker);
    }
    ////////////////////////////////////////////////////
    /// Creating map for a bype.
    /// number of ones,  00 coded for 0 in  a additive model, 11 for 2, 10 or 01 for 1
    /// This step assumes no missingness. (10 for missing.)
    ///    mapper : 65536 = 2 ^ (2*8 bits)
    size_t const sizeOfMap = (size_t) (pow( 2 , (sizeof(dataType) * 8) ));
    // invCountOnes is a "function (aa)  aa = ~aa; bitset<8> (aa).count"
    // mapper2   is a "function (aa) aa = ~aa; aa & (aa<<1) & 0xaa" 
    // mapper3   is a "function (aa)  aa & (aa<<1) & 0xaa" 
    // countOnes is a "function (aa)  bitset<8> (aa).count"
    // mapper5   is a "function (aa)  aa = ~aa ; ( ((aa ^ aa <<1) & 0xaa ) >>1 ) *3"
    uint      invCountOnes  [sizeOfMap ] = {0};
    uint      mapper2[sizeOfMap ]        = {0}; 
    uint      mapper3[sizeOfMap ]        = {0}; 
    uint      countOnes[sizeOfMap ]      = {0}; 
    dataType  mapper5[sizeOfMap ]        = {0}; 
    dataType  markMissing[sizeOfMap ]        = {0}; 
#pragma omp parallel for 
    for (unsigned int i = 0; i < sizeOfMap ; i ++)
        countOnes[i] = (std::bitset<sizeof(dataType)*8> (i)).count();
    dataType screener = dataType (0xaaaaaaaaaaaaaaaa); /// the 64bit value automatically truncated to the dataType size.
#pragma omp parallel for 
    for (unsigned int i = 0; i < sizeOfMap ; i ++)
    {
        dataType aa = (dataType)(i);
        invCountOnes [i] = countOnes[aa];
        mapper2[i]       = countOnes[(dataType)(aa & (aa << 1 ) & screener )];
        mapper3[i]       = countOnes[(dataType)((dataType)(i) & ((dataType)(i) << 1 ) & screener) ];
        markMissing[i]   =  ((aa ^ aa <<1) & screener & ~aa ) ;
        markMissing[i]  = ~( markMissing[i] | (markMissing[i] >> 1) );
    }


    // ******************************************************
    ///     1. Calculating the var(SNP_i), E(SNP_i), E(SNP_i^2)
    //      2. Avert the genotype if specified.
    // ******************************************************
    unsigned long loadSize = perMakerSizeOrig * sizeof(uchar) * (theMarkIdx[sizeOfMarkIdx-1] - theMarkIdx[0] +1) ;
    uchar* bufferAllMarkers = new uchar [loadSize ];

    printf("[info] Buffer size is %d Mb. \n", int(loadSize/1e6));
    
    fseek (bedFileReader , perMakerSizeOrig * sizeof(uchar) * (theMarkIdx[0]) + nThrowAway, SEEK_SET );
    fread (bufferAllMarkers, 1, loadSize, bedFileReader);
    dataType* bufferMaker = NULL;  // This is pointer to the memory of a particular marker.
    
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    
 
    // *****************************************************************
    //
    //        Calculating cov(SNP_i, SNP_j),  LD(SNP_i, SNP_j)
    //
    //  When there is no missingness no-missingness, 
    //     cov(SNP_i, SNP_j) =  SNP_i,k & SNP_j,k 
    //                            + (SNP_i,k==2) + (SNP_j,k ==2) 
    //                            - (SNP_i,k == 2 & SNP_i,j == 0)
    //    Notably, - padding bits (when nSample %mod% 4 != 0) are set to 0.
    // *****************************************************************
    //
    //


    dataType** GENO = new dataType* [sizeOfMarkIdx] ;
 #pragma omp parallel for
    for (unsigned long long i = 0; i < sizeOfMarkIdx; i ++)
    {
        GENO[i] = (dataType*) (perMakerSizeOrig * (theMarkIdx[i] - theMarkIdx[0]) + bufferAllMarkers);
    }




    std::cout << "jump" << *jump << std::endl;

 #pragma omp parallel for
    for (unsigned long long i = 0; i < sizeOfMarkIdx; i ++)
    {
        dataType* GENO_i = GENO[i] ;


        int starting = (*jump) > (i+1)? (*jump): (i +1);
        //for (unsigned int j =i+1; j< sizeOfMarkIdx; j ++)
        for (unsigned long long j = starting; j< sizeOfMarkIdx; j ++)
        {
            int sign = 1;
            if(toAvert[i] != toAvert[j]) 
                sign = -1;
            uint nMissing = 0;
            double   E_i = 0;
            double   E_j = 0;
            double   sum11_i = 0;
            double   sum11_j = 0;
            double   sum_i = 0;
            double   sum_j = 0;
            double   E_i_sq = 0;
            double   E_j_sq = 0;
            double   var_i = 0 ;
            double   var_j = 0 ;
            long double   sum_XY = 0;
            double   sum_0011 = 0; // 00 OR 11
            dataType marker = 0  ;

            dataType* GENO_j = GENO[j] ;

            if (j > i + cutoff)  break;
            for (unsigned long long k =0; k < perMakerSize; k ++)
            {

                marker = markMissing[ GENO_i[k]  ]  & markMissing[GENO_j[k]]   ;
                nMissing += countOnes[(dataType)(~marker)];
                sum_i    += countOnes[GENO_i [k] & marker];
                sum_j    += countOnes[GENO_j [k] & marker];
                sum11_i  += mapper2[(dataType)(GENO_i [k] & marker)];
                sum11_j  += mapper2[(dataType)(GENO_j [k] & marker)];
                // dataType aa = GENO_i[k] &  marker;
                // dataType bb = GENO_j[k] &  marker;
                sum_XY   +=  countOnes  [ (dataType) (GENO_i[k]  & GENO_j[k])  &  marker ];
                sum_0011 +=  mapper3[ (dataType) (GENO_i[k]  ^ GENO_j[k])  &  marker ];
            }

            nMissing = nMissing /2;
            E_i         = double(sum_i  ) / (nKeptSample - nMissing);
            E_j         = double(sum_j  ) / (nKeptSample - nMissing);
            E_i_sq     = double(sum_i   + sum11_i*2 )/ (nKeptSample -nMissing);
            E_j_sq     = double(sum_j   + sum11_j*2 )/ (nKeptSample -nMissing);
            var_i =  E_i_sq - E_i * E_i ;
            var_j =  E_j_sq - E_j * E_j ;



            sum_XY += sum11_i + sum11_j - sum_0011;
            //double cov_XY = (sum_XY - sum_i * sum_j) / (nSample - nMissing);
            double cov_XY = sum_XY / (nKeptSample - nMissing) - E_i * E_j;
            //// double LD = cov_XY / sqrt ( (var_i * var_j) ) * sign;
            /// This is for rare cases when SNP i or j has 0 variance after removing the NAs genotype individuals.
            double LD = 0.001;
            if(!(var_i <= 0 || var_j <= 0 || var_i * var_j <= 0))
                LD = cov_XY / sqrt ( (var_i * var_j) ) * sign;
            else
            {
                printf("[error] LD was set to 0, because var_i = %f, var_j  %f nMissing = %d E_i_sq (%f) - E_i * E_i\n", var_i, var_j, nMissing, E_i_sq);
                printf( "[error] , when parse the markerIdx %d-%d.\n", theMarkIdx[i],  theMarkIdx[j]);
                printf("[error] %d - %d \n", theMarkIdx[sizeOfMarkIdx-1] , theMarkIdx[0] );
                exit(-1);
            }

            if(sizeof(*result) == 1) // char
            {
                char aa =  char(LD*250);
                printf("(%f, %f)\n",LD, (0L |255 & aa)/250.0  );
                result[sizeOfMarkIdx * i + j] = char(LD * 250);
                result[sizeOfMarkIdx * j + i] = char(LD * 250);
            }
            else if(sizeof(*result) == 2) // short
            {
                result[sizeOfMarkIdx * i + j] = short(LD * 1000);
                result[sizeOfMarkIdx * j + i] = short(LD * 1000);
            }
            else
            {
                result[sizeOfMarkIdx * i + j] = LD ;
                result[sizeOfMarkIdx * j + i] = LD ;
            }

        }
    }



#pragma omp parallel for                                       
    for (unsigned int i = 0; i < sizeOfMarkIdx; i ++)          
    {                                                          
         if(sizeof(*result) == 1) // char
             result[sizeOfMarkIdx * i + i] = char(1 * 250);
         else if(sizeof(*result) == 2) // short
             result[sizeOfMarkIdx * i + i] = short(1000); 
         else
             result[sizeOfMarkIdx * i + i] = 1; 
    }                                            
  

    delete[] GENO;
    delete[] bufferAllMarkers ;
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    D(printf("[info] Step2:  Time elapsed is %f. \n", elapsed ););
    fclose (bedFileReader);
    return 0;

}



// ***********************************************************
//              
//    The   typedef  dataType short
//          or typedef  dataType uchar determines the processing of data
//          in one byte or two byte at one time
//    The template <class T> sets the LD return in which datatype, T.
//
/// Caveats: 
//    a) For speed, the function uses the number of sample dividable by
//        4 (when dateType is uchar) or 8 (when dataType is short).
//        e.g.   when nSample = 3642, 
//                the nKeptSample = int(nSample / 4.0) *4  = 3640
//
// ***********************************************************

template <class T>
int calcLDFromBfile_quicker_nomissing (std::string bedFile, uint nSample, long nMarker, uint* theMarkIdx, unsigned long long sizeOfMarkIdx, uint* toAvert, int cutoff, int* ncpus, T* result, int* jump)
{
    // **************************************************************
    ///                    set timer
    // **************************************************************
    struct timespec start, finish;
    double elapsed;
    
    // **************************************************************
    ///                    set number of cpus
    // **************************************************************
    int nProcessors = *ncpus;
    omp_set_num_threads( nProcessors );
    D(printf("[info] Calc LD  based on %d cpus \n", nProcessors););
    // **************************************************************
    //                     Variables
    // **************************************************************
    uint processSamplePerRound = sizeof(dataType)*8 /2;
    const int individualsPerByte = 4;
    uint perMakerSizeOrig = ceil ( (nSample) / 1.0 / individualsPerByte );
    uint perMakerSize = uint( (nSample) / 1.0 / processSamplePerRound );
    uint nBlanks   = ( processSamplePerRound - (nSample) % processSamplePerRound  ) % processSamplePerRound; 
    long lSize =0;
    /// headerCoding is 9 bytes in total for plink 1.9 bedfile format, 3 bytes in older version.
    int const nByteHeader = 9;
    uchar headerCoding[nByteHeader +1 ] = { 0x6c, 0x1b, 0x01, 0xdc, 0x0f, 0xe7, 0x0f, 0x6b, 0x01};
    int const nByteHeader_older = 3;
    std::string formatVersion = "";
    int nThrowAway = 0;
    uchar  headerBuf[nByteHeader ] ;
    uint nKeptSample = perMakerSize * processSamplePerRound;
    std::vector<double>         GENO_VAR (sizeOfMarkIdx, -1 );
    std::vector<double>         GENO_Ex  (sizeOfMarkIdx, -1 );
    std::vector<double>         GENO_sum11  (sizeOfMarkIdx, -1 ); // sum of sample with genotype 11
    if(cutoff > sizeOfMarkIdx  ) cutoff = sizeOfMarkIdx;


    // **************************************************************
    //                1. validation of the bfile version
    // 2. validation of the bed file size, given the number of markers and samples
    // **************************************************************
    FILE* bedFileReader = fopen ( bedFile.c_str()  , "rb" );
    if (bedFileReader ==NULL) {fputs ("File not found error",stderr); }
    fseek (bedFileReader, 0 , SEEK_END);
    lSize = ftell (bedFileReader);
    rewind (bedFileReader);
    /// This checks the version of bed file. Examine if there is damage to
    /// bedfile, in which case, the estimated bed file size would be inconsistent with the
    /// acutal size.
    fread (&headerBuf,1, nByteHeader, bedFileReader);
    if(!memcmp(headerCoding, headerBuf, nByteHeader)) {
        D(printf("[info]This bed file is plink 1.9 bedfile format. (Newer) \n"); );
        formatVersion="1.9"; nThrowAway = nByteHeader;
    };
    if(!memcmp(headerCoding, headerBuf, nByteHeader_older)) {
        D(printf("[info]This bed file is plink 1.0 bedfile format. (Older)\n"););
        formatVersion="1.0"; nThrowAway = nByteHeader_older;
    };
    if(lSize  != long(perMakerSizeOrig * nMarker + nThrowAway) )
    {
        printf("[error] The size of bedFile %ld is inconsistenty with the estimated %u basd on %u samples and %d markers. \n",
                lSize, perMakerSizeOrig * nMarker + nThrowAway, perMakerSizeOrig, nMarker);
    }
    ////////////////////////////////////////////////////
    /// Creating map for a bype.
    /// number of ones,  00 coded for 0 in  a additive model, 11 for 2, 10 or 01 for 1
    /// This step assumes no missingness. (10 for missing.)
    ///    mapper : 65536 = 2 ^ (2*8 bits)
    /// invCountOnes is a "function (aa)  aa = ~aa; bitset<8> (aa).count"
    /// mapper2      is a "function (aa) aa = ~aa; aa & (aa<<1) & 0xaa" 
    /// mapper3      is a "function (aa)  aa & (aa<<1) & 0xaa" 
    /// countOnes    is a "function (aa)  bitset<8> (aa).count"
    size_t const sizeOfMap = (size_t) (pow( 2 , (sizeof(dataType) * 8) ));
    uint      invCountOnes  [sizeOfMap ]     = {0};
    uint      mapper2[sizeOfMap ]            = {0}; 
    uint      mapper3[sizeOfMap ]            = {0}; 
    uint      countOnes[sizeOfMap ]          = {0}; 
    dataType  markMissing[sizeOfMap ]        = {0}; 
#pragma omp parallel for 
    for (unsigned int i = 0; i < sizeOfMap ; i ++)
        countOnes[i] = (std::bitset<sizeof(dataType)*8> (i)).count();
    dataType screener = dataType (0xaaaaaaaaaaaaaaaa); /// the 64bit value automatically truncated to the dataType size.
#pragma omp parallel for 
    for (unsigned int i = 0; i < sizeOfMap ; i ++)
    {
        dataType aa = (dataType)(i);
        invCountOnes [i] = countOnes[aa];
        mapper2[i]       = countOnes[(dataType)(aa & (aa << 1 ) & screener )];
        mapper3[i]       = countOnes[(dataType)((dataType)(i) & ((dataType)(i) << 1 ) & screener) ];
        markMissing[i]   =  ((aa ^ aa <<1) & screener & ~aa ) ;
        markMissing[i]  = ~( markMissing[i] | (markMissing[i] >> 1) );
    }


    // ******************************************************
    ///     1. Calculating the var(SNP_i), E(SNP_i), E(SNP_i^2)
    //      2. Avert the genotype if specified.
    // ******************************************************
    unsigned long long loadSize = perMakerSizeOrig * sizeof(uchar) * (theMarkIdx[sizeOfMarkIdx-1] - theMarkIdx[0] +1) ;
    uchar* bufferAllMarkers = new uchar [loadSize ];

    D(printf("[info] Buffer size is %d Mb. \n", int(loadSize/1e6)););
    printf("[info] LD matrix size is %d Mb. \n", int(sizeOfMarkIdx * 1.0 * sizeOfMarkIdx/1e6 * 8)); 
    D(printf("[info] Map size is %d Mb. \n", int(sizeOfMap * 1.0 * sizeof(uint)* 4 + sizeOfMap * 1.0 * sizeof(dataType) * 1)););
    
    fseek (bedFileReader , perMakerSizeOrig * sizeof(uchar) * (theMarkIdx[0]) + nThrowAway, SEEK_SET );
    fread (bufferAllMarkers, 1, loadSize, bedFileReader);
    dataType* bufferMaker = NULL;  // This is pointer to the memory of a particular marker.
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    // *****************************************************************
    //
    //        Calculating cov(SNP_i, SNP_j),  LD(SNP_i, SNP_j)
    //
    //  When there is no missingness no-missingness, 
    //     cov(SNP_i, SNP_j) =  SNP_i,k & SNP_j,k 
    //                            + (SNP_i,k==2) + (SNP_j,k ==2) 
    //                            - (SNP_i,k == 2 & SNP_i,j == 0)
    //    Notably, - padding bits (when nSample %mod% 4 != 0) are set to 0.
    // *****************************************************************
    //
    //


    dataType** GENO = new dataType* [sizeOfMarkIdx] ;
#pragma omp parallel for
    for (unsigned long long i = 0; i < sizeOfMarkIdx; i ++)
    {
        GENO[i] = (dataType*) (perMakerSizeOrig * (long long)(theMarkIdx[i] - theMarkIdx[0]) + bufferAllMarkers);
    }


    std::vector<double>  E (sizeOfMarkIdx);
    std::vector<double>  E_sq (sizeOfMarkIdx);
    std::vector<double>  sum11 (sizeOfMarkIdx);
    std::vector<double>  VAR(sizeOfMarkIdx);
#pragma omp parallel for
   for (unsigned long long i = 0; i < sizeOfMarkIdx; i ++)
   {
      double sum_i     = 0;
      double sum11_i   = 0;
      for (unsigned long long k =0; k < perMakerSize; k ++) {
          sum_i    += countOnes[GENO[i][k]];
          sum11_i  += mapper2[(dataType)(GENO[i][k])];

      }
      sum11 [i] = sum11_i;
      E     [i] = double(sum_i  ) / (nKeptSample);
      E_sq  [i] = double(sum_i   + sum11_i*2 )/ (nKeptSample);
      VAR   [i] =  E_sq[i] - E[i] * E[i] ;
   }


    D(std::cout << "jump" << *jump << std::endl;);


#pragma omp parallel for
     for (unsigned long long i = 0; i < sizeOfMarkIdx; i ++)
     {
         dataType* GENO_i = GENO[i] ;
         unsigned long long starting = (*jump) > (i+1)? (*jump): (i +1);
         for (unsigned long long j = starting; j< sizeOfMarkIdx; j ++)
         {
             long double   sum_XY = 0;
             double   sum_0011 = 0; // 00 OR 11
             dataType* GENO_j = GENO[j] ;
             int sign = 1;
             if(toAvert[i] != toAvert[j]) 
                 sign = -1;
             if (j > i + cutoff)  break;
             for (unsigned long long k =0; k < perMakerSize; k ++)
             {
                 sum_XY   +=  countOnes  [ (dataType) (GENO_i[k]  & GENO_j[k]) ];
                 sum_0011 +=  mapper3[ (dataType) (GENO_i[k]  ^ GENO_j[k])     ];
             }
             sum_XY += sum11[i] + sum11[j] - sum_0011;
             double cov_XY = sum_XY / (nKeptSample ) - E[i] * E[j];
             double LD = 0;
             if(!(VAR[i] == 0 || VAR[j] == 0))
                 LD = cov_XY / sqrt ( (VAR[i] * VAR[j]) ) * sign;
             else
             {
                 
                 printf("[error] LD was set to 0, because var_[%d] = %f, var_[%d]  %f  E_i_sq (%f) - E_i * E_i\n", i, VAR[i], j, VAR[j],  E_sq[i]);
                 printf( "[error] , when parse the markerIdx %d-%d.\n", theMarkIdx[i],  theMarkIdx[j]);
                 printf("[error] %d - %d \n", theMarkIdx[sizeOfMarkIdx-1] , theMarkIdx[0] );
                 exit(-1);
             }

             if(sizeof(*result) == 1) // char
             {
                 char aa =  char(LD*250);
                 result[sizeOfMarkIdx * i + j] = char(LD * 250);
                 result[sizeOfMarkIdx * j + i] = char(LD * 250);
             }
             else if(sizeof(*result) == 2) // short
             {
                 result[sizeOfMarkIdx * i + j] = short(LD * 10000);
                 result[sizeOfMarkIdx * j + i] = short(LD * 10000);
             }
             else if (std::is_same<T, int>::value) // int
             {
                 result[sizeOfMarkIdx * i + j] = int(LD * 1e6);
                 result[sizeOfMarkIdx * j + i] = int(LD * 1e6);
             }
             else
             {
                 result[sizeOfMarkIdx * i + j] = LD ;
                 result[sizeOfMarkIdx * j + i] = LD ;
             }
         }
     }



#pragma omp parallel for                                       
    for (unsigned long long i = 0; i < sizeOfMarkIdx; i ++)          
    {                                                          
        if(sizeof(*result) == 1) // char
            result[sizeOfMarkIdx * i + i] = char(1 * 250);
        else if(sizeof(*result) == 2) // short
            result[sizeOfMarkIdx * i + i] = short(10000); 
         else if (std::is_same<T, int>::value) // int
            result[sizeOfMarkIdx * i + i] = T(1000000); 
        else
            result[sizeOfMarkIdx * i + i] = 1; 
    }                                            
  

    delete[] GENO;
    delete[] bufferAllMarkers ;
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    D(printf("[info] Step2:  Time elapsed is %f. \n", elapsed ););
    fclose (bedFileReader);
    return 0;

}




















extern "C"
{

    int _LDFromBfile_old(char** bedFileCstr, uint* nMarkers, uint* nSamples, uint* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, double* result)
    {
        std::string bedFile( *bedFileCstr) ;
        int nToAvert = 0;
#pragma omp parallel for reduction(+:nToAvert)
        for (unsigned int i = 0; i < *arrSize; i ++)
            nToAvert += toAvert[i];
        printf("Parameter settting:  \n ");
        printf("\t bedFile = %s \n"           , *bedFileCstr);
        printf("\t nMarkers = %d \n"          , *nMarkers);
        printf("\t nSamples = %d \n"          , *nSamples);
        printf("\t nTargetMarkers = %d \n"    , *arrSize);
        printf("\t nGenoToBeAverted = %d \n"  , nToAvert );
        printf("\t ldCalcCutoff = %d \n"      , *cutoff);
        printf("\t ncpus = %d \n"             , *ncpus);
        
        calcLDFromBfile_old <double> (bedFile, *nSamples, *nMarkers, theMarkIdx, *arrSize, toAvert, *cutoff, ncpus, result);
        return 0;
    }



    long getLineNum(std::string filename)
    {
        long number_of_lines = 0;
        std::string line;
        std::ifstream myfile(filename.c_str());
        while (std::getline(myfile, line))
            ++number_of_lines;
        return number_of_lines;
    }

//uint* nMarkers, uint* nSamples, uint* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, double* result
//
//
//

    int max(int a, int b)
    {
        if(a<b) return b;
        return a;
    }
   int min(int a, int b)
    {
        if(a>b) return b;
        return a;
    }


    int _getLDscore(char** bfile1Cstr, uint* theMarkIdx1,  uint* arrSize, uint* toAvert,int* cutoff,  int* ncpus, double* ldscore)
    {
        std::string bfile1( *bfile1Cstr) ;
        long nSamples1 = getLineNum(bfile1 + ".fam");
        long nMarkers1 = getLineNum(bfile1 + ".bim");
        int nToAvert = 0;
        //*ncpus = 24;
#pragma omp parallel for reduction(+:nToAvert)
        for (unsigned int i = 0; i < *arrSize; i ++)
            nToAvert += toAvert[i];
        printf("Parameter settting:  \n ");
        printf("\t bedFile = %s \n"          , bfile1.c_str());
        printf("\t nMarkers = %d \n"         , nMarkers1);
        printf("\t nSamples = %d \n"         , nSamples1);
        printf("\t ldCalcCutoff = %d \n"      , *cutoff);
        printf("\t nGenoToBeAverted = %d \n"  , nToAvert );
        printf("\t ncpus = %d \n"             , *ncpus);
        printf("\t nMarkerToCheck = %d \n"    , *arrSize);
        long winSize = min(*arrSize, 240000);
        //long winSize = min(*arrSize, 250);
        typedef  short v_size;
        v_size* result1  = new v_size[winSize * winSize]();
        uint start = 0;
        uint offset = 0;
        bool ifEnds = false;
        int jump = 0;
        for(uint k = 0; k < *arrSize; )
        {
            k =  offset + *cutoff * (offset != 0);
            //printf("\t offset : %d \n", offset);
            if(offset + winSize >= *arrSize)
            {
                winSize = *arrSize - offset;
                ifEnds = true;
            }

            calcLDFromBfile <v_size> (bfile1 + ".bed", nSamples1, nMarkers1, theMarkIdx1 + offset,  winSize , toAvert + offset, *cutoff, ncpus, result1, &jump);
            for (unsigned int i = start; i < winSize-  ((*cutoff) * (!ifEnds)) ; i ++)
            {
                ldscore[k]  = 0;
                for (unsigned int j = max(0, i - (*cutoff)); j <= i + (*cutoff); j ++)
                {
                    if(j + offset >= *arrSize) break;
                    double rsq_tmp = result1[i * winSize + j] * result1[i * winSize + j]/1000000.0;
                    ldscore[k] +=  rsq_tmp - (1-rsq_tmp)/nSamples1;
                }

            //printf("\t k = %d lds = %f \n"    , k, ldscore[k]);
                k++;
            }

            if(ifEnds) break;

            offset +=  winSize - 2 * (*cutoff);
            start  = *cutoff;

        }
        delete[] result1;
        return 0;
    }


    void getBoundIdx(double* pos, int* l_idx, int* r_idx, int* arrSize, double* cutoff)
    {
        int leftMost_idx = 0, rightMost_idx = *arrSize -1;
        for (int i = 0; i < *arrSize; i ++)
        {
            while ( pos[leftMost_idx] + *cutoff < pos[i]  )
                leftMost_idx ++;
            l_idx[i] = leftMost_idx;
        }

        for (int i = *arrSize -1; i >=0; i --)
        {
            while ( pos[rightMost_idx] - *cutoff > pos[i]  )
                rightMost_idx --;
            r_idx[i] = rightMost_idx;
        }

        
    }



    int _getLDscore2(char** bfile1Cstr, uint* theMarkIdx1, int* markerLeftIdx, int* markerRightIdx,   uint* arrSize, uint* toAvert, uint* cutoff,  int* ncpus, double* ldscore)
    {
        std::string bfile1( *bfile1Cstr) ;
        long nSamples1 = getLineNum(bfile1 + ".fam");
        long nMarkers1 = getLineNum(bfile1 + ".bim");
        int nToAvert = 0;
        //*ncpus = 24;
#pragma omp parallel for reduction(+:nToAvert)
        for (unsigned int i = 0; i < *arrSize; i ++)
            nToAvert += toAvert[i];
        printf("Parameter settting:  \n ");
        printf("\t bedFile = %s \n"          , bfile1.c_str());
        printf("\t nMarkers = %d \n"         , nMarkers1);
        printf("\t nSamples = %d \n"         , nSamples1);
        printf("\t ldCalcCutoff = %d \n"      , *cutoff);
        printf("\t nGenoToBeAverted = %d \n"  , nToAvert );
        printf("\t ncpus = %d \n"             , *ncpus);
        printf("\t nMarkerToCheck = %d \n"    , *arrSize);
        typedef  short v_size;
        // typedef  short v_size;
        long winSize = 240000;
        uint start = 0;
        uint offset = 0; winSize = min(winSize, *arrSize-offset);
        v_size* result1  = new v_size[winSize * winSize]();
        
        bool ifEnds = false;
        int jump = 0;

        calcLDFromBfile <v_size> (bfile1 + ".bed", nSamples1, nMarkers1, theMarkIdx1 + offset,  winSize , toAvert + offset, *cutoff, ncpus, result1, &jump);


        for(uint k = 0; k < *arrSize; k++)
        {


                ldscore[k]  = 0;
            if( markerRightIdx[k] - offset >= winSize ) // load the next block when hits the block boundary
            {
                offset =  markerLeftIdx[k] ;
                winSize = min(winSize, *arrSize-offset);
                calcLDFromBfile <v_size> (bfile1 + ".bed", nSamples1, nMarkers1, theMarkIdx1 + offset,  winSize , toAvert + offset, *cutoff, ncpus, result1, &jump);
                printf("[info] reading next block %d : %d + %d\n", offset, offset, winSize);
            }

            //if(k == 2881)
            //{
            //    printf("-------------------------------------------- \n");
            //    printf("-------------------------------------------- \n");
            //    printf("%d  - %d \n", markerLeftIdx[k], markerRightIdx[k]);
            //}
            for (unsigned int j = markerLeftIdx[k] -offset; j <= markerRightIdx[k] - offset ; j ++)
            {
                int i = k -offset;
                double rsq_tmp = result1[i * winSize + j] * result1[i * winSize + j]/1000000.0;
                // double rsq_tmp = result1[i * winSize + j] * result1[i * winSize + j]; /// /1000000.0;
                ldscore[k] +=  rsq_tmp - (1-rsq_tmp)/nSamples1;

            }


        }


  
        delete[] result1;
        return 0;
    }






    int _LDDiffBfile_new(char** bfile1Cstr, char** bfile2Cstr,
                         uint* theMarkIdx1, uint* theMarkIdx2,
                         int* markerLeftIdx, int* markerRightIdx,
                         uint* arrSize, uint* toAvert,
                         uint* cutoff,  int* ncpus, double* ldscore, double* ldscore2, double* chisq, int* numInconsist)
    {
        typedef  short v_size;
        std::string bfile1( *bfile1Cstr) ;
        std::string bfile2( *bfile2Cstr) ;
        long nSamples1 = getLineNum(bfile1 + ".fam");
        long nMarkers1 = getLineNum(bfile1 + ".bim");

        long nSamples2 = getLineNum(bfile2 + ".fam");
        long nMarkers2 = getLineNum(bfile2 + ".bim");

        int nToAvert = 0;
        int jump = 0;
        //*ncpus = 24;
#pragma omp parallel for reduction(+:nToAvert)
        for (unsigned int i = 0; i < *arrSize; i ++)
            nToAvert += toAvert[i];
        printf("Parameter settting:  \n ");
        printf("\t bedFile = %s \n"          , bfile1.c_str());
        printf("\t nMarkers = %d \n"         , nMarkers1);
        printf("\t nSamples = %d \n"         , nSamples1);
        printf("\t bedFile = %s \n"          , bfile2.c_str());
        printf("\t nMarkers = %d \n"         , nMarkers2);
        printf("\t nSamples = %d \n"         , nSamples2);
        printf("\t ldCalcCutoff = %d \n"      , *cutoff);
        printf("\t nGenoToBeAverted = %d \n"  , nToAvert );
        printf("\t ncpus = %d \n"             , *ncpus);
        printf("\t nMarkerToCheck = %d \n"    , *arrSize);
        long winSize = 240000;
        uint start = 0;
        uint offset = 0; winSize = min(winSize, *arrSize-offset);
        v_size* result1  = new v_size[winSize * winSize]();
        v_size* result2  = new v_size[winSize * winSize]();

        double varRdiff = (1.0 / nSamples1  + 1.0/nSamples2) ;
        double varRdiff_scaled =   (1.0 / nSamples1  + 1.0/nSamples2) * 1000000L;
        double significance = 32;
        bool ifEnds = false;

        calcLDFromBfile <v_size> (bfile1 + ".bed", nSamples1, nMarkers1, theMarkIdx1 + offset,  winSize , toAvert + offset, *cutoff, ncpus, result1, &jump);
        calcLDFromBfile <v_size> (bfile2 + ".bed", nSamples2, nMarkers2, theMarkIdx2 + offset,  winSize , toAvert + offset, *cutoff, ncpus, result2, &jump);

        for(uint k = 0; k < *arrSize; k++)
        {


            ldscore[k]      = 0;
            ldscore2[k]     = 0;
            chisq[k]        = 0;
            numInconsist[k] = 0;
            if( markerRightIdx[k] - offset >= winSize ) // load the next block when hits the block boundary
            {
                offset =  markerLeftIdx[k] ;
                winSize = min(winSize, *arrSize-offset);
                calcLDFromBfile <v_size> (bfile1 + ".bed", nSamples1, nMarkers1, theMarkIdx1 + offset,  winSize , toAvert + offset, *cutoff, ncpus, result1, &jump);
                calcLDFromBfile <v_size> (bfile2 + ".bed", nSamples2, nMarkers2, theMarkIdx2 + offset,  winSize , toAvert + offset, *cutoff, ncpus, result2, &jump);
                printf("[info] reading next block %d : %d + %d\n", offset, offset, winSize);
            }


            int count = 0;
            for (unsigned int j = markerLeftIdx[k] -offset; j <= markerRightIdx[k] - offset ; j ++)
            {
                int i = k -offset;
                double rsq_tmp = result1[i * winSize + j] * result1[i * winSize + j]/1000000.0;
                double rsq_tmp2 = result2[i * winSize + j] * result2[i * winSize + j]/1000000.0;
                ldscore[k]  +=  rsq_tmp  - (1 -  rsq_tmp)/nSamples1;
                ldscore2[k] +=  rsq_tmp2 - (1 - rsq_tmp2)/nSamples2;
                double diff  =  result1[i * winSize + j]  -  result2[i * winSize + j] ;

                chisq[k] += diff * diff / varRdiff_scaled;
                count ++;
                if(diff * diff / varRdiff_scaled > significance)  // chisq, df = 1
                    numInconsist[k] ++;
                
            }

            if  (count >0) chisq[k] /= count;


        }


  
        delete[] result1;
        delete[] result2;
        return 0;
    }



    int _LDDiffBfile(char** bfile1Cstr, char** bfile2Cstr, uint* theMarkIdx1, uint* theMarkIdx2, uint* theMarkerPos, uint* arrSize, uint* toAvert,int* cutoff,  int* ncpus, int* countOutliers)
    {
        std::string bfile1( *bfile1Cstr) ;
        std::string bfile2( *bfile2Cstr) ;
        long nSamples1 = getLineNum(bfile1 + ".fam");
        long nMarkers1 = getLineNum(bfile1 + ".bim");
        long nSamples2 = getLineNum(bfile2 + ".fam");
        long nMarkers2 = getLineNum(bfile2 + ".bim");
        int jump = 0;
        int nToAvert = 0;
        *ncpus = 24;
        
#pragma omp parallel for reduction(+:nToAvert)
        for (unsigned int i = 0; i < *arrSize; i ++)
            nToAvert += toAvert[i];

        printf("Parameters:  \n ");
        printf("\t bedFile1 = %s \n"          , bfile1.c_str());
        printf("\t bedFile2 = %s \n"          , bfile2.c_str());
        printf("\t nMarkers1 = %d \n"         , nMarkers1);
        printf("\t nSamples1 = %d \n"         , nSamples1);
        printf("\t nMarkers2 = %d \n"         , nMarkers2);
        printf("\t nSamples2 = %d \n"         , nSamples2);
        printf("\t ldCalcCutoff = %d \n"      , *cutoff);
        printf("\t nGenoToBeAverted = %d \n"  , nToAvert );
        printf("\t ldCalcCutoff = %d \n"      , *cutoff);
        printf("\t ncpus = %d \n"             , *ncpus);

        /// For examining only.
        /// *arrSize = 130;
        /// *cutoff = 30;
        

        long winSize =  6L * (*cutoff) ;

        typedef  short v_size;
        v_size* result1  = new v_size[winSize * winSize]();
        v_size* result2  = new v_size[winSize * winSize]();
        uint*   toAvert2 = new uint[winSize]();


        //uint*  countOutliers = new uint[*arrSize]();

        calcLDFromBfile <v_size> (bfile1 + ".bed", nSamples1, nMarkers1, theMarkIdx1 ,  winSize, toAvert, *cutoff, ncpus, result1 , &jump);
        calcLDFromBfile <v_size> (bfile2 + ".bed", nSamples2, nMarkers2, theMarkIdx2 ,  winSize, toAvert2, *cutoff, ncpus, result2, &jump);

//        for (unsigned int i = 0; i < winSize  ; i ++)
//            printf("%d:%d  ", toAvert[i], theMarkIdx1[i]);
//        printf("\n");
//        for (unsigned int i = 0; i < winSize  ; i ++)
//            printf("%d:%d  ", toAvert[i], theMarkIdx2[i]);
//        printf("\n");


        double varRdiff = 1.0 / nSamples1  + 1.0/nSamples2;
        double threshold = 32;
        uint start = 0;
        uint offset = 0;
        //for(uint k = 0; k < 10; )
        for(uint k = 0; k < *arrSize; )
        {
            for (unsigned int i = start; i < winSize - (*cutoff) ; i ++)
            {
                //countOutliers [k]  = 0;
                for (unsigned int j = max(0, i - (*cutoff)); j <= i + (*cutoff); j ++)
                {
                    //if(j + offset >= *arrSize) break;
                    long diff = result1[i * winSize + j] - result2[i * winSize + j];
                    //if(abs(diff) > 200 )  // chisq, df = 1
                    if(diff * diff / varRdiff > threshold * 1000000L)  // chisq, df = 1
                        countOutliers [k] ++;
                    //result1[i * winSize + j] = result2[i * winSize + j] =0;
                }
                k++;
            }

            offset +=  winSize - 2 * (*cutoff);
            start  = k - offset;
            if(offset + winSize >= *arrSize)
            {
                winSize = *arrSize - offset;
                //printf("[info] window size reduced to %ld.\n\n\n", winSize);
            }



            calcLDFromBfile <v_size> (bfile1 + ".bed", nSamples1, nMarkers1, theMarkIdx1 + offset,  winSize , toAvert + offset, *cutoff, ncpus, result1, &jump);
            calcLDFromBfile <v_size> (bfile2 + ".bed", nSamples2, nMarkers2, theMarkIdx2 + offset,  winSize, toAvert2 + offset, *cutoff, ncpus, result2, &jump);


            if(offset + winSize >= *arrSize)
            {
                for (unsigned int i = start; i < winSize ; i ++)
                {
                   for (unsigned int j = max(0, i - (*cutoff)); j <= i + (*cutoff); j ++)
                   {
                       if(j + offset >= *arrSize) break;
                       long diff = result1[i * winSize + j] - result2[i * winSize + j];
                       if(diff * diff / varRdiff > threshold * 1000000L)  // chisq, df = 1
                           countOutliers [k] ++;
                   }
                   k++;

                   if(k == *arrSize -1 ) break;
                }
                break;

             }
        }
        delete[] toAvert2;
        delete[] result1;
        delete[] result2;

       
        return 0;
    }


}
















template <class T>
int _LDFromBfile(char** bedFileCstr, uint* nMarkers, uint* nSamples, uint* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, T* result, int* jump, int* withNA)
    {
        std::string bedFile( *bedFileCstr) ;
        int nToAvert = 0;
#pragma omp parallel for reduction(+:nToAvert)
        for (unsigned int i = 0; i < *arrSize; i ++)
            nToAvert += toAvert[i];
        D(printf("Parameter settting:  \n ");                    );
        D(printf("\t bedFile = %s \n"           , *bedFileCstr); );
        D(printf("\t nMarkers = %d \n"          , *nMarkers);    );
        D(printf("\t nSamples = %d \n"          , *nSamples);    );
        D(printf("\t nTargetMarkers = %d \n"    , *arrSize);     );
        D(printf("\t nGenoToBeAverted = %d \n"  , nToAvert );    );
        D(printf("\t ldCalcCutoff = %d \n"      , *cutoff);      );
        D(printf("\t ncpus = %d \n"             , *ncpus);       );
        D(printf("\t jump = %d \n"             , *jump);         );
        D(printf("[info] Calculating LD matrix \n "); );
        if(*withNA)
        {
            calcLDFromBfile <T> (bedFile, *nSamples, *nMarkers, theMarkIdx, *arrSize, toAvert, *cutoff, ncpus, result, jump);
        }
        else
        {
            calcLDFromBfile_quicker_nomissing <T> (bedFile, *nSamples, *nMarkers, theMarkIdx, *arrSize, toAvert, *cutoff, ncpus, result, jump);
        }
        return 0;
    }


template int _LDFromBfile <short> (char** bedFileCstr, uint* nMarkers, uint* nSamples, uint* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, short* result, int* jump, int* withNA);
template int _LDFromBfile <int> (char** bedFileCstr, uint* nMarkers, uint* nSamples, uint* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, int* result, int* jump, int* withNA);
template int _LDFromBfile <float> (char** bedFileCstr, uint* nMarkers, uint* nSamples, uint* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, float* result, int* jump, int* withNA);
template int _LDFromBfile <double> (char** bedFileCstr, uint* nMarkers, uint* nSamples, uint* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus,double* result, int* jump, int* withNA);

