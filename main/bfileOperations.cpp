#include "dataTypeDef.h"
#include <bitset>


// ******************************************************************************
/// This calcualting the LD between SNPi and SNPj(missing genotypes) the GCTA way.
//   The var(SNPi) = var(SNPi^nonmissing)
//       mean (SNPi) = mean(SNPi^nonmissing)
//   The cov(SNPi, SNPj) were calcualted at nonmissing (SNPi & SNPj).
//
//   Note the the nonmissing (SNPi & SNPj) is a subset of individuals of 
//     nonmissing (SNPi), or   nonmissing (SNPj).
// ******************************************************************************
template <class T>
int calcLDFromBfile_gcta(std::string bedFile, int64 nSample, int64 nMarker, int64* theMarkIdx,
        int64 sizeOfMarkIdx, int64* toAvert, int cutoff, int ncpus, T* result, int64 jump)
{

    // **************************************************************
    ///                    set timer
    // **************************************************************
    struct timespec start, finish;
    double elapsed;
    
    // **************************************************************
    ///                    set number of cpus
    // **************************************************************
    int nProcessors = ncpus;
    omp_set_num_threads( nProcessors );
        //omp_set_num_threads( 1);
    D(printf("[info] Calc LD  based on %d cpus \n", nProcessors););
    // **************************************************************
    //                     Variables
    // **************************************************************
    int64 processSamplePerRound = sizeof(dataType)*8 /2;
    //uint64 perMakerSize = ceil ( (nSample) / 1.0 / processSamplePerRound );
    int64 perMakerSizeOrig = ceil ( (nSample) / 1.0 / 4);
    int64 perMakerSize = int( (nSample) / 1.0 / processSamplePerRound );
    int64 nBlanks   = ( processSamplePerRound - (nSample) % processSamplePerRound  ) % processSamplePerRound; 
    int64 lSize =0;
    /// headerCoding is 9 bytes in total for plink 1.9 bedfile format, 3 bytes in older version.
    int const nByteHeader = 9;
    uchar headerCoding[nByteHeader +1 ] = { 0x6c, 0x1b, 0x01, 0xdc, 0x0f, 0xe7, 0x0f, 0x6b, 0x01};
    int const nByteHeader_older = 3;
    std::string formatVersion = "";
    int nThrowAway = 0;
    uchar  headerBuf[nByteHeader ] ;
    //uint64 nKeptSample = nSample;
    int64 nKeptSample = perMakerSize * processSamplePerRound;
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
    size_t readSize = fread (&headerBuf,1, nByteHeader, bedFileReader);
    if(!memcmp(headerCoding, headerBuf, nByteHeader)) {D(printf("[info]This bed file is plink 1.9 bedfile format. (Newer) \n")); formatVersion="1.9"; nThrowAway = nByteHeader;};
    if(!memcmp(headerCoding, headerBuf, nByteHeader_older)) {D(printf("[info]This bed file is plink 1.0 bedfile format. (Older)\n")); formatVersion="1.0"; nThrowAway = nByteHeader_older;};
    if(lSize  != long(perMakerSizeOrig * nMarker + nThrowAway) )
    {
        printf("[error] The size of bedFile %lld is inconsistenty with the estimated %lld basd on %lld samples and %lld markers. \n", lSize, perMakerSizeOrig * nMarker + nThrowAway, perMakerSizeOrig, nMarker);
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
    unsigned long long loadSize = perMakerSizeOrig * sizeof(uchar) * (theMarkIdx[sizeOfMarkIdx-1] - theMarkIdx[0] +1) ;
    uchar* bufferAllMarkers = new uchar [loadSize ];
    D(printf("[info] Buffer size is %d Mb. \n", int(loadSize/1e6)););
    fseek (bedFileReader , perMakerSizeOrig * sizeof(uchar) * (theMarkIdx[0]) + nThrowAway, SEEK_SET );
    readSize = fread (bufferAllMarkers, 1, loadSize, bedFileReader);
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

    std::vector<double>  E (sizeOfMarkIdx);
    std::vector<double>  E_sq (sizeOfMarkIdx);
    std::vector<double>  VAR(sizeOfMarkIdx);

   for (int64 i = 0; i < sizeOfMarkIdx; i ++)
   {

      dataType marker = 0  ;
      double sum_i     = 0;
      double sum11_i   = 0;
      double nMissing  = 0;
      dataType* GENO_i = GENO[i] ;
      for (int64  k =0; k < perMakerSize; k ++) {
          marker = markMissing[ GENO_i[k]  ]   ;
          nMissing += countOnes[(dataType)(~marker)];
          sum_i    += countOnes[GENO_i [k] & marker];
          sum11_i  += mapper2[(dataType)(GENO_i [k] & marker)];
      }
      nMissing = nMissing /2;
      E     [i] = double(sum_i  ) / (nKeptSample - nMissing);
      E_sq  [i] = double(sum_i   + sum11_i*2 )/ (nKeptSample -nMissing);
      VAR   [i] =  E_sq[i] - E[i] * E[i] ;
      if(VAR[i] == 0)
      {
          throw(i) ;
          return -1;
      }
   }


 #pragma omp parallel for
    for (int64 i = 0; i < sizeOfMarkIdx; i ++)
    {
        dataType* GENO_i = GENO[i] ;


        int starting = (jump) > (i+1)? (jump): (i +1);
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
            double E_i2         = double(sum_i  ) / (nKeptSample );
            double E_j2         = double(sum_j  ) / (nKeptSample );


            nMissing = nMissing /2;
            E_i         = E[i];
            E_j         = E[j];
            E_i_sq      = E_sq[i];
            E_j_sq      = E_sq[j];
            var_i =  VAR[i] ;
            var_j =  VAR[j] ;

            sum_XY += sum11_i + sum11_j - sum_0011;
            //double cov_XY = (sum_XY - sum_i * sum_j) / (nSample - nMissing);
            //double cov_XY = sum_XY / (nKeptSample ) - E_i * E_j;
            double cov_XY = sum_XY / (nKeptSample ) + E_i * E_j * (nKeptSample -nMissing) / nKeptSample
                                    - E_i * E_j2 - E_i2 * E_j;
            //// double LD = cov_XY / sqrt ( (var_i * var_j) ) * sign;
            /// This is for rare cases when SNP i or j has 0 variance after removing the NAs genotype individuals.
            double LD = 0.001;
            if(!(var_i <= 0 || var_j <= 0 || var_i * var_j <= 0))
                LD = cov_XY / sqrt ( (var_i * var_j) ) * sign;
            else
            {
                printf("[error] LD was set to 0, because var_[%lld] = %f, var_[%lld]  %f  E_i_sq (%f) - E_i * E_i\n",
                        i, VAR[i], j, VAR[j],  E_sq[i]);
                printf("[error] , when parse the markerIdx %lld-%lld.\n", theMarkIdx[i],  theMarkIdx[j]);
                printf("[error] %lld - %lld \n", theMarkIdx[sizeOfMarkIdx-1] , theMarkIdx[0] );
                exit(-1);

            }

            saveData<T>(LD, i, j, result, sizeOfMarkIdx);

           }
     }

     setDiag<T>(result, sizeOfMarkIdx);

 

    delete[] GENO;
    delete[] bufferAllMarkers ;
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
    D(printf("[info] Step2:  Time elapsed is %f. \n", elapsed ););
    fclose (bedFileReader);
    return 0;

}



///  This tried to take all the missing values out, before calculating 
///   correlation between a SNPi and SNPj.
///  This tricky part is when the removed individuals with missing genotypes 
//    at SNPi or SNPj, the var(SNP_i^nomissing) == 0, which leading to cov/0 problem.
//   

template <class T>
int calcLDFromBfile       (std::string bedFile, int64 nSample, int64 nMarker, uint64* theMarkIdx, uint64 sizeOfMarkIdx, uint64* toAvert, int cutoff, int ncpus, T* result, int jump)
{
    // **************************************************************
    ///                    set timer
    // **************************************************************
    struct timespec start, finish;
    double elapsed;
    
    // **************************************************************
    ///                    set number of cpus
    // **************************************************************
    int nProcessors = ncpus;
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
    int64 lSize =0;
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
    size_t readSize = fread (&headerBuf,1, nByteHeader, bedFileReader);
    if(!memcmp(headerCoding, headerBuf, nByteHeader)) {
        D(printf("[info]This bed file is plink 1.9 bedfile format.  \n"););
        formatVersion="1.9"; nThrowAway = nByteHeader;
    };
    if(!memcmp(headerCoding, headerBuf, nByteHeader_older)) {
        D(printf("[info]This bed file is plink 1.0 bedfile format. \n"););
        formatVersion="1.0"; nThrowAway = nByteHeader_older;
    };
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
    unsigned long long loadSize = perMakerSizeOrig * sizeof(uchar) * (theMarkIdx[sizeOfMarkIdx-1] - theMarkIdx[0] +1) ;
    uchar* bufferAllMarkers = new uchar [loadSize ];

    D(printf("[info] Buffer size is %d Mb. \n", int(loadSize/1e6)););
    
    fseek (bedFileReader , perMakerSizeOrig * sizeof(uchar) * (theMarkIdx[0]) + nThrowAway, SEEK_SET );
    readSize = fread (bufferAllMarkers, 1, loadSize, bedFileReader);
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



 #pragma omp parallel for
    for (unsigned long long i = 0; i < sizeOfMarkIdx; i ++)
    {
        dataType* GENO_i = GENO[i] ;


        int starting = (jump) > (i+1)? (jump): (i +1);
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
            saveData<T> (LD, i, j, result, sizeOfMarkIdx);
        }
    }



///#pragma omp parallel for                                       
///    for (unsigned long long i = 0; i < sizeOfMarkIdx; i ++)          
///    {                                                          
///         if(sizeof(*result) == 1) // char
///             result[sizeOfMarkIdx * i + i] = T(1 * 250);
///         else if(sizeof(*result) == 2) // short
///             result[sizeOfMarkIdx * i + i] = T(1000); 
///         else
///             result[sizeOfMarkIdx * i + i] = 1; 
///    }                                            
  

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
int calcLDFromBfile_quicker_nomissing (std::string bedFile, int64 nSample, int64 nMarker, int64* theMarkIdx, uint64 sizeOfMarkIdx, int64* toAvert, int cutoff, int ncpus, T* result, int jump)
{
    // **************************************************************
    ///                    set timer
    // **************************************************************
    struct timespec start, finish;
    double elapsed;
    
    // **************************************************************
    ///                    set number of cpus
    // **************************************************************
    int nProcessors = ncpus;
    omp_set_num_threads( nProcessors );
    D(printf("[info] Calc LD  based on %d cpus \n", nProcessors););
    // **************************************************************
    //                     Variables
    // **************************************************************
    const int individualsPerByte = 4;
    int64    processSamplePerRound = sizeof(dataType)*8 /2;
    int64    perMakerSizeOrig = ceil ( (nSample) / 1.0 / individualsPerByte );
    int64    perMakerSize = uint( (nSample) / 1.0 / processSamplePerRound );
    int64    nBlanks   = ( processSamplePerRound - (nSample) % processSamplePerRound  ) % processSamplePerRound; 
    int64 lSize =0;
    /// headerCoding is 9 bytes in total for plink 1.9 bedfile format, 3 bytes in older version.
    int const nByteHeader = 9;
    uchar headerCoding[nByteHeader +1 ] = { 0x6c, 0x1b, 0x01, 0xdc, 0x0f, 0xe7, 0x0f, 0x6b, 0x01};
    int const nByteHeader_older = 3;
    std::string formatVersion = "";
    int nThrowAway = 0;
    uchar  headerBuf[nByteHeader ] ;
    uint64 nKeptSample = perMakerSize * processSamplePerRound;
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
    size_t readSize = fread (&headerBuf,1, nByteHeader, bedFileReader);
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
        printf("[error] The size of bedFile %lld is inconsistenty with the estimated %lld basd on %lld samples and %lld markers. \n",
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
    readSize = fread (bufferAllMarkers, 1, loadSize, bedFileReader);
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
    for (unsigned long long i = 0; i < sizeOfMarkIdx; i ++)
    {
        GENO[i] = (dataType*) (perMakerSizeOrig * (theMarkIdx[i] - theMarkIdx[0]) + bufferAllMarkers);
    }


    std::vector<double>  E (sizeOfMarkIdx);
    std::vector<double>  E_sq (sizeOfMarkIdx);
    std::vector<double>  sum11 (sizeOfMarkIdx);
    std::vector<double>  VAR(sizeOfMarkIdx);
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
       if(VAR[i] == 0)
       {
           throw(i) ;
           return -1;
       }
    }

#pragma omp parallel for
     for (unsigned long long i = 0; i < sizeOfMarkIdx; i ++)
     {
         dataType* GENO_i = GENO[i] ;
         unsigned long long starting = (jump) > (i+1)? (jump): (i +1);
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
                 printf("[error] LD was set to 0, because var_[%lld] = %f, var_[%lld]  %f  E_i_sq (%f) - E_i * E_i\n",
                         i, VAR[i], j, VAR[j],  E_sq[i]);
                 printf("[error] , when parse the markerIdx %lld-%lld.\n", theMarkIdx[i],  theMarkIdx[j]);
                 printf("[error] %lld - %lld \n", theMarkIdx[sizeOfMarkIdx-1] , theMarkIdx[0] );
                 exit(-1);
             }
             saveData<T> (LD, i, j, result, sizeOfMarkIdx);
         }
     }

     setDiag<T>(result, sizeOfMarkIdx);


/// #pragma omp parallel for                                       
///     for (unsigned long long i = 0; i < sizeOfMarkIdx; i ++)          
///     {                                                          
///         if(sizeof(*result) == 1) // char
///             result[sizeOfMarkIdx * i + i] = char(1 * 250);
///         else if(sizeof(*result) == 2) // short
///             result[sizeOfMarkIdx * i + i] = short(10000); 
///          else if (std::is_same<T, int>::value) // int
///             result[sizeOfMarkIdx * i + i] = T(1000000); 
///         else
///             result[sizeOfMarkIdx * i + i] = 1; 
///     }                                            
  

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

    
    long getLineNum(std::string filename)
    {
        long number_of_lines = 0;
        std::string line;
        std::ifstream myfile(filename.c_str());
        while (std::getline(myfile, line))
            ++number_of_lines;
        return number_of_lines;
    }
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
}




template <class T>
int _LDFromBfile(char** bedFileCstr, uint* nMarkers, uint* nSamples, int64* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, T* result, int* jump, int* withNA)
{
    std::string bedFile( *bedFileCstr) ;
    int nToAvert = 0;
#pragma omp parallel for reduction(+:nToAvert)
    for (unsigned int i = 0; i < *arrSize; i ++)
        nToAvert += toAvert[i];
    int64*  theMarkIdx_tmp = new int64 [*arrSize];
    int64*     toAvert_tmp = new int64 [*arrSize];
    for(uint i = 0; i < *arrSize ; i ++)
    {
        theMarkIdx_tmp[i] = theMarkIdx[i];
        toAvert_tmp[i]    = toAvert[i];
    }

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
        calcLDFromBfile_gcta <T> (bedFile, int64(*nSamples), int64(*nMarkers), theMarkIdx_tmp, uint64(*arrSize), toAvert_tmp, *cutoff, *ncpus, result, *jump);
    }
    else
    {

        calcLDFromBfile_quicker_nomissing<T> (bedFile, int64(*nSamples), int64(*nMarkers), theMarkIdx_tmp, uint64(*arrSize), toAvert_tmp, *cutoff, *ncpus, result, *jump);
    }

    delete [] theMarkIdx_tmp;
    delete [] toAvert_tmp;
    return 0;
}


template int _LDFromBfile <int>        (char** bedFileCstr, uint* nMarkers, uint* nSamples, int64* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, int* result, int* jump, int* withNA);
template int _LDFromBfile <short>        (char** bedFileCstr, uint* nMarkers, uint* nSamples, int64* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus,short* result, int* jump, int* withNA);
template int _LDFromBfile <char>        (char** bedFileCstr, uint* nMarkers, uint* nSamples, int64* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, char* result, int* jump, int* withNA);
template int _LDFromBfile <float>      (char** bedFileCstr, uint* nMarkers, uint* nSamples, int64* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, float* result, int* jump, int* withNA);
template int _LDFromBfile <double>     (char** bedFileCstr, uint* nMarkers, uint* nSamples, int64* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus,double* result, int* jump, int* withNA);
template int _LDFromBfile <smatrix_d > (char** bedFileCstr, uint* nMarkers, uint* nSamples, int64* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, smatrix_d* result, int* jump, int* withNA);
template int _LDFromBfile <smatrix_f > (char** bedFileCstr, uint* nMarkers, uint* nSamples, int64* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, smatrix_f* result, int* jump, int* withNA);
template int _LDFromBfile <smatrix_i > (char** bedFileCstr, uint* nMarkers, uint* nSamples, int64* theMarkIdx, uint* arrSize, uint* toAvert, int* cutoff, int* ncpus, smatrix_i* result, int* jump, int* withNA);




