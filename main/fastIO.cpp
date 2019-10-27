#include "dataTypeDef.h"



// ****************************************************************************
/// This is for faster io (than read.table in R) of bim, frq, summary Stat.
//  
//

void readBim(char** prefix, char** chrID, char** rsID,  int* M)
{
    std::string prefixStr = std::string(*prefix);


    FILE* bimReader = NULL;
    if ((bimReader = fopen ( (prefixStr  + ".bim").c_str(), "r"))==NULL) {
        printf("\t[error] Cannot open file %s\n", (prefixStr  + ".bim").c_str());
    }

//    for(uint i = 0; i < M; i ++)
//        fscanf(bimReader,"%s %s", chrID,rsID);
    fclose(bimReader);
}


////   template <class T>
////   int calcLDFromBfile (std::string bedFile, uint nSample, long nMarker, uint* theMarkIdx, long sizeOfMarkIdx, uint* toAvert, int cutoff, int* ncpus, T* result)
////   {
////       // **************************************************************
////       ///                    set timer
////       // **************************************************************
////       struct timespec start, finish;
////       double elapsed;
////       
////       // **************************************************************
////       ///                    set number of cpus
////       // **************************************************************
////       int nProcessors = omp_get_max_threads();
////       if(*ncpus < nProcessors) nProcessors = *ncpus;
////       omp_set_num_threads( nProcessors );
////           //omp_set_num_threads( 1);
////       printf("[info] Calc LD  based on %d cpus \n", nProcessors);
////       // **************************************************************
////       //                     Variables
////       // **************************************************************
////       typedef unsigned short dataType;
////       uint processSamplePerRound = sizeof(dataType)*8 /2;
////       //uint perMakerSize = ceil ( (nSample) / 1.0 / processSamplePerRound );
////       uint perMakerSizeOrig = ceil ( (nSample) / 1.0 / 4);
////       uint perMakerSize = int( (nSample) / 1.0 / processSamplePerRound );
////       uint nBlanks   = ( processSamplePerRound - (nSample) % processSamplePerRound  ) % processSamplePerRound; //
////       long lSize =0;
////       /// headerCoding is 9 bytes in total for plink 1.9 bedfile format, 3 bytes in older version.
////       int const nByteHeader = 9;
////       uchar headerCoding[nByteHeader +1 ] = { 0x6c, 0x1b, 0x01, 0xdc, 0x0f, 0xe7, 0x0f, 0x6b, 0x01};
////       int const nByteHeader_older = 3;
////       std::string formatVersion = "";
////       int nThrowAway = 0;
////       uchar  headerBuf[nByteHeader ] ;
////       //uint nKeptSample = nSample;
////       uint nKeptSample = perMakerSize * processSamplePerRound;
////       std::vector<std::vector<dataType> > GENO      (sizeOfMarkIdx, std::vector<dataType> ( nSample, 0) );
////       std::vector<double>         GENO_VAR (sizeOfMarkIdx, -1 );
////       std::vector<double>         GENO_Ex  (sizeOfMarkIdx, -1 );
////       std::vector<double>         GENO_sum11  (sizeOfMarkIdx, -1 ); // sum of sample with genotype 11
////       if(cutoff > sizeOfMarkIdx  ) cutoff = sizeOfMarkIdx;
////   
////       // **************************************************************
////       //                1. validation of the bfile version
////       // 2. validation of the bed file size, given the number of markers and samples
////       // **************************************************************
////       FILE* bedFileReader = fopen ( bedFile.c_str()  , "rb" );
////       if (bedFileReader ==NULL) {fputs ("File not found error",stderr); error("");}
////       fseek (bedFileReader, 0 , SEEK_END);
////       lSize = ftell (bedFileReader);
////       rewind (bedFileReader);
////       /// This checks the version of bed file. Examine if there is damage to
////       /// bedfile, in which case, the estimated bed file size would be inconsistent with the
////       /// acutal size.
////       fread (&headerBuf,1, nByteHeader, bedFileReader);
////       if(!memcmp(headerCoding, headerBuf, nByteHeader)) {printf("[info] This bed file is plink 1.9 bedfile format. (Newer) \n"); formatVersion="1.9"; nThrowAway = nByteHeader;};
////       if(!memcmp(headerCoding, headerBuf, nByteHeader_older)) {printf("[info] This bed file is plink 1.0 bedfile format. (Older)\n"); formatVersion="1.0"; nThrowAway = nByteHeader_older;};
////       if(lSize  != long(perMakerSizeOrig * nMarker + nThrowAway) )
////       {
////           printf("[error] The size of bedFile %ld is inconsistenty with the estimated %u basd on %u samples and %d markers. \n", lSize, perMakerSizeOrig * nMarker + nThrowAway, perMakerSizeOrig, nMarker);
////           error("");
////       }
////       ////////////////////////////////////////////////////
////       /// Creating map for a bype.
////       /// number of ones,  00 coded for 0 in  a additive model, 11 for 2, 10 or 01 for 1
////       /// This step assumes no missingness. (10 for missing.)
////       ///    mapper : 65536 = 2 ^ (2*8 bits)
////       size_t const sizeOfMap = (size_t) (pow( 2 , (sizeof(dataType) * 8) ));
////       // invCountOnes is a "function (aa)  aa = ~aa; bitset<8> (aa).count"
////       // mapper2   is a "function (aa) aa = ~aa; aa & (aa<<1) & 0xaa" 
////       // mapper3   is a "function (aa)  aa & (aa<<1) & 0xaa" 
////       // countOnes is a "function (aa)  bitset<8> (aa).count"
////       // mapper5   is a "function (aa)  aa = ~aa ; ( ((aa ^ aa <<1) & 0xaa ) >>1 ) *3"
////       uint      invCountOnes  [sizeOfMap ] = {0};
////       uint      mapper2[sizeOfMap ]        = {0}; 
////       uint      mapper3[sizeOfMap ]        = {0}; 
////       uint      countOnes[sizeOfMap ]      = {0}; 
////       dataType  mapper5[sizeOfMap ]        = {0}; 
////       dataType  markMissing[sizeOfMap ]        = {0}; 
////   #pragma omp parallel for 
////       for (unsigned int i = 0; i < sizeOfMap ; i ++)
////           countOnes[i] = (std::bitset<sizeof(dataType)*8> (i)).count();
////       dataType screener = dataType (0xaaaaaaaaaaaaaaaa); /// the 64bit value automatically truncated to the dataType size.
////   #pragma omp parallel for 
////       for (unsigned int i = 0; i < sizeOfMap ; i ++)
////       {
////           dataType aa = (dataType)(i);
////           invCountOnes [i] = countOnes[aa];
////           mapper2[i]       = countOnes[(dataType)(aa & (aa << 1 ) & screener )];
////           mapper3[i]       = countOnes[(dataType)((dataType)(i) & ((dataType)(i) << 1 ) & screener) ];
////           markMissing[i]   =  ((aa ^ aa <<1) & screener & ~aa ) ;
////           markMissing[i]  = ~( markMissing[i] | (markMissing[i] >> 1) );
////       }
////   
////   
////       // ******************************************************
////       ///     1. Calculating the var(SNP_i), E(SNP_i), E(SNP_i^2)
////       //      2. Avert the genotype if specified.
////       // ******************************************************
////       unsigned long loadSize = perMakerSizeOrig * sizeof(uchar) * (theMarkIdx[sizeOfMarkIdx-1] - theMarkIdx[0] +1) ;
////       uchar* bufferAllMarkers = new uchar [loadSize ];
////   
////       printf("[info] Buffer size is %d Mb. \n", int(loadSize/1e6));
////       
////       fseek (bedFileReader , perMakerSizeOrig * sizeof(uchar) * (theMarkIdx[0]) + nThrowAway, SEEK_SET );
////       fread (bufferAllMarkers, 1, loadSize, bedFileReader);
////       dataType* bufferMaker = NULL;  // This is pointer to the memory of a particular marker.
////       for (unsigned int i = 0; i < sizeOfMarkIdx; i ++)
////       {
////           bufferMaker = reinterpret_cast<dataType*> (perMakerSizeOrig * (theMarkIdx[i] - theMarkIdx[0]) + bufferAllMarkers);
////           for (unsigned int j =0; j<perMakerSize; j ++)
////               GENO[i][j] =  bufferMaker[j];
////   
////       }
////   
////       
////       clock_gettime(CLOCK_MONOTONIC, &start);
////   
////       // *****************************************************************
////       //
////       //        Calculating cov(SNP_i, SNP_j),  LD(SNP_i, SNP_j)
////       //
////       //  When there is no missingness no-missingness, 
////       //     cov(SNP_i, SNP_j) =  SNP_i,k & SNP_j,k 
////       //                            + (SNP_i,k==2) + (SNP_j,k ==2) 
////       //                            - (SNP_i,k == 2 & SNP_i,j == 0)
////       //    Notably, - padding bits (when nSample %mod% 4 != 0) are set to 0.
////       // *****************************************************************
////       //
////       //
////   
////   
////   
////       ///// modified upto here
////   
////   
////   
////   #pragma omp parallel for
////       for (unsigned int i = 0; i < sizeOfMarkIdx; i ++)
////       {
////           for (unsigned int j =i+1; j< sizeOfMarkIdx; j ++)
////           {
////               int sign = 1;
////               if(toAvert[i] != toAvert[j]) 
////                   sign = -1;
////               uint nMissing = 0;
////               double   E_i = 0;
////               double   E_j = 0;
////               double   sum11_i = 0;
////               double   sum11_j = 0;
////               double   sum_i = 0;
////               double   sum_j = 0;
////               double   E_i_sq = 0;
////               double   E_j_sq = 0;
////               double   var_i = 0 ;
////               double   var_j = 0 ;
////               long double   sum_XY = 0;
////               double   sum_0011 = 0; // 00 OR 11
////               dataType marker = 0  ;
////               
////   
////               if (j > i + cutoff)  break;
////               for (unsigned int k =0; k < perMakerSize; k ++)
////               {
////   
////                   marker = markMissing[GENO [i][k]]  & markMissing[GENO [j][k]]   ;
////                   nMissing += countOnes[(dataType)(~marker)];
////                   sum_i    += countOnes[GENO[i][k] & marker];
////                   sum_j    += countOnes[GENO[j][k] & marker];
////                   sum11_i  += mapper2[(dataType)(GENO[i][k] & marker)];
////                   sum11_j  += mapper2[(dataType)(GENO[j][k] & marker)];
////                   GENO [i][k] &=  marker;
////                   GENO [j][k] &=  marker;
////                   sum_XY   +=  countOnes  [ (dataType) (GENO [i][k] & GENO [j][k]) ];
////                   sum_0011 +=  mapper3[ (dataType) (GENO [i][k] ^ GENO [j][k]) ];
////               }
////   
////               nMissing = nMissing /2;
////               E_i         = double(sum_i  ) / (nKeptSample - nMissing);
////               E_j         = double(sum_j  ) / (nKeptSample - nMissing);
////               E_i_sq     = double(sum_i   + sum11_i*2 )/ (nKeptSample -nMissing);
////               E_j_sq     = double(sum_j   + sum11_j*2 )/ (nKeptSample -nMissing);
////               var_i =  E_i_sq - E_i * E_i ;
////               var_j =  E_j_sq - E_j * E_j ;
////   
////   
////   
////               sum_XY += sum11_i + sum11_j - sum_0011;
////               //double cov_XY = (sum_XY - sum_i * sum_j) / (nSample - nMissing);
////               double cov_XY = sum_XY / (nKeptSample - nMissing) - E_i * E_j;
////               double LD = cov_XY / sqrt ( (var_i * var_j) ) * sign;
////   
////   
////   
////               if(sizeof(*result) == 1) // char
////               {
////                   char aa =  char(LD*250);
////                   printf("(%f, %f)\n",LD, (0L |255 & aa)/250.0  );
////                   result[sizeOfMarkIdx * i + j] = char(LD * 250);
////                   result[sizeOfMarkIdx * j + i] = char(LD * 250);
////               }
////               else if(sizeof(*result) == 2) // short
////               {
////                   result[sizeOfMarkIdx * i + j] = short(LD * 1000);
////                   result[sizeOfMarkIdx * j + i] = short(LD * 1000);
////               }
////               else
////               {
////                   result[sizeOfMarkIdx * i + j] = LD ;
////                   result[sizeOfMarkIdx * j + i] = LD ;
////               }
////           }
////       }
////   
////       delete[] bufferAllMarkers ;
////       clock_gettime(CLOCK_MONOTONIC, &finish);
////       elapsed = (finish.tv_sec - start.tv_sec);
////       elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
////       printf("[info] Step2:  Time elapsed is %f. \n", elapsed );
////       fclose (bedFileReader);
////       return 0;
////   
////   }
////   
////   
////   
////   
////   
