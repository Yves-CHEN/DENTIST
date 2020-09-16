#include "dataTypeDef.h"



int getGeno (std::string bedFile, uint nSample, long nMarker, int64* theMarkIdx, long sizeOfMarkIdx,  int* result)
{

    const int individualsPerByte = 4;

    int64 perMakerSize = ceil ( (nSample) / 1.0 / individualsPerByte  );
    uint nBlanks   = ( 4 - (nSample) % individualsPerByte  ) % individualsPerByte; //
    int64 lSize =0;
    /// headerCoding is 9 bytes in total for plink 1.9 bedfile format, 3 bytes in older version.
    int const nByteHeader = 9;
    uchar headerCoding[nByteHeader +1 ] = { 0x6c, 0x1b, 0x01, 0xdc, 0x0f, 0xe7, 0x0f, 0x6b, 0x01};
    int const nByteHeader_older = 3;
    std::string formatVersion = "";
    int nThrowAway = 0;
    uchar  headerBuf[nByteHeader ] ;
    uint nKeptSample = nSample;
    std::vector<std::vector<uchar> > GENO      (sizeOfMarkIdx, std::vector<uchar> ( nSample, 0) );

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
    size_t ret = fread (&headerBuf,1, nByteHeader, bedFileReader);
    if(!memcmp(headerCoding, headerBuf, nByteHeader)) {printf("[info] This bed file is plink 1.9 bedfile format. (Newer) \n"); formatVersion="1.9"; nThrowAway = nByteHeader;};
    if(!memcmp(headerCoding, headerBuf, nByteHeader_older)) {printf("[info] This bed file is plink 1.0 bedfile format. (Older)\n"); formatVersion="1.0"; nThrowAway = nByteHeader_older;};
    if(lSize  != long(perMakerSize * nMarker + nThrowAway) )
    {
        printf("[error] The size of bedFile %lld is inconsistenty with the estimated %lld based on %lld samples and %ld markers. \n", lSize, perMakerSize * nMarker + nThrowAway, perMakerSize, nMarker);
    }


    unsigned long loadSize = perMakerSize * sizeof(char) * (theMarkIdx[sizeOfMarkIdx-1] - theMarkIdx[0] +1) ;
    uchar* bufferAllMarkers = new unsigned char [loadSize ];
    printf("[info] Buffer size is %d Mb. \n", int(loadSize/1e6));
    fseeko64 (bedFileReader , perMakerSize * sizeof(char) * (theMarkIdx[0]) + nThrowAway, SEEK_SET );
    ret = fread (bufferAllMarkers, 1, loadSize, bedFileReader);


    uchar* bufferMaker = NULL;  // This is pointer to the memory of a particular marker.
    for (unsigned int i = 0; i < sizeOfMarkIdx; i ++)
    {
        double sum = 0, sum2 = 0;
        double  E_x=0, E_x_square=0;
        bufferMaker = perMakerSize * (theMarkIdx[i] - theMarkIdx[0]) + bufferAllMarkers;

        // avert the genotypes. For 00 and 11, it would ~ geno. For 01 and 10, it is more difficult.
        for (unsigned int j =0; j<perMakerSize; j ++)
        {
            //GENO[i][j] =  
            //
            uchar aa = bufferMaker[j];
            for (uint k = 0; k < 4; k ++)
            {
                result[i * perMakerSize*4 + j * 4 +k ] = ceil(uchar(aa & uchar(3)) /2.0);
                aa = aa >> 2;
            }
               
        }

    }

    return 0;
}



extern "C"
{

    int getGenoInvoker(char** bedFileCstr, uint* nMarkers, uint* nSamples, int64* theMarkIdx, uint* arrSize, int* result)
    {
        std::string bedFile( *bedFileCstr) ;

        printf("Parameter settting:  \n ");
        printf("\t bedFile = %s \n"           , *bedFileCstr);
        printf("\t nMarkers = %d \n"          , *nMarkers);
        printf("\t nSamples = %d \n"          , *nSamples);
        printf("\t nTargetMarkers = %d \n"    , *arrSize);
        

        getGeno (bedFile, *nSamples, *nMarkers, theMarkIdx, *arrSize,  result);
        return 0;
    }
}


