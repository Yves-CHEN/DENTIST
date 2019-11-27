
#ifndef __BFILE__
#define __BFILE__

#include "headers.h"
#include <bitset>

typedef  unsigned char  uchar ;

class BedFile
{
public:
    vector<string> A1;
    vector<string> A2;
    vector<string> chrID;
    vector<string> rs;
    vector<long int> seqNo;
    vector<int>    bp;
    vector<uint>   include;
    vector<double> maf;
    vector<double> geneticDist;
    bool           moganMissing;
    string         bfileStr = "";
    long int M;
    long int N;
    inline long int size(){return (rs.size());};
    vector<double> calcMaf (string bfileName, long int N, long int M, uint ncpus);
    BedFile(string bfileName, float maf, uint ncpus);
    BedFile(string bfileName);



};
///************************************************************************* 
//  Todo:
//     1. warning about the presence of mogan distance
//     2. Check if per-chr bed file is used, than genomewide.
//        Also the bp are ordered
///*************************************************************************
BedFile::BedFile(string  bfileName ) 
    :bfileStr(bfileName)
{

    string bimFile = bfileName + ".bim";
    ifstream Bim(bimFile.c_str());

    if(!Bim.good()) printf( "[error] [%s] not found.\n", bimFile.c_str());
    int ibuf = 0;
    string cbuf = "0";
    double dbuf = 0.0;
    string str_buf;
    long int rowNum = 0;
    int ncpus = 1;

    moganMissing = false;
    int p_bp = -1;
    while (Bim) 
    {
        Bim >> str_buf;
        if (Bim.eof()) break;
        this->chrID.push_back( str_buf);
        Bim >> str_buf;
        this->rs.push_back(str_buf);
        Bim >> dbuf;
        this->geneticDist.push_back(dbuf);
        Bim >> ibuf;
        this->bp.push_back(ibuf);
        if(ibuf >= p_bp) // ordered bp is expected
            p_bp = ibuf;
        else
            break;
        Bim >> cbuf;
        to_upper(cbuf);
        this->A1.push_back(cbuf.c_str());
        Bim >> cbuf;
        to_upper(cbuf);
        this->A2.push_back(cbuf.c_str());
        this->seqNo.push_back(rowNum);
        rowNum  ++;
    }
    Bim.close();
    if(ibuf < p_bp)
    {
        cout << "[error] The bed file should ordered by bp and it should be per-chromosome file than genomewide." << endl;
        assert( ibuf >= p_bp );
    }
    uint dist = this->bp[ this->bp.size() -1] - this->bp[ 0];
    if(dist> 5e-6 && this->geneticDist[this->bp.size() -1] == 0 ) // warning mogan missing, if last - first > 5Mb, but the mogan is 0
    {
        moganMissing = true;
        cout << "[warning] mogan is like to be missing, but it will only be a problem when calling for it." << endl;
    }

    M = this -> rs.size();
    string famFile = bfileName + ".fam";
    ifstream fin(famFile.c_str());
    int N = 0;
    string line;
    while (getline(fin, line)) N++;
    fin.close();
    this->N = N;
    for (uint i = 0; i < this->maf.size(); i ++)
        this->include.push_back(i);

    this->maf = calcMaf (bfileName, N, M, ncpus);


};

BedFile::BedFile(string  bfileName, float minMaf, uint ncpus )
:bfileStr(bfileName)
{
    string bimFile = bfileName + ".bim";
    ifstream Bim(bimFile.c_str());
    int ibuf = 0;
    string cbuf = "0";
    double dbuf = 0.0;
    string str_buf;
    long int rowNum = 0;
    moganMissing = false;
    int p_bp =-1;
    while (Bim) 
    {
        Bim >> str_buf;
        if (Bim.eof()) break;
        this->chrID.push_back( str_buf);
        Bim >> str_buf;
        this->rs.push_back(str_buf);
        Bim >> dbuf;
        this->geneticDist.push_back(dbuf);
        Bim >> ibuf;
        this->bp.push_back(ibuf);
        if(ibuf >= p_bp) // ordered bp is expected
            p_bp = ibuf;
        else
            break;

        Bim >> cbuf;
        to_upper(cbuf);
        this->A1.push_back(cbuf.c_str());
        Bim >> cbuf;
        to_upper(cbuf);
        this->A2.push_back(cbuf.c_str());
        this->seqNo.push_back(rowNum);
        rowNum  ++;
    }
    Bim.close();
    if(ibuf < p_bp)
    {
        cout << "[error] The bed file should ordered by bp and it should be per-chromosome file than genomewide." << endl;
        exit(-1);
    }

    uint dist = this->bp[ this->bp.size() -1] - this->bp[ 0];
    if(dist> 5e-6 && this->geneticDist[this->bp.size() -1] == 0 ) // warning mogan missing, if last - first > 10Mb, but the mogan is 0
    {
        cout << "[warning] mogan is like to be missing, but it will only be a problem when calling for it." << endl;
        moganMissing = true;
    }


    M = this -> rs.size();
    string famFile = bfileName + ".fam";
    ifstream fin(famFile.c_str());
    int N = 0;
    string line;
    while (getline(fin, line)) N++;
    fin.close();
    this->N = N;
    this->maf = calcMaf (bfileName, N, M, ncpus);

    for (uint i = 0; i < this->maf.size(); i ++)
        if( (this->maf[i] <= 0.5 && this->maf[i] > minMaf) ||  (this->maf[i] > 0.5 && 1- this->maf[i] > minMaf)  )
            this->include.push_back(i);


}


vector<double>   BedFile::calcMaf (string bfileName, long int N, long int M, uint ncpus)
{
    long int nSample = N;
    long int nMarker = M;
    long int sizeOfMarkIdx = M;
    string bedFile = bfileName + ".bed";
    string bimFile = bfileName + ".bim";
    // **************************************************************
    ///                    set timer
    // **************************************************************
    struct timespec start, finish;
    double elapsed;
    // **************************************************************
    ///                    set number of cpus
    // **************************************************************
    int nProcessors = omp_get_max_threads();
    if(ncpus < nProcessors) nProcessors = ncpus;
    omp_set_num_threads( nProcessors );
        //omp_set_num_threads( 1);
    printf("[info] Calculating frequencies with %d cpus \n", nProcessors);
    // **************************************************************
    //                     Variables
    // **************************************************************
    typedef unsigned short dataType;
    uint processSamplePerRound = sizeof(dataType)*8 /2;
    //uint perMakerSize = ceil ( (nSample) / 1.0 / processSamplePerRound );
    uint perMakerSizeOrig = ceil ( (nSample) / 1.0 / 4);
    uint perMakerSize = int( (nSample) / 1.0 / processSamplePerRound );
    uint nBlanks   = ( processSamplePerRound - (nSample) %
                processSamplePerRound  ) % processSamplePerRound; //
    long lSize =0;
    /// headerCoding is 9 bytes in total for plink 1.9 bedfile format, 3 bytes
    //in older version.
    int const nByteHeader = 9;
    uchar headerCoding[nByteHeader +1 ] = { 0x6c, 0x1b, 0x01, 0xdc, 0x0f, 0xe7,
                                            0x0f, 0x6b, 0x01};
    int const nByteHeader_older = 3;
    std::string formatVersion = "";
    int nThrowAway = 0;
    uchar  headerBuf[nByteHeader ] ;
    //uint nKeptSample = nSample;
    uint nKeptSample = perMakerSize * processSamplePerRound;
    std::vector<double>   GENO_VAR (sizeOfMarkIdx, -1 );
    std::vector<double>   GENO_Ex  (sizeOfMarkIdx, -1 );
    std::vector<double>   GENO_sum11  (sizeOfMarkIdx, -1 ); // sum of sample with genotype 11

    // **************************************************************
    //                1. validation of the bfile version
    // 2. validation of the bed file size, given the number of markers and  samples
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

    if(!memcmp(headerCoding, headerBuf, nByteHeader)) 
        {printf("[info] This bed file is plink 1.9 bedfile format. (Newer) \n"); 
        formatVersion="1.9"; nThrowAway = nByteHeader;};
    if(!memcmp(headerCoding, headerBuf, nByteHeader_older)) 
        {printf("[info] This bed file is plink 1.0 bedfile format. (Older)\n");
         formatVersion="1.0"; nThrowAway = nByteHeader_older;};
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
    ///
    size_t const sizeOfMap = (size_t) (pow( 2 , (sizeof(dataType) * 8) ));
    // invCountOnes is a "function (aa)  aa = ~aa; bitset<8> (aa).count"
    // mapper2   is a "function (aa) aa = ~aa; aa & (aa<<1) & 0xaa" 
    // mapper3   is a "function (aa)  aa & (aa<<1) & 0xaa" 
    // countOnes is a "function (aa)  bitset<8> (aa).count"
    // mapper5   is a "function (aa)  aa = ~aa ; ( ((aa ^ aa <<1) & 0xaa ) >>1
    // ) *3"
    uint      invCountOnes  [sizeOfMap ] = {0};
    uint      mapper2[sizeOfMap ]        = {0};
    uint      mapper3[sizeOfMap ]        = {0};
    uint      countOnes[sizeOfMap ]      = {0};
    dataType  mapper5[sizeOfMap ]        = {0};
    dataType  markMissing[sizeOfMap ]    = {0};
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



    const uint maxBlockSize = 80000000;  // 80M
    vector <uint> startingIdx; 
    vector <uint> readLen; 
    assert(maxBlockSize>nSample *2);
    uint oneUnit = uint( maxBlockSize / nSample);
    int counter = 0;
    while( sizeOfMarkIdx >  (oneUnit) ) 
    {
        startingIdx.push_back (counter);
        readLen.push_back(oneUnit);
        counter  +=  oneUnit;
        sizeOfMarkIdx -= oneUnit;
    }
    if(sizeOfMarkIdx > 0)
    {
        startingIdx.push_back (counter);
        readLen.push_back(sizeOfMarkIdx);
    }


//    for (int i = 0; i < readLen.size(); i ++ )
//        cout << readLen[i] << "\t";
//    cout << endl;
//
//    for (int i = 0; i < readLen.size(); i ++ )
//        cout << startingIdx[i] << "\t";
//    cout << endl;

    vector<double>   maf (nMarker, -1);
    for (uint block_i =0; block_i < readLen.size(); block_i ++)
    {

        unsigned long loadSize = perMakerSizeOrig * sizeof(uchar) * readLen[block_i] ;
        uchar* bufferAllMarkers = new uchar [loadSize ];
        fseek (bedFileReader , perMakerSizeOrig * sizeof(uchar) * (startingIdx[block_i]) + nThrowAway, SEEK_SET );
        fread (bufferAllMarkers, 1, loadSize, bedFileReader);
        dataType*  bufferMaker = NULL;  // This is pointer to the memory of a particular marker.
        dataType** GENO = new dataType* [readLen[block_i]] ;
#pragma omp parallel for
        for (unsigned int i = 0; i < readLen[block_i]; i ++)
        {
            GENO[i] = (dataType*) (perMakerSizeOrig * i + bufferAllMarkers);
        }
#pragma omp parallel for
        for (unsigned int i = 0; i < readLen[block_i] ; i ++)
        {
            dataType* GENO_i = GENO[i] ;
            uint      nMissing = 0;
            double    E_i = 0;
            dataType  marker = 0  ;
            double    sum_i = 0;
            for (unsigned int k =0; k < perMakerSize; k ++)
            {
                marker    = markMissing[ GENO_i[k]  ]  & markMissing[GENO_i[k]]   ;
                nMissing += countOnes[(dataType)(~marker)];
                sum_i    += countOnes[GENO_i [k] & marker];
            }
            E_i         = double(sum_i  ) / (nKeptSample - nMissing);
            maf[startingIdx[block_i] + i] = 1- E_i / 2;
        }



        delete [] bufferAllMarkers;
        delete [] GENO;

    }

    return maf;

};



#endif  // __BFILE__

