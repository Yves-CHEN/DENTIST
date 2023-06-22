
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
    vector<int64>  seqNo;
    vector<int>    bp;
    vector<uint>   include;
    vector<double> maf;
    vector<double> geneticDist;
    bool           moganMissing;
    string         bfileStr = "";
    long int       M;
    long int       N;
    inline long int size(){return (rs.size());};
    vector<double> calcMaf (string bfileName, long int N, long int M, uint ncpus);
    BedFile(string bfileName, float maf,  uint ncpus);
    BedFile(string bfileName);
    inline BedFile() {};
    inline  string print(uint i) const 
    {
        assert(i<include.size());
        int j = include[i];
        string output = (chrID[j]) + "\t" +(rs[j])  + "\t" + to_string(geneticDist[j]) + "\t"
            + to_string(bp[j])  + "\t" + (A1[j])  + "\t" +(A2[j])  ;
        return output;

    };



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
    int    ibuf = 0;
    string cbuf = "0";
    double dbuf = 0.0;
    string str_buf;
    int64  rowNum = 0;
    int    ncpus = 1;

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
    this->maf = calcMaf (bfileName, N, M, ncpus);
    for (size_t i = 0; i < this->bp.size(); i ++)
        this->include.push_back(i);


    

};



BedFile::BedFile(string  bfileName, float minMaf,  uint ncpus)
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
        p_bp = ibuf;
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
    uint dist = this->bp[ this->bp.size() -1] - this->bp[ 0];
    if(dist> 5e-6 && this->geneticDist[this->bp.size() -1] == 0 ) // warning mogan missing, if last - first > 10Mb, but the mogan is 0
    {
        cout << "[warning] mogan is like to be missing, but it will only be a problem when calling for it." << endl;
        moganMissing = true;
    }
    //this->M = this -> rs.size();
    this->M = rowNum;
    string famFile = bfileName + ".fam";
    ifstream fin(famFile.c_str());
    int N = 0;
    string line;
    while (getline(fin, line)) N++;
    fin.close();
    this->N = N;
    this->maf = calcMaf (bfileName, N,  this -> rs.size(), ncpus);
    decltype(this->maf) tmpMaf;
    for (size_t i = 0; i < this->maf.size(); i ++) 
    {
        if( (this->maf[i] <= 0.5 && this->maf[i] > minMaf) ||  (this->maf[i] > 0.5 && 1- this->maf[i] > minMaf)  )
        {
            this->include.push_back(i);
            tmpMaf.push_back(this->maf[i]);
        }
    }

    if(std::find(tmpMaf.begin(), tmpMaf.end(), 0) != tmpMaf.end())
        stop("[error]  A SNP(s) with maf of 0 is found.\n");
}




vector<double>   BedFile::calcMaf (string bfileName, long int N, long int M, uint ncpus)
{
    long int nSample = N;
    long int nMarker = this->M;
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
    typedef unsigned char dataType;

    uint processSamplePerRound = sizeof(dataType)*8 /2;
    //uint perMakerSize = ceil ( (nSample) / 1.0 / processSamplePerRound );
    int64 perMakerSizeOrig = ceil ( (nSample) / 1.0 / 4);
    //int64 perMakerSize = ceil( (nSample) / 1.0 / processSamplePerRound );
    int64 perMakerSize = perMakerSizeOrig;
    uint nBlanks   = ( processSamplePerRound - (nSample) %
                processSamplePerRound  ) % processSamplePerRound; //
    int64 lSize =0;
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
    if (bedFileReader== NULL)
        stop("[error] cannot write to [%s]", bedFile.c_str());
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

    if(lSize  != int64(perMakerSizeOrig) * nMarker + nThrowAway )
    {
        cout << "[error] The size of bedFile " << lSize << " is inconsistenty with the estimated " << perMakerSizeOrig * nMarker + nThrowAway 
            << " basd on " << perMakerSizeOrig << " samples and " << nMarker << "markers. \n";
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

    // segmenting
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

    assert(nBlanks == 0);
    vector<double>   maf (nMarker, -1);
    for (size_t block_i =0; block_i < readLen.size(); block_i ++)
    {

        int64 loadSize = perMakerSizeOrig * sizeof(uchar) * readLen[block_i] ;
        uchar* bufferAllMarkers = new uchar [loadSize ];

        auto orderIdx = seqNo[startingIdx[block_i]];
        fseek (bedFileReader , perMakerSizeOrig * sizeof(uchar) * (orderIdx) + nThrowAway, SEEK_SET );
        fread (bufferAllMarkers, 1, loadSize, bedFileReader);
        dataType** GENO = new dataType* [readLen[block_i]] ;
#pragma omp parallel for
        for (size_t i = 0; i < readLen[block_i]; i ++) {
            GENO[i] = (dataType*) (perMakerSizeOrig * i + bufferAllMarkers);
        }
        
#pragma omp parallel for
        for (size_t i = 0; i < readLen[block_i] ; i ++)
        {
            dataType* GENO_i = GENO[i] ;
            uint      nMissing = 0;
            double    E_i = 0;
            dataType  marker = 0  ;
            double    sum_i = 0;
            for (unsigned int k =0; k < perMakerSize; k ++)
            {
                marker    = markMissing[ GENO_i[k]  ]    ;
                nMissing += countOnes[(dataType)(~marker)]/2;
                sum_i    += countOnes[GENO_i [k] & marker];
            }
            // sum of alternative alleles / (2*sample size)
            maf[startingIdx[block_i] + i]  = (double(sum_i  ) / 2 ) / (nKeptSample - nMissing ); 
        }
        delete [] bufferAllMarkers;
        delete [] GENO;
    }
    ofstream fout (bfileName + ".DENTIST.maf.txt");
    for (size_t i = 0; i < maf.size(); i ++)
        fout <<  maf[i] << endl;
    fout.close();
    return maf;

};


class BLDFILE : public BedFile
{
    uint LDwindSize;
    uint elementInBypes;
public:
    BLDFILE (string filename);

};

BLDFILE::BLDFILE (string filename)
        :BedFile() 
{

    string idxFile = filename+ ".ridx";
    ifstream Bim(idxFile.c_str());

    if(!Bim.good()) printf( "[error] [%s] not found.\n", idxFile.c_str());
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

        int throwaway;
        Bim >> throwaway >> throwaway >> throwaway;

        // disregard the idx information for now
        //cin.ignore (std::numeric_limits<std::streamsize>::max(), '\n');
    }
    Bim.close();
    if(ibuf < p_bp)
    {
        cout << "[error] The BLD file should ordered by bp and it should be per-chromosome file than genomewide." << endl;
        assert( ibuf >= p_bp );
    }
    uint dist = this->bp[ this->bp.size() -1] - this->bp[ 0];
    if(dist> 5e-6 && this->geneticDist[this->bp.size() -1] == 0 ) // warning mogan missing, if last - first > 5Mb, but the mogan is 0
    {
        moganMissing = true;
        cout << "[warning] mogan is like to be missing, but it will only be a problem when calling for it." << endl;
    }

    M = this -> rs.size();
    // string famFile = bfileName + ".fam";
    // ifstream fin(famFile.c_str());
    // int N = 0;
    // string line;
    // while (getline(fin, line)) N++;
    // fin.close();
    // this->N = N;
    // this->maf = calcMaf (bfileName, N, M, ncpus);
    cout << "parsing BLD" << endl;

    for (size_t i = 0; i < this->bp.size(); i ++)
    {
        this->include.push_back(i);
    }


    FILE*  datFile      = fopen((filename+".bld").c_str(), "r");
    if (datFile == NULL)
        stop("[error] cannot write to [%s]", (filename+".bld").c_str());
    // ----------------------------------------------------------------------------
    // read header save in integer (4bytes)
    int headerInfo [5] = {}; // fileSubtype, N, M, LDwindSize, LD in number of Bytes
    fread (headerInfo,  sizeof(int), 5 ,datFile);

    fclose(datFile);
    this -> N = headerInfo[1];
    LDwindSize = headerInfo[3];
    elementInBypes = headerInfo[4];


};
void setChr(BedFile ref, string targetChrID)
{
    if(targetChrID == "")
    {
        printf("[info] Guessing the chrID.\n");
        targetChrID =  ref.chrID[0];
        for (size_t i = 0 ; i < ref.chrID.size(); i ++)
            if(ref.chrID[i] != targetChrID)
                stop("[error] More than one chromosome is found. Please specify --chrID");
    }
    cout <<  "[info] chrID == " << targetChrID << endl;
    BedFile tmpBed ;
    decltype(ref.include) tmpInclude ;
    for (size_t i = 0 ; i < ref.chrID.size(); i ++)
    {
        if(ref.chrID[i] == targetChrID)
            tmpInclude.push_back(ref.include[i]);
    }
    ref.include = tmpInclude;

}
#endif  // __BFILE__

