#ifndef __DENTIST__

#define __DENTIST__

#define MAX_LINE_SIZE 0x10000

#include "invoker.h"
#include "headers.h"
//#include "gzstream/gzstream.h"
#include <bitset>
#include <numeric>

#include <boost/math/distributions/inverse_chi_squared.hpp>
#include <assert.h>
#include "utils.h"
#include <zlib.h>



using namespace std;


//#void runQCForEQTL(string bfile, SMRWK& smrwk, map<string,long>& seqNoMap, int thread_num, SMRWK& smrwk_tmp, double Degree_QC);


class GWAS
{
public:
    vector<string> A1;
    vector<string> A2;
    vector<double> b;
    vector<double> se;
    vector<double> zscore;
    vector<string> rs;
    vector<double> maf;
    vector<int>    _include;
    vector<int>    splSize;
    vector<double> pValue;
    vector<long int>  aligned;
    inline GWAS (){};
    GWAS (string summmaryFile );
    GWAS (string summmaryFile, bool ifGZ );
    long int M;
    long int N;
    inline long int size(){return (rs.size());};
};

GWAS::GWAS (string summmaryFile, bool ifGZ)
{
    bool warnnullfreq=false;
    //igzstream gwasFile (summmaryFile.c_str());
    cout << "Reading GWAS summary data from [" + summmaryFile  + "]." << endl;
    int lineNum = 0;
    //string line;
    //getline(gwasFile, line);
    //const char* buf = line.c_str();


    gzFile  file = gzopen(summmaryFile.c_str(), "r");
    if (file == Z_NULL) {
        fprintf(stderr, "gzopen error: not a gzfile. \n");
        exit(EXIT_FAILURE); 
    }
    const uint maxSize = maxSummaryRowSize;
    char buf[maxSummaryRowSize ] = "";
    printf("Reserving %d M memory for reading. \n", maxSummaryRowSize  / 1000000);
    // assumming there is a header and the header can be discarded.
    if(gzgets(file, buf, maxSize) == NULL)
    {
        printf("ERROR: the first row of the file %s is empty.\n",summmaryFile.c_str());
        exit(EXIT_FAILURE);
    }



    if(buf[0]=='\0')
    {
        printf("ERROR: the first row of the file %s is empty.\n", summmaryFile.c_str());
        exit(EXIT_FAILURE);
    }
    vector<string> vs_buf;
    split_string(buf, vs_buf, " \t\n");
    to_upper(vs_buf[0]);
    if(vs_buf[0]!="SNP") {
        printf("ERROR: %s should have headers that start with \"SNP\" rather than %s.\n", summmaryFile.c_str(), vs_buf);
        exit(EXIT_FAILURE);
    }


    while(gzgets(file, buf, maxSize))
    {
        if(buf[0]!='\0'){
            vs_buf.clear();
            int col_num = split_string(buf, vs_buf, " \t\n");
            if(col_num!=8) {
                printf("ERROR: column number is not correct in row %d!\n", lineNum+2);
                exit(EXIT_FAILURE);
            }
            if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
                printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
                exit(EXIT_FAILURE);
            }
            this->rs.push_back(vs_buf[0]);
            if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
                printf("ERROR: allele1 is \'NA\' in row %d.\n", lineNum+2);
                exit(EXIT_FAILURE);
            }
            to_upper(vs_buf[1]);
            this->A1.push_back(vs_buf[1]);
            
            if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
                printf("ERROR: allele2 is \'NA\' in row %d.\n", lineNum+2);
                exit(EXIT_FAILURE);
            }
            to_upper(vs_buf[2]);
            this->A2.push_back(vs_buf[2]);
            
            if(vs_buf[3]=="NA" || vs_buf[3]=="na")
            {
                if(!warnnullfreq){
                    warnnullfreq=true;
                    printf("WARNING: frequency is \'NA\' in one or more rows.\n");
                }
                this->maf.push_back(-9);
                
            }
            else {
                this->maf.push_back(atof(vs_buf[3].c_str()));
            }
           
            
            if(vs_buf[4]=="NA" || vs_buf[4]=="na"){
                printf("WARNING: effect size is \'NA\' in row %d.\n", lineNum+2);
                this->b.push_back(0);
            } else {
                this->b.push_back(atof(vs_buf[4].c_str()));
            }
            if(vs_buf[5]=="NA" || vs_buf[5]=="na"){
                printf("WARNING: standard error is \'NA\' in row %d.\n", lineNum+2);
                this->se.push_back(-9);
            } else {
                this->se.push_back(atof(vs_buf[5].c_str()));
            }

            
            this->zscore.push_back( this->b [b.size() -1]/ this->se [se.size() -1] );
            this->pValue.push_back(atof(vs_buf[6].c_str()));
            this->splSize.push_back(atoi(vs_buf[7].c_str()));
            this->_include.push_back(lineNum);
            lineNum++;
        }
    }
    this->M = this->_include.size();
    cout <<"GWAS summary data of "<< this->M << " SNPs to be included from [" + string( summmaryFile) + "]." << endl;
    //gwasFile.close();


    gzclose(file);

};
/// GWAS::GWAS (string summmaryFile)
/// {
///     bool warnnullfreq=false;
///     ifstream gwasFile (summmaryFile);
///     if (!gwasFile.is_open())
///     {
///         //fprintf (stderr, "%s: Couldn't open file %s\n", summmaryFile, strerror (errno));
///         exit (EXIT_FAILURE);
///     }
///     cout << "Reading GWAS summary data from [" + summmaryFile  + "]." << endl;
///     int lineNum = 0;
///     char buf[MAX_LINE_SIZE];
///     gwasFile.getline(buf,MAX_LINE_SIZE);// the header
/// 
///    
/// 
/// 
///     if(buf[0]=='\0')
///     {
///         printf("ERROR: the first row of the file %s is empty.\n", summmaryFile.c_str());
///         exit(EXIT_FAILURE);
///     }
///     vector<string> vs_buf;
///     split_string(buf, vs_buf, " \t\n");
///     to_upper(vs_buf[0]);
///     if(vs_buf[0]!="SNP") {
///         printf("ERROR: %s should have headers that start with \"SNP\" rather than %s.\n", summmaryFile.c_str(), vs_buf);
///         exit(EXIT_FAILURE);
///     }
/// 
/// 
///     while(!gwasFile.eof())
///     {
///         gwasFile.getline(buf,MAX_LINE_SIZE);
///         
///         if(buf[0]!='\0'){
///             vs_buf.clear();
///             int col_num = split_string(buf, vs_buf, " \t\n");
///             if(col_num!=8) {
///                 printf("ERROR: column number is not correct in row %d!\n", lineNum+2);
///                 exit(EXIT_FAILURE);
///             }
///             if(vs_buf[0]=="NA" || vs_buf[0]=="na"){
///                 printf("ERROR: the SNP name is \'NA\' in row %d.\n", lineNum+2);
///                 exit(EXIT_FAILURE);
///             }
///             this->rs.push_back(vs_buf[0]);
///             if(vs_buf[1]=="NA" || vs_buf[1]=="na"){
///                 printf("ERROR: allele1 is \'NA\' in row %d.\n", lineNum+2);
///                 exit(EXIT_FAILURE);
///             }
///             to_upper(vs_buf[1]);
///             this->A1.push_back(vs_buf[1]);
///             
///             if(vs_buf[2]=="NA" || vs_buf[2]=="na"){
///                 printf("ERROR: allele2 is \'NA\' in row %d.\n", lineNum+2);
///                 exit(EXIT_FAILURE);
///             }
///             to_upper(vs_buf[2]);
///             this->A2.push_back(vs_buf[2]);
///             
///             if(vs_buf[3]=="NA" || vs_buf[3]=="na")
///             {
///                 if(!warnnullfreq){
///                     warnnullfreq=true;
///                     printf("WARNING: frequency is \'NA\' in one or more rows.\n");
///                 }
///                 this->maf.push_back(-9);
///                 
///             }
///             else {
///                 this->maf.push_back(atof(vs_buf[3].c_str()));
///             }
///            
///             
///             if(vs_buf[4]=="NA" || vs_buf[4]=="na"){
///                 printf("WARNING: effect size is \'NA\' in row %d.\n", lineNum+2);
///                 this->b.push_back(0);
///             } else {
///                 this->b.push_back(atof(vs_buf[4].c_str()));
///             }
///             if(vs_buf[5]=="NA" || vs_buf[5]=="na"){
///                 printf("WARNING: standard error is \'NA\' in row %d.\n", lineNum+2);
///                 this->se.push_back(-9);
///             } else {
///                 this->se.push_back(atof(vs_buf[5].c_str()));
///             }
/// 
///             
///             this->zscore.push_back( this->b [b.size() -1]/ this->se [se.size() -1] );
///             this->pValue.push_back(atof(vs_buf[6].c_str()));
///             this->splSize.push_back(atoi(vs_buf[7].c_str()));
///             this->_include.push_back(lineNum);
///             lineNum++;
///         }
///     }
///     this->M = this->_include.size();
///     cout <<"GWAS summary data of "<< this->M << " SNPs to be included from [" + string( summmaryFile) + "]." << endl;
///     gwasFile.close();
/// };


class Options
{
public:

    string summmaryFile;
    string bfileName   ;
    string outPrefix   ;
    int    maxDim      ;
    int    minDim      ;
    int    maxDist;
    int    thread_num  ;
    double Degree_QC   ;
    double lambda      ;
    string targetSNP   ;
    string extractFile;
    bool   ignoreWarnings ;
    int    withNA;
    //const char *flgs[] ;
    vector<string> flags;
    float mafThresh     ;
    int   debugMode;
    /// to be implemented
    //  double dupThresh ;  the thresh of LD r^2 between two SNPs to be considered duplicates.
    //  int nIterations; the number interations to be performed. At least one.
    //  int qcMethod;    1) use zscore diff 2) zscore_diff/se
    //  int distMethod;  1) by number of snps, 2) by bp  3)by morgan distance.
    //
    //
    //
    static inline bool not_in_flags(vector<string> &flags, string str) { 
        return find(flags.begin(),flags.end(),str) == flags.end(); 
    };
    static inline bool has_prefix(const string &str, const string &prefix)
    {
        return str.size() >= prefix.size() &&
                str.compare(0, prefix.size(), prefix) == 0;
    };


    static inline void FileExist(string filename)
    {
        ifstream ifile(filename.c_str());
        if(!ifile) throw("Error: can not open the file ["+filename+"] to read.");
    };



    static inline void bool_FLAG_VALID_CK(string str, const char* flag)
    {
        if( !(flag==NULL || has_prefix(flag, "--"))  )
        {
            // fprintf (stderr, "Please verify the flag %s!: \n",
            //          str.c_str());
            exit (EXIT_FAILURE);
        }


    }
    static inline void FLAG_VALID_CK(string str, const char* flag)
    {
        if(flag==NULL || has_prefix(flag, "--"))
        {
            // fprintf (stderr, "Please verify the flag %s!: \n",
            //          str.c_str());
            exit (EXIT_FAILURE);
        }
    };



    void FLAGS_VALID_CK(int option_num, char* option_str[]);
    


    inline Options()
    {
        summmaryFile  = ""; 
        flags.push_back("--gwas-summary"); //summmaryFile
        bfileName     = "";
        flags.push_back("--bfile");  // bed file
        outPrefix     = "out";  // default output prefix
        flags.push_back("--out");
        thread_num    = 1;      // number of threads for QC
        flags.push_back("--thread-num");
        Degree_QC     = 0.3;    // percentage of probes to be filtered
        // flags.push_back("--degree-of-QC");
        lambda        = 0.1;    // the lambda for Ridge regression 
        // flags.push_back("--lambda");
        maxDim        = 25000;   // default max window size for imputation
        // flags.push_back("--min-wind");
        mafThresh     = -1;   // default -1 for no restrictions on maf.
        flags.push_back("--maf");
        extractFile = "";
        flags.push_back("--extract");
        targetSNP = "";
        flags.push_back("--target");
        withNA = 0;
        flags.push_back("--with-NA-geno");
        maxDist = -1;
        flags.push_back("--wind-dist");
        minDim        = 2000;   // default min window size for imputation
        flags.push_back("--wind");
        debugMode  = 0;
        flags.push_back("--debug");
        ignoreWarnings = false;
        flags.push_back("--ignore-warnings");




// flags.push_back("--distance");
// flags.push_back("--use-pvalue");

    }


    void parseOptions(int nArgs, char* option_str[]);



};


void Options::FLAGS_VALID_CK(int option_num, char* option_str[])
{
    if(option_num<3)
    {
        cout<<"Flags include:"<<endl;
        int cur_mark=0;
        for(int i=0;i<flags.size();i++)
        {
            int tmp=i>>2;
            if(tmp>cur_mark)
            {
                cout<<endl;
                cur_mark=tmp;
            }
            cout<<flags[i]<<",";
        }
        cout<<endl;
        exit (EXIT_FAILURE);
    }
    for(int i=0;i<option_num;i++)
    {
        if(has_prefix(option_str[i],"--"))
            if(not_in_flags(flags, option_str[i]))
            {
                //fprintf (stderr, "%s: Invalid option\n", option_str[i]);
                exit (EXIT_FAILURE);
            }
    }

}

void Options::parseOptions(int nArgs, char* option_str[])
{
    FLAGS_VALID_CK(nArgs, option_str);
    for(int i = 0; i < nArgs; i ++)
    {
        if(strcmp(option_str[i],"--gwas-summary")==0)
        {
            summmaryFile = option_str[++i];
            Options::FLAG_VALID_CK(string("--gwas-summmary"), summmaryFile.c_str());
            cout<< option_str[i-1] << " " <<  summmaryFile <<endl;
            FileExist(summmaryFile);
        }


        if(strcmp(option_str[i], "--bfile") == 0)
        {
            bfileName = option_str[++i];
            Options::FLAG_VALID_CK(string("--bfile"), bfileName.c_str());
            cout<< option_str[i-1] << " " <<  bfileName <<endl;
            FileExist(bfileName + ".bed");
            FileExist(bfileName + ".bim");
            FileExist(bfileName + ".fam");
        }
        if(strcmp(option_str[i], "--out") == 0)
        {
            outPrefix = option_str[++i];
            Options::FLAG_VALID_CK(string("--out"), outPrefix.c_str());
            cout<< option_str[i-1] << " "<< outPrefix <<endl;
            // CommFunc::FileExist(oproblstName);
        }
        if(strcmp(option_str[i], "--thread-num") == 0){
            thread_num = atoi(option_str[++i]);
            cout << option_str[i-1] << " " << thread_num<< endl;
            if(thread_num < 0)
            {
                fprintf (stderr, "Error: --thread-num should be over 0.\n");
                exit (EXIT_FAILURE);
            }
        }

        if(strcmp(option_str[i], "--degree-of-QC") == 0){
            Degree_QC = atof(option_str[++i]);
            cout << option_str[i-1] << " " << Degree_QC  << endl;
            if(Degree_QC  < 0 || Degree_QC >= 1) 
            {
                fprintf (stderr, "Error: --degree-of-QC should be between 0 and 1.\n");
                exit (EXIT_FAILURE);
            }
        }
        if(strcmp(option_str[i], "--lambda") == 0){
            lambda = atof(option_str[++i]);
            cout << option_str[i-1] << " " << lambda << endl;
            if(Degree_QC  < 0 || Degree_QC >= 0.1) 
            {
                fprintf (stderr, "Error: --lambda should be between 0 and 0.1.\n");
                exit (EXIT_FAILURE);
            }
        }
        if(strcmp(option_str[i], "--maf") == 0)
        {
            mafThresh = atof ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--maf"), bfileName.c_str());
            cout<< option_str[i-1] << " " <<  mafThresh <<endl;
        }


        if(strcmp(option_str[i], "--wind-dist") == 0)
        {
            maxDist = atof ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--wind-dist"), bfileName.c_str());
            cout<< option_str[i-1] << " " <<  maxDist <<endl;
        }

        if(strcmp(option_str[i], "--wind") == 0)
        {
            maxDim = atof ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--wind"), bfileName.c_str());
            cout<< option_str[i-1] << " " <<  maxDim <<endl;
        }

        if(strcmp(option_str[i], "--target") == 0)
        {
            targetSNP = ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--target"), targetSNP.c_str());
            cout<< option_str[i-1] << " " <<  targetSNP <<endl;
        }
        if(strcmp(option_str[i], "--extract") == 0)
        {
            extractFile = ( option_str[++i] );
            Options::FLAG_VALID_CK(string("--extract"), extractFile.c_str());
            cout<< option_str[i-1] << " " <<  extractFile <<endl;
        }

        if(strcmp(option_str[i], "--ignore-warnings") == 0)
        {
            ignoreWarnings = true;
            if(i+1 < nArgs)
            {
                Options::bool_FLAG_VALID_CK(string("--ignore-warnings"), option_str[i+1]);
            }
            cout<< option_str[i] << " " <<  " TRUE" <<endl;
        }
        if(strcmp(option_str[i], "--with-NA-geno") == 0)
        {
            withNA  = 1;
            if(i+1 < nArgs)
            {
                Options::bool_FLAG_VALID_CK(string("--with-NA-geno"), option_str[i+1]);
            }
            cout<< option_str[i] << " " <<  " TRUE" <<endl;
        }

        if(strcmp(option_str[i], "--debug") == 0)
        {
            debugMode = 1;
            if(i+1 < nArgs)
            {
                Options::bool_FLAG_VALID_CK(string("--debug"), option_str[i+1]);
            }
            cout<< option_str[i] << " " <<  " TRUE" <<endl;
        }



            

    }

}




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
    long int M;
    long int N;
    inline long int size(){return (rs.size());};
    vector<double> calcMaf (string bfileName, long int N, long int M, uint ncpus);
    BedFile(string bfileName, float maf, uint ncpus);


};


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



    uint maxBlockSize = 80000000;  // 80M
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
            maf[startingIdx[block_i] + i] = E_i / 2;
        }



        delete [] bufferAllMarkers;
        delete [] GENO;

    }

    return maf;






};


BedFile::BedFile(string  bfileName, float minMaf, uint ncpus )
{
    string bimFile = bfileName + ".bim";
    ifstream Bim(bimFile.c_str());
    int ibuf = 0;
    string cbuf = "0";
    double dbuf = 0.0;
    string str_buf;
    long int rowNum = 0;
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
    {
        if( (this->maf[i] <= 0.5 && this->maf[i] > minMaf) ||  (this->maf[i] > 0.5 && 1- this->maf[i] > minMaf)  )
        {
            this->include.push_back(i);
        }

    }



}
uint moveKeep(double* LD, uint arrSize, uint currentDim, uint keepFromIdx)
{
    uint m = 0;
    for (uint i = keepFromIdx; i < arrSize ; i ++, m ++ )
        for (uint j = keepFromIdx, n = 0; j < arrSize ; j ++, n ++ )
            LD[m *  currentDim + n ] = LD[i * arrSize + j];
    return m;
}



// this different from moveKeep by creating a tmp matrix store the previous dat
//  before coping dat to LD mat. This avoid reading and writing from the same
//  LD mat, which can lead to problems when arrSize < currentDim.
template<class T>
uint moveKeepProtect(T* LD, uint arrSize, uint currentDim, uint keepFromIdx)
{
    if(arrSize  - keepFromIdx) return 0;
    uint tmp_dim = arrSize  - keepFromIdx;
    double* LD_tmp = new double [tmp_dim * tmp_dim];
    for (uint i = keepFromIdx, m=0; i < arrSize ; i ++, m ++ )
        for (uint j = keepFromIdx, n = 0; j < arrSize ; j ++, n ++ )
            LD_tmp[m *  tmp_dim + n ] = LD[i * arrSize + j];



    uint m = 0;
    for (m = 0; m < tmp_dim; m ++ )
        for (uint n = 0; n < tmp_dim;  n ++ )
            LD[m *  currentDim + n ] = LD_tmp[m * tmp_dim + n];


    delete[] LD_tmp;
    return m;
}
double logPvalueChisq(double stat)
{
    boost::math::inverse_chi_squared_distribution<double> mydist(1);
    double p = boost::math::cdf(mydist,1/(stat));
    return ( -log10(p) ) ;
}

void segmentedQCed_dist (string bfileName, string qcFile, long int nSamples, long int nMarkers,
        vector<string>& rsIDs,vector<uint>& bp, vector<double>& zScores, vector<long int>& seqNos,
        vector<bool>& flipped, const Options& opt)
{

    
    int    maxDim      = opt.maxDim; 
    int    minDim      = opt.minDim;
    int    thread_num  = opt.thread_num;
    double Degree_QC   = opt.Degree_QC;
    int    doDebug     = opt.debugMode;
    double lambda      = opt.lambda;
    string bedFile = bfileName + ".bed";
    vector<int> toKeep;
    int cutoff = maxDim;
    int maxBlockSize  = opt.maxDim;
    int minBlockSize  = 2000;
    cutoff = opt.maxDist;

    // The computation should be on large enough region.
    assert(zScores.size()  - minDim > 0);
    
    vector<double>  imputed (zScores.size(), 0);
    vector<double>  rsq     (zScores.size(), 0);
    vector<double>  zscore_e(zScores.size(), 0);
    vector<bool>    ifDup   (zScores.size(), false);


    // The most distant neighbor for i, given the distance cutoff.
    vector<int> nextIdx;
    for (uint i =1; i < bp.size(); i ++ )
        if(bp[i] - bp[0] >= cutoff) {nextIdx.push_back(i); break;}
    for (uint i =1; i < bp.size() ; i ++ ) {
        uint j = nextIdx[nextIdx.size() -1];
        while( bp[j] - bp[i] < cutoff) {
            if(j >= bp.size()-1 ) break;
            j ++;
        }
        nextIdx.push_back(j);
    }
    // A quater of distance cutoff away from i.
    vector<int> quaterIdx;
    for (uint i =1; i < bp.size(); i ++ )
        if(bp[i] - bp[0] >= cutoff/4) { quaterIdx.push_back(i); break;}
    for (uint i =1; i < bp.size() ; i ++ ) {
        uint j = quaterIdx[quaterIdx.size() -1];
        while( bp[j] - bp[i] < cutoff/4) {
            if(j >= bp.size()-1 ) break;
            j ++;
        }
        quaterIdx.push_back(j);
    }



    
    // ************************************************************************
    // Finding all the gaps greater than 1 Mb in size including the centromere
    // The starting and ending of a chromosome are counted as two gaps.
    //-------------------------------------------
    vector<long int> diff;
    diff.resize(bp.size());
    adjacent_difference (bp.begin(), bp.end(), diff.begin());

    vector<uint> allGaps;
    allGaps.push_back(0);
    uint gapSizeThresh =  1000000;
    if(gapSizeThresh > cutoff) gapSizeThresh = cutoff;
    for (uint i = 1; i < diff.size(); i ++)  
        if(diff[i] > gapSizeThresh) allGaps.push_back(i);
    allGaps.push_back(diff.size());
    // **********************************************************************
    // Divide into ranges
    vector<uint> startList;
    vector<uint> endList;


    vector<uint> fillStartList;
    vector<uint> fillEndList;

    for (uint k = 0; k < allGaps.size() -1; k ++)
    {
        uint rangeSize = allGaps[k+1];
        if(rangeSize < minBlockSize /2) continue;
        if(rangeSize  - minDim < 0) continue;
        int  startIdx = allGaps[k], endIdx = allGaps[k+1]; 
        endIdx   =  quaterIdx[quaterIdx[quaterIdx[quaterIdx[startIdx] ]] ] ;
        //endIdx   =  nextIdx[startIdx] ;
        int notStartInterval = 0, notLastInterval = 1;
        do
        {
            startList.push_back(startIdx), endList.push_back(endIdx);
            if(rangeSize <= endIdx ) notLastInterval = 0;
            //startIdx =  ceil( (endIdx + startIdx)/2.0 );
            fillStartList.push_back(quaterIdx[startIdx] ); // one quater of dist cutoff away
            fillEndList.push_back(quaterIdx[quaterIdx[quaterIdx[startIdx] ]] ); // three quaters of dist cutoff away
            D(printf("%d %d %d %d fill:%d\n", rangeSize , 
                    startIdx, endIdx, fillStartList[fillStartList.size()-1], fillEndList[fillEndList.size()-1]););
            startIdx =  quaterIdx[quaterIdx[startIdx] ];
            //endIdx   =  nextIdx[startIdx] ;
            endIdx   =  quaterIdx[quaterIdx[quaterIdx[quaterIdx[startIdx] ]] ] ;
            //if(endIdx - startIdx > maxBlockSize) endIdx = startIdx + maxBlockSize;
            endIdx   =  rangeSize -  minDim > endIdx ? endIdx: rangeSize;  // This is for the last region.
            notStartInterval  =1;
        }
        while ( notLastInterval );
    }
    fillStartList [0] = startList[0];
    fillEndList [fillEndList.size()-1] = endList[endList.size() -1];
    vector<double>  zScores_tmp;
    vector<long>     seqNos_tmp;
    vector<string>    rsIDs_tmp;
    vector<string>       A1_tmp;
    vector<int>         toAvert;
    uint pre_end   = 0;
    uint pre_start = 0;
    uint preDim    = 0;


    LDType* LD = new LDType[ (opt.maxDim + opt.minDim) *  (opt.maxDim + opt.minDim) ]();

    for (uint k = 0; k < startList.size(); k ++ )
    {
        uint startIdx = startList[k];
        uint endIdx   =   endList[k];
        uint fillStartIdx = fillStartList[k];
        uint fillEndIdx   = fillEndList[k];

        uint rangeSize = endIdx - startIdx;
        printf("..%.1f%%", k*100.0 / startList.size());
        int nKept = moveKeepProtect<LDType>( LD, preDim, rangeSize, startIdx - pre_start); // reUse LD part

        //int nKept = 0;
        zScores_tmp.resize(rangeSize);
        seqNos_tmp.resize(rangeSize);
        rsIDs_tmp.resize(rangeSize);
        A1_tmp.resize(rangeSize);
        toAvert.resize(rangeSize);
#pragma omp parallel for
        for (uint i = startIdx; i < endIdx; i ++) {
            zScores_tmp[i -startIdx] = (zScores[i]);
            seqNos_tmp[i - startIdx] = (seqNos[i]);
            rsIDs_tmp[i - startIdx] = ( rsIDs[i]);
            toAvert[i - startIdx] = (1-flipped[i]);
        }
        bool performed = true;
        if((endIdx - startIdx) > minDim/5 )
        {
            testMethods(bedFile, rsIDs_tmp, seqNos_tmp, toAvert, zScores_tmp, nMarkers, 
                    nSamples, endIdx - startIdx +1,lambda, qcFile, thread_num, toKeep,
                    Degree_QC, imputed, rsq, zscore_e, ifDup, 
                    startIdx,
                    fillStartIdx, fillEndIdx, LD, nKept, opt.withNA );
                    //(startIdx < pre_end) *  ((endIdx - startIdx)/4) + startIdx,   endIdx, LD, nKept );
        }
        else
            performed = false;

        pre_end = endIdx, pre_start = startIdx, preDim = seqNos_tmp.size();
        if(!performed) preDim = 0;
    }
    cout << endl;

            
    delete[] LD; 

    ofstream qout (qcFile);

    if(doDebug)
    {
        for (uint i =0; i < zScores.size(); i ++ )
            qout << rsIDs[i] << "\t" << zScores[i] << "\t" << imputed[i] << "\t" << rsq[i] << "\t" << ifDup[i] << endl;
    }
    else
    {
        for (uint i =0; i < zScores.size(); i ++ )
        {
            double stat = pow(zScores[i]-imputed[i], 2) / (1-rsq[i]);
            qout << rsIDs[i] << "\t" << stat << "\t" << logPvalueChisq(stat) << "\t" << ifDup[i] << endl;
        }
    }
    qout.close();
   
}
void segmentedQCed (string bfileName, string qcFile, long int nSamples, long int nMarkers, vector<string>& rsIDs,vector<uint>& bp, vector<double>& zScores, vector<long int>& seqNos, vector<bool>& flipped, const Options& opt)
{

    
    int    maxDim      = opt.maxDim; 
    int    minDim      = opt.minDim;
    int    thread_num  = opt.thread_num;
    double Degree_QC   = opt.Degree_QC;
    double lambda      = opt.lambda;
    string bedFile = bfileName + ".bed";
    vector<int> toKeep;
    int cutoff = maxDim;

    
    vector<double>  imputed (zScores.size(), 0);
    vector<double>  rsq     (zScores.size(), 0);
    vector<double>  zscore_e(zScores.size(), 0);
    vector<bool>    ifDup   (zScores.size(), false);

    //-------------------------------------------
    //   identify the gap (centromemere)
    //-------------------------------------------
    vector<long int> diff;
    diff.resize(bp.size());
    adjacent_difference (bp.begin(), bp.end(), diff.begin());
    int gapIdx = distance(diff.begin(), max_element(diff.begin()+1, diff.end())) ;
    assert(gapIdx > 0);
    

    vector<uint> allGaps;
    allGaps.push_back(0);
    if(diff[gapIdx] > 1e6  && gapIdx != diff.size()-1) // this takes care of centromeric region
    {
        allGaps.push_back(gapIdx);
        printf( "A gap (%d) is found at %d - %d. \n", diff[gapIdx], bp[gapIdx], bp[gapIdx-1]);
    }

    LDType* LD = new LDType[ (opt.maxDim + opt.minDim) *  (opt.maxDim + opt.minDim) ]();

    allGaps.push_back(diff.size());
    assert(zScores.size()  - minDim > 0);
    for (uint k = 0; k < allGaps.size() -1; k ++)
    {
        cout << "k = " << k << endl;
        int rangeSize = allGaps[k+1];
        if(rangeSize  - minDim < 0) continue;
        int  startIdx = allGaps[k]; int endIdx = allGaps[k+1]; 
        endIdx   =  rangeSize -  minDim > startIdx + maxDim ? startIdx + maxDim: rangeSize;   
        int notStartInterval = 0;
        int notLastInterval = 1;
        int pre_start = 0;
        int preDim = 0;
        do
        {
            vector<double>  zScores_tmp;
            vector<long>     seqNos_tmp;
            vector<string>    rsIDs_tmp;
            vector<string>       A1_tmp;
            vector<int>       toAvert;

            //uint nKept = 0;
            //
            uint newDim= endIdx - startIdx;

            // uint nKept = moveKeep( LD, preDim, newDim, startIdx - pre_start); // reUse LD part
            uint nKept = 0;

            if(rangeSize == endIdx )
                notLastInterval = 0;
            //int cutoff = endIdx - startIdx +1;
            for (uint i = startIdx; i < endIdx; i ++)
            {
                zScores_tmp.push_back(zScores[i]);
                seqNos_tmp.push_back(seqNos[i]);
                rsIDs_tmp.push_back( rsIDs[i]);
                toAvert.push_back(1-flipped[i]);
            }
            cout << "start : " << rsIDs_tmp[0] << " end: " << rsIDs_tmp[rsIDs_tmp.size()-1] << endl;

            testMethods(bedFile, rsIDs_tmp, seqNos_tmp, toAvert, zScores_tmp, nMarkers, 
                            nSamples, endIdx - startIdx +1,lambda, qcFile, thread_num, toKeep,
                            Degree_QC, imputed, rsq, zscore_e, ifDup, 
                            startIdx,
                            notStartInterval * (cutoff/4) + startIdx,   endIdx -  (cutoff/4) * notLastInterval, LD, nKept, opt.withNA);
                            //startIdx,   endIdx);
            pre_start = startIdx;
            preDim   = seqNos_tmp.size();
            startIdx = startIdx + cutoff/2;
            endIdx   =  rangeSize - minDim >  startIdx + maxDim ? startIdx + maxDim: rangeSize;   
            notStartInterval  =1;
            printf("%d %d %d \n", zScores.size(), startIdx, endIdx);
        }
        while ( notLastInterval );
    }

    delete[] LD;
    ofstream qout (qcFile);

    for (uint i =0; i < zScores.size(); i ++ )
        qout << rsIDs[i] << "\t" << zScores[i] << "\t" << imputed[i] << "\t" << rsq[i] << "\t" << ifDup[i] << endl;
    qout.close();
   
}




void alignGWAS (GWAS& gtab, BedFile& btab,  vector<double>& zScore, vector<long int>& seqNo, vector<bool>& toFlip, vector<string>& rsID,vector<uint>& bp )
{
    printf("[info] Aligning GWAS to bedfile assumming the bfile SNPs are ordered by BP.\n");
    vector<long int> alignToWhich (gtab.size(), -1);
    vector<bool    > haveFliped   (gtab.size(), false);
    int ncpus= omp_get_num_threads();
    map<string, long int> id_map;
    map<string, long int>::iterator iter;
    for (long int i =0 ; i < gtab.size(); i ++)
        id_map [ gtab.rs[i] ] = i;
    //for (long int j =0 ; j < btab.rs.size(); j ++)
    for (long int k =0 ; k < btab.include.size(); k ++)
    {
        uint j = btab.include[k];

        iter = id_map.find(btab.rs[j]);
        if(iter != id_map.end() )
        {
            long int i = iter->second;
            if(gtab.A1[i] == btab.A1[j] && gtab.A2[i] == btab.A2[j]  )
            {
                zScore.push_back(gtab.zscore[i]);
                seqNo.push_back(btab.seqNo[j]);
                rsID.push_back(btab.rs[j]);
                bp.push_back(btab.bp[j]);
                toFlip.push_back(false);
                //alignToWhich[i] = btab.seqNo[j];
            }
            else if(gtab.A1[i] == btab.A2[j] && gtab.A2[i] == btab.A1[j]  )
            {
                zScore.push_back(gtab.zscore[i]);
                seqNo.push_back(btab.seqNo[j]);
                rsID.push_back(btab.rs[j]);
                bp.push_back(btab.bp[j]);
                toFlip.push_back(true);
                //alignToWhich[i] = btab.seqNo[j];
                //haveFliped[i] = true;
            } 
        }
    }
    printf("[info]%d SNPs (rsID) were shared between the summary and reference data. \n", rsID.size());



}





//void runQCForEQTL(string bfile, bInfo* bdata, eqtlInfo* esdata, int probechr, vector<string>& allSNPs, vector<string>& toFlip, char* outFileName, int thread_num, double Degree_QC )
void runQC(const Options& opt)
{
    string summmaryFile = opt.summmaryFile;
    string bfileName = opt.bfileName;
    string outPrefix = opt.outPrefix;


    string qcFile = outPrefix + ".qc.txt";
    // read summary
    // gzopen() can be used to read a file which is not in gzip format; in this case gzread() will directly read from the file without decompression
    // Therefore, there is no need to judge if it is plain file or gz file.
    GWAS gwasDat;
    gwasDat = GWAS  (summmaryFile, true);
    // read bedfile
    BedFile  ref (bfileName, opt.mafThresh, opt.thread_num);

    if(opt.targetSNP != "")
    {
        D(cout << "Extracting SNPs at the target SNP : " << opt.targetSNP << endl;);
        uint foundAt = ~0;
        // find the SNP
        for (uint kk =0; kk < ref.bp.size(); kk ++)
            if(ref.rs[kk] == opt.targetSNP)
            {
                foundAt = kk; 
                break;
            }

        if(foundAt != ~0)
        {
            printf("[info] Target %s is found at %d.\n", opt.targetSNP.c_str(),   foundAt);
            vector<uint> updatedInclude;
            for (uint kk =0; kk < ref.include.size(); kk ++)
            {
                uint  i = ref.include[kk];
                if(fabs( ref.bp[i] -  ref.bp[foundAt] ) <= 10e6)
                    updatedInclude.push_back(i);
            }
            ref.include = updatedInclude;
        }
        else
        {
            printf("[Warning] The target SNP [%s] is not found.\n", opt.targetSNP.c_str() );
            if(opt.ignoreWarnings == false) exit(-1);
        }
        printf("[info] %d were exracted.\n", ref.include.size());
    }

    // apply extraction
    if(opt.extractFile != "")
    {
        D(cout << "[info] Reading SNPs to be extracted (kept) from : " << opt.extractFile << endl;);
        ifstream eFin (opt.extractFile.c_str());
        if (!eFin.is_open())
        {
            cout << "Fail to read from " << opt.extractFile << endl;
            exit (EXIT_FAILURE);
        }
        string tmpStr ;
        map<string, bool> snpList;
        while(getline(eFin, tmpStr))
        {
            snpList [tmpStr] = true;
        }

        vector<uint> updatedInclude;
        for (uint kk =0; kk < ref.include.size(); kk ++)
        {
            uint  i = ref.include[kk];
            if(snpList.find(ref.rs[i] ) !=snpList.end())
            {
                updatedInclude.push_back(i);
            }
        }
        ref.include = updatedInclude;
        eFin.close();
        printf("[info] %d SNPs remained after --extract", updatedInclude.size());
    }

    // Align bfile and summary data.
    vector<double> zScore;
    vector<long int> seqNo;
    vector<bool> toFlip;
    vector<string> rsID;
    vector<uint> bp;
    alignGWAS (gwasDat, ref,   zScore,  seqNo,  toFlip, rsID, bp);

    

    if(opt.maxDist != -1) {
        printf("[Notice] --wind-dist is set, which would overwrite --wind option.\n" );
        segmentedQCed_dist (bfileName, qcFile, ref.N, ref.M,  rsID, bp, zScore, seqNo, toFlip, opt);
    } 
    else {
        segmentedQCed (bfileName, qcFile, ref.N, ref.M,  rsID, bp, zScore, seqNo, toFlip, opt);
    }


}

#endif  //  __DENTIST__

