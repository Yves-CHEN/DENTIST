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
#include <zlib.h>
#include  "bedfile.h"
#include  "bld.io.h"
#include  "options.h"



using namespace std;


//#void runQCForEQTL(string bfile, SMRWK& smrwk, map<string,long>& seqNoMap, int thread_num, SMRWK& smrwk_tmp, double Degree_QC);






class GWAS
{
public:
    vector<string> rs;
    vector<string> A1;
    vector<string> A2;
    vector<double> b;
    vector<double> se;
    vector<double> zscore;
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

map<string, double> createMap(const vector<string>& rsID, const vector<double>& maf, const vector<string>& A1, const vector<string>& A2)
{

    map<string, double> mm;
    for(uint i =0; i< rsID.size(); i ++)
    {
        string key = "";
        double ff = maf[i];
        if(A1[i].compare(A2[i]) > 0)
            key = rsID[i] + A1[i] + A2[i];
        else
        {
            key = rsID[i] + A2[i] + A1[i];
            ff = 1 -ff;
        }

        if(mm.find(key) == mm.end())
            mm [key]   = ff;
        else
            stop("[error] duplicated key [%s]. \n", key);
    }
    return mm;
}

void deltaMAF(GWAS&   gwas, BedFile& ref)
{
    auto m1 = createMap (gwas.rs, gwas.maf, gwas.A1, gwas.A2);
    auto m2 = createMap (ref.rs, ref.maf, ref.A1, ref.A2);

    float threshold = 0.1;
    vector<uint> updatedInclude;

    int before = 0;
    for (uint kk =0; kk < ref.include.size(); kk ++)
    {
        uint  i = ref.include[kk];
        string key = "";
        if(ref.A1[i].compare(ref.A2[i]) > 0)
            key = ref.rs[i] + ref.A1[i] + ref.A2[i];
        else
            key = ref.rs[i] + ref.A2[i] + ref.A1[i];
        if(m1.find(key) != m1.end() && m2.find(key) != m2.end()    ) before ++;
    }


    for (uint kk =0; kk < ref.include.size(); kk ++)
    {
        uint  i = ref.include[kk];
        string key = "";
        if(ref.A1[i].compare(ref.A2[i]) > 0)
            key = ref.rs[i] + ref.A1[i] + ref.A2[i];
        else
            key = ref.rs[i] + ref.A2[i] + ref.A1[i];

        if(m1.find(key) != m1.end() && m2.find(key) != m2.end()   
                && fabs(m1[key] - m2[key]) < threshold)
        {
            updatedInclude.push_back(i);
        }
    }
    ref.include = updatedInclude;
    printf("[info] %d (%.1f%%) SNPs were filtered when applying deltaMaf %f. \n",  before - ref.include.size(),  (before - ref.include.size() ) * 1.0 * 100/ before, threshold);



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
    if(arrSize  - keepFromIdx <= 0) return 0;
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

    bool  readLD = opt.loadLD;
    if(readLD  && Options::FileExist2(bfileName + ".bld") == -1 )
        stop("Cannot find the LD matrix in BLD file at [%s]", (bfileName + ".bld").c_str());

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
    // Finding all the gaps greater than k Mb in size including the centromere
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
            printf("%d %d %d %d fill:%d\n", rangeSize , 
                    startIdx, endIdx, fillStartList[fillStartList.size()-1], fillEndList[fillEndList.size()-1]);
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

    uint theMaxDim = (opt.maxDim + opt.minDim);
    LDType* LD = new LDType[  theMaxDim *  theMaxDim]();
    //readLD = false;
    //
    //

    for (uint k = 0; k < startList.size(); k ++ )
    {
        uint startIdx = startList[k];
        uint endIdx   =   endList[k];
        uint fillStartIdx = fillStartList[k];
        uint fillEndIdx   = fillEndList[k];

        cout << endIdx << ", "<< startIdx<< endl;
        //cout << endIdx - startIdx<< endl;
        //cout << seqNos.size()<< endl;
        //cout << flipped.size()<< endl;


        if(readLD)
        {
            // if(endIdx == seqNos.size())
            // {
            //     seqNos.push_back(seqNos[seqNos.size() -1]+1);
            // }
            bool ifPrint = false;
            string bldLDFile = bfileName;
            int dim = seqNos[endIdx-1] - seqNos[startIdx] +1;
            assert (theMaxDim > dim +1);
            float* LDFromFile = readLDFromFile_FromTo(bldLDFile,
                    dim , seqNos[startIdx], seqNos[endIdx-1] +1, ifPrint);
            vector<long> rtSeqNos; // target seqNos in bed, relative the 0th seqNo.
            vector<long> toAvert; // target seqNos in bed, relative the 0th seqNo.
            toAvert.resize(endIdx - startIdx);
            rtSeqNos.resize(endIdx - startIdx);
            for (uint i = startIdx; i < endIdx; i ++) {
                rtSeqNos[i - startIdx] = (seqNos[i] - seqNos[startIdx] );
                if(flipped[i] == 0)
                    toAvert [i - startIdx] = 1;
                else 
                    toAvert [i - startIdx] = -1;
            }

            int zeros = 0;
            int greaterThanOnes = 0;
            double diagSum = 0;
            for (uint i =0; i < rtSeqNos.size(); i ++)
            {
                for (uint j =0; j < rtSeqNos.size(); j ++)
                {
                    int sign = toAvert[i] * toAvert[j];
                    LD [i * rtSeqNos.size() + j]
                        = LDFromFile [rtSeqNos[i] * dim + rtSeqNos[j] ] *sign ;

                    if(LD [i * rtSeqNos.size() + j] == 0) 
                        zeros ++;

                    
                    if(fabs(LD [i * rtSeqNos.size() + j]) > 1) 
                        greaterThanOnes ++;
                               
                }
            }
            cout << greaterThanOnes << ", " << zeros << endl;
            cout << diagSum/ rtSeqNos.size() << endl;
            delete[] LDFromFile;


        }



        cout << startIdx << endl;
        cout << endIdx   << endl;
        //readLDFromFile_FromTo("test2", 2e6,startIdx,endIdx);  
        uint rangeSize = endIdx - startIdx;
        printf("..%.1f%%", k*100.0 / startList.size());
        int nKept = 0;
        if(!opt.loadLD )
            nKept = moveKeepProtect<LDType>( LD, preDim, rangeSize, startIdx - pre_start); // reUse LD part

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
                    fillStartIdx, fillEndIdx, LD, nKept, opt.withNA, readLD);
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

            bool readLD = true;
            testMethods(bedFile, rsIDs_tmp, seqNos_tmp, toAvert, zScores_tmp, nMarkers, 
                            nSamples, endIdx - startIdx +1,lambda, qcFile, thread_num, toKeep,
                            Degree_QC, imputed, rsq, zscore_e, ifDup, 
                            startIdx,
                            notStartInterval * (cutoff/4) + startIdx,   endIdx -  (cutoff/4) * notLastInterval, LD, nKept, opt.withNA, readLD);
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
    int sum = 0;
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
            sum ++;
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
            sum ++;
            } 
        }
    }
    printf("[info] %d SNPs (rsID) were shared between the summary and reference data. \n", sum);



}





//void runQCForEQTL(string bfile, bInfo* bdata, eqtlInfo* esdata, int probechr, vector<string>& allSNPs, vector<string>& toFlip, char* outFileName, int thread_num, double Degree_QC )
void runQC(const Options& opt)
{
    string summmaryFile = opt.summmaryFile;
    string bfileName = opt.bfileName;
    string outPrefix = opt.outPrefix;


    string qcFile = outPrefix + ".qc.txt";
    // read summary
    // gzopen() can be used to read a file which is not in gzip format;
    //   in this case gzread() will directly read from the file without decompression
    // Therefore, there is no need to judge if summmaryFile is plain file or gz file.
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


    deltaMAF  (gwasDat, ref);
    

    if(opt.maxDist != -1) {
        segmentedQCed_dist (bfileName, qcFile, ref.N, ref.M,  rsID, bp, zScore, seqNo, toFlip, opt);
    } 
    else {
        segmentedQCed (bfileName, qcFile, ref.N, ref.M,  rsID, bp, zScore, seqNo, toFlip, opt);
    }


}

#endif  //  __DENTIST__

