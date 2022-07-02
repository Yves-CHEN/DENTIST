#ifndef __DENTIST__
#define __DENTIST__


#include "invoker.h"
#include "headers.h"
//#include "gzstream/gzstream.h"
#include <bitset>
#include <numeric>

//#include <boost/math/distributions/inverse_chi_squared.hpp>
#include <boost/math/distributions/chi_squared.hpp> 
#include <assert.h>
#include <zlib.h>
#include  "bedfile.h"
#include  "bld.io.h"
#include  "options.h"



using namespace std;




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
    vector<uint>    splSize;
    vector<double> pValue;
    vector<long int>  aligned;
    inline GWAS (){};
    GWAS (string summmaryFile, bool ifGZ );
    // this created an empty GWAS from bedfile
    GWAS(BedFile& btab);  
    bool fillGWAS (GWAS&);
    long int M;
    long int N;
    inline long int size(){return (rs.size());};
};

GWAS::GWAS(BedFile& btab)
{
    for (size_t i = 0; i < btab.include.size(); i ++)
    {
        uint j = btab.include[i];
        rs.push_back( btab.rs[j]);
        A1.push_back( btab.A1[j]);
        A2.push_back( btab.A2[j]);
        maf.push_back( btab.maf[j]);
    }
    b      = vector<double> (rs.size(),  0);
    se     = vector<double> (rs.size(), -1);
    zscore = vector<double> (rs.size(), HUGE_VAL  );
};

bool GWAS::fillGWAS (GWAS& gtab)
{
    map<string, long int> id_map;
    map<string, long int>::iterator iter;
    for (long int i =0 ; i < gtab.rs.size(); i ++)
        id_map [ gtab.rs[i] ] = i;


    for (size_t j =0 ; j < this->rs.size(); j ++)
    {
        iter = id_map.find(rs[j]);
        if(iter != id_map.end() )
        {
            long int i = iter->second;
            if(gtab.A1[i] == A1[j] && gtab.A2[i] == A2[j]  )
            {
                b[j]      = gtab.b     [i];
                se[j]     = gtab.se    [i];
                zscore[j] = gtab.zscore[i];
                splSize[j]= gtab.splSize[i];
                maf[j]    = gtab.maf[i];

            }
            else if(gtab.A1[i] == A2[j] && gtab.A2[i] == A1[j]  )
            {
                b[j]      = gtab.b     [i];
                se[j]     = gtab.se    [i];
                zscore[j] = -gtab.zscore[i];
                splSize[j]= gtab.splSize[i];
                maf[j]    = 1 - gtab.maf[i];
            } 
        }
    }

    return 0;
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
            stop("[error] duplicated key [%s]. \n", key.c_str());
    }
    return mm;
}

void deltaMAF(GWAS&   gwas, BedFile& ref, double threshold, vector<bool>& toFlip,
        vector<string>& rsID,vector<uint>& bp,  vector<int64>& seqNo,
        vector<double>& zScores, vector<uint>& perSNP_N, string outFilePrefix )
{
    auto m1 = createMap (gwas.rs, gwas.maf, gwas.A1, gwas.A2);
    auto m2 = createMap (ref.rs, ref.maf, ref.A1, ref.A2);
    vector<bool>   toFlip_tmp ;
    vector<string> rsID_tmp   ;
    vector<uint>   bp_tmp     ;
    vector<int64>   seqNo_tmp;
    vector<double>   zScore_tmp;
    vector<uint>   perSNP_N_tmp;

    vector<uint> updatedInclude;
    int before = 0;
    for (size_t kk =0; kk < ref.include.size(); kk ++)
    {
        uint  i = ref.include[kk];
        string key = "";
        if(ref.A1[i].compare(ref.A2[i]) > 0)
            key = ref.rs[i] + ref.A1[i] + ref.A2[i];
        else
            key = ref.rs[i] + ref.A2[i] + ref.A1[i];
        if(m1.find(key) != m1.end() && m2.find(key) != m2.end()    ) before ++;
    }


    ofstream exclOut (outFilePrefix + ".DENTIST.excluded.txt",  std::fstream::app);
    vector<string> msg;
    msg.push_back("notFoundInGWAS");
    msg.push_back("AF_QC");

    for (size_t kk =0; kk < ref.include.size(); kk ++)
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
            toFlip_tmp.push_back(toFlip[kk]);
            rsID_tmp.push_back(rsID[kk]);
            bp_tmp.push_back(bp[kk]);
            seqNo_tmp.push_back(seqNo[kk]);
            zScore_tmp.push_back(zScores[kk]);
            perSNP_N_tmp.push_back(perSNP_N[kk]);

        }
        else
        {
            exclOut << key << "\t" <<
                  msg [ (fabs(m1[key] - m2[key]) < threshold ) + 1] << " gwas: "<< m1[key]  << " v.s. ref: " << m2[key]  << endl;
        }
    }

    exclOut.close();
    toFlip = toFlip_tmp;
    rsID   = rsID_tmp  ;
    bp     = bp_tmp    ;
    seqNo  = seqNo_tmp    ;
    zScores = zScore_tmp    ;
    perSNP_N = perSNP_N_tmp    ;

    ref.include = updatedInclude;
    printf("[info] %d (%.1f%%) SNPs were filtered when applying deltaMaf %f. %d retained. \n",  before - ref.include.size(),  (before - ref.include.size() ) * 1.0 * 100/ before, threshold, ref.include.size());



}
/// 
/// // this different from moveKeep by creating a tmp matrix store the previous dat
/// //  before coping dat to LD mat. This avoid reading and writing from the same
/// //  LD mat, which can lead to problems when arrSize < currentDim.
/// template<class T>
/// uint moveKeepProtect(T* LD, uint arrSize, uint currentDim, uint keepFromIdx)
/// {
///     if(long(arrSize)  - long(keepFromIdx) <= 0) return 0;
///     uint tmp_dim = arrSize  - keepFromIdx;
///     double* LD_tmp = new double [tmp_dim * tmp_dim];
///     for (uint i = keepFromIdx, m=0; i < arrSize ; i ++, m ++ )
///         for (uint j = keepFromIdx, n = 0; j < arrSize ; j ++, n ++ )
///             LD_tmp[m *  tmp_dim + n ] = LD[i * arrSize + j];
///     uint m = 0;
///     for (m = 0; m < tmp_dim; m ++ )
///         for (uint n = 0; n < tmp_dim;  n ++ )
///             LD[m *  currentDim + n ] = LD_tmp[m * tmp_dim + n];
///     delete[] LD_tmp;
///     return m;
/// }
/// 
/// 


uint moveKeep(double* LD, size_t arrSize, size_t currentDim, size_t keepFromIdx)
{
    uint m = 0;
    for (size_t i = keepFromIdx; i < arrSize ; i ++, m ++ )
        for (size_t j = keepFromIdx, n = 0; j < arrSize ; j ++, n ++ )
            LD[m *  currentDim + n ] = LD[i * arrSize + j];
    return m;
}




double minusLogPvalueChisq(double stat)
{
    boost::math::chi_squared_distribution<double> mydist(1);
    double p = boost::math::cdf( boost::math::complement(mydist, stat) );
    return ( -log10(p) ) ;
}


void segmentingByDist_new(vector<uint>& bp, vector<uint>& startList, vector<uint>& endList, vector<uint>& fillStartList, vector<uint>& fillEndList, const Options& opt)
{

}



void segmentingByDist(vector<uint>& bp, vector<uint>& startList, vector<uint>& endList, vector<uint>& fillStartList, vector<uint>& fillEndList, const Options& opt)
{
    int    minDim      = opt.minDim;
    int    cutoff      = opt.maxDist;
    int    maxBlockSize= opt.maxDim;
    const int minBlockSize  = 2000;

    // The most distant neighbor for i, given the distance cutoff.
    vector<int> nextIdx;
    for (size_t i =1; i < bp.size(); i ++ )
        if(bp[i] - bp[0] >= cutoff || i == bp.size()-1 ) {nextIdx.push_back(i); break;}
    for (size_t i =1; i < bp.size() ; i ++ ) {
        auto j = nextIdx[nextIdx.size() -1];
        while( bp[j] - bp[i] < cutoff) {
            if(j >= bp.size()-1 ) break;
            j ++;
        }
        nextIdx.push_back(j);
    }

    // A quater of distance cutoff away from i.
    vector<int> quaterIdx;
    for (size_t i =1; ; i ++ )
        if(i == bp.size() || bp[i] - bp[0] >= cutoff/4 ) { quaterIdx.push_back(i-1); break;}
    for (size_t i =1; i < bp.size() ; i ++ ) {
        auto j = quaterIdx[quaterIdx.size() -1];
        while(bp[j] < cutoff/4 + long(bp[i]) && j < bp.size() ) {
            j ++;
        }
        quaterIdx.push_back(j-1);
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
    uint gapSizeThresh =  cutoff/4 ;
    if(gapSizeThresh > cutoff) gapSizeThresh = cutoff;
    for (size_t i = 1; i < diff.size(); i ++)  
    {
        if(diff[i] > gapSizeThresh) 
        {
            allGaps.push_back(i);
        }
    }
    allGaps.push_back(diff.size());

    int times = 0;
    for (size_t k = 0; k < allGaps.size() -1; k ++)
    {
        auto rangeSize = allGaps[k+1];
        if(rangeSize < minBlockSize /2) continue;
        if(rangeSize  - minDim < 0) continue;
        int  startIdx = allGaps[k], endIdx = allGaps[k+1]; 
        endIdx   =  quaterIdx[quaterIdx[quaterIdx[quaterIdx[startIdx] ]] ] +1;
        int  notLastInterval = 1;
        auto old_startIdx = startIdx;
        do
        {
            if(times ++ > 400) stop("p");
            //  save
            fillStartList.push_back(quaterIdx[startIdx] ); // one quater of dist cutoff away
            fillEndList.push_back(quaterIdx[quaterIdx[quaterIdx[startIdx] ]] ); // three quaters of dist cutoff away
            if(rangeSize <= endIdx ){
                notLastInterval = 0;
                //if last interval is small, startIdx goes back one steps.
                // but this does affect the region to fill
                if(long(bp[endIdx-1]) - long(bp[ quaterIdx[old_startIdx] ]) < cutoff )
                    startIdx = quaterIdx[old_startIdx];
            }
            startList.push_back(startIdx), endList.push_back(endIdx);
            // update to next
            old_startIdx = startIdx;
            startIdx =  quaterIdx[quaterIdx[startIdx] ];
            endIdx   =  quaterIdx[quaterIdx[quaterIdx[quaterIdx[startIdx] ]] ] +1 ;
            //if(endIdx - startIdx > maxBlockSize) endIdx = startIdx + maxBlockSize;
            //endIdx   =  rangeSize -  minDim > endIdx ? endIdx: rangeSize;  // This is for the last region.
        }
        while ( notLastInterval );
    }

    if (startList.size()<=0) stop ( "[error] No. of interval is 0\n");
    fillStartList [0] = startList[0];
    fillEndList [fillEndList.size()-1] = endList[endList.size() -1];
    for (size_t i = 0; i < startList.size(); i ++)
        cout << startList[i] << "-" << endList[i] << ", " << bp[endList[i]-1] - bp[startList[i]] << ", " << allGaps[0] <<std::endl;

}

class ImputeOperator
{ 
public:
    int64 *   seqNos    ;
    double*   zScores   ;
    string*   rsIDs     ;
    uint  *   toAvert   ;
    double*   imputed   ;
    double*   rsq       ;
    double*   zScores_e ;
    bool  *   ifDup     ;
    double*   sigma_ii_sq;
    int*      iterID;      // this notes the z-score was calcualted at which iter.
    
    inline ImputeOperator ( vector<int64>&   seqNos , vector<double>&  zScores,    
            vector<string>&   rsIDs    , vector<bool  >&   toAvert  ,
            vector<double>&   imputed  , vector<double>&   rsq      ,
            vector<double>&   zScores_e, vector<bool  >&   ifDup     )
    {
        this->seqNos   = new int64 [zScores.size()]();
        this->zScores  = new double[zScores.size()]();
        this->rsIDs    = new string[zScores.size()]();
        this->toAvert  = new uint  [zScores.size()]();
        this->imputed  = new double[zScores.size()]();
        this->rsq      = new double[zScores.size()]();
        this->zScores_e= new double[zScores.size()]();
        this->ifDup    = new bool  [zScores.size()]();
        this->sigma_ii_sq = new double[zScores.size()]();
        this->iterID   = new int[zScores.size()]();
        std::copy(seqNos.begin(),  seqNos.end(), this->seqNos);
        std::copy(zScores.begin(), zScores.end(),this->zScores);
        std::copy(rsIDs.begin(),   rsIDs.end(),  this->rsIDs);
        std::copy(toAvert.begin(), toAvert.end(),this->toAvert);
        std::copy(imputed.begin(), imputed.end(),this->imputed   );
        std::copy(rsq.begin(),     rsq.end(),    this->rsq       );
        std::copy(zScores_e.begin(), zScores_e.end(),this->zScores_e );
        std::copy(ifDup.begin(),   ifDup.end(),  this->ifDup     );
    };
    ~ImputeOperator()
    {
        delete[] seqNos    ;
        delete[] zScores   ;
        delete[] rsIDs     ;
        delete[] toAvert   ;
        delete[] imputed   ;
        delete[] rsq       ;
        delete[] zScores_e ;        
        delete[] ifDup     ;    
        delete[] sigma_ii_sq;
        delete[] iterID;
    };

};

template<class T>
void loadLDFromBLD (string bldLDFile, int64* seqNos, vector<bool>& flipped, uint startIdx, uint endIdx, T* LD)
{
    cout << "bldLDFile" << bldLDFile << endl;

    bool ifPrint = false;
    int dim = (seqNos[endIdx-1]) -(seqNos[startIdx]) +1;
    assert (theMaxDim > dim +1);
    double* LDFromFile = readLDFromFile_FromTo(bldLDFile,
            dim , seqNos[startIdx], seqNos[endIdx-1] +1, ifPrint);
    vector<long> rtSeqNos;// target seqNos in bed, relative the 0th seqNo.
    vector<long> toAvert; // target seqNos in bed, relative the 0th seqNo.
    toAvert.resize(endIdx - startIdx);
    rtSeqNos.resize(endIdx - startIdx);
    for (auto i = startIdx; i < endIdx; i ++) {
        rtSeqNos[i - startIdx] = (seqNos[i] - seqNos[startIdx] );
        if(flipped[i] == 0)
            toAvert [i - startIdx] = 1;
        else 
            toAvert [i - startIdx] = -1;
    }
    for (size_t i =0; i < rtSeqNos.size(); i ++)
    {
        for (size_t j =0; j < rtSeqNos.size(); j ++)
        {
            int sign = toAvert[i] * toAvert[j];
            saveData<T>(LDFromFile [rtSeqNos[i] * dim + rtSeqNos[j] ] *sign, i, j,LD, rtSeqNos.size());
        }
    }

    setDiag<T>(LD, rtSeqNos.size());
    delete[] LDFromFile;




}
void segmentedQCed_dist (string bfileName, string qcFile, uint nSamples,
        uint nMarkers, vector<string>& rsIDs,vector<uint>& bp,
        vector<double>& zScores, vector<uint> perSNP_N, vector<int64>& seqNos,
        vector<bool>& flipped, const Options& opt)
{
    // **********************************************************************
    // Divide into ranges
    vector<uint> startList;
    vector<uint> endList;
    vector<uint> fillStartList;
    vector<uint> fillEndList;
    segmentingByDist(bp, startList, endList, fillStartList, fillEndList, opt);

    int    minDim      = opt.minDim;
    int    thread_num  = opt.thread_num;
    int    doDebug     = opt.debugMode;
    string bedFile = bfileName + ".bed";

    char bedFileCstr[1000] = "";
    std::strcpy ( bedFileCstr, bedFile.c_str());
    char*   head       = bedFileCstr;
    //vector<int> toKeep;
    bool  readLD = opt.loadLD;
    if(readLD  && Options::FileExist2(bfileName + ".bld") == -1 )
        stop("Cannot find the LD matrix in BLD file at [%s]",
                (bfileName + ".bld").c_str());
    // The computation should be on large enough region.
    assert(zScores.size()  - minDim > 0);
    vector<double>  imputed (zScores.size(), 0);
    vector<double>  rsq     (zScores.size(), 0);
    vector<double>  zscore_e(zScores.size(), 0);
    vector<bool>    ifDup   (zScores.size(), false);

    ImputeOperator impOp( seqNos,  zScores, rsIDs, flipped, 
                                        imputed, rsq, zscore_e, ifDup   ); 
    uint pre_end   = 0;
    uint pre_start = 0;
    uint preDim    = 0;
    int nKept = 0;
    uint theMaxDim = (opt.maxDim + opt.minDim);
    if(bp.size() < theMaxDim) theMaxDim = bp.size() + 1;
       cout << "[info]  At least " << theMaxDim * theMaxDim /2 *8/1000000 << " Mb of memory is required." << endl;




    //LDType* LD = new LDType[  theMaxDim *  theMaxDim]();
    LDType2* LD = createStorage<LDType2>(long(theMaxDim));
    for (size_t k = 0; k < startList.size(); k ++ )
    {
        uint startIdx = startList[k], endIdx   =   endList[k];
        uint fillStartIdx = fillStartList[k], fillEndIdx   = fillEndList[k];
        if((endIdx - startIdx) <= minDim/5 ) 
        {
            preDim = 0;
            continue;
        }
        cout <<   startIdx << "th" << " - "<< endIdx << "th" << endl;
        printf("..%.1f%%", k*100.0 / startList.size());
        if(!opt.loadLD )
            nKept = moveKeepProtect<LDType2>( LD, preDim, endIdx - startIdx,
                    startIdx - pre_start); // reUse LD part
        bool performed = true;
        if((endIdx - startIdx) > minDim/5 )
        {
            int withNA = opt.withNA;
            int cutoff = endIdx - startIdx +1;
            uint arrSize = endIdx - startIdx;
            resize(LD,  arrSize);
            if(!readLD)
            {
                try {
                    _LDFromBfile <LDType2>(&head, &nMarkers, &nSamples, 
                       impOp.seqNos + startIdx, &arrSize, 
                       impOp.toAvert + startIdx, &cutoff,
                       &thread_num, LD, &nKept, &(withNA));
                } 
                catch (long long whichError) {
                    cout << "[error] SNP [" << impOp.rsIDs[startIdx + whichError] << "] has variance of 0" << endl;
                    exit(-1);
                }
            }

            else 
                loadLDFromBLD<LDType2>(bfileName, impOp.seqNos, flipped, startIdx, endIdx, LD);
            

            if(opt.doImpute)
                runImpute<LDType2>( nSamples,
                     impOp.zScores + startIdx,
                     cutoff, thread_num,
                     impOp.imputed, impOp.rsq, impOp.zScores_e, impOp.ifDup,
                     startIdx, fillStartIdx, fillEndIdx, LD,  opt);
            if(opt.doQC)
            {
                runDENTIST<LDType2>( nSamples,
                     impOp.zScores + startIdx,
                     cutoff, thread_num,
                     impOp.imputed, impOp.rsq, impOp.zScores_e, impOp.ifDup, impOp.iterID,
                     startIdx, fillStartIdx, fillEndIdx, LD,  opt);
            }
        }
        else
            performed = false;
        pre_end = endIdx, pre_start = startIdx, preDim = endIdx - startIdx;
        if(!performed) preDim = 0;
    }
    cout << endl;

            
    //delete[] LD; 
    deleteStorage<LDType2>(LD);

    ofstream qout (qcFile+".DENTIST.full.txt");
    ofstream outLierout (qcFile+".DENTIST.short.txt");
    if(doDebug)
    {
        double lambda  = 1;
        if(opt.gcControl)
        {
            vector<double> Td;
            for (size_t i =0; i < zScores.size(); i ++ )
            {
                double Td_i = pow(impOp.zScores[i]-impOp.imputed[i], 2) /(1-impOp.rsq[i]);
                Td.push_back(Td_i);
            }
            std::nth_element(Td.begin(), Td.begin() + Td.size()/2, Td.end());
            lambda =  Td [Td.size()/2]/ 0.46;
        }

        for (size_t i =0; i < zScores.size(); i ++ )
        {
            double stat = pow(impOp.zScores[i]-impOp.imputed[i], 2) /(1-impOp.rsq[i]);
            qout << impOp.rsIDs[i] << "\t" << impOp.zScores[i] << "\t" 
                << impOp.imputed[i] << "\t" 
                << impOp.rsq[i] << "\t" << impOp.ifDup[i] <<"\t" << stat/lambda << endl;
            if( minusLogPvalueChisq(stat/lambda) > -log10(5e-8))
            {
                outLierout  << impOp.rsIDs[i]  << endl;
            }
        }
    }
    else
    {
        double lambda  = 1;
        if(opt.gcControl)
        {
            vector<double> Td;
            for (size_t i =0; i < zScores.size(); i ++ )
            {
                double Td_i = pow(impOp.zScores[i]-impOp.imputed[i], 2) /(1-impOp.rsq[i]);
                Td.push_back(Td_i);
            }
            std::nth_element(Td.begin(), Td.begin() + Td.size()/2, Td.end());
            lambda =  Td [Td.size()/2]/ 0.46;
        }


                 
        for (size_t i =0; i < zScores.size(); i ++ )
        {
            double stat = pow(impOp.zScores[i]-impOp.imputed[i], 2) /(1-impOp.rsq[i]);
            qout << impOp.rsIDs[i] << "\t" << stat /lambda << "\t" 
                <<  minusLogPvalueChisq(stat/lambda) << "\t" << impOp.ifDup[i] << endl;
            if( minusLogPvalueChisq(stat/lambda) > -log10(5e-8))
                outLierout  << impOp.rsIDs[i]  << endl;
        }

    }
    qout.close();
    outLierout.close();
   
}
void segmentedQCed (string bfileName, string qcFile, long int nSamples, long int nMarkers, vector<string>& rsIDs,vector<uint>& bp, vector<double>& zScores, vector<int64>& seqNos, vector<bool>& flipped, const Options& opt)
{

    
    int    maxDim      = opt.maxDim; 
    int    minDim      = opt.minDim;
    int    thread_num  = opt.thread_num;
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
    for (size_t k = 0; k < allGaps.size() -1; k ++)
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
            vector<int64>    seqNos_tmp;
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
            for (auto i = startIdx; i < endIdx; i ++)
            {
                zScores_tmp.push_back(zScores[i]);
                seqNos_tmp.push_back(seqNos[i]);
                rsIDs_tmp.push_back( rsIDs[i]);
                toAvert.push_back(1-flipped[i]);
            }
            cout << "start : " << rsIDs_tmp[0] << " end: " << rsIDs_tmp[rsIDs_tmp.size()-1] << endl;

            bool readLD = true;
            testMethods(bedFile, rsIDs_tmp, seqNos_tmp, toAvert, zScores_tmp, nMarkers, 
                            nSamples, endIdx - startIdx +1, qcFile, thread_num, toKeep,
                             imputed, rsq, zscore_e, ifDup, 
                            startIdx,
                            notStartInterval * (cutoff/4) + startIdx,   endIdx -  (cutoff/4) * notLastInterval, LD, nKept, opt, readLD);
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

    for (size_t i =0; i < zScores.size(); i ++ )
        qout << rsIDs[i] << "\t" << zScores[i] << "\t" << imputed[i] << "\t" << rsq[i] << "\t" << ifDup[i] << endl;
    qout.close();
   
}




void alignGWAS (GWAS& gtab, BedFile& btab,  vector<double>& zScore, vector<uint>& perSNPsampleSize, vector<int64>& seqNo, vector<bool>& toFlip, vector<string>& rsID,vector<uint>& bp, string outFilePrefix )
{
    printf("[info] Aligning GWAS to the reference sample assumming both files are ordered.\n");
    vector<long int> alignToWhich (gtab.size(), -1);
    vector<bool    > haveFliped   (gtab.size(), false);
    int ncpus= omp_get_num_threads();
    map<string, long int> id_map;
    map<string, long int>::iterator iter;
    for (long int i =0 ; i < gtab.size(); i ++)
        id_map [ gtab.rs[i] ] = i;
    //for (long int j =0 ; j < btab.rs.size(); j ++)
    int sum = 0;
    auto include_tmp = btab.include;
    ofstream exclOut (outFilePrefix + ".DENTIST.ignored.txt"); 

    vector<string> msg;
    msg.push_back("notFoundInGWAS");
    msg.push_back("unmatchedAlleles");
    include_tmp.resize(0);
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
                perSNPsampleSize.push_back(gtab.splSize [i]);
                seqNo.push_back(btab.seqNo[j]);
                rsID.push_back(btab.rs[j]);
                bp.push_back(btab.bp[j]);
                toFlip.push_back(false);
                //alignToWhich[i] = btab.seqNo[j];
                include_tmp.push_back(j);
            sum ++;
            }
            else if(gtab.A1[i] == btab.A2[j] && gtab.A2[i] == btab.A1[j]  )
            {
                zScore.push_back(gtab.zscore[i]);
                perSNPsampleSize.push_back(gtab.splSize [i]);
                seqNo.push_back(btab.seqNo[j]);
                rsID.push_back(btab.rs[j]);
                bp.push_back(btab.bp[j]);
                toFlip.push_back(true);
                //alignToWhich[i] = btab.seqNo[j];
                //haveFliped[i] = true;
                include_tmp.push_back(j);
            sum ++;
            }  else
                exclOut << btab.rs[j] << "\t" << btab.A1[j] << "\t" << btab.A2[j]   << "\t" << msg [1] << endl;
        }
        else
            exclOut << btab.rs[j] << "\t" << msg [0] << endl;

    }
    exclOut.close();
    printf("[info] Performing DENTIST at %d SNPs shared between the summary and reference data. \n", sum);
    btab.include = include_tmp;
}




void runSummaryImpute(const Options& opt)
{
    string summmaryFile = opt.summmaryFile;
    string bfileName = opt.bfileName;
    string outPrefix = opt.outPrefix;
    string impFile = outPrefix  ;
    // read summary
    // gzopen() can be used to read a file which is not in gzip format;
    //   in this case gzread() will directly read from the file without decompression
    // Therefore, there is no need to judge if summmaryFile is plain file or gz file.
    GWAS gwasDat;
    gwasDat = GWAS  (summmaryFile, true);
    // read bedfile
    BedFile  ref (bfileName, opt.mafThresh, opt.thread_num);

    setChr(ref, opt.chrID);
    
    if(opt.targetSNP != "")
    {
        D(cout << "Extracting SNPs at the target SNP : " << opt.targetSNP << endl;);
        uint foundAt = ~0;
        // find the SNP
        for (size_t ii =0; ii < ref.include.size(); ii ++)
        {
            auto kk = ref.include[ii];
            if(ref.rs[kk] == opt.targetSNP)
            {
                foundAt = kk; 
                break;
            }
        }

        if(foundAt != ~0)
        {
            printf("[info] Target %s is found at %d.\n", opt.targetSNP.c_str(),   foundAt);
            vector<uint> updatedInclude;
            for (size_t kk =0; kk < ref.include.size(); kk ++)
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
        printf("[info] %d SNPs are extracted at the target region.\n", ref.include.size());
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
        for (size_t kk =0; kk < ref.include.size(); kk ++)
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

    GWAS  impTab(ref);
    impTab.fillGWAS(gwasDat);

    // Align bfile and summary data.
    vector<double> zScore = impTab.zscore;
    vector<bool>   toFlip (impTab.zscore.size(), false);
    vector<uint>   perSNPsampleSize = impTab.splSize;
    vector<string> rsID   = impTab.rs;
    vector<int64> seqNo;
    vector<uint> bp;

    for (size_t i = 0; i < ref.include.size(); i ++)
    {
        uint j = ref.include[i];
        seqNo.push_back(ref.seqNo[j]);
        bp.push_back(ref.bp[j]);
    }

    if(opt.deltaMAF != -1)
        deltaMAF  (gwasDat, ref,  opt.deltaMAF, toFlip, rsID, bp, seqNo, zScore, perSNPsampleSize, impFile);
    
        if(opt.maxDist != -1) {
        //stop("%d %d %d %d", rsID.size(), bp.size(), zScore.size(), seqNo.size() );
        segmentedQCed_dist (bfileName, impFile, ref.N, ref.M, rsID, bp, zScore, perSNPsampleSize, seqNo, toFlip, opt);
    } 
    else {
        segmentedQCed (bfileName, impFile, ref.N, ref.M,  rsID, bp, zScore, seqNo, toFlip, opt);
    }


}

void getFreq(const Options& opt)
{
    string summmaryFile = opt.summmaryFile;
    string bfileName = opt.bfileName;
    string outPrefix = opt.outPrefix;
    string qcFile = outPrefix;
    // read summary
    // gzopen() can be used to read a file which is not in gzip format;
    //   in this case gzread() will directly read from the file without decompression
    // Therefore, there is no need to judge if summmaryFile is plain file or gz file.
    GWAS gwasDat;
    gwasDat = GWAS  (summmaryFile, true);
    BedFile  ref;
    if(opt.loadLD)
    {
        ref  = BLDFILE(opt.bldLDFile); // read BLD file
        bfileName = opt.bldLDFile;
    }
    else
        ref  = BedFile(bfileName, opt.mafThresh, opt.thread_num); // read bedfile

    setChr(ref, opt.chrID);



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
        for (size_t kk =0; kk < ref.include.size(); kk ++)
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



    ofstream mafOut ( opt.outPrefix + ".frq");
    for (size_t i = 0; i <ref.include.size(); i ++)
        mafOut << ref.maf[ref.include[i]] << endl;
    mafOut.close();

}


//
// Require: 1. summary file, 2. plink bed file / BLD file
void runQC(const Options& opt)
{
    string summmaryFile = opt.summmaryFile;
    string bfileName = opt.bfileName;
    string outPrefix = opt.outPrefix;
    string qcFile = outPrefix;
    // read summary
    // gzopen() can be used to read a file which is not in gzip format;
    //   in this case gzread() will directly read from the file without decompression
    // Therefore, there is no need to judge if summmaryFile is plain file or gz file.
    GWAS gwasDat;
    gwasDat = GWAS  (summmaryFile, true);
    BedFile  ref;
    if(opt.loadLD)
    {
        ref  = BLDFILE(opt.bldLDFile); // read BLD file
        bfileName = opt.bldLDFile;
    }
    else
        ref  = BedFile(bfileName, opt.mafThresh, opt.thread_num); // read bedfile

    if(!opt.loadLD)
        setChr(ref,opt.chrID);


    if(opt.targetSNP != "" && opt.targetBP != -1)
    {
        stop("[error] both --target and --target-bp were set, when only one is needed. \n");
    }
    if(opt.targetBP != -1)
    {
        D(cout << "Extracting SNPs at the target SNP : " << opt.targetBP << endl;);
        vector<uint> updatedInclude;
        for (size_t kk =0; kk < ref.include.size(); kk ++)
        {
            uint  i = ref.include[kk];
            if(fabs( ref.bp[i] -  opt.targetBP ) <= opt.radius)
                updatedInclude.push_back(i);
        }
        ref.include = updatedInclude;
        printf("[info] %d SNPs are extracted at the target region.\n", ref.include.size());
        if(ref.include.size() <10)
            stop("[error] less than 10 SNPs has left at the target region.\n");

    }

    
    if(opt.targetSNP != "")
    {
        D(cout << "Extracting SNPs at the target SNP : " << opt.targetSNP << endl;);
        uint foundAt = ~0;
        // find the SNP
        for (size_t kk =0; kk < ref.bp.size(); kk ++)
            if(ref.rs[kk] == opt.targetSNP)
            {
                foundAt = kk; 
                break;
            }

        if(foundAt != ~0)
        {
            printf("[info] Target %s is found at %d.\n", opt.targetSNP.c_str(),   foundAt);
            vector<uint> updatedInclude;
            for (size_t kk =0; kk < ref.include.size(); kk ++)
            {
                uint  i = ref.include[kk];
                if(fabs( ref.bp[i] -  ref.bp[foundAt] ) <= opt.radius)
                    updatedInclude.push_back(i);
            }
            ref.include = updatedInclude;
        }
        else
        {
            printf("[Warning] The target SNP [%s] is not found.\n", opt.targetSNP.c_str() );
            if(opt.ignoreWarnings == false) exit(-1);
        }
        printf("[info] %d SNPs remain after --target.\n", ref.include.size());

        if(ref.include.size() <10)
            stop("[error] less than 10 SNPs has left at the target region.\n");


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
        for (size_t kk =0; kk < ref.include.size(); kk ++)
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
    vector<uint> perSNPsampleSize;
    vector<int64> seqNo;
    vector<bool> toFlip;
    vector<string> rsID;
    vector<uint> bp;
    alignGWAS (gwasDat, ref,   zScore, perSNPsampleSize,  seqNo,  toFlip, rsID, bp, outPrefix);

    if(opt.deltaMAF != -1)
        deltaMAF  (gwasDat, ref, opt.deltaMAF, toFlip, rsID, bp, seqNo, zScore, perSNPsampleSize, outPrefix);
    
    
    if(opt.maxDist != -1) {
        // stop("%d %d %d %d", rsID.size(), bp.size(), zScore.size(), seqNo.size() );
        segmentedQCed_dist (bfileName, qcFile, ref.N, ref.M, rsID, bp, zScore, perSNPsampleSize, seqNo, toFlip, opt);
    } 
    else {
        segmentedQCed (bfileName, qcFile, ref.N, ref.M,  rsID, bp, zScore, seqNo, toFlip, opt);
    }


}


#endif  //  __DENTIST__

