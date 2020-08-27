#ifndef BLD_IO_H
#define BLD_IO_H
//#include <vector>
#include "bedfile.h"
//#include <cstring>
//#include <iostream>
//#include <fstream>
//#include <string>
//#include <cstdlib>
//#include <stdio.h>
//using namespace std;

#include "headers.h"

template<class T>
uint moveKeepProtect(T* LD, uint arrSize, uint currentDim, uint keepFromIdx);

template<class T>
void saveLD (Options& opt);
//
//template<class T>
//void saveLD (string bedfileName, const char* outFilePrefix, uint LDwinSize, int ncpus);


// This loads in a block of LD start from  the [fromIdx], 
//    the size of the block is defined the window size of the bld data.
template <class T>
T* readLDFromFile_at(string ldDatFilePrefix, int windowSize, int fromIdx);
double* readLDFromFile_FromTo(string ldDatFilePrefix, int windowSize, int fromIdx,
        int toIdx);

//void writeLD2File_wind(FILE* outfile, LDType LD[], int dim, int*  rightMostIdx, int fromIdx,  int toIdx);
//int writeLD2File_wind(FILE* outfile, LDType LD[], int dim, int* bp, int ldWind);

//template<class T>
//int writeLD2File_wind(FILE* outfile, ofstream& idxfile, T* LD, int dim, const BedFile& bp, int ldWind, uint startIdx, bool fully);

template<class T>
int _LDFromBfile(char** bedFileCstr, uint* nMarkers, uint* nSamples, uint* theMarkIdx,
        uint* arrSize, uint* toAvert, int* cutoff, int* ncpus,T* result, int* jump, int* withNA);

//typedef short LDType ;
int const RESERVEDUNITS  = 10000;
void findRight(const int* bp, int dim, vector<int>& right, int wind)
{
    int j = 0;
    for (int i =1; i < dim; i ++)
        while (bp[i] > bp[j] + wind )
            right[j++] = i; 
    while(j != dim )  // fills the last bits at end of the vector
        right[j++] = dim ;
}

void trimingBed(BedFile& ref, const Options& opt)
{
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
        if(updatedInclude.size()<10)
            stop("[error] <10 SNPs left after --extract.");
    }

}


// returns where it stops saving, when the LD are not full calculated
//   mode1: write the LD matrix fully
//   mode2: write only those satisfying ldWind > ? Mb
template <class T>
int writeLD2File_wind(FILE* outfile, ofstream& idxfile, T* LD, int dim, const BedFile& ref, int ldWind, uint startIdx, bool fully)
{
    stop("[error] Not implemented.");
    return 0;
}

// int writeLD2File_wind<int>(FILE* outfile, ofstream& idxfile, int* LD, int dim,
//         const BedFile& ref, int ldWind, uint startIdx, bool fully)


//template int writeLD2File_wind <int>(FILE* outfile, ofstream& idxfile, int* LD, int dim, const BedFile& ref, int ldWind, uint startIdx, bool fully);
template <>
int writeLD2File_wind <int>(FILE* outfile, ofstream& idxfile, int* LD, int dim, const BedFile& ref, int ldWind, uint startIdx, bool fully)
{
    vector<int> right(dim, 0);
    auto bp = ref.bp;
    bp.resize(ref.include.size());
    for (uint i =0; i < ref.include.size(); i ++)
        bp [i]    = ref.bp   [ ref.include[i] ];
    findRight(bp.data()+startIdx, dim, right, ldWind);
    uint atIdx = 0 ;
    if(!fully)
    {
        for (atIdx =0; atIdx < dim; atIdx ++)
            if(right[atIdx] == dim) 
                break;
    }
    else
        atIdx = dim;
    /// write index
    for(int i =0; i < atIdx;i ++)
        idxfile <<  ref.print(startIdx +i) << "\t" << i + startIdx << "\t" << right[i] 
            << "\t" << right[i] - i << endl;
    /// write dat
    for(int i =0; i < atIdx;i ++)
    {
        int j = right[i];
        fwrite(LD + i*dim+i +1, sizeof(int), j -i -1 ,  outfile);  
    }
    return atIdx;
}



// Require: bed file
// Accept: --extract --maf
// Reject: --target
//
template<class T>
void saveLD (Options& opt)
{
    stop("Not implemented yet.\n");
}

template<> void saveLD <int> (Options& opt)
{
    string bedfileName = opt.bfileName; 
    const char* outFilePrefix= opt.outPrefix.c_str();
    uint LDwinSize =   opt.maxDist; int ncpus = opt.thread_num;
    uint extendBy     = LDwinSize / 5;//slighly extended LDwind which ease the saving of LD
    int cutoff          = LDwinSize + extendBy;
    string outPrefix (outFilePrefix);
    string rightIdxFile = outPrefix + ".ridx";
    string bldFile      = outPrefix + ".bld";
    FILE*  bldWriter    = fopen(bldFile.c_str(), "w");
    ofstream  idxfile2 (rightIdxFile.c_str());
    const     uint maxDim = 200000;
    // read bedfile
    BedFile  ref (bedfileName, opt.mafThresh, opt.thread_num); // apply MAF thresh
    // --maf --extract
    trimingBed(ref, opt);
    auto  bp    = ref.bp;
    auto  seqNo = ref.seqNo;
    bp.resize(ref.include.size());
    seqNo.resize(ref.include.size());
    for (uint i =0; i < ref.include.size(); i ++)
    {
        bp [i]    = ref.bp   [ ref.include[i] ];
        seqNo [i] = ref.seqNo[ ref.include[i] ];
    }

    string bedFile = bedfileName + ".bed";
    vector<int> right(bp.size(), 0);
    findRight(bp.data(),bp.size(), right, cutoff);

    // write header to bld file.
    //    header includes 4bype * RESERVEDUNITS
    assert (RESERVEDUNITS >=5 );
    int reserved [RESERVEDUNITS] = {};
    reserved[0]  = 1;
    reserved[1]  = ref.N;
    reserved[2]  = ref.include.size();
    reserved[3]  = LDwinSize;
    reserved[4]  = sizeof(int); // number bytes
    for(int i = 5; i < RESERVEDUNITS; i ++) reserved[i]=-9;
    fwrite(&reserved[0],sizeof(int), RESERVEDUNITS, bldWriter);

    int* LD = new int[ maxDim *  maxDim ]();
    //T* LD = createStorage(maxDim);
    int    jump   = 0;
    int    withNA = 0;
    uint startIdx = 0;
    assert (bp.size() >1);
    string genotypeFile = ref.bfileStr + ".bed";
    std::cout << genotypeFile << std::endl;
    char bedFileCstr[10000] = "";
    std::strcpy ( bedFileCstr, genotypeFile.c_str());
    char* head = bedFileCstr;

    uint M      = ref.M;
    uint N      = ref.N;
    uint toMove = 0;
    uint rangeSize =0;
    uint endIdx = right[startIdx];
    do
    {
        endIdx = right[startIdx];
        int nKept = moveKeepProtect<int>( LD, rangeSize, endIdx - startIdx, toMove); 
        jump = nKept;
        rangeSize = endIdx - startIdx;
        assert(rangeSize < maxDim);
        printf("..%.1f%%", startIdx *1.0 / bp.size () * 100);
        //  get LD at the range from startIdx to endIdx
        uint*   theMarkIdx = new uint[rangeSize]  ;
        uint*   toAvert    = new uint[rangeSize] () ;
#pragma omp parallel for
        for (uint i = startIdx; i < endIdx; i ++) 
             theMarkIdx[i-startIdx] = seqNo[i];

        _LDFromBfile <int>(&head, &M, &N, theMarkIdx, &rangeSize,
                toAvert, &cutoff,  &ncpus, LD, &jump, &withNA);
        delete[] theMarkIdx;
        delete[] toAvert;
        //  Save LD
        //int* bp_tmp= new int[rangeSize]();
        //for(uint i =0; i < rangeSize; i ++)
        //    bp_tmp[i] = ref.bp[startIdx + i];
        //int* bp_tmp = ref.bp.data() + startIdx;
        toMove = writeLD2File_wind<int>(bldWriter, idxfile2, LD, rangeSize, ref,
                LDwinSize, startIdx, endIdx == bp.size());
        //delete[] bp_tmp;
        assert(toMove !=0);
        startIdx += toMove;

        cout << "startIdx: " << startIdx << ", bp: "<<  bp [startIdx] << endl;
    }
    while (endIdx !=  bp.size());
    delete[] LD;
    //deleteStorage( LD);
    idxfile2.close();
    fclose(bldWriter);
}












// assume LD (r) is needed.
// assume per-chr calculation
// assume bed is ordered by bp
// from which SNPi to SNPj
// given the right most element
//  assume dim >1
//  assume length ( rightMostIdx) >= dim -1
//void writeLD2File_wind(FILE* outfile, LDType LD[], int dim, int*  rightMostIdx, int fromIdx,  int toIdx)
//{
//    
//    int* right = rightMostIdx;
//    for(int i = fromIdx; i < toIdx;i ++)
//    {
//        int j = right[i];
//        fwrite(LD + i*dim+i +1, sizeof(LDType), j -i -1 ,  outfile);  
//    }
//}







/// // assume LD (r) is needed.
/// // assume per-chr calculation
/// // assume bed is ordered by bp
/// template  <class T>
/// void writeLD2File(T LD[], int dim, const int*  bp,  int ldWind, string outFileName)
/// {
///     
///     string bldname      = string(outFileName)+".bld";
///     string rightIdxFile = string(outFileName)+".ridx";
///     FILE*  outfile     = fopen(bldname.c_str(), "w");
///     ofstream  idxfile (rightIdxFile.c_str());
///     vector<int> reserved(RESERVEDUNITS);
///     int sampleSize = 1000;
///     reserved[1] = sampleSize;
///     reserved[2] = dim;
///     reserved[3] = ldWind;
///     reserved[0] = 0;
///     /// Need to save the LDType as well?
///     for(int i = 4; i < RESERVEDUNITS; i ++) reserved[i]=-9;
///     ///fwrite(&reserved[0],sizeof(int), RESERVEDUNITS, outfile);
///     vector<int> right(dim, 0);
///     findRight(bp, dim,right, ldWind);
///     for (int i =0; i < dim; i ++)
///         idxfile << right[i] << "\n";
/// 
///     // if(sizeof(LDType ) == 2)
///     // {
///     //     // convert LD
///     // }
/// 
///     for(int i = 0; i < dim-1;i ++)
///     {
///         int j = right[i];
///         fwrite(LD + i*dim+i +1, sizeof(LDType), j -i -1 ,  outfile);  
///     }
///     fclose(outfile);
///     idxfile.close();
///     printf("LD information is saved in the binary file %s.\n",bldname.c_str());
/// }
/// 
// another function to load indexing and dim first
//  int   readIdxFile (string bldFile) : return dim and indexing


//   Todos: 
//    1. generalizing it to handle different return types, Template <class T>
//
//    Note: the LD is saved in certain Type (int, float or double), 
//       the caller of this function want other Types, specified in T.
double* readLDFromFile_FromTo(string ldDatFilePrefix, int windowSize, int fromIdx,
        int toIdx, bool ifPrint)
{
    cout <<  "Expected to read LD in a window of " 
        << windowSize  << " bp" << endl;
    cout << "Read from " << fromIdx << ", "  << toIdx << endl;
    string outFilePrefix( ldDatFilePrefix);
    string rightIdxFile = string(outFilePrefix)+".ridx";
    string bldname      = string(outFilePrefix)+".bld";
    FILE*  datFile      = fopen(bldname.c_str(), "r");
    ifstream  idxfile (rightIdxFile.c_str());
    const uint maximalStorage = 1e10; /// 10G
    uint dim = toIdx - fromIdx;
    assert(dim>0);
    if(dim > sqrt(maximalStorage/4) ) 
        stop("[error] Too many variants [%d] to load into memory.\n", dim);
    double* LD_typeT = new double[dim * dim](); /// where LD matrix will be loaded.
    // ----------------------------------------------------------------------------
    // read header save in integer (4bytes)
    int headerInfo [5] = {}; //fileSubtype, N, M, LDwindSize, LD in number of Bytes
    fread (headerInfo,  sizeof(int), 5 ,datFile);
    if(windowSize > headerInfo[3]) // sanity check
        stop("[error] expected WindSize[%d] is larger than the data [%d].\n", 
                    windowSize, headerInfo[3]);
    int SaveInNumByte = headerInfo[4];
    assert (toIdx - fromIdx <= windowSize);
    rewind(datFile);
    cout << "SaveInNumByte " << SaveInNumByte << endl;

    // --------------------------------------------------------------------------
    // load the right most element Idx. 
    //    Notice!! had a problem with int for loc, the int overflowing problem
    vector<long long> rowLen;
    long long sum  = RESERVEDUNITS * sizeof(int);//jump over bytes reserved for header
    vector<long long> loc ;
    if (idxfile.is_open()) 
    {
        uint ridx   = 0, snpIdx = 0, dist   = 0; 
        int ibuf = 0;
        string cbuf = "0";
        double dbuf = 0.0;
        string str_buf;
        //while(idxfile >> str_buf >> str_buf >> dbuf >> ibuf >> cbuf >> cbuf >> snpIdx  >> ridx >> dist)
        while (idxfile >> snpIdx  >> ridx >> dist)  // expect 3-column file format
        {
            loc.push_back(sum ); 
            // *sizeof(LDTtype) :  multiply number of bytes to jump to next_i
            sum += (dist -1) * SaveInNumByte; // dist-1 as diagnal value is not saved
            rowLen.push_back(dist);
        }
        if(sum> pow(2, 63)) stop( "[error] long long overflowing\n");
    }
    else
        stop("Unable to open file [%s]", rightIdxFile.c_str() );
    assert (rowLen.size() > (toIdx - fromIdx));

    /// --------------------------------
    /// Reading LD matrix.
    if(SaveInNumByte == 2)
    {
        short* LD = new short[dim * dim](); /// new and set zeros
        //for (uint i = 0; i < dim; i ++)
        //    LD[i*dim + i] =10000;
        int SaveUpTo = dim -1;
        for(int i = fromIdx, k = 0; i < toIdx -1;i ++, k ++)
        {
            // fwind
            fseek ( datFile,  loc[i] , SEEK_SET );
            fread(&LD [ k * dim + k +1], SaveInNumByte, SaveUpTo ,datFile);  
            SaveUpTo  --;
        }
        for (uint i = 0; i < dim-1; i ++)
            for (uint j = i+1; j < dim; j ++)
                LD[j*dim + i] = LD[i*dim + j]  ;
        for (uint i = 0; i < dim; i ++)
            for (uint j = 0; j < dim; j ++)
                LD_typeT[i*dim + j] = double(LD[i*dim + j] * 1.0 / 10000)  ;
        delete[] LD;
    }
    if(SaveInNumByte == 4)
    {
        int* LD = new int[dim * dim](); /// new and set zeros
        uint scaler  = 1e8;
        for (uint i = 0; i < dim; i ++)
            LD[i*dim + i] = scaler;
        int SaveUpTo = dim -1;
        for(int i = fromIdx, k = 0; i < toIdx -1;i ++, k ++)
        {
            // fwind
            fseek ( datFile,  loc[i] , SEEK_SET );
            fread(&LD [ k * dim + k +1], SaveInNumByte, SaveUpTo ,datFile);  
            SaveUpTo  --;
        }
        for (uint i = 0; i < dim-1; i ++)
            for (uint j = i+1; j < dim; j ++)
                LD[j*dim + i] = LD[i*dim + j]  ;

        for (uint i = 0; i < dim; i ++)
            for (uint j = 0; j < dim; j ++)
                LD_typeT[i*dim + j] = double(LD[i*dim + j] * 1.0 / scaler)  ;
        delete[] LD;
    }
    if(ifPrint)
        for (uint i = 0; i < dim; i ++)
        {
            for (uint j = 0; j < dim; j ++)
                cout << LD_typeT[j*dim + i] << "\t"  ;
            cout << endl;
        }
    fclose(datFile);
    idxfile.close();
    return LD_typeT;

}

template <class T>
T* readLDFromFile_at(string ldDatFilePrefix, int windowSize, int fromIdx)
{
    string outFilePrefix( ldDatFilePrefix);
    string rightIdxFile = string(outFilePrefix)+".ridx";
    string bldname      = string(outFilePrefix)+".bld";
    FILE*  datFile      = fopen(bldname.c_str(), "r");
    ifstream  idxfile (rightIdxFile.c_str());
    //const int tmpReadBufferSize = 2000000; // read buffer for a row of LD
    // ----------------------------------------------------------------------------
    // read header save in integer (4bytes)
    int headerInfo [4] = {}; // N, M, LDwindSize, LD in number of Bytes
    fread (headerInfo,  sizeof(int), 4 ,datFile);
    if(windowSize > headerInfo[2]) // sanity check
    {
        printf("[error] expected WindSize[%d] is larger than the data [%d].", 
                    windowSize, headerInfo[2]);
        exit(-1);
    }
    rewind(datFile);
    int SaveInNumByte = headerInfo[3];
    //if(headerInfo[3] == 2) // LD saved in short
    //{
    //} else if(headerInfo[3] == sizeof(float) )
    //{
    //} else if(headerInfo[3] == sizeof(double) )
    //{
    //}else
    //{
    //    printf("Unrecognized LD data type, when it should be short float or double",);
    //    exit(-1);
    //}

    // --------------------------------------------------------------------------
    // load the right most element Idx. 
    vector<int> right;
    uint sum  = RESERVEDUNITS * sizeof(int); // jump over bytes reserved for header
    vector<uint> loc ;
    if (idxfile.is_open()) 
    {
        int ridx =0; 
        int snpIdx = 0;
        int dist  = 0; // in number of SNPs
        // expect 3-column file format:
        while (idxfile >> snpIdx  >> ridx >> dist) 
        {
            loc.push_back(sum ); 
            // *sizeof(LDTtype) :  multiply number of bytes to jump to next_i
            //  -1  :  diagnal value is not saved
            sum += (dist -1)* SaveInNumByte; 
            right.push_back(ridx);
        }
    }
    else
        cout << "Unable to open file [" << rightIdxFile << "]"   <<endl;

    int lookAt = fromIdx;
    uint dim = right[lookAt] - lookAt;
    assert(dim>0);
    char* LD_untyped = new char[dim * SaveInNumByte * dim *SaveInNumByte]();;


    int SaveUpTo = dim -1;
    for(int i = lookAt, k = 0; i < right[lookAt] -1;i ++, k ++)
    {
        // fwind
        fseek ( datFile,  loc[i] , SEEK_SET );
        fread(LD_untyped + (k * dim + k +1) * SaveInNumByte, SaveInNumByte, SaveUpTo * SaveInNumByte ,datFile);  
        SaveUpTo  --;
    }

    fclose(datFile);
    idxfile.close();
    cout << "loaded" << endl;

    cout <<  SaveInNumByte << endl;
        cout << sizeof(T) << endl;
    T* LD_typeT =  NULL;

    if(sizeof(T) == SaveInNumByte)
    {
        LD_typeT = (T*)LD_untyped;
    }
    else
    {
        T* LD_typeT = new T[dim * dim](); /// new and set zeros
        if(SaveInNumByte ==2) // saved in short, but expecting float or double
        {

        cout << "tt2" << endl;
            short* LD = (short*) LD_untyped;
            for (uint i = 0; i < dim-1; i ++)
                for (uint j = i+1; j < dim; j ++)
                    LD_typeT[j*dim + i] = float(LD[i*dim + j] * 1.0 / 10000)  ;

        } else if (is_floating_point<T>::value || std::is_same<T, double>::value ) 
        {

            if (SaveInNumByte == sizeof(float)) {
                float* LD = (float*)LD_untyped;
                for (uint i = 0; i < dim-1; i ++)
                    for (uint j = i+1; j < dim; j ++)
                        LD_typeT[j*dim + i] = LD[i*dim + j] ;
            }
            if (SaveInNumByte == sizeof(double)) {
                double* LD = (double*)LD_untyped;
                for (uint i = 0; i < dim-1; i ++)
                    for (uint j = i+1; j < dim; j ++)
                        LD_typeT[j*dim + i] = LD[i*dim + j] ;
            }
        } else if (is_same<T,short>::value) /// saved in float/double, but expecting short
        {

            T* LD = (T*)LD_untyped;
            for (uint i = 0; i < dim-1; i ++)
                for (uint j = i+1; j < dim; j ++)
                    LD_typeT[j*dim + i] = short(LD[i*dim + j] * 10000) ;

        } else
        {
            printf("[error] Unexpected types \n");
            exit(-1);
        }
        delete[]  LD_untyped;
    }
    cout << "recreawted" << endl;
    // for (uint i = 0; i < dim-1; i ++)
    //    for (uint j = i+1; j < dim; j ++)
    //        LD_typeT[j*dim + i] = LD_typeT[i*dim + j]  ;
    cout << LD_typeT << endl;

    for (uint i = 0; i < dim-1; i ++)
        cout << LD_typeT[i*dim + dim -1] << endl;
    return LD_typeT;
}




#endif  ///BLD_IO_H


