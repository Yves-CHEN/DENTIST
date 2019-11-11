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


// This loads in a block of LD start from  the [fromIdx], 
//    the size of the block is defined the window size of the bld data.
template <class T>
T* readLDFromFile_at(string ldDatFilePrefix, int windowSize, int fromIdx);
float* readLDFromFile_FromTo(string ldDatFilePrefix, int windowSize, int fromIdx,
        int toIdx);

void writeLD2File_wind(FILE* outfile, LDType LD[], int dim, int*  rightMostIdx, int fromIdx,  int toIdx);
//int writeLD2File_wind(FILE* outfile, LDType LD[], int dim, int* bp, int ldWind);

int writeLD2File_wind(FILE* outfile, ofstream& idxfile, LDType LD[], int dim, int* bp, int ldWind, uint startIdx);
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
    while(j != dim -1)  // fills the last bit at end of the vector
        right[j++] = dim ;
}





void saveLD (string bedfileName, const char* outFilePrefix, uint LDwinSize, int ncpus)
{
    int offset =  200000; // generate offset to LDwind, which ease the saving of LD.
    int cutoff =  LDwinSize + offset;
    string outPrefix (outFilePrefix);
    string rightIdxFile = outPrefix + ".ridx";
    string bldFile      = outPrefix + ".bld";
    FILE*  bldWriter   = fopen(bldFile.c_str(), "w");
    ofstream  idxfile2 (rightIdxFile.c_str());
    uint maxDim = 10000;
    // read bedfile
    BedFile  ref (bedfileName);
    string bedFile = bedfileName + ".bed";
    auto& bp   = ref.bp;

    vector<int> right(bp.size(), 0);
    findRight(&bp[0],bp.size(), right, cutoff);

    // write header to bld file.
    //    header includes 4bype * RESERVEDUNITS
    assert (RESERVEDUNITS >=4 );
    int reserved [RESERVEDUNITS] = {};
    reserved[0]  = ref.N;
    reserved[1]  = ref.M;
    reserved[2]  = LDwinSize;
    reserved[3]  = sizeof(LDType); // number bytes
    for(int i = 4; i < RESERVEDUNITS; i ++) reserved[i]=-9;
    fwrite(&reserved[0],sizeof(int), RESERVEDUNITS, bldWriter);

    //ofstream  idxfile (rightIdxFile.c_str());
    //for (int i =0; i < bp.size(); i ++)
    //    idxfile << ref.rs[i] << "\t"
    //            << ref.bp[i] << "\t"
    //            << ref.A1[i] << "\t"
    //            << ref.A2[i] << "\t"
    //            << ref.maf[i] << "\t"
    //            << nextIdx[i] << "\n";
    //idxfile.close();

    LDType* LD = new LDType[ maxDim *  maxDim ]();
    int    jump  = 0;
    auto&  seqNo = ref.seqNo;
    int withNA = 0;
    uint startIdx = 0;
    assert (bp.size() >1);
    string genotypeFile = ref.bfileStr + ".bed";
    char bedFileCstr[10000] = "";
    std::strcpy ( bedFileCstr, genotypeFile.c_str());
    char* head = bedFileCstr;

    uint M  = ref.M;
    uint N  = ref.N;
    int  toMove = 0;
    uint rangeSize =0;
    do
    {
        uint     endIdx = right[startIdx];
        int nKept = 
            moveKeepProtect<LDType>( LD, rangeSize, endIdx - startIdx, toMove); 
        jump = nKept;
        rangeSize = endIdx - startIdx;
        printf("..%.1f%%", startIdx *1.0 / bp.size () * 100);
        //  get LD at the range from startIdx to endIdx
        uint*   theMarkIdx = new uint[rangeSize]  ;
        uint*   toAvert    = new uint[rangeSize] () ;
#pragma omp parallel for
        for (uint i = startIdx; i < endIdx; i ++) 
             theMarkIdx[i-startIdx] = seqNo[i];
        _LDFromBfile <LDType>(&head, &M, &N, theMarkIdx, &rangeSize,
                toAvert, &cutoff,  &ncpus, LD, &jump, &withNA);
        delete[] theMarkIdx;
        delete[] toAvert;
       //  Save LD
        int* bp_tmp= new int[rangeSize]();
        for(uint i =0; i < rangeSize; i ++)
            bp_tmp[i] = ref.bp[startIdx + i];
        toMove= writeLD2File_wind(bldWriter, idxfile2, LD, rangeSize, bp_tmp,
                LDwinSize, startIdx);
        assert(toMove !=0);
        startIdx += toMove;

    }
    while (right[startIdx] !=  bp.size());

    delete[] LD;
    idxfile2.close();
    fclose(bldWriter);
}



void saveLD_2 (string bedfileName, int ncpus)
{

    //string summmaryFile = opt.summmaryFile;
    //string bfileName = opt.bfileName;
    //string outPrefix = opt.outPrefix;
    //string  idxFile  = outPrefix + ".qc.txt";
    //
    //int    thread_num  = opt.thread_num;
    //int cutoff = opt.maxDist;
    //
    int offset =  200000;
    int cutoff = 2000000 ;
    int    thread_num  = 8;
    string  bfileName = "./hrs_hm2_par1_chr22";
    string outPrefix =  "test2";
    string rightIdxFile = outPrefix + ".ridx";
    string bldFile      = outPrefix + ".bld";
    FILE*  bldWriter   = fopen(bldFile.c_str(), "w");
    ofstream  idxfile2 (rightIdxFile.c_str());
    uint maxDim = 10000;
    
    // read bedfile
    BedFile  ref (bfileName);
    string bedFile = bfileName + ".bed";
    auto& bp   = ref.bp;



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

    int minBlockSize = 200;
    int minDim = 200;

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

    for (uint i =0; i < startList.size(); i ++)
        cout << startList [i] << "\t" << endList[i] << endl;



  

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
void writeLD2File_wind(FILE* outfile, LDType LD[], int dim, int*  rightMostIdx, int fromIdx,  int toIdx)
{
    
    int* right = rightMostIdx;
    for(int i = fromIdx; i < toIdx;i ++)
    {
        int j = right[i];
        fwrite(LD + i*dim+i +1, sizeof(LDType), j -i -1 ,  outfile);  
    }
}

// returns where it stops saving, when the LD are not full calculated
int writeLD2File_wind(FILE* outfile, ofstream& idxfile, LDType LD[], int dim, int* bp, int ldWind, uint startIdx)
{

    vector<int> right(dim, 0);
    findRight(bp, dim, right, ldWind);
    int atIdx =0;
    for (uint i =0; i < dim; i ++)
        if(right[atIdx++] == dim) break;

    for(int i =0; i < atIdx;i ++)
        idxfile << i + startIdx  << "\t" << right[i] << "\t"
            << right[i] - i << endl;
    //int* right = rightMostIdx;
    for(int i =0; i < atIdx;i ++)
    {
        int j = right[i];
        //cout << LD[i * dim + i +1] << endl;
        fwrite(LD + i*dim+i +1, sizeof(LDType), j -i -1 ,  outfile);  
    }
    return atIdx;
}

// assume LD (r) is needed.
// assume per-chr calculation
// assume bed is ordered by bp
template  <class T>
void writeLD2File(T LD[], int dim, const int*  bp,  int ldWind, string outFileName)
{
    
    string bldname      = string(outFileName)+".bld";
    string rightIdxFile = string(outFileName)+".ridx";
    FILE*  outfile     = fopen(bldname.c_str(), "w");
    ofstream  idxfile (rightIdxFile.c_str());
    vector<int> reserved(RESERVEDUNITS);
    int sampleSize = 1000;
    reserved[1] = sampleSize;
    reserved[2] = dim;
    reserved[3] = ldWind;
    reserved[0] = 0;
    /// Need to save the LDType as well?
    for(int i = 4; i < RESERVEDUNITS; i ++) reserved[i]=-9;
    ///fwrite(&reserved[0],sizeof(int), RESERVEDUNITS, outfile);
    vector<int> right(dim, 0);
    findRight(bp, dim,right, ldWind);
    for (int i =0; i < dim; i ++)
        idxfile << right[i] << "\n";

    // if(sizeof(LDType ) == 2)
    // {
    //     // convert LD
    // }

    for(int i = 0; i < dim-1;i ++)
    {
        int j = right[i];
        fwrite(LD + i*dim+i +1, sizeof(LDType), j -i -1 ,  outfile);  
    }
    fclose(outfile);
    idxfile.close();
    printf("LD information is saved in the binary file %s.\n",bldname.c_str());
}

// another function to load indexing and dim first
//  int   readIdxFile (string bldFile) : return dim and indexing


//   Todos: 
//    1. generalizing it to handle different return types, Template <class T>
//
//    Note: the LD is saved in certain Type (int, float or double), 
//       the caller of this function want other Types, specified in T.
float* readLDFromFile_FromTo(string ldDatFilePrefix, int windowSize, int fromIdx,
        int toIdx, bool ifPrint)
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
    int SaveInNumByte = headerInfo[3];
    assert (toIdx - fromIdx <= windowSize);
    rewind(datFile);

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
            sum += (dist -1) * SaveInNumByte; 
            right.push_back(ridx);
        }
    }
    else
        cout << "Unable to open file [" << rightIdxFile << "]"   <<endl;
    uint dim = toIdx - fromIdx;
    assert(dim>0);
    
    short* LD = new short[dim * dim](); /// new and set zeros
    for (uint i = 0; i < dim; i ++)
        LD[i*dim + i] =10000;

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

    float* LD_typeT = new float[dim * dim](); /// new and set zeros
    for (uint i = 0; i < dim; i ++)
        for (uint j = 0; j < dim; j ++)
            LD_typeT[i*dim + j] = float(LD[i*dim + j] * 1.0 / 10000)  ;

    if(ifPrint)
    {
        cout << "dim : " << dim << endl;
        for (uint i = 0; i < dim; i ++)
        {
            for (uint j = 0; j < dim; j ++)
                cout << LD_typeT [i*dim + j]  << '\t';
            cout << endl;
        }
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
        cout << "tt" << endl;
         LD_typeT = (T*)LD_untyped;
    }
    else
    {
        cout << "tt2" << endl;
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




// 1. check what window size was used?
// 2. window size for loading
// 3. loading LD for which markers
// 4. change the signs of correlation specified by toAvert
LDType * readLDFromFile(string outFilePrefix, int& dim)
{
    string rightIdxFile = string(outFilePrefix)+".ridx";
    string bldname      = string(outFilePrefix)+".bld";
    FILE*  datFile      = fopen(bldname.c_str(), "r");
    ifstream  idxfile (rightIdxFile.c_str());
    // load the right most element Idx. 
    vector<int> right;
    if (idxfile.is_open()) 
    {
        int ridx =0; 
        while (idxfile >> ridx) 
            right.push_back(ridx);
        dim = right.size();
    }
    else
        cout << "Unable to open file [" << rightIdxFile << "]"   <<endl;
    

    int const tmpReadBufferSize = 100000;
    LDType buffer [tmpReadBufferSize];
    cout << "buffer size: " 
        << tmpReadBufferSize * sizeof(LDType)*1.0 /1e6 << " Mb" << endl;
    LDType* LD = new LDType[dim * dim](); /// new and set zeros
    for(int i = 0; i < dim-1;i ++)
    {
        fread(buffer, sizeof(LDType), right[i] -i -1 ,datFile);  
        for (int j = i+1; j <  right [i] ; j ++) 
        {
            LD[i*dim  + j] = buffer [j - (i+1)];
        }
    }
    fclose(datFile);
    idxfile.close();
    return LD;

}


void test1()
{
    uint nMarkers = 27664;
    uint nSamples = 8652;
    uint arrSize  = 1000;
    uint theMarkIdx [10000]    = {};
    uint     toAvert[10000]    = {};

    int cutoff = 1000;
    int ncpus  = 4;
    int jump = 0;
    int withNA = 0;
    double *  LD =  new double[arrSize *arrSize]();
    string bedfilePrefix = 
        "/home/uqwche11/30days/simulation.UK10K/Data/hrs.hm2/hrs_hm2_chr22";

    char bedFileCstr[1000] = "";
    std::strcpy ( bedFileCstr, (bedfilePrefix + ".bed").c_str());
    BedFile  refBedFile (bedfilePrefix);
    cout << "Bim file : " << refBedFile.bp.size() << endl;
    char* head =  bedFileCstr;
    for(uint i =0; i < arrSize ; i ++)
        theMarkIdx[i] = i;

    _LDFromBfile<double> (&head, &nMarkers, &nSamples, &theMarkIdx[0],
        &arrSize, &toAvert[0], &cutoff, &ncpus, LD, &jump, &withNA);


    short*    LD_twoBytes =  new short[arrSize *arrSize]();
    for (uint i = 0; i < arrSize; i ++)
        for (uint j = 0; j < arrSize; j ++)
            LD_twoBytes[i*arrSize + j] = short( LD[i*arrSize + j] * 10000);

    uint windSize = 10000;
    writeLD2File<short>(LD_twoBytes, arrSize, &refBedFile.bp[0], windSize, "tt2");
    int testDim = 0;
//    short* LD_loaded = readLDFromFile("tt2",  testDim);
//    for (uint i =0; i < arrSize; i ++)
//        LD_loaded[i*arrSize + i]  = 1 * 10000;
//
//
//    double diff = 0;
//    double sum  = 0;
//    for (uint i = 0; i < arrSize-1; i ++)
//        for (uint j = i+1; j < arrSize; j ++)
//        {
//            if(refBedFile.bp[j] - refBedFile.bp[i] <= windSize )
//                diff += fabs(LD_twoBytes[i*arrSize + j] - LD_loaded[i*arrSize + j]);
//        }
//    cout << diff << endl;
//    cout << LD_loaded [1] << endl;
//    cout << LD_twoBytes [1] << endl;
        
    //delete[] LD_loaded;
    delete[] LD;
    delete[] LD_twoBytes;


}


#endif  ///BLD_IO_H


