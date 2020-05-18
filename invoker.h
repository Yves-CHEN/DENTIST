#ifndef __INVOKER__
#define __INVOKER__

#include "headers.h"
#include <options.h> 

//#include <template>
using namespace std;
typedef unsigned int uint;
typedef unsigned char uchar;

template<class T>  int _LDFromBfile(char** bedFileCstr, uint* nMarkers, uint* nSamples, uint* theMarkIdx,
        uint* arrSize, uint* toAvert, int* cutoff, int* ncpus,T* result, int* jump, int* withNA);
/// 
/// 
template<class T> void DENTIST(T* LDmat, uint* markerSize, uint* nSample, double* zScore,
        double* imputedZ, double* rsq, double* zScore_e, double pValueThresh,  int* interested, float propSVD, bool gcControl, int* ncpus);

template<class T> void impute(T* LDmat, uint* markerSize, uint* nSample,
        double* zScore, double* imputedZ, double* rsq, double* zScore_e,
        double pValueThresh, float propSVD, int* ncpus);






// void generateData (uint theMarkIdx[], uint toAvert[], double zScores[], uint M)
// {
// 
//     for (int i =0; i < M; i++)
//         theMarkIdx[i] = i +1;
//     for (int i =0; i < M; i++)
//         toAvert[i] = 0;
//     for (int i =0; i < M; i++)
//         zScores[i] = i +2;
// 
// }

void findDup2(double* r,  double rThreshold, vector<int>& dupBearer, vector<int>& sign)
{
    uint dim  = dupBearer.size();
    vector <uint> ifDup ( dim, 0);
    for (int i =0; i < dim; i ++)
    {
        if(ifDup [i] != 0 ) continue;

        for (int j =i+1; j < dim ; j ++)
        {
            ifDup[j] += fabs(r[i*dim + j]) > rThreshold;
            if((i != j ) && fabs(r[i * dim + j ]) > rThreshold) 
            {
                if(r[i * dim + j ] < 0) sign[j] = -1;
                dupBearer[j] = i;
            }
        }
    }

}

template<class T>
void runDENTIST( uint nSamples,
        double* the_zScores,  
        int cutoff,   int ncpus,  double* imputedZ,
        double* rsq, double* zScore_e, bool* ifDup,
        int thePos, int startIdx, int endIdx, T* result,
        const Options& opt)
{

    uint    arrSize    = cutoff -1;
    double* zScores    = the_zScores  ;
    double rThreshold = round(sqrt(opt.dupThresh ) * 1000 )/1000.0;
    vector<int>    dupBearer( arrSize, -1);
    vector<double> corABS( arrSize, -1);
    vector<int>    sign( arrSize, 1);
    findDup<T>(result, rThreshold,  dupBearer, corABS, sign);
    uint sum =0;
    for (int i =0; i < dupBearer.size(); i ++)
        sum += (dupBearer[i] !=-1);

    int interested = -1;
    double*   zScores_tmp = new double[ arrSize]() ;
    double*  imputedZ_tmp = new double[ arrSize]() ;
    double*       rsq_tmp = new double[ arrSize]() ;
    double*  zScore_e_tmp = new double[ arrSize]() ;
    uint      arrSize_tmp = arrSize - sum;
    //LDType* resultNoDup = new LDType[ arrSize * arrSize]() ;
    T* resultNoDup = createStorage<T>(long(arrSize_tmp));

    int count =0;
    for (uint i = 0; i < dupBearer.size(); i ++)
    {
        if(dupBearer[i] == -1)
        {
           zScores_tmp[count] = zScores[i];
           int count2 = 0;
           for (uint j =0; j < dupBearer.size(); j ++)
                if(dupBearer[j] == -1)
                {
                    saveData(readMatrix(result, arrSize, i,j ), count, count2, resultNoDup, arrSize_tmp);
                    //resultNoDup[count * arrSize_tmp + count2] = result[i*arrSize + j];
                    count2 ++;
                }
           count ++ ;
        }
    }
    DENTIST<T>(resultNoDup, &arrSize_tmp,  &nSamples, zScores_tmp,
            imputedZ_tmp, rsq_tmp, zScore_e_tmp, opt.pValueThreshold,
            &interested, opt.propPCtrunc, opt.gcControl, &ncpus);
    deleteStorage(resultNoDup);
    count =0;
    vector<uint> assignIdx (dupBearer.size() , 0);
    for (uint i =0; i <dupBearer.size(); i ++)
    {
        if(dupBearer[i ] == -1)
            assignIdx[i] = count ++;
        else
            assignIdx[i] = dupBearer[i];
    }
    for (uint i =startIdx; i < endIdx; i ++)
    {
        if(i - thePos > dupBearer.size()) stop("[error] function runDENTIST"); 
        imputedZ[i]  = imputedZ_tmp[assignIdx[i - thePos]] * sign[i-thePos];
        //rsq[i]       = rsq_tmp     [assignIdx[i - thePos]] ;
        uint j       = assignIdx[i - thePos];
        rsq[i]       =  rsq_tmp [j] ;
        zScore_e[i]       =  zScore_e_tmp[j] ;
        //rsq[i]       = result[arrSize * i + i] - rsq_tmp [j] ;
        ifDup[i]     = dupBearer[i - thePos] != -1 ;
    }

    delete[]  zScores_tmp;
    delete[] imputedZ_tmp;
    delete[]      rsq_tmp;
    delete[] zScore_e_tmp;
}

template<class T>
void runImpute(uint nSamples, double* zScores,  
        int cutoff,   int ncpus,  double* imputedZ,
        double* rsq, double* zScore_e, bool* ifDup,
        int thePos, int startIdx, int endIdx, T* result,
        const Options& opt)
{
    uint    arrSize    = cutoff -1;
    double rThreshold = round(sqrt(opt.dupThresh ) * 1000 )/1000.0;
    vector<int>    dupBearer( arrSize, -1);
    vector<double> corABS( arrSize, -1);
    vector<int>    sign( arrSize, 1);
    findDup<T>(result, rThreshold,  dupBearer, corABS, sign);
    uint sum =0;
    for (int i =0; i < dupBearer.size(); i ++)
        sum += (dupBearer[i] !=-1);
    double*   zScores_tmp = new double[ arrSize]() ;
    double*  imputedZ_tmp = new double[ arrSize]() ;
    double*       rsq_tmp = new double[ arrSize]() ;
    double*  zScore_e_tmp = new double[ arrSize]() ;
    uint      arrSize_tmp = arrSize - sum;
    //LDType*   resultNoDup = new LDType[ arrSize * arrSize]() ;
    T* resultNoDup = createStorage<T>(long(arrSize_tmp));
    for (uint i = 0, new_i =0; i < dupBearer.size(); i ++)
    {
        if(dupBearer[i] == -1)
        {
            zScores_tmp[new_i] = zScores[i];
            for (uint j =0, new_j= 0; j < dupBearer.size(); j ++)
                if(dupBearer[j] == -1)
                {
                    //resultNoDup[new_i* arrSize_tmp + new_j] = result[i*arrSize + j];
                    saveData(readMatrix(result, arrSize, i,j ), new_i, new_j, resultNoDup, arrSize_tmp);
                    new_j ++;
                }
            new_i ++ ;
        }
    }
    impute<T>(resultNoDup, &arrSize_tmp,  &nSamples, zScores_tmp,
            imputedZ_tmp, rsq_tmp, zScore_e_tmp, opt.pValueThreshold, 
            opt.propPCtrunc, &ncpus);
    deleteStorage(resultNoDup);
    int count =0;
    vector<uint> assignIdx (dupBearer.size() , 0);
    for (uint i =0; i <dupBearer.size(); i ++)
    {
        if(dupBearer[i ] == -1)
            assignIdx[i] = count ++;
        else
            assignIdx[i] = dupBearer[i];
    }
    for (uint i =startIdx; i < endIdx; i ++)
    {
        if(i - thePos > dupBearer.size()) stop("[error] function runImpute"); 
        imputedZ[i]  = imputedZ_tmp[assignIdx[i - thePos]] * sign[i-thePos];
        uint j = assignIdx[i - thePos];
        rsq[i]       = readMatrix(resultNoDup,arrSize_tmp,  j, j) - rsq_tmp     [assignIdx[i - thePos]] ;
        ifDup[i]     = dupBearer[i - thePos] != -1 ;
    }
    delete[]  zScores_tmp;
    delete[] imputedZ_tmp;
    delete[]      rsq_tmp;
    delete[] zScore_e_tmp;
}



void testMethods(string bedFile,vector<string>& rsIDs,  vector<long>& seqNos, vector<int>& toAvert_Null, vector<double>& the_zScores, uint nMarkers, uint nSamples, 
        int cutoff,  string outFileName, int ncpus, vector<int>& toKeep,  
        vector<double>& imputedZ, vector<double>& rsq, vector<double>& zScore_e, vector<bool>& ifDup, int thePos, int startIdx, int endIdx, LDType* result, uint nKept, const Options& opt, bool preCalculated)
{

    std::cout << "testMethods" << std::endl;
    char bedFileCstr[1000] = "";
    std::strcpy ( bedFileCstr, bedFile.c_str());
    char*   head       = bedFileCstr;
    uint    arrSize    = seqNos.size();
    uint*   theMarkIdx = new unsigned int[arrSize]  ;
    uint*   toAvert    = new unsigned int[arrSize] ;
    double* zScores    = new double[arrSize] ;

    for (uint i =0; i < arrSize; i ++)
    {
         theMarkIdx[i] = seqNos[i];
         toAvert[i]    = toAvert_Null[i];
         zScores[i]    = the_zScores[i];
    }


    int jump = nKept ;
    int withNA = opt.withNA;
    //int jump = 0;
    //_LDFromBfile (&head, &nMarkers, &nSamples, theMarkIdx, &arrSize, toAvert, &cutoff,  &ncpus, result, &jump, &withNA);
    
    /// ofstream log("t1.txt");
    /// for (uint i =0; i < arrSize; i ++)
    /// {
    ///     for (uint j =0; j < arrSize; j ++)
    ///         log << result[ i *arrSize +j ] << "\t";
    ///     log << endl;
    /// }
    /// log.close();

    if(!preCalculated)
       _LDFromBfile <LDType>(&head, &nMarkers, &nSamples, theMarkIdx,
                &arrSize, toAvert, &cutoff,  &ncpus, result, &jump, &withNA);

    //for (uint i =0; i < arrSize; i ++)
    //    for (uint j =0; j < arrSize; j ++)
    //        result[ i *arrSize +j ]  = int(result[ i *arrSize +j ]  * 1e5) * 1.0 /1e5;
    //
    //
    /// LDType* tmpLD = new LDType[arrSize * arrSize]();
    /// _LDFromBfile <LDType>(&head, &nMarkers, &nSamples, theMarkIdx,
    ///             &arrSize, toAvert, &cutoff,  &ncpus, tmpLD, &jump, &withNA);
    /// double largestDiff = 0;
    /// uint  whichi = 0;
    /// uint  whichj = 0;
    /// bool found = false;
    /// for (uint i =0; i < arrSize; i ++)
    /// {
    ///     if(found) break;
    ///     for (uint j =0; j < arrSize; j ++)
    ///     {
    ///     double tmpDiff = fabs(result[i*arrSize + j] - tmpLD[i*arrSize + j]);
    ///         if( tmpDiff > largestDiff)
    ///         {
    ///             largestDiff = tmpDiff;
    ///             whichi = i;
    ///             whichj = j;
    ///             found = true;
    ///             break;
    ///         }
    ///     }
    /// }
    /// cout << rsIDs[whichi] << endl;
    /// cout << rsIDs[whichj] << endl;
    /// cout << whichi  << endl;
    /// cout << whichj  << endl;
    /// cout <<seqNos[whichi] << endl;
    /// cout <<seqNos[whichj] << endl;
    /// cout << result[whichi * arrSize + whichj] << endl;
    /// cout << tmpLD[whichi * arrSize + whichj] << endl;
    /// cout << "largestDiff " << largestDiff << endl;
    /// delete[] tmpLD;

   


   

    // omp_set_num_threads(ncpus);
    // Eigen::setNbThreads(ncpus);
    // regularizeLD(result, arrSize);
   


    //double rThreshold = 0.95;
    // double rThreshold = 0.97;
    //double rThreshold = 0.99;
    //double rThreshold = 0.995;
    double rThreshold = round(sqrt(opt.dupThresh ) * 1000 )/1000.0;
    vector<int> dupBearer( arrSize, -1);
    vector<double> corABS( arrSize, -1);
    vector<int> sign( arrSize, 1);
    findDup<LDType>(result, rThreshold,  dupBearer, corABS, sign);

    
    uint sum =0;
    for (int i =0; i < dupBearer.size(); i ++)
        sum += (dupBearer[i] !=-1);

         

    int interested = -1;
    double*   zScores_tmp = new double[ arrSize]() ;
    double*  imputedZ_tmp = new double[ arrSize]() ;
    double*       rsq_tmp = new double[ arrSize]() ;
    double*  zScore_e_tmp = new double[ arrSize]() ;
    uint      arrSize_tmp = arrSize - sum;
    LDType* resultNoDup = new LDType[ arrSize * arrSize]() ;
    vector<string> tmp_rsIDs;

    int count =0;
    for (uint i = 0; i < dupBearer.size(); i ++)
    {
        if(dupBearer[i] == -1)
        {
           zScores_tmp[count] = zScores[i];
           int count2 = 0;
           for (uint j =0; j < dupBearer.size(); j ++)
                if(dupBearer[j] == -1)
                {
                    resultNoDup[count * arrSize_tmp + count2] = result[i*arrSize + j];
                    count2 ++;
                }
           count ++ ;
           tmp_rsIDs.push_back(rsIDs[i]);
        }
    }
    
    impute<LDType>(resultNoDup, &arrSize_tmp,  &nSamples, zScores_tmp,
            imputedZ_tmp, rsq_tmp, zScore_e_tmp, opt.pValueThreshold,
            opt.propPCtrunc, &ncpus);

    count =0;
    vector<double> imputedZ_tmp_unfold(dupBearer.size(),-1);
    vector<double> rsq_tmp_unfold(dupBearer.size(), -1);
    vector<bool> ifDup_unfold(dupBearer.size(), false);
    for (uint i =0; i <dupBearer.size(); i ++)
    {
        if(dupBearer[i ] == -1)
        {
            imputedZ_tmp_unfold[i]  = imputedZ_tmp[count];
            rsq_tmp_unfold[i]       = rsq_tmp [count];
            count ++;
        }
        else
        {
            imputedZ_tmp_unfold[i]  = imputedZ_tmp[dupBearer[i]] * sign[i];
            rsq_tmp_unfold[i]       = rsq_tmp [dupBearer[i]];
            if(rsq_tmp[dupBearer[i]] < 0)
            {
                cout << dupBearer[i] << " " << rsq_tmp[dupBearer[i]] << endl;
                exit(-1);
            }
            //imputedZ_tmp_unfold[i]  = zScores_tmp[dupBearer[i]] * corABS[i] * sign[i];
            //rsq_tmp_unfold[i]       = corABS[i] * corABS[i] - 0.0001;
            ifDup_unfold[i]         = true;
        }
    }
    delete[] resultNoDup;
    for (uint i =startIdx; i < endIdx; i ++)
    {

        if(i - thePos > imputedZ_tmp_unfold.size()) {cout << "[error] i - thePos > imputedZ_tmp_unfold" << endl; exit(-1);}
        imputedZ[i]  = imputedZ_tmp_unfold[i - thePos] ;
        rsq[i]       = rsq_tmp_unfold[i - thePos]      ;
        ifDup[i]     = ifDup_unfold[i - thePos]        ;
    }

    delete[]  zScores_tmp;
    delete[] imputedZ_tmp;
    delete[]      rsq_tmp;
    delete[] zScore_e_tmp;
    delete[] theMarkIdx;
    delete[] toAvert;
    delete[] zScores;

}
#endif


