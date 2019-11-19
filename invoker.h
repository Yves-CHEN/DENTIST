#ifndef __INVOKER__
#define __INVOKER__

#include "headers.h"
//#include <template>
using namespace std;
typedef unsigned int uint;
typedef unsigned char uchar;

template<class T>  int _LDFromBfile(char** bedFileCstr, uint* nMarkers, uint* nSamples, uint* theMarkIdx,
        uint* arrSize, uint* toAvert, int* cutoff, int* ncpus,T* result, int* jump, int* withNA);
/// 
/// 
template<class T> void DENTIST(T* LDmat, uint* markerSize, uint* nSample, double* zScore,
        double* imputedZ, double* rsq, double* zScore_e,  double* lambda, int* interested, int* ncpus);
/// 
//
// int _LDFromBfile(char** bedFileCstr, uint* nMarkers, uint* nSamples, uint* theMarkIdx,
//         uint* arrSize, uint* toAvert, int* cutoff, int* ncpus,double* result, int* jump, int* withNA);



// extern void DENTIST(double* LDmat, uint* size, double* zScore, double* imputedZ, double* rsq, double* zScore_e,  double* lambda, int* interested, int* ncpus);

// extern void DENTIST(double* LDmat, uint* size, uint* nSample, double* zScore, double* imputedZ, double* rsq, double* zScore_e,  double* lambda, int* interested, int* ncpus);
//
//
// extern void  DENTIST(double*, unsigned int*, double*, double*, double*, double*, double*, int*, int*);


void  regularizeLD (double* LDmat, int markerSize)
{   
    Eigen::MatrixXd VV (markerSize,  markerSize );


    for (uint i = 0; i < markerSize; i ++)
      for (uint j = 0; j < markerSize; j ++)
          VV(i,j) =   LDmat[i * (markerSize) + j];
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(VV);
    int nRank = es.eigenvectors().rows();
    // int nRank = es.rank();
    Eigen::MatrixXd  wi = Eigen::MatrixXd::Identity(nRank, nRank);
    int nZeros = 0;
    for (int j=0; j<nRank; j++) 
    {
        wi(j,j) = es.eigenvalues()(j);
        if ( es.eigenvalues()(j) < 0.0001)
        {
            nZeros ++;
        }
    }
    nRank = nRank - nZeros;


  
    for (int j=0; j< int(nRank * 3/4); j++) 
        wi(j,j) = 0;

    Eigen::MatrixXd ui = es.eigenvectors();
    Eigen::MatrixXd aa = ui * wi * ui.transpose();
    for (uint i = 0; i < markerSize; i ++)
      for (uint j = 0; j < markerSize; j ++)
      {
             LDmat[i * (markerSize) + j] = aa(i,j) ;
             if(i == j )  LDmat[i * (markerSize) + j]  = 1.0001;
      }
}

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

template <class T>
void findDup(T result[], double rThreshold, vector<int>& dupBearer, vector<double>& corABS, vector<int>& sign)
{
    int count = 0;
    unsigned int matSize = dupBearer.size();
    double minValue = 1;
    for(int i =0 ; i < matSize; i ++)
    {
        if(dupBearer[i] != -1)  continue; // otherwise i is not a dup.
        for(int j =i+1 ; j < matSize ; j ++)
        {
            if((dupBearer[j]  == -1 ) && fabs(result[i*matSize + j ]) > rThreshold) 
            {
                if(result[i*matSize + j ] < 0) sign[j] = -1;
                corABS[j] = fabs(result[i*matSize + j ]);
                dupBearer[j] = count;
            }
            if(fabs(result[i*matSize + j ]) < minValue)
                minValue =  fabs(result[i*matSize + j ]);
        }
        count ++;
    }


}

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






//void testMethods(string bedFile,vector<string> rsIDs, vector<string> A1,  vector<long> seqNos, vector<double> the_zScores, uint nMarkers, uint nSamples, int cutoff,double lambda, string outFileName , int ncpus, vector<int>& toKeep, double Degree_QC)
void testMethods(string bedFile,vector<string>& rsIDs,  vector<long>& seqNos, vector<int>& toAvert_Null, vector<double>& the_zScores, uint nMarkers, uint nSamples, 
        int cutoff, double lambda, string outFileName, int ncpus, vector<int>& toKeep, double Degree_QC, 
        vector<double>& imputedZ, vector<double>& rsq, vector<double>& zScore_e, vector<bool>& ifDup, int thePos, int startIdx, int endIdx, LDType* result, uint nKept, int withNA, bool preCalculated)
{

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

    ///         if( tmpDiff > 0.001)
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
    double rThreshold = 0.995;
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
    
    // DENTIST ( resultNoDup, &arrSize_tmp,  &nSamples, zScores_tmp, imputedZ_tmp, rsq_tmp, zScore_e_tmp,  &lambda, &interested, &ncpus);
    DENTIST<LDType> ( resultNoDup, &arrSize_tmp,  &nSamples, zScores_tmp, imputedZ_tmp, rsq_tmp, zScore_e_tmp,  &lambda, &interested, &ncpus);

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


