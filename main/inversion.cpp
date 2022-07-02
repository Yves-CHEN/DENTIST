#include <functional>
#include "dataTypeDef.h"
#include "utils.h"
#include "stats/FDR.h"
#include <Eigen/Dense> 
// #include <boost/math/distributions/inverse_chi_squared.hpp>
#include <boost/math/distributions/chi_squared.hpp>




template <class T>
T getQuantile(const std::vector<T>& dat, double whichQuantile)
{
     std::vector<T> diff2 = dat;
     std::sort (diff2.begin(), diff2.end());
     T threshold = diff2[ ceil(diff2.size()  * ( whichQuantile))  -1 ] ;
     return (threshold);
}

template <class T>
T getQuantile2(const std::vector<T>& dat, std::vector<uint> grouping, double whichQuantile)
{

    int sum = std::accumulate(grouping.begin(), grouping.end(), 0);

    if(sum <50)
    {
        return 0;
    }

     std::vector<T> diff2;
     for (uint i =0; i < dat.size(); i ++)
     {
         if(grouping[i] ==1) 
             diff2.push_back(dat[i]);
     }
     std::sort (diff2.begin(), diff2.end());
     T threshold = diff2[ ceil(diff2.size()  * ( whichQuantile))  -1 ] ;
     return threshold;
}



template <class T>
std::vector<T>& operator!(std::vector<T>& logicalDat)
{
    for (uint i =0 ; i < logicalDat.size(); i ++)
        logicalDat[i] = 1- logicalDat[i] ;
    return logicalDat;
}

// 
// 
// extern "C"
// {
//     bool calInversion(double* Vi_mkl, double* b, int* size, int* msg, int* ncpus)
//     {
//         int nProcessors = omp_get_max_threads();
//         if(nProcessors > *ncpus) nProcessors = *ncpus;
//         omp_set_num_threads(nProcessors );
//         printf("[info] matrix inversion based on %d cpus \n", nProcessors);
//         int  N =*size;
//         int ipiv[N];
//         struct timespec start, finish;
//         double elapsed;
//         clock_gettime(CLOCK_MONOTONIC, &start);
//         dgesv( &N, &N, Vi_mkl, &N, ipiv, b, &N, msg);
//         clock_gettime(CLOCK_MONOTONIC, &finish);
//         elapsed = (finish.tv_sec - start.tv_sec);
//         elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
//         D(printf("Time taken %f seconds.\n", elapsed););
// 
//         if( *msg > 0 ) {
//                 printf( "[error] The diagonal element of the triangular factor of A,\n" );
//                 printf( "        U(%i,%i) is zero, so that A is singular;\n", *msg, *msg);
//                 printf( "        the solution could not be computed.\n" );
//                 exit(-1);
//         }
//         return 0;
//     }
// }
// 
// 




// bool imputeOneZscore(double* Vi_mkl, double* imputeToTest, double* zscore, double* bb, int* size, double* res, int* cpus)
// {
//     int N = *size;
//     const float lambda = 0.001;
// 
//     double backup[N*N];
// 
// #pragma omp parallel for 
//         for (unsigned i =0; i < N; i ++)
//             for (unsigned j =0; j < N; j ++)
//                 backup[i*N + j] = Vi_mkl[i*N + j] ;
// 
//     for (unsigned i =0; i < N; i ++) Vi_mkl [i*N+ i] += 0.001;
//     int msg = 0;
//     calInversion(Vi_mkl,  bb, size, &msg, cpus);
//     double imputed_se = 0;
//     double beta [N] = {};
//     double imputed = 0;
// #pragma omp parallel for reduction(+:imputed)
//     for (unsigned i =0; i < N; i ++)
//     {
//         for (unsigned j =0; j < N; j ++)
//             beta[i] += bb[i*N + j] * imputeToTest[j];
//         imputed +=  beta[i] * zscore[i];
//     }
//            
// #pragma omp parallel for reduction(+:imputed_se)
//     for (unsigned i =0; i < N; i ++)
//     {
//         imputed_se  +=  beta[i] * imputeToTest[i];
//     }
// 
//     if(imputed_se < 1e-6)
//     {
//         printf("zScore %f, var %f \n", imputed,  imputed_se);
//         exit(-1);
//     }
//     res[0] = imputed ;
//     res[1] = sqrt(imputed_se);
//     return 0;
// }


///
/// 
// bool impVecOfZscore(double* LDmat, double* zScore, int* markerSize, int* winSize, double* resVecImp, double* seImp, int* ncpus)
// {
// 
//     struct timespec start, finish;
//     double elapsed;
//     clock_gettime(CLOCK_MONOTONIC, &start);
//     int halfWinSize = int((*winSize -1) / 2.0);
//     int startIdx    = halfWinSize ;
//     int endIdx      = *markerSize - halfWinSize ;
//     int MM          = halfWinSize * 2 ;
//     // checking winSize and markerSize
//     if(startIdx >= endIdx)
//         printf("[error] The markerSize (%d) is not smaller than winSize (%d)", *markerSize, *winSize ), exit(-1);
//     int nProcessors = omp_get_max_threads();
//     //if(nProcessors > *ncpus) nProcessors = *ncpus;
//     nProcessors = *ncpus;
//     mkl_set_num_threads( *ncpus );
//     printf("[info]  DENTIST based on %d cpus \n", nProcessors);
//     // create subset of matrix
//     double* subLDmat     = new double [ MM * MM] ; // variable array size is only allowed by g++
//     double* subZscore    = new double [ MM ] ; // variable array size is only allowed by g++
//     double* imputeToTest = new double [ MM ] ; // variable array size is only allowed by g++
//     double* bb = (double*) calloc(MM * MM, sizeof(double)); // Diag mat
//     for (unsigned i =0; i < MM; i ++) bb [i*MM+ i] =1;
// 
//     for (unsigned i = startIdx; i < endIdx ; i ++)
//     {
//         for (unsigned j = i-halfWinSize, m =0; j <= i + halfWinSize; j ++, m ++)
//         {
//             if (j == i) { m --; continue; }
//             subZscore    [m] = zScore[j];
//             imputeToTest [m] = LDmat[i*(*markerSize) + j];
//             for (unsigned k = i -halfWinSize, n = 0; k <= i + halfWinSize; k ++, n ++)
//             {
//                 if(k == i) {n --; continue;}
//                 subLDmat[m*MM + n ]  = LDmat[j*(*markerSize) + k ] ;
//             }
//         }
//         int tmpSize = halfWinSize *2;
//         double res[2] = {0, 0};
//         imputeOneZscore(subLDmat, imputeToTest, subZscore, bb, &tmpSize,  res, ncpus);
//         resVecImp[i] = res[0];
//         seImp[i] = res[1];
// 
// #pragma omp parallel for 
//         for (unsigned j =0; j < MM; j ++) 
//             for (unsigned k =0; k < MM; k ++) 
//             {
//                 bb [j*MM+ k] = 0;
//                 if (j == k ) bb  [j*MM+ k] = 1;
//             }
// 
//     }
//     delete[] subLDmat;
//     delete[] subZscore;
//     delete[] imputeToTest;
//     free(bb);
//     clock_gettime(CLOCK_MONOTONIC, &finish);
//     elapsed = (finish.tv_sec - start.tv_sec);
//     elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
//     D(printf("Time taken %f seconds", elapsed););
// 
//     return (0);
// }

bool printMat (double* mat, int nrow, int ncol)
{
    printf("----------------------------------------\n");
    printf("dim: %d\t%d.\n", nrow, ncol);
    for (uint i =0; i < nrow; i ++)
    {
        for (uint k =0; k < ncol; k ++) printf("%f(%d,%d)\t", mat[i*ncol + k], i, k);
        printf("\n");
    }
    printf("----------------------------------------\n");
    return 0;
}


bool matMultiply(double* mat1, int* dimMat1, double* mat2, int* dimMat2, double* res, int* ncpus)
{
    int nProcessors = omp_get_max_threads();
    if(nProcessors > *ncpus) nProcessors = *ncpus;
    omp_set_num_threads(*ncpus);
    printf("[info] matrix multiply based on %d cpus \n", *ncpus);
    int I = dimMat1[0], J = dimMat1[1], K=dimMat2[1];
#pragma omp parallel for
    for (uint i =0; i < I; i ++)
        for (uint k =0; k < K; k ++)
        {
            double sum = 0;
            for (uint j =0; j < J; j ++)
                sum += mat1[i * J + j] * mat2[k  + j* K];
            res[i*K + k] = sum; 

        }
    return 0;
}






bool matMultiplyDiag(double* mat1, int* dimMat1, double* mat2, int* dimMat2, double* res, int* ncpus)
{
    int nProcessors = omp_get_max_threads();
    if(nProcessors > *ncpus) nProcessors = *ncpus;
    omp_set_num_threads(*ncpus);
    printf("[info] matrix multiply based on %d cpus \n", *ncpus);
    int I = dimMat1[0], J = dimMat1[1], K=dimMat2[1];
#pragma omp parallel for
    for (uint i =0; i < I; i ++)
    {
        uint k = i;
        double sum = 0;
        for (uint j =0; j < J; j ++)
            sum += mat1[i * J + j] * mat2[k* K  + j];
        res[i] = sum; 
    }
    return 0;
}



bool multiRegress(double* matXY, int* dimMat1, double* matV, int* dimMat2, double* matX, int* dimMat3, double* prediction, int* ncpus)
{
    
    int nProcessors = omp_get_max_threads();
    double* beta = new double[dimMat1[0] * dimMat2[1]];
    int beta_dim[2] = {dimMat1[0], dimMat2[1]};
    double* r_sq = new double [dimMat1[0] ];

    printf("[info] start .....\n \n \n");
    if(nProcessors > *ncpus) nProcessors = *ncpus;
    omp_set_num_threads(*ncpus);
    printf("[info] matrix multiply based on %d cpus \n", *ncpus);
    matMultiply(matXY, dimMat1,  matV,  dimMat2,  beta,        ncpus);
    matMultiply(beta,  beta_dim, matX,  dimMat3,  prediction,  ncpus);
    matMultiplyDiag(beta,  beta_dim, matXY,  dimMat1,  r_sq,  ncpus);
#pragma omp parallel for
    for(uint i = 0; i < dimMat1[0]; i ++)
        prediction[i] = prediction[i] / sqrt(r_sq[i]);
//
    //printMat (res, I, K);
    delete[] beta; 
    delete[] r_sq; 
    return 0;
}






template <class T>
void oneIteration (T* LDmat, uint* matSize, double* zScore, std::vector<uint>& idx, std::vector<uint>& idx2,
        double* imputedZ, double* rsqList, double* zScore_e, uint nSample, float probSVD, int* ncpus)
{
    omp_set_num_threads(*ncpus);
    Eigen::setNbThreads(*ncpus);
    Eigen::initParallel();
    uint K = std::min(uint(idx.size()),nSample) * probSVD;
    //double* LDmatSubset = new double[idx.size() * idx.size()];
    Eigen::MatrixXd LD_it (idx2.size(),idx.size() );
    Eigen::VectorXd  zScore_eigen (idx.size() );
    Eigen::MatrixXd VV (idx.size(),idx.size() );

#pragma omp parallel for
    for (uint i = 0; i < idx2.size(); i ++)
        for (uint k = 0; k < idx.size(); k ++)
            LD_it(i,k) = readMatrix<T>(LDmat, *matSize, idx2[i], idx[k]) ;
            //LD_it(i,k) = LDmat[idx2[i] * (*matSize) + idx[k] ] ;

#pragma omp parallel for
    for (uint i = 0; i < idx.size(); i ++)
        zScore_eigen(i)  = zScore[idx[i] ];

#pragma omp parallel for
    for (uint i = 0; i < idx.size(); i ++)
        for (uint j = 0; j < idx.size(); j ++)
            VV(i,j) = readMatrix<T>(LDmat, *matSize, idx[i], idx[j]) ;
            //VV(i,j) =   LDmat[idx[i] * (*matSize) + idx[j]];

    struct timespec start, finish;
    double elapsed;
    clock_gettime(CLOCK_MONOTONIC, &start);

 

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(VV);
    D(printf("Eigen decomposition is done." ););
    clock_gettime(CLOCK_MONOTONIC, &finish);
    elapsed = (finish.tv_sec - start.tv_sec);
    elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

    D(printf("[info] Time elapsed is %f. \n", elapsed ););

    int nRank = es.eigenvectors().rows();
    int nZeros = 0;
    for (int j=0; j<nRank; j++) 
        if ( es.eigenvalues()(j) < 0.0001) nZeros ++;
    nRank = nRank - nZeros;
    if(K>nRank) K = nRank ;
    if (K<=1) { printf("[error] Rank of eigen matrix <=1\n"); exit(-1); }
    Eigen::MatrixXd  ui = Eigen::MatrixXd::Identity(es.eigenvectors().rows(), K);
    Eigen::MatrixXd  wi = Eigen::MatrixXd::Identity(K, K);

#pragma omp parallel for
    for (int m=0; m<K; m++) {
        int j = es.eigenvectors().rows() - m -1;
        for (int i=0; i< es.eigenvectors().rows(); i++) {
            ui(i,m)  = es.eigenvectors() (i,j);
        }
        wi(m,m) = 1/es.eigenvalues()(j);
    }

    Eigen::MatrixXd  beta = LD_it * ui  * wi ;
    Eigen::VectorXd  zScore_eigen_imp = beta * (ui.transpose() * zScore_eigen);
    Eigen::VectorXd  rsq_eigen  = (beta * (ui.transpose() * LD_it.transpose()) ).diagonal();


#pragma omp parallel for
    for (uint i = 0; i < idx2.size(); i ++)
    {
        imputedZ[idx2[i] ] = zScore_eigen_imp (i);
        rsqList [idx2[i]] = rsq_eigen(i);
        if(rsq_eigen(i) >= 1)
        {
            printf("[error] Divividing zero : Rsq = %f \n", rsq_eigen(i));
            exit(-1);
        }
        uint j = idx2[i];
        zScore_e[j] = (zScore[j ] - imputedZ[j] ) /  sqrt( readMatrix<T>(LDmat, *matSize, j, j) - rsqList [j] );
        //zScore_e[j] = (zScore[j ] - imputedZ[j] ) /  sqrt(  LDmat[*matSize * j + j] - rsqList [j] );
        //if(rsqList [idx2[i]]  >0.7) imputedZ[idx2[i]] /= sqrt(rsqList [idx2[i]] );
    }




    //delete[] LDmatSubset;
}



void testSymetry(double* LDmat, uint* markerSize)
{

    uint len = *markerSize;
    for (uint i = 0; i < *markerSize; i ++)
        for (uint j = i+1; j < *markerSize; j ++)
        {

            if(LDmat[i *len +j ] !=LDmat[j *len +i ] )
            {
                printf("unsymetric\n");
                exit(-1);
            }

            
        }
}

//double logPvalueChisq1(double stat)
//{
//    boost::math::inverse_chi_squared_distribution<double> mydist(1);
//    double p = boost::math::cdf(mydist,1/(stat));
//    return ( -log10(p) ) ;
//}
//


double minusLogPvalueChisq2(double stat)
{
    boost::math::chi_squared_distribution<double> mydist(1);
    double p = boost::math::cdf( boost::math::complement(mydist, stat) );
    return ( -log10(p) ) ;
}




// void  regularizeLD (double* LDmat, int markerSize)
// {   
//     Eigen::MatrixXd VV (markerSize,  markerSize );
// 
// 
//     for (uint i = 0; i < markerSize; i ++)
//       for (uint j = 0; j < markerSize; j ++)
//           VV(i,j) =   LDmat[i * (markerSize) + j];
//     Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(VV);
//     int nRank = es.eigenvectors().rows();
//     // int nRank = es.rank();
//     Eigen::MatrixXd  wi = Eigen::MatrixXd::Identity(nRank, nRank);
//     int nZeros = 0;
//     for (int j=0; j<nRank; j++) 
//     {
//         wi(j,j) = es.eigenvalues()(j);
//         if ( es.eigenvalues()(j) < 0.0001)
//         {
//             nZeros ++;
//         }
//     }
//     nRank = nRank - nZeros;
// 
//     printf("Rank %d \n", nRank);
// 
//   
//     for (int j=0; j< int(nRank/2); j++) 
//         wi(j,j) = 0;
// 
//     Eigen::MatrixXd ui = es.eigenvectors();
//     Eigen::MatrixXd aa = ui * wi * ui.transpose();
//     for (uint i = 0; i < markerSize; i ++)
//       for (uint j = 0; j < markerSize; j ++)
//              LDmat[i * (markerSize) + j] = aa(i,j) ;
// }   
///  template <class T>
///  void ImpG(T* LDmat, uint* markerSize, uint* nSample, double* zScore,
///          double* imputedZ, double* rsq, double* zScore_e, double pValueThreshold,  
///          int* interested, int* ncpus)
///  {
///   
///      //double internalControl = 5.45;
///      D(printf("\n------------------------\n");                 );
///      D(printf("--------- DENTIST  ---------\n");               );
///      D(printf("ncpus = %d \n", *ncpus);                        );
///      D(printf("size  = %d \n", *markerSize);                   );
///      D(printf("P-value threshold = %e \n", pValueThreshold);   );
///      int nProcessors = omp_get_max_threads();
///      if(*ncpus < nProcessors) nProcessors = *ncpus;
///          omp_set_num_threads( nProcessors );
///  
///      // init
///      std::vector<size_t> randOrder = generateSetOfNumbers(*markerSize, 10);
///      std::vector<uint> idx;
///      std::vector<uint> idx2;
///      std::vector<size_t> fullIdx = randOrder;
///      for (uint i = 0; i < *markerSize; i ++)
///      {
///          if(randOrder[i] > (*markerSize)/2)
///              idx.push_back(i);
///          else
///              idx2.push_back(i);
///      }
///      std::vector<double> diff;
///  
///      //for (uint t =0; t < 1; t ++)
///      for (uint t =0; t < 8; t ++)
///      {
///          std::vector<uint> idx2_QCed;
///          std::vector<size_t> fullIdx_tmp;
///          double threshold = 0;
///          oneIteration<T>(LDmat, markerSize, zScore, idx, idx2, imputedZ, rsq, zScore_e, *nSample,  ncpus);
///          diff.resize(idx2.size());
///          // for(uint i = 0; i < diff.size(); i ++) diff[i] = fabs(zScore[idx2[i]] - imputedZ[idx2[i]]);
///          for(uint i = 0; i < diff.size(); i ++) diff[i] = fabs(zScore_e[idx2[i]]);
///          threshold = getQuantile <double> (diff , (99/100.0)) ;
///          D(printf("thresh %f \n", threshold););
///          for(uint i = 0; i < diff.size(); i ++)
///          {
///              if( (diff[i] < threshold) ) 
///                  idx2_QCed.push_back(idx2[i]);
///              // else
///              //     idx.push_back(idx2[i]);
///          }
///          oneIteration<T>(LDmat, markerSize, zScore, idx2_QCed, idx, imputedZ, rsq, zScore_e, *nSample, ncpus);
///          D(printf("%d, %d\n", idx.size(), idx2.size()););
///          diff.resize(fullIdx.size());
///  
///          for(uint i = 0; i < diff.size(); i ++) diff[i] = fabs(zScore_e[fullIdx[i]]);
///          threshold = getQuantile <double> (diff , (99.5/100.0)) ;
///          std::vector<double> chisq;
///          for(uint i = 0; i < diff.size(); i ++)
///          {
///   
///              if( !(diff[i] > threshold && logPvalueChisq1(  diff[i] * diff[i] ) >  -log10(pValueThreshold)) ) 
///              {
///                  fullIdx_tmp.push_back(fullIdx[i]);
///                  chisq.push_back( zScore_e[fullIdx[i]] * zScore_e[fullIdx[i]] ); 
///              }
///          }
///          D(cout << "max Chisq = " << *std::max_element(chisq.begin(), chisq.end()) << endl;);
///          assert(chisq.size()>2);
///          double sum = std::accumulate(std::begin(chisq), std::end(chisq), 0.0);
///          double inflationFactor =  sum / chisq.size();
///          D(printf("[Notice] Inflation factor %f \n", inflationFactor);)
///          fullIdx = fullIdx_tmp;
///          randOrder = generateSetOfNumbers(fullIdx.size() , 20000 + t*20000);
///          idx.resize(0);
///          idx2.resize(0);
///          for (uint i = 0; i < fullIdx.size(); i ++)
///          {
///              if(randOrder[i] > (fullIdx.size())/2)
///                  idx.push_back(fullIdx[i]);
///              else
///                  idx2.push_back(fullIdx[i]);
///          }
///  
///          D(printf("%d, %d\n", idx.size(), idx2.size());)
///      }
///  
///      // //  Rescue
///      // diff.resize(fullIdx.size());
///      // for(uint i = 0; i < diff.size(); i ++) diff[i] = fabs(rsq[fullIdx[i]]);
///      // double threshold = getQuantile <double> (diff , (20/100.0)) ;
///      // for(uint i = 0; i < diff.size(); i ++)
///      // {
///      //     if( (diff[i] < threshold) )
///      //         idx.push_back(fullIdx[i]);
///      //     else
///      //         idx2.push_back(fullIdx[i]);
///      // }
///      // oneIteration<T>(LDmat, markerSize, zScore, idx2, idx, imputedZ, rsq, zScore_e,  ncpus);
///  
///  
///  }

template <class T>
void impute(T* LDmat, uint* markerSize, uint* nSample, double* zScore,
        double* imputedZ, double* rsq, double* zScore_e, double pValueThreshold,  
        float propSVD, int* ncpus)
{
    D(printf("\n------------------------\n");                 );
    D(printf("--------- DENTIST  ---------\n");               );
    D(printf("ncpus = %d \n", *ncpus);                        );
    D(printf("size  = %d \n", *markerSize);                   );
    D(printf("P-value threshold = %e \n", pValueThreshold);   );
    omp_set_num_threads( *ncpus );

    // init
    std::vector<size_t> randOrder = generateSetOfNumbers(*markerSize, 10);
    std::vector<uint> idx;
    std::vector<uint> idx2;
    std::vector<size_t> fullIdx = randOrder;
    for (uint i = 0; i < *markerSize; i ++)
    {
        if(zScore[i] == HUGE_VAL)
            idx2.push_back(i);
        else
            idx.push_back(i);
    }
    oneIteration<T>(LDmat, markerSize, zScore, idx, idx2, imputedZ, rsq,
            zScore_e, *nSample, propSVD,  ncpus);

}

template <class T>
void DENTIST(T* LDmat, uint* markerSize, uint* nSample, double* zScore,
        double* imputedZ, double* rsq, double* zScore_e, int* iterID,
        double pValueThreshold,  
        int* interested, float propSVD, bool gcControl, int nIter,
        double groupingPvalue_thresh, int* ncpus)
{
 
    D(printf("\n------------------------\n");                 );
    D(printf("--------- DENTIST  ---------\n");               );
    D(printf("ncpus = %d \n", *ncpus);                        );
    D(printf("size  = %d \n", *markerSize);                   );
    D(printf("P-value threshold = %e \n", pValueThreshold);   );
    int nProcessors = omp_get_max_threads();
    if(*ncpus < nProcessors) nProcessors = *ncpus;
        omp_set_num_threads( nProcessors );

    // init
    int seed = 10;
    std::vector<size_t> randOrder = generateSetOfNumbers(*markerSize, seed);
    std::vector<uint> idx;
    std::vector<uint> idx2;
    std::vector<size_t> fullIdx = randOrder;
    for (uint i = 0; i < *markerSize; i ++)
    {
        if(randOrder[i] > (*markerSize)/2)
            idx.push_back(i);
        else
            idx2.push_back(i);
    }


    std::vector<uint> groupingGWAS (*markerSize,0);
    double groupThresh =  groupingPvalue_thresh;
    for (uint i =0; i < *markerSize; i ++)
    {
        if( minusLogPvalueChisq2(zScore[i]*zScore[i]) > -log10(groupThresh) )
            groupingGWAS [i] =1;
    }




    std::vector<double> diff;
    std::vector<uint> grouping_tmp;

    for (uint t =0; t < nIter; t ++)
    {
        std::vector<uint> idx2_QCed;
        std::vector<size_t> fullIdx_tmp;
        double threshold = 0;
        oneIteration<T>(LDmat, markerSize, zScore, idx, idx2, imputedZ, rsq, zScore_e, *nSample, propSVD,  ncpus);



        diff.resize(idx2.size());
        grouping_tmp.resize(idx2.size());
        for(uint i = 0; i < diff.size(); i ++) { 
            diff[i] = fabs(zScore_e[idx2[i]]);
            grouping_tmp[i] = groupingGWAS[idx2[i]];

        }
        // threshold1 for large Zscores, threshold2 for smaller
        threshold = getQuantile <double> (diff,  (99.5/100.0)) ; 
        double threshold1 = getQuantile2 <double> (diff,grouping_tmp ,  (99.5/100.0)) ; 
        double threshold0 = getQuantile2 <double> (diff,!grouping_tmp , (99.5/100.0)) ;
        if(threshold1 == 0){ threshold1 = threshold; threshold0 = threshold; }

        if(t > nIter -2 )
        {
            threshold0 =  threshold;
            threshold1 =  threshold;
        }


        for(uint i = 0; i < diff.size(); i ++)
        {
            if( grouping_tmp[i]==1 && (diff[i] <= threshold1)   ) 
                idx2_QCed.push_back(idx2[i]);
            if( grouping_tmp[i]==0 && (diff[i] <= threshold0)   ) 
                idx2_QCed.push_back(idx2[i]);

        }
        oneIteration<T>(LDmat, markerSize, zScore, idx2_QCed, idx, imputedZ, rsq, zScore_e, *nSample, propSVD, ncpus);
        diff.resize(fullIdx.size());
        grouping_tmp.resize(fullIdx.size());

        for(uint i = 0; i < diff.size(); i ++) {
            diff[i] = fabs(zScore_e[fullIdx[i]]);
            grouping_tmp[i] = groupingGWAS[fullIdx[i]];
        }
        threshold = getQuantile <double> (diff , (99.5/100.0)) ;
        threshold1 = getQuantile2 <double> (diff,grouping_tmp ,  (99.5/100.0)) ; 
        threshold0 = getQuantile2 <double> (diff,!grouping_tmp , (99.5/100.0)) ;

        if(t > nIter -2 )
        {
            threshold0 =  threshold;
            threshold1 =  threshold;
        }
        if(threshold1 == 0){ threshold1 = threshold; threshold0 = threshold; }
        std::vector<double> chisq;
        for(uint i = 0; i < diff.size(); i ++)
        {
            chisq.push_back( zScore_e[fullIdx[i]] * zScore_e[fullIdx[i]] ); 
        }
        D(cout << "max Chisq = " << *std::max_element(chisq.begin(), chisq.end()) << endl;);
        assert(chisq.size()>2);
        std::nth_element(chisq.begin(), chisq.begin() + chisq.size()/2, chisq.end());
        double inflationFactor = chisq[chisq.size()/2] / 0.46;
        // double sum = std::accumulate(std::begin(chisq), std::end(chisq), 0.0);
        // double inflationFactor =  sum / chisq.size();
        if(gcControl == true)
        {
            D(printf("[Notice] Correction for inflation factor %f \n", inflationFactor);)
        }
        for(uint i = 0; i < diff.size(); i ++)
        {
            if(gcControl == true)
            {
                if( !(diff[i] > threshold && minusLogPvalueChisq2(  diff[i] * diff[i] /inflationFactor  ) >  -log10(pValueThreshold)) ) 
                    fullIdx_tmp.push_back(fullIdx[i]);
            } 
            else
            {
                //if( !(diff[i] > threshold && minusLogPvalueChisq2(  diff[i] * diff[i] ) >  -log10(pValueThreshold)) ) 
                if(  minusLogPvalueChisq2(  diff[i] * diff[i] ) <  -log10(pValueThreshold) ) 
                {
                
                    if(grouping_tmp[i] ==1 && diff[i] <= threshold1)
                    {
                        fullIdx_tmp.push_back(fullIdx[i]);
                        iterID [ fullIdx[i] ] ++;
                    }
                    else if(grouping_tmp[i] ==0 && diff[i] <= threshold0)
                    {
                        fullIdx_tmp.push_back(fullIdx[i]);
                        iterID [ fullIdx[i] ] ++;
                    }

                }
            }
        }
        fullIdx = fullIdx_tmp;
        randOrder = generateSetOfNumbers(fullIdx.size() , 20000 + t*20000);
        idx.resize(0);
        idx2.resize(0);
        for (uint i = 0; i < fullIdx.size(); i ++)
        {
            if(randOrder[i] > (fullIdx.size())/2)
                idx.push_back(fullIdx[i]);
            else
                idx2.push_back(fullIdx[i]);
        }
    }

}
// instantiate T into float
template void DENTIST <float>(float* LDmat, uint* markerSize, uint* nSample, double* zScore,
        double* imputedZ, double* rsq, double* zScore_e, int* iterID,  double pValueThreshold,
        int* interested, float propSVD, bool gcControl, int nIter, double, int* ncpus);
template void  DENTIST<double >(double* LDmat, uint* markerSize, uint* nSample, double* zScore,
        double* imputedZ, double* rsq, double* zScore_e, int* iterID, double pValueThreshold,
        int* interested, float propSVD, bool gcControl, int nIter, double, int* ncpus);

template void DENTIST <smatrix_d >(smatrix_d* LDmat, uint* markerSize, uint* nSample, double* zScore,
        double* imputedZ, double* rsq, double* zScore_e,  int* iterID, double pValueThreshold,
        int* interested, float propSVD, bool gcControl, int nIter, double, int* ncpus);
template void DENTIST <smatrix_f >(smatrix_f* LDmat, uint* markerSize, uint* nSample, double* zScore,
        double* imputedZ, double* rsq, double* zScore_e,  int* iterID, double pValueThreshold,
        int* interested, float propSVD, bool gcControl, int nIter, double, int* ncpus);
template void DENTIST <smatrix_i >(smatrix_i* LDmat, uint* markerSize, uint* nSample, double* zScore,
        double* imputedZ, double* rsq, double* zScore_e,  int* iterID, double pValueThreshold,
        int* interested, float propSVD, bool gcControl, int nIter, double, int* ncpus);




template void oneIteration <float> (float* LDmat, uint* matSize, double* zScore, std::vector<uint>& idx, std::vector<uint>& idx2, double* imputedZ, double* rsqList, double* zScore_e, uint nSample, float propSVD,  int* ncpus);

template void impute<double >(double* LDmat, uint* markerSize, uint* nSample,
        double* zScore, double* imputedZ, double* rsq, double* zScore_e,
        double pValueThreshold, float propSVD, int* ncpus);

template void impute<float>(float* LDmat, uint* markerSize, uint* nSample,
        double* zScore, double* imputedZ, double* rsq, double* zScore_e,
        double pValueThreshold, float propSVD, int* ncpus);

template void impute<smatrix_i>(smatrix_i* LDmat, uint* markerSize, uint* nSample,
        double* zScore, double* imputedZ, double* rsq, double* zScore_e,
        double pValueThreshold, float propSVD, int* ncpus);

template void impute<smatrix_f>(smatrix_f* LDmat, uint* markerSize, uint* nSample,
        double* zScore, double* imputedZ, double* rsq, double* zScore_e,
        double pValueThreshold, float propSVD, int* ncpus);

template void impute<smatrix_d>(smatrix_d* LDmat, uint* markerSize, uint* nSample,
        double* zScore, double* imputedZ, double* rsq, double* zScore_e,
        double pValueThreshold, float propSVD, int* ncpus);






template void oneIteration <double> (double* LDmat, uint* matSize, double* zScore, std::vector<uint>& idx, std::vector<uint>& idx2, double* imputedZ, double* rsqList, double* zScore_e, uint nSample, float propSVD,  int* ncpus);




