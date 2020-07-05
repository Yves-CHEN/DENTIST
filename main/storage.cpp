#ifndef STORAGE_H
#define STORAGE_H
#include "storage.h"
#include <vector>
#include <iostream>

template<class T>
void saveData (double ld_r, long i, long j, T* storage, long dim) {
    storage[dim* i + j] = ld_r ;
    storage[dim* j + i] = ld_r ;
}
// accutate at the 2th decimal place
template<>
void saveData <char>(double ld_r, long i, long j, char* storage, long dim) {
    storage[dim* i + j] = char(ld_r * 250);
    storage[dim* j + i] = char(ld_r * 250);
}
// accutate at the 4th decimal place
template<> void saveData <short>(double ld_r, long i, long j, short* storage, long dim) {
    storage[dim* i + j] = short(ld_r * 10000);
    storage[dim* j + i] = short(ld_r * 10000);
}
// accuate at the 8th decimal place
template<> void saveData <int>(double ld_r, long i, long j, int* storage, long dim) {
    storage[dim* i + j] = int(ld_r * 1e8);
    storage[dim* j + i] = int(ld_r * 1e8);
}
// symmetric_matrix <double> saves a half of the storage
template<> 
void saveData < smatrix_d >(double ld_r, long i, long j, smatrix_d* storage, long dim) {
    (*storage) (i,j) =  ld_r;
}
// symmetric_matrix <int> saves a half*half of the storage
template<>
void saveData < smatrix_i >(double ld_r, long i, long j, smatrix_i* storage, long dim) {
    (*storage) (i,j) = int(ld_r * 1e8);  
}

template<> 
void saveData < smatrix_f >(double ld_r, long i, long j, smatrix_f* storage, long dim) {
    (*storage) (i,j) =  ld_r;
}

// ***************************************************************************
//            readMatrix Functions
// ***************************************************************************
template<class T> double readMatrix(T* storage,  long dim, long i, long j ) {
    return double(storage[dim* i + j] );
}
// accutate at the 2th decimal place
template<> double readMatrix<char>(char* storage,  long dim, long i, long j ) {
    return double(storage[dim *i + j ]) / 250;
}
// accutate at the 4th decimal place
template<> double readMatrix<short>(short* storage,  long dim, long i, long j ) {
    return double(storage[dim *i + j ]) / 1e4;
}
// accutate at the 8th decimal place
template<> double readMatrix<int>(int* storage,  long dim, long i, long j ) {
    return double(storage[dim *i + j ]) / 1e8;
}
template<> double readMatrix<smatrix_d>(smatrix_d* storage,  long dim, long i, long j ) {
    return (*storage) (i,j);
}
template<> double readMatrix<smatrix_f >(smatrix_f* storage,  long dim, long i, long j ) {
    return (*storage) (i,j);
}
template<> double readMatrix<smatrix_i >(smatrix_i* storage,  long dim, long i, long j )
{
    return double( (*storage) (i,j) / 1e8 );
}


template<class T>  T* createStorage(long  maxDim){ T* LD = new T[  maxDim*maxDim](); return LD; }
template<> 
smatrix_f* createStorage <smatrix_f >(long maxDim){smatrix_f* LD = new smatrix_f(maxDim); return LD; }
template<> 
smatrix_d* createStorage <smatrix_d > (long maxDim){ smatrix_d* LD = new smatrix_d(maxDim); return LD; }
template<> 
smatrix_i* createStorage <smatrix_i > (long maxDim){
    smatrix_i* LD = new smatrix_i(maxDim);
    return LD;
}

template<class T>  void deleteStorage (T* data) { delete[] data; }
template<> void deleteStorage<smatrix_f > (smatrix_f* data) { delete data; data = NULL; }
template<> void deleteStorage< smatrix_d> (smatrix_d* data) { delete data; data = NULL; }
template<> void deleteStorage< smatrix_i> (smatrix_i* data) { delete data; data = NULL; }


template<class T> void resize(T* data, long newSize) { }

template<> void resize <smatrix_i > (smatrix_i * data, long newSize) { data->resize(newSize); }
template<> void resize <smatrix_f > (smatrix_f * data, long newSize) { data->resize(newSize); }
template<> void resize <smatrix_d > (smatrix_d * data, long newSize) { data->resize(newSize); }

template<class T> 
void setDiag(T* data, long dim) { for(uint i =0; i < dim; i++) saveData(1, i, i, data,dim); }
template<> void setDiag<int>(int* data, long dim) { 
    for(uint i =0; i < dim; i++)
    saveData(1e8, i, i, data,dim);
}
template<>
void setDiag<smatrix_i>(smatrix_i* data, long dim) {
    for(uint i =0; i < dim; i++)
        saveData(1e8, i, i, data,dim);
}


// this different from moveKeep by creating a tmp matrix store the previous dat
//  before coping dat to LD mat. This avoid reading and writing from the same
//  LD mat, which can lead to problems when arrSize < currentDim.
template<class T>
uint moveKeepProtect(T* LD, uint arrSize, uint currentDim, uint keepFromIdx)
{
    if(long(arrSize)  - long(keepFromIdx) <= 0) return 0;
    uint tmp_dim = arrSize  - keepFromIdx;
    T* LD_tmp = createStorage<T> (tmp_dim);
    for (uint i = keepFromIdx, m=0; i < arrSize ; i ++, m ++ )
        for (uint j = keepFromIdx, n = 0; j < arrSize ; j ++, n ++ )
            saveData<T>(readMatrix<T> (LD, arrSize, i, j ), m, n, LD_tmp,tmp_dim);
            //LD_tmp[m *  tmp_dim + n ] = readMatrix (LD, arrSize, i, j );
    uint m = 0;
    resize(LD, currentDim);
    for (m = 0; m < tmp_dim; m ++ )
        for (uint n = 0; n < tmp_dim;  n ++ )
            saveData<T>(readMatrix<T> (LD_tmp,tmp_dim, m, n ), m, n, LD, currentDim);
            // LD[m *  currentDim + n ] = LD_tmp[m * tmp_dim + n];
    deleteStorage<T>( LD_tmp );
    return m;
}

template <class T>
void findDup(T* result, double rThreshold, std::vector<int>& dupBearer, std::vector<double>& corABS, std::vector<int>& sign)
{
    int count = 0;
    unsigned int matSize = dupBearer.size();
    double minValue = 1;
    for(int i =0 ; i < matSize; i ++)
    {
        if(dupBearer[i] != -1)  continue; // otherwise i is not a dup.
        for(int j =i+1 ; j < matSize ; j ++)
        {
            
            if((dupBearer[j]  == -1 ) && fabs(readMatrix(result, matSize, i ,j)) > rThreshold) 
            {
                if(readMatrix(result, matSize, i ,j) < 0) sign[j] = -1;
                corABS[j] = fabs(readMatrix(result, matSize, i ,j));
                dupBearer[j] = count;
            }
            if(fabs(readMatrix(result, matSize, i ,j)) < minValue)
                minValue =  fabs(readMatrix(result, matSize, i ,j));
        }
        count ++;
    }


}






template smatrix_d* createStorage <smatrix_d > ( long maxDim);
template smatrix_f*  createStorage <smatrix_f > ( long maxDim);
template smatrix_i*    createStorage <smatrix_i > ( long maxDim);

template double* createStorage <double>( long maxDim);
template float*  createStorage <float> ( long maxDim);
template int*    createStorage <int>   ( long maxDim);
template void    deleteStorage <int>   (int* );
template void    deleteStorage <float> (float* );
template void    deleteStorage <double>(double* );
template void resize<double>(double* data, long newSize);
template void resize<int>(int* data, long newSize);
template void resize<float>(float* data, long newSize);
template void resize<smatrix_i>(smatrix_i* data, long newSize);
template void resize<smatrix_d>(smatrix_d* data, long newSize);
template void resize<smatrix_f>(smatrix_f* data, long newSize);

template void setDiag<double>(double* data, long);
template void setDiag<int>(int* data, long);
template void setDiag<float>(float* data, long);
template void setDiag<smatrix_i>(smatrix_i* data, long);
template void setDiag<smatrix_d>(smatrix_d* data, long);
template void setDiag<smatrix_f>(smatrix_f* data, long);






template void saveData <double> (double ld_r, long i, long j, double* storage, long dim);
template void saveData <float>  (double ld_r, long i, long j,  float* storage, long dim);
template void saveData <int>    (double ld_r, long i, long j,    int* storage, long dim);


template double readMatrix<double>(double* storage,  long dim, long i, long j );
template double readMatrix<float> (float* storage,  long dim, long i, long j );
template double readMatrix<int>   (int* storage,  long dim, long i, long j );

template uint    moveKeepProtect<double> (double* LD, uint arrSize, uint currentDim, uint keepFromIdx);
template uint    moveKeepProtect<float> (float* LD, uint arrSize, uint currentDim, uint keepFromIdx);
template uint    moveKeepProtect<int> (int* LD, uint arrSize, uint currentDim, uint keepFromIdx);
template uint    moveKeepProtect<smatrix_i> (smatrix_i* LD, uint arrSize, uint currentDim, uint keepFromIdx);
template uint    moveKeepProtect<smatrix_f> (smatrix_f* LD, uint arrSize, uint currentDim, uint keepFromIdx);
template uint    moveKeepProtect<smatrix_d> (smatrix_d* LD, uint arrSize, uint currentDim, uint keepFromIdx);



template void findDup<double>(double* result, double rThreshold, std::vector<int>& dupBearer, std::vector<double>& corABS, std::vector<int>& sign);
template void findDup<float>(float* result, double rThreshold, std::vector<int>& dupBearer, std::vector<double>& corABS, std::vector<int>& sign);
template void findDup<int>(int* result, double rThreshold, std::vector<int>& dupBearer, std::vector<double>& corABS, std::vector<int>& sign);
template void findDup<smatrix_d>(smatrix_d* result, double rThreshold, std::vector<int>& dupBearer, std::vector<double>& corABS, std::vector<int>& sign);
template void findDup<smatrix_f> (smatrix_f* result, double rThreshold, std::vector<int>& dupBearer, std::vector<double>& corABS, std::vector<int>& sign);
template void findDup<smatrix_i>   (smatrix_i* result, double rThreshold, std::vector<int>& dupBearer, std::vector<double>& corABS, std::vector<int>& sign);

#endif //STORAGE_H

