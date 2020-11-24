#include <boost/numeric/ublas/symmetric.hpp>
#include  <boost/numeric/ublas/io.hpp>
//using namespace boost::numeric::ublas;

typedef boost::numeric::ublas::symmetric_matrix<double, boost::numeric::ublas::upper>  smatrix_d;
typedef boost::numeric::ublas::symmetric_matrix<float, boost::numeric::ublas::upper>   smatrix_f;
typedef boost::numeric::ublas::symmetric_matrix<int, boost::numeric::ublas::upper>     smatrix_i;
typedef boost::numeric::ublas::symmetric_matrix<short, boost::numeric::ublas::upper>     smatrix_s;
typedef boost::numeric::ublas::symmetric_matrix<char, boost::numeric::ublas::upper>     smatrix_c;

template<class T> void saveData (double ld_r, long i, long j, T* storage, long dim);
template<class T> double readMatrix(T* storage,  long dim, long i, long j );
template<class T> T* createStorage(long  maxDim);
template<class T> void deleteStorage(T* data);
template<class T> void resize (T* data, long newSize);
template<class T> void setDiag(T* data, long dim);
template<class T> uint moveKeepProtect(T* LD, uint arrSize, uint currentDim, uint keepFromIdx);
template<class T> 
void findDup(T* result, double rThreshold, std::vector<int>& dupBearer, std::vector<double>& corABS, std::vector<int>& sign);




