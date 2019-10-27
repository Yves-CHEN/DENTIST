#ifndef FDR_H
#define FDR_H

#include "cdflib.h"
#include <vector>
using namespace std;

template <class T>
T MIN (T aa, T bb)
{
    if(aa < bb) return aa;
    return bb;
}


double pchisq(double x, double df) {
    if (x < 0) return -9;

    double p, q;
    int st = 0; // error variable
    int w = 1; // function variable
    double bnd = 1; // boundary function
    // NCP is set to 0
    cdfchi(&w, &p, &q, &x, &df, &st, &bnd);
    // Check status
    if (st != 0) return -9;
    // Return p-value
    return q;
}


// vector<double> ControlFDR_BH(double  p_value[], int n) 
// {
//    vector<double> fdr(n);
//    vector<pair<double, int>> pval_buf(n);
//    vector<pair<int, int>> indx_buf(n);
//    double min_val = 1.0;
// 
//    for(uint i = 0; i < n; i++ ) pval_buf[i] = make_pair(p_value[i], i);
//    stable_sort(pval_buf.begin(), pval_buf.end(),
//            [](const pair<double,int> a, const pair<double,int> b) {return a.first > b.first; });
//    for(uint i = 0; i < n; i++ ) indx_buf[i] = make_pair(pval_buf[i].second, i);
//    stable_sort(indx_buf.begin(), indx_buf.end());
//    for(uint i = 0; i < n; i++ ) {
//        double c = 0.0; 
//        c = (double) n / (double) (n-i) * pval_buf[i].first;
//        if(c < min_val) min_val = c;
//        fdr[indx_buf[i].second] = MIN(1.0, min_val) ;
//    }
//    return (fdr);
// }

vector<double> ControlFDR_BH(double  p_value[], int n) {
   int i = 0;    double c = 0.0, min_val = 1.0;
   vector<double> fdr(n);
   vector<pair<double, int>> pval_buf(n);
   for( i = 0; i < n; i++ ) pval_buf[i] = make_pair(p_value[i], i);
   stable_sort(pval_buf.begin(), pval_buf.end(), [](const pair<double,int> a, const pair<double,int> b) {return a.first > b.first; });
   for( i = 0; i < n; i++ ) {
       c = (double) n / (double) (n-i) * pval_buf[i].first;
       if(c < min_val) min_val = c;
       fdr[pval_buf[i].second] = MIN(1.0, min_val);
   }
   return (fdr);
}



vector<bool> ControlFDR_BH_marking(double  p_value[], int n) {
  int i = 0;    double c = 0.0, min_val = 1.0;
   vector<bool> fdr(n);
   vector<pair<double, int>> pval_buf(n);
   for( i = 0; i < n; i++ ) pval_buf[i] = make_pair(p_value[i], i);
   stable_sort(pval_buf.begin(), pval_buf.end(), [](const pair<double,int> a, const pair<double,int> b) {return a.first > b.first; });
   for( i = 0; i < n; i++ ) {
       c = (double) n / (double) (n-i) * pval_buf[i].first;
       if(c < min_val) min_val = c;
       fdr[pval_buf[i].second] = MIN(1.0, min_val) >= 0.05;
   }
   return (fdr);
}



void maskSpuriousCor(double LDmat[], int markerSize, int nSample)
{
    int df =1;
    double* pV = new double [markerSize * markerSize];
    for (uint i = 0; i< markerSize * markerSize ; i ++)
        pV [i] = pchisq(LDmat[i] *  LDmat[i] * nSample, df);
    vector<bool>  marked = ControlFDR_BH_marking(pV, markerSize * markerSize);
    vector<int> allMarked;
    double maxV = 0;
    for (uint i = 0; i< markerSize  ; i ++)
    {
        int count = 0;
        for (uint j = 0; j< markerSize  ; j ++)
        {
            if(marked[i*markerSize + j ]){
                if(maxV < fabs(LDmat[i * markerSize +j])) maxV = fabs(LDmat[i * markerSize +j]) ;
                count ++;  LDmat[i*markerSize + j] = 0;
            }
        }
    }


    delete []  pV;


}




#endif

