//#include <cstdio>
#include <omp.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>

class Genotypes
{
public:
    std::vector<std::vector<int> > genoMatrix;
};


//******************************************
//  For parsing the 1D vector into table,
//    given the data, geno and size,
//******************************************
class GenoTable
{
public:
    int* geno; // This is 1D vector converted from 2D table, by  row or by  col
    int* size; // nrow * ncol
    std::vector<std::vector<int> >  genoTab;
    inline GenoTable(    int* geno, int* size, std::string method)
    {
        this->geno    = geno;
        this->size    = size;
        this->genoTab =  std::vector<std::vector<int> > (size[0], std::vector<int> ( size[1], 0));
        if(method == "byCol")
        {
#pragma omp parallel for 
            for(int i = 0; i < genoTab.size(); i ++)
                for(int j = 0; j < genoTab[0].size(); j ++)
                    this->genoTab[i][j] = geno[i + j * size[0]];
        }
        else if(method == "byRow")
            printf("Read genotype byRow is not implemented.");
    };
    inline void print()
    {
        printf("[info] nrow = %d, ncol = %d, \n", genoTab.size(), genoTab[0].size());
        for(int i = 0; i < genoTab.size(); i ++)
            for(int j = 0; j < genoTab[0].size(); j ++)
                printf("%d\t", genoTab[i][j]);
            printf("\n");

    };

};


double var(std::vector<int>& geno)
{
    //E( (x - E(x))^2 )
    double E_x = 0;
    double var_x = 0;

#pragma omp parallel for reduction(+:E_x)
    for (int i =0; i < geno.size(); i ++)
        E_x  += geno[i];
    E_x /= geno.size();


#pragma omp parallel for reduction(+:var_x)
    for (int i =0; i < geno.size(); i ++)
        var_x  += (geno[i] - E_x) * (geno[i] - E_x);
    var_x /= (geno.size()-1);

    return var_x;
}


double covar(std::vector<int>& geno1,  std::vector<int>& geno2)
{
    //(x1 -E(x1)* (x2 - E(x2))
    double E_x1 = 0, E_x2 = 0, E_x1_x2 = 0;

#pragma omp parallel for reduction(+:E_x1,E_x2)
    for (int i =0; i < geno1.size(); i ++)
        E_x1  += geno1[i] , E_x2  += geno2[i] ;
    E_x1 /= geno1.size();
    E_x2 /= geno2.size();


#pragma omp parallel for reduction(+:E_x1_x2)
    for (int i =0; i < geno1.size(); i ++)
        E_x1_x2 += (geno1[i] - E_x1) *  (geno2[i] - E_x2);
    E_x1_x2 /= (geno1.size()-1);
    return E_x1_x2;
}

void dim(std::vector<std::vector<double> >  genoTab)
{
    printf("nrow = %d, ncol = %d, \n", genoTab.size(), genoTab[0].size());
}


// std::vector<double> operator- (std::vector<int>& geno1,  std::vector<int>& geno2)
// {
//     //(x1 -E(x1)* (x2 - E(x2))
//     vector<double> res (geno1.size(), 0);
//     for (int i =0; i < geno1.size(); i ++)
//         res[i] = geno1[i] - geno2[i];
//     return res;
// }
// 
// 
// std::vector<double> operator * (double beta,  std::vector<int>& geno1)
// {
//     //(x1 -E(x1)* (x2 - E(x2))
//     vector<double> res (geno1.size(), 0);
//     for (int i =0; i < geno1.size(); i ++)
//         res[i] = beta * geno2[i];
//     return res;
// }
//







void  calcLD (GenoTable& gg1s, GenoTable& gg2s, double* r_mat, int nrow, int ncol)
{
    // std::vector<std::vector<double> > r_mat = std::vector< std::vector<double> >( gg1s.genoTab.size(), std::vector<double>(gg2s.genoTab.size(), 0));
   
    printf("[info] LD dim: nrow = %d, ncol = %d\n", nrow, ncol);
    double* varMat  = new double[nrow];


#pragma omp parallel for
    for ( int i =0; i < nrow; i ++)
        varMat[i] = var(gg1s.genoTab[i]);

#pragma omp parallel for
    for ( int i =0; i < nrow; i ++)
        for ( int j =i; (j< ncol) ; j ++)
        {
            if (j > i + 1000) continue;
            double covg1g2 = covar(gg1s.genoTab[i], gg2s.genoTab[j]) ;
            double scaler  =(varMat[i])*(varMat[j]);
            r_mat[i*ncol + j]     = covg1g2 ==0 ? 0:covg1g2 / sqrt(scaler);
        }

    delete[] varMat;
            
    return;
}

void  calcLD (std::vector<std::vector<int> >& ggs, double* r_mat, int nrow, int ncol)
{
    // std::vector<std::vector<double> > r_mat = std::vector< std::vector<double> >( gg1s.genoTab.size(), std::vector<double>(gg2s.genoTab.size(), 0));
   
#pragma omp parallel for
    for ( int i =0; i < nrow; i ++)
        for ( int j =i; j < ncol; j ++)
        {
            double covg1g2 = covar(ggs[i], ggs[j]) ;
            double scaler  =(var(ggs[i])* var(ggs[j]));
            r_mat[i*ncol + j]     = covg1g2 ==0 ? 0:covg1g2 / sqrt(scaler);
        }
    return;
}






void  calcBeta (GenoTable& gg1s, GenoTable& gg2s, double* r_mat, int nrow, int ncol)
{
    // std::vector<std::vector<double> > r_mat = std::vector< std::vector<double> >( gg1s.genoTab.size(), std::vector<double>(gg2s.genoTab.size(), 0));
   
    for ( int i =0; i < nrow; i ++)
        for ( int j =0; j < ncol; j ++)
        {
            double covg1g2 = covar(gg1s.genoTab[i], gg2s.genoTab[j]) ;
            double scaler  = var(gg2s.genoTab[j]);
            r_mat[i*ncol + j]     = covg1g2 /(scaler);
        }
    return;
}


// standard error of Beta (regression coefficient)
//   
// beta (estimated_from_least_square) =  cov(x,y) / (n * var(x))
//
// var(beta) = ( var(y) - var(beta * x) ) / (n * var(x))
//
//
void  calcBetaSe (GenoTable& gg1s, GenoTable& gg2s, double* r_mat, int nrow, int ncol)
{
    // std::vector<std::vector<double> > r_mat = std::vector< std::vector<double> >( gg1s.genoTab.size(), std::vector<double>(gg2s.genoTab.size(), 0));
    int nn = gg2s.genoTab[0].size();
    for ( int i =0; i < nrow; i ++)
    {
        double var_y = var(gg1s.genoTab[i]);
        for ( int j =0; j < ncol; j ++)
        {
            double var_x = var(gg2s.genoTab[j]);
            double beta    =  covar(gg1s.genoTab[i],      gg2s.genoTab[j]) / var_x;
            r_mat[i*ncol + j]     = (var_y - beta * beta * var_x) / (nn * var_x);
        }
    }
    return;
}


// This load the output (plink.traw) from plink2 --recode A-transpose
std::vector<std::vector<int> >  readPlinkLines(char genoFile[], int N, int M, int* snpToAvert)
{
    std::vector<std::vector<int> > markers  (M, std::vector<int>(N, 0));
    std::vector<int> aMarker (N, 0);
    std::ifstream   fin (genoFile);
    std::string header;
    getline (fin,  header);
    std::string chrID, rsID, cM, refAllele, altAllele;
    long int pos;

    int geno;
    for (unsigned i =0; i < M; i ++)
    {
        fin >> chrID >> rsID >> cM >> pos >> refAllele >> altAllele;
        int sum = 0;
        int sign = 1;
        if ( snpToAvert[i] !=0 ) {sign = -1;}
        for (unsigned j =0; j < N; j ++)
        {
            fin >> geno;
            sum +=geno * sign + snpToAvert[i];
            markers[i][j] = geno * sign + snpToAvert[i];
        }
    }
    fin.close();
    return markers;

}

extern "C" 
{
    // -------------------------------------------------------------------------
    // LD between two sets of markers, given the genoSet1 and genoSet2
    // matrices.
    int calcLD_geno (int* geno1, int* geno2, int* size1, int* size2, int* ncpus, double* r_mat)
    {
        //int nProcessors = omp_get_max_threads();
        int nProcessors = omp_get_num_procs();
        if(nProcessors > *ncpus) nProcessors = *ncpus;
        omp_set_num_threads( nProcessors );
        printf("[info] Calculating geno based on %d cpus \n", nProcessors);
        GenoTable gg1s = GenoTable(geno1, size1, "byCol");
        //gg1s.print();
        GenoTable gg2s = GenoTable(geno2, size2, "byCol");
        calcLD (gg1s, gg2s, r_mat, gg1s.genoTab.size(), gg2s.genoTab.size());
        return 0;
    }

    // ---------------------------------------------------------------------------
    // LD within one set of markers, given  the file that stores the genotypes
    // (The file is generated from plink2 --recode transpose-A).
    int calcLD_geno2 (char* genoFile[], int* N, int* M,int* genoToAvert,  int* ncpus, double* r_mat)
    {
        //int nProcessors = omp_get_max_threads();
        int nProcessors = omp_get_num_procs();
        if(nProcessors > *ncpus) nProcessors = *ncpus;
        omp_set_num_threads( nProcessors );
        printf("[info] Calculating geno based on %d cpus \n", nProcessors);
        std::vector<std::vector<int> >  ggs = readPlinkLines(*genoFile, *N, *M, genoToAvert);
        calcLD (ggs, r_mat, ggs.size(), ggs.size());
        return 0;
    }

    int calcBeta_geno (int* geno1, int* geno2, int* size1, int* size2, double* r_mat)
    {
        GenoTable gg1s = GenoTable(geno1, size1, "byCol");
        //gg1s.print();
        GenoTable gg2s = GenoTable(geno2, size2, "byCol");
        calcBeta (gg1s, gg2s, r_mat, gg1s.genoTab.size(), gg2s.genoTab.size());
        return 0;
    }
    int calcBetaSE_geno (int* geno1, int* geno2, int* size1, int* size2, double* r_mat)
    {
        GenoTable gg1s = GenoTable(geno1, size1, "byCol");
        //gg1s.print();
        GenoTable gg2s = GenoTable(geno2, size2, "byCol");
        calcBetaSe (gg1s, gg2s, r_mat, gg1s.genoTab.size(), gg2s.genoTab.size());
        return 0;
    }

    void lookup(char** ids, int* pos, char** value, int* m, int* n, int* res)
    {

        std::map <std::string, int> dict;
        for (int j=0; j < *n; j ++)
            dict[ std::string(value[j]) ] = pos[j];
        for (int j=0; j < *m; j ++)
            res[j] = dict[ std::string(ids[j]) ] ;
        return ;
    }


    void lookup2(int* posToLook, int* pos, char** value, int* m, int* n, int* res)
    {

        std::map <int, int> dict;
        for (int j=0; j < *n; j ++)
            dict[ pos[j] ] = j+1;
        for (int j=0; j < *m; j ++)
            res[j] = dict[ posToLook[j] ] ;
        return ;
    }

};



