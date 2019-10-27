#include<dataTypeDef.h>



/* Auxiliary routine: printing a matrix */
void print_matrix(const char desc[], int m, int n, double* a, int lda ) {
        int i, j;
        printf( "\n %s\n", desc );
        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
                printf( "\n" );
        }
}


extern "C"
{


    bool calInversion_test1(double* Vi_mkl, int* size, int* msg)
    {
        // MKL's Cholesky decomposition
        int info = 0, int_n = *size;
        int overall = *size * (*size);
        char uplo = 'L';
        // Compute the LU factorization of Vi_mkl
        dpotrf(&uplo, &int_n, Vi_mkl, &int_n, &info);
        for(int i =0; i < *size; i++)
        {
            for(int j =0; j < *size; j++)
                printf("%f \t ", Vi_mkl[i*(*size) + j]);
            printf("\n");
        }
        printf("info : %d \n", info);
    //    *msg  = info;
    //    if (info < 0) throw ("Error: Cholesky decomposition failed. Invalid values found in the matrix.\n");
    //    else if (info > 0) return false;
    //    else 
    //    {
          // Calcualte V inverse
            dpotri(&uplo, &int_n, Vi_mkl, &int_n, &info);

       for(int i =0; i < *size; i++)
        {
            for(int j =0; j < *size; j++)
                printf("%f \t ", Vi_mkl[i*(*size) + j]);
            printf("\n");
        }
        printf("info : %d \n", info);


        printf("info : %d \n", info);
    //        *msg  = info;
    //        if (info < 0) throw ("Error: invalid values found in the varaince-covaraince (V) matrix.\n");
    //        if (info > 0) return false;
    //    }
        
        return true;

    }

    bool calInversion_test2(double* Vi_mkl, int* size, int* msg)
    {

        int N = *size, LDA = *size, info, lwork;
        double wkopt;
        //double* work;
        double w[N];
        /* Executable statements */
        /* Query and allocate the optimal workspace */
        lwork = -1;
        dsyev( "Vectors", "Upper", &N, Vi_mkl, &LDA, w, &wkopt, &lwork, &info );
        lwork = (int)wkopt;
        double* work  = (double*)malloc( lwork*sizeof(double) );
        /* Solve eigenproblem */
        dsyev( "Vectors", "Upper", &N, Vi_mkl, &LDA, w, work, &lwork , &info );

        /* Print eigenvalues */
        //print_matrix( "Eigenvalues", 1, N, w, 1 );
        /* Print eigenvectors */
        //print_matrix( "Eigenvectors (stored columnwise)", N, N, Vi_mkl, LDA );


        //
        // aa_inv = evec%*%diag(eval)%*%t(evec)
        double tmp_row[N];
        for(unsigned int i =0; i < N; i ++)
        {
             
            for(unsigned int k =i; k < N; k ++)
            {
                tmp_row[k] = 0;

                for(unsigned int j =0; j < N; j ++)
                {
                   tmp_row[k] += Vi_mkl[i + j*N] / w[j] *  Vi_mkl[k + j*N];
                }
            }
            for(unsigned int k =0; k < N; k ++)
                Vi_mkl[i + k*N] = tmp_row[k];

        }




        /* Free workspace */
        free( (void*)work );

        return true;

    }
}

