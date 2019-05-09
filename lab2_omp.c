#include <malloc.h>
#include <omp.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define TOLERANCE 1e-3

#pragma GCC optimize("Ofast")
#pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx,tune=native")

float matrix_similarity(float *M_1, int m, int n, float *M_2)
{
    float l2_diff=0.0;
    for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; j++)
        {
            l2_diff+=(M_1[i*n+j]-M_2[i*n+j])*(M_1[i*n+j]-M_2[i*n+j]);
        }
    }
    l2_diff = sqrtf(l2_diff);
    printf("L2-diff b/w D_T's: %f\n", l2_diff);
    return l2_diff;
}

void transpose(float *M, int m, int n, float *M_T)
{
    int i, j;
#pragma omp parallel for num_threads(1) private(i, j) collapse(2) schedule(static) 
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            M_T[j * m + i] = M[i * n + j];
        }
    }
}

void multiply(float *M_1, int m1, int n1, float *M_2, int m2, int n2, float *result)
{
    assert(n1 == m2);
    float sum = 0.0;
    //compute M_2_T:
    float *M_2_T = (float *) malloc(sizeof(float)*n2*m2);
    transpose(M_2, m2, n2, M_2_T);
    int i, j, k, temp1, temp2;
#pragma omp parallel for private(i, j, k, sum, temp1, temp2) schedule(static)
    for (i = 0; i < m1; i++)
    {
        temp1 = i*n1; 
        for (j = 0; j < n2; j++)
        {
            sum = 0.0;
            temp2 = j*m2;
            for (k = 0; k < n1; k++)
            {
                sum += M_1[temp1 + k] * M_2_T[temp2 + k];
            }
            result[i * n2 + j] = sum;
        }
    }
    free(M_2_T);
}


float* initialize_identity(int size)
{
    float *I = (float *)calloc(size * size, sizeof(float));
    for (int i = 0; i < size; i++)
    {
        I[i * size + i] = 1.0;
    }
    return I;
}

float l2_norm(float *v_col, int length)
{
    float norm, sum_sq = 0.0;
    for (int i = 0; i < length; i++)
    {
        sum_sq += v_col[i] * v_col[i];
    }
    return norm = sqrtf(sum_sq);
}

float l2_norm_diagonal_diff(float *A_next, float *A_current, int P, float *E_next, float *E_current)
{
    int i,j;
    float norm, sum_sq = 0.0;
    for (i = 0; i < P; i++)
    {
        //sum_sq += fabs(A_next[i * P + i] - A_current[i * P + i]);
        sum_sq += (A_next[i * P + i] - A_current[i * P + i]) * (A_next[i * P + i] - A_current[i * P + i]);
    }
    if (sum_sq>TOLERANCE)
        return norm = sqrtf(sum_sq);

#pragma omp parallel for private(i, j) reduction(+: sum_sq)
    for (i=0; i<P; i++)
    {
        for (j=0; j<P; j++)
        {
            //sum_sq+=fabs(E_next[i*P+j]-E_current[i*P+j]);
            sum_sq+=(E_next[i*P+j]-E_current[i*P+j])*(E_next[i*P+j]-E_current[i*P+j]);
        }
    }
    return norm = sqrtf(sum_sq);
}

void print_matrix(float *A, int M, int N, bool console)
{
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (!console)
                fprintf(stderr, "%f ", A[i * N + j]);
            else
                printf("%f ", A[i * N + j]);
        }
        if (!console)
            fprintf(stderr, "\n");
        else
            printf("\n");
    }
}

void classicalGS(float *A_current, float *A_T, int P, float *Q_current, float *R_current)
{
    /*QR-factorisation of A_current (|=PXP)*/
    float *v_col = (float *)malloc(sizeof(float) * P);
    int col, row, row_;
    float result;
    for (col = 0; col < P; col++)
    {
        memcpy(v_col, A_T + col * P, sizeof(float) * P);
        //v_col = a_col
        for (row = 0; row < col; row++)
        {
            //r[row][col] computation:
            result = 0.0;
            for (row_ = 0; row_ < P; row_++)
            {
                result += (Q_current[row_ * P + row] * (A_T[col * P + row_]));
                //result += (Q_current[row_ * P + row] * (v_col[row_]));
            }
            R_current[row * P + col] = result;
            //v[col] computation:
            for (row_ = 0; row_ < P; row_++)
            {
                v_col[row_] -= R_current[row * P + col] * Q_current[row_ * P + row];
            }
        }
        R_current[col * P + col] = l2_norm(v_col, P);
        for (row = 0; row < P; row++)
        {
            Q_current[row * P + col] = v_col[row] / R_current[col * P + col];
        }
    }
    free(v_col);
}

void modifiedGS(float *A_current, int P, float *Q_current, float *R_current)
{
    /*QR-factorisation of A_current (|=PXP)*/
    float *V = (float *)malloc(sizeof(float) * P*P);
    memcpy(V, A_current, sizeof(float)*P*P);
    int i, j, k;
    float l2_norm=0.0, inner_product=0.0;
    for (i=0; i<P; i++)
    {
        l2_norm=0.0;
        for (j=0; j<P; j++)
        {
            l2_norm+=V[j*P+i]*V[j*P+i];
        }
        l2_norm = sqrtf(l2_norm);
        R_current[i*P+i] = l2_norm;
        for (j=0; j<P; j++)
        {
            Q_current[j*P+i] = V[j*P+i]/l2_norm;
        }
        for (j=i+1; j<P; j++)
        {
            inner_product=0.0;
            for (k=0; k<P; k++)
            {
               inner_product+=Q_current[k*P+i]*V[k*P+j]; 
            }
            R_current[i*P+j] = inner_product;
            for (k=0; k<P; k++)
            {
                V[k*P+j]-=R_current[i*P+j]*Q_current[k*P+i];
            }
        }
    }
    free(V);
}
void compute_V(float **SIGMA, float *D_T, float **U, float **V_T, int N, int P)
{
    //V_T = INV-SIGMA * U_T * M
    float *INV_SIGMA = (float *)calloc(N * P, sizeof(float)); //|=NXP
    for (int i = 0; i < P; i++)
    {
        INV_SIGMA[i * P + i] = 1.0 / ((*SIGMA)[i]);
    }
    float *U_T = (float *)malloc(sizeof(float) * P * P);
    transpose(*U, P, P, U_T);
    //first, multiply INV-SIGMA X U_T |=(NXP)
    float *product = (float *)malloc(sizeof(float) * N * P);
    multiply(INV_SIGMA, N, P, U_T, P, P, product);
    //now, multiply product X D_T |=(NXN)
    multiply(product, N, P, D_T, P, N, *V_T);
    free(INV_SIGMA);
    free(U_T);
    free(product);
}

void SVD(int N, int P, float *D, float **U, float **SIGMA, float **V_T)
{
    /* 1.Perform SVD for D_T */
    // Get eigen-values & eigen-vectors for D_T*D
    float *D_T = (float *)malloc(sizeof(float) * P * N);
    transpose(D, N, P, D_T);
    int i, j;
#pragma omp parallel for private(i, j) collapse(2) schedule(static) 
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < P; j++)
        {
            D_T[j * N + i] = D[i * P + j];
        }
    }

    float *A = (float *)calloc(P * P, sizeof(float));   //A=D_T*D|(PxP)
    float *A_T = (float *)calloc(P * P, sizeof(float)); //A_T|(PXP)
    multiply(D_T, P, N, D, N, P, A);
   
    //begin QR-algorithm for A:
    float *A_current = (float *)malloc(sizeof(float) * P * P);
    memcpy(A_current, A, sizeof(float) * P * P); //PxP; initialised with A_0
    float *E_current = initialize_identity(P);   //PXP; initialised to E_0
    float *Q_ = (float *)calloc(P * P, sizeof(float));
    float *R_ = (float *)calloc(P * P, sizeof(float));
    float diff_norm;
    int iter = 0;
    do //QR-algorithm loop
    {
        /*Both Classical Gram Schmidt & Modified Gram Schmidt procedures have been implemented*/
        //transpose(A_current, P, P, A_T);
        //classicalGS(A_current, A_T, P, Q_, R_);
        modifiedGS(A_current, P, Q_, R_);
        float *A_next = (float *)calloc(P * P, sizeof(float));
        multiply(R_, P, P, Q_, P, P, A_next);
        float *E_next = (float *)calloc(P * P, sizeof(float));
        multiply(E_current, P, P, Q_, P, P, E_next);
        diff_norm = l2_norm_diagonal_diff(A_next, A_current, P, E_next, E_current);
        //printf("diff_norm: %f\n", diff_norm);
        free(A_current);
        free(E_current);
        A_current = A_next;
        E_current = E_next;
        iter++;
    } 
    while(diff_norm > TOLERANCE && iter<6000);
    //printf("num of iters:%d\n", iter);
    
    //eigenvalues are diagonals of A_current
    float temp = FLT_MAX;
    for (int i = 0; i < P; i++)
    {
        (*SIGMA)[i] = sqrtf(A_current[i * P + i]);
        if ((*SIGMA)[i] > temp)
        {
            printf("EXCEPTION!\n");
            exit(0);
        }
        temp = (*SIGMA)[i];
    }

    //eigenvectors matrix (U for D_T*D) is E_current
    for (int i = 0; i < P; i++)
    {
        for (int j = 0; j < P; j++)
        {
            (*U)[i * P + j] = E_current[i * P + j];
        }
    }
   
    // float *temp_sigma = (float *)calloc(P * N, sizeof(float));
    // for (int i = 0; i < P; i++)
    // {
    //     temp_sigma[i * N + i] = (*SIGMA)[i];
    // }
    
    //compute V_T
    compute_V(SIGMA, D_T, U, V_T, N, P);
  
    /*SVD verification*/ //COMMENTED CODE ONLY FOR SVD VERIFICATION
    //D_T == U*SIGMA*V_T
    // float *product_one = (float *)malloc(sizeof(float) * P * N);
    // multiply(*U, P, P, temp_sigma, P, N, product_one); //U*SIGMA
    // float *product_two = (float *)malloc(sizeof(float) * P * N);
    // multiply(product_one, P, N, *V_T, N, N, product_two); //(U*SIGMA)*V_T [==A_0]
    // //printf("\nORIGINAL D_T:\n");
    // //print_matrix(D_T, P, N, 0);
    // //printf("\nORIGINAL D:\n");
    // //print_matrix(D, N, P, 0);
    // //printf("\nVERIFIED D_T:\n");
    // //print_matrix(product_two, P, N, 0);
    // //printf("\n A0 = D_TXD: \n");
    // //print_matrix(A, P, P, 0);
    // matrix_similarity(D_T, P, N, product_two);
    // free(product_one);
    // free(temp_sigma);
    // free(product_two);
}

void PCA(int retention, int N, int P, float *D, float *U, float *SIGMA, float **D_HAT, int *K)
{
    float sum_eigenvalues = 0.0;
    int i,j,k, temp1, temp2;
    float sum;
    for (i = 0; i < P; i++)
    {
        sum_eigenvalues += SIGMA[i]*SIGMA[i];
    }
    *K = 0;
    float retention_ = 0.0;
    i = 0;
    while ((retention_ < retention) && (i < P))
    {
        retention_ += (SIGMA[i]*SIGMA[i] / sum_eigenvalues) * 100;
        (*K)++;
        i++;
    }
    //fprintf(stdout, "K: %d, retention_: %f\n", *K, retention_);
    *D_HAT = (float *)malloc(sizeof(float) * N * (*K));

#pragma omp parallel for private(i, j, k, sum, temp1, temp2) schedule(static)
    for (i=0; i<N; i++)
    {
        temp1 = i*P;
        for (j=0; j<(*K); j++)
        {
            sum = 0.0;
            for (k=0; k<P; k++)
            {
                sum += D[temp1+k]*U[k*P+j];
            }
            (*D_HAT)[i*(*K)+j] = sum;
        }
    }
    //printf("PRINTING D_HAT:\n");
    //print_matrix(*D_HAT, N, *K, 1);
}
