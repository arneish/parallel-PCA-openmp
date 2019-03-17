#include <malloc.h>
#include <omp.h>
#include <assert.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#define TOLERANCE 1e-3
// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */

// int compare_eigen(const void *pa, const void *pb)
// {
//     const int *a = pa;
//     const int *b = pb;
//     if (a[0]==b[0])
//         return a[1]-b[1];
//     else
//         return a[0]-b[0];
//     // float fa = eigen_[*(const int*)a];
//     // float fb = eigen_[*(const int*)b];
//     // return (fa<fb)-(fa>fb);
// }
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
    //printf("L2-diff b/w D_T's: %f\n", l2_diff);
    return l2_diff;
}
void transpose(float *M, int m, int n, float *M_T)
{
    int i, j;
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
    int i, j, k;
    for (i = 0; i < m1; i++)
    {
        for (j = 0; j < n2; j++)
        {
            for (k = 0; k < n1; k++)
            {
                sum += M_1[i * n1 + k] * M_2_T[j * m2 + k];
            }
            result[i * n2 + j] = sum;
            sum = 0.0;
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

float l2_norm_diagonal_diff(float *A_next, float *A_current, int P)
{
    float norm, sum_sq = 0.0;
    for (int i = 0; i < P; i++)
    {
        sum_sq += (A_next[i * P + i] - A_current[i * P + i]) * (A_next[i * P + i] - A_current[i * P + i]);
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
    fprintf(stderr, "classical GS over:\n");
    fprintf(stderr, "Printing Q_current:\n");
    print_matrix(Q_current, P, P, 0);
    fprintf(stderr, "Printing R_current:\n");
    print_matrix(R_current, P, P, 0);
}

void compute_V(float **SIGMA, float *D_T, float **U, float **V_T, int N, int P)
{
    //V_T = INV-SIGMA * U_T * M
    float *INV_SIGMA = (float *)calloc(N * P, sizeof(float)); //|=NXP
    for (int i = 0; i < P; i++)
    {
        INV_SIGMA[i * P + i] = 1.0 / ((*SIGMA)[i]);
    }
    printf("\n inv-sigma:\n");
    print_matrix(INV_SIGMA, N, P, 0);
    float *U_T = (float *)malloc(sizeof(float) * P * P);
    transpose(*U, P, P, U_T);
    //first, multiply INV-SIGMA X U_T |=(NXP)
    float *product = (float *)malloc(sizeof(float) * N * P);
    multiply(INV_SIGMA, N, P, U_T, P, P, product);
    //now, multiply product X D_T |=(NXN)
    multiply(product, N, P, D_T, P, N, *V_T);

    printf("\n compute_V:\n");
    print_matrix(*V_T, N, N, 0);
    free(INV_SIGMA);
    free(U_T);
    free(product);
    /*float *U_T = (float *) malloc (sizeof(float)*P*P);
    transpose(*U, P, P, U_T);
    float *U_col = (float *) malloc (sizeof(float)*P);
    float *V_col = (float *) malloc (sizeof(float)*P);
    for (int cols=0; cols<P; cols++) //iterating over cols of U
    {
        memcpy(U_col, U_T+cols*P, sizeof(float)*P);
        multiply(D, P, P, U_col, P, 1, V_col);
        for (int i=0; i<P; i++)
        {
            V_col[i]/=(sqrtf((*SIGMA)[i]));
            //write V_col[i] in V_T's row id:cols
            (*V_T)[cols*P+i] = V_col[i]; 
        }
    }
    free(U_T);
    free(U_col);
    free(V_col);
    */
}

void SVD(int N, int P, float *D, float **U, float **SIGMA, float **V_T)
{
    /* 1.Perform SVD for D_T */
    // Get eigen-values & eigen-vectors for D_T*D
    printf("Printing Matrix D:\n");
    print_matrix(D, N, P, 0);
    float *D_T = (float *)malloc(sizeof(float) * P * N);
    transpose(D, N, P, D_T);
    printf("Printing Matrix D_T:\n");
    print_matrix(D_T, P, N, 0);
    float *A = (float *)calloc(P * P, sizeof(float));   //A=D_T*D|(PxP)
    float *A_T = (float *)calloc(P * P, sizeof(float)); //A_T|(PXP)
    multiply(D_T, P, N, D, N, P, A);
    printf("Printing Matrix A=D_T*D|(PXP)\n");
    print_matrix(A, P, P, 0);

    //begin QR-algorithm for A:
    float *A_current = (float *)malloc(sizeof(float) * P * P);
    memcpy(A_current, A, sizeof(float) * P * P); //PxP; initialised with A_0
    float *E_current = initialize_identity(P);   //PXP; initialised to E_0
    printf("Printing Matrix E_0|(PXP)\n");
    print_matrix(E_current, P, P, 0);

    float *Q_ = (float *)malloc(sizeof(float) * P * P);
    float *R_ = (float *)malloc(sizeof(float) * P * P);
    float diff_norm;
    printf("\n");
    int iter = 0;
    do //convergence condition for QR-algorithm
    {
        printf("iter:%d\n", ++iter);
        transpose(A_current, P, P, A_T);
        classicalGS(A_current, A_T, P, Q_, R_);
        float *A_next = (float *)malloc(sizeof(float) * P * P);
        multiply(R_, P, P, Q_, P, P, A_next);
        float *E_next = (float *)malloc(sizeof(float) * P * P);
        multiply(E_current, P, P, Q_, P, P, E_next);
        fprintf(stderr, "A_current:\n");
        print_matrix(A_current, P, P, 0);
        fprintf(stderr, "A_next:\n");
        print_matrix(A_next, P, P, 0);
        fprintf(stderr, "E_current:\n");
        print_matrix(E_current, P, P, 0);
        fprintf(stderr, "E_next:\n");
        print_matrix(E_next, P, P, 0);
        diff_norm = l2_norm_diagonal_diff(A_next, A_current, P);
        free(A_current);
        free(E_current);
        A_current = A_next;
        E_current = E_next;
        printf("diff_norm: %f, tol:%f\n", diff_norm, TOLERANCE);
    } while (diff_norm > TOLERANCE);
    
    //eigenvalues are diagonals of A_current
    float temp = FLT_MAX;
    printf("\nPrinting singular-values: ");
    for (int i = 0; i < P; i++)
    {
        (*SIGMA)[i] = sqrtf(A_current[i * P + i]);
        if ((*SIGMA)[i] > temp)
        {
            printf("EXCEPTION!\n");
            exit(0);
        }
        temp = (*SIGMA)[i];
        printf("%f ", (*SIGMA)[i]);
    }
    printf("\n");
    printf("\nE: ");
    print_matrix(E_current, P, P, 0);

    //qsort(eigen_, P, sizeof(float), compare);

    //eigenvectors matrix (U for D_T*D) is E_current
    for (int i = 0; i < P; i++)
    {
        for (int j = 0; j < P; j++)
        {
            (*U)[i * P + j] = E_current[i * P + j];
        }
    }
    printf("\n U:\n");
    print_matrix(*U, P, P, 0);
    float *temp_sigma = (float *)calloc(P * N, sizeof(float));
    for (int i = 0; i < P; i++)
    {
        temp_sigma[i * N + i] = (*SIGMA)[i];
    }
    printf("\n SIGMA:\n");
    print_matrix(temp_sigma, P, N, 0);

    //compute V_T
    compute_V(SIGMA, D_T, U, V_T, N, P);
    printf("\n V_T:\n");
    print_matrix(*V_T, N, N, 0);

    /*SVD verification*/
    //D_T == U*SIGMA*V_T
    float *product_one = (float *)malloc(sizeof(float) * P * N);
    multiply(*U, P, P, temp_sigma, P, N, product_one); //U*SIGMA
    float *product_two = (float *)malloc(sizeof(float) * P * N);
    multiply(product_one, P, N, *V_T, N, N, product_two); //(U*SIGMA)*V_T [==A_0]
    free(product_one);
    printf("\nORIGINAL D_T:\n");
    print_matrix(D_T, P, N, 0);
    printf("\nORIGINAL D:\n");
    print_matrix(D, N, P, 0);
    printf("\nVERIFIED D_T:\n");
    print_matrix(product_two, P, N, 0);
    printf("\n A0 = D_TXD: \n");
    print_matrix(A, P, P, 0);
    matrix_similarity(D_T, P, N, product_two);
    free(temp_sigma);
    free(product_two);
}

// /*
// 	*****************************************************
// 		TODO -- You must implement this function
// 	*****************************************************
// */
void PCA(int retention, int N, int P, float *D, float *U, float *SIGMA, float **D_HAT, int *K)
{
    float sum_eigenvalues = 0.0;
    int i;
    for (i = 0; i < P; i++)
    {
        sum_eigenvalues += SIGMA[i];
    }
    *K = 0;
    float retention_ = 0.0;
    i = 0;
    while ((retention_ < retention) && (i < P))
    {
        printf("adding to retention: %f\n", SIGMA[i] / sum_eigenvalues);
        retention_ += SIGMA[i] / sum_eigenvalues;
        (*K)++;
        i++;
    }
    fprintf(stderr, "K: %d, retention_: %f\n", *K, retention_);
    *D_HAT = (float *)malloc(sizeof(float) * N * (*K));
    multiply(D, N, P, U, P, *K, *D_HAT);
    printf("PRINTING D_HAT:\n");
    print_matrix(*D_HAT, N, *K, 0);
}
